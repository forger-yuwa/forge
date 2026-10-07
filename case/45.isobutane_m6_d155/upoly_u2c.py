"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c: `integral_bl` の実効 `rtol` だけを変える診断 (CFD 0 step、forge は起動しない)。

事前登録 (commit 1f272036、2026-10-07、結果を見る前。諮問 notes/reviews/2026-10-07-upoly-sizing-noise-floor-diagnose.md):
  - 変えるのは `integral_bl` の実効 rtol だけ: A = 1e-6 (今)、B = 1e-10。MOC・k_f・熱条件・平滑化・atol・max_step・逆算の許容差 (1e-9 m)・
    反復の上限 (30 回) は固定。rtol は診断用の引数 (`integral_delta_r(rtol=)`・`solve_rt(integral_rtol=)`) で注入し、
    solve_ivp に実際に渡った値を記録する (YAML には書かない)。
  - 測るもの: r_t = 0.07665396532029795 m を中心に U2 の `noise_floor` と同じ 11 点で、生の δ_r と平滑化後の δ_r、
    物理スロートの半径と出口半径の局所の散らばり (1 次の傾向を除いた max − min)。同じ入力での再評価の一致。
    出口の目標 0.775 m を両腕で固定点反復 (最大 30 回) で解き、返った寸法から同じ腕の精度で生産の壁を作り直して往復の誤差を測る。
  - 判定: B で両半径の局所の散らばり ≤ 1e-10 m かつ出口の往復の誤差 ≤ 1e-9 m → 「今の積分の精度が障害だった」説を支持。
    B でもどちらかを満たさない → 「精度を上げるだけで直る」案を棄却。途中まで改善しただけの結果を PASS にしない。

固定入力は U2・U2b と同じ (`_band_ab/upoly/inputs/u3_A_poly_legacy.problem.yaml`、poly 明示 + legacy)。
11 点の定義は `upoly_verify.noise_floor` と同じ式 (U2b の記録 [同じ中心] の drt_m と一致することを実行時に確かめる)。
往復は `upoly_verify.roundtrip_production` (target exit) と同じ手順に rtol の注入だけを足したもの。

usage: python upoly_u2c.py → _band_ab/upoly/U2c.json (+ U2c_profiles.npz: 11 点の δ_r の分布)
"""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
DESIGN = C.parents[1] / "design"
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, str(C))
import upoly_verify as UV  # noqa: E402  (problem_of・UP・U2_TOL_M)

import forge_design.feedback.deltastar_integral as DI  # noqa: E402
from forge_design.evaluate.runner_axismach import build_physical_wall, design_chain, integral_delta_r, load_problem  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

REG_COMMIT = "1f272036"
RT_C = 0.07665396532029795          # U2b の目標 0.775 の最後の反復の r_t (U2b の雑音の床の中心)
R_EXIT = 0.775
ARMS = {"A": 1e-6, "B": 1e-10}
TOL_SPREAD_M = 1e-10                # 判定: 両半径の局所の散らばり
TOL_RT_M = UV.U2_TOL_M              # 判定: 出口の往復の誤差 (1e-9 m)
PROB = UV.problem_of("u3_A_poly_legacy")
OUT = UV.UP / "U2c.json"
OUT_NPZ = UV.UP / "U2c_profiles.npz"
CHANGED_SRC = ("forge_design/feedback/deltastar_integral.py", "forge_design/evaluate/runner_axismach.py",
               "forge_design/feedback/deltastar_loop.py")

# solve_ivp に実際に渡った rtol の記録 (コード側の自己申告とは別に、呼び出しの引数そのものを見る)
_REAL_SOLVE_IVP = DI.solve_ivp
IVP_CALLS: list = []


def _recording_solve_ivp(*a, **kw):
    IVP_CALLS.append(kw.get("rtol", "<not passed>"))
    return _REAL_SOLVE_IVP(*a, **kw)


DI.solve_ivp = _recording_solve_ivp


def sha_file(f: Path) -> str:
    return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def drt_points(rt: float) -> np.ndarray:
    """`upoly_verify.noise_floor` の 11 点 (同じ式)。"""
    return np.r_[-1e-10, -1e-11, -3e-12, -1e-12, -3 * np.spacing(rt), 0.0, 3 * np.spacing(rt), 1e-12, 3e-12, 1e-11, 1e-10]


def detrend(drt, vals) -> dict:
    """`upoly_verify.noise_floor` と同じ: 1 次の傾向を除いた散らばり (max − min) と rms。"""
    vals = np.asarray(vals, dtype=float)
    slope = np.polyfit(drt, vals, 1)
    res = vals - np.polyval(slope, drt)
    return {"slope": float(slope[0]), "detrended_spread": float(res.max() - res.min()), "detrended_rms": float(np.sqrt(np.mean(res ** 2)))}


def detrend_profile(drt, Y) -> dict:
    """x ごとに 1 次の傾向を除いた散らばり (Y: 点 × x)。最大とその位置。"""
    coef = np.polyfit(drt, Y, 1)
    res = Y - (np.outer(drt, coef[0]) + coef[1])
    sp = res.max(axis=0) - res.min(axis=0)
    i = int(np.argmax(sp))
    return {"max_spread": float(sp[i]), "i_at_max": i, "spread_per_x": sp}


def evaluate(p, d, rt: float, rtol) -> dict:
    """生産と同じ経路 (integral_delta_r → build_physical_wall) で 1 点を評価する。rtol=None は注入しない (既定)。"""
    res, drx, info = integral_delta_r(p, d, p.raw["deltastar_initializer"], scale=rt, rtol=rtol)
    W = build_physical_wall(p, d, rt, delta_r_x=drx, offset="radial")
    x_e = float(W.x_e)
    return {"rt": rt, "raw": np.asarray(res["delta_r_raw_integral"], dtype=float), "smooth": np.asarray(res["delta_r"], dtype=float),
            "x": np.asarray(res["x"], dtype=float), "throat_m": rt * float(W.r_throat), "exit_m": rt * float(W.r(np.r_[x_e])[0]),
            "x_throat_rt": float(W.x_throat), "x_e": x_e, "solve_ivp": dict(res["solve_ivp"]),
            "delta_r_x0_smooth_fn": float(drx(np.r_[0.0])[0]), "delta_r_xe_smooth_fn": float(drx(np.r_[x_e])[0])}


def compare_eval(a: dict, b: dict) -> dict:
    out = {}
    for k in ("raw", "smooth"):
        out[k] = {"bit_identical": bool(np.array_equal(a[k], b[k])), "max_abs_diff_rt": float(np.max(np.abs(a[k] - b[k])))}
    for k in ("throat_m", "exit_m", "x_throat_rt"):
        out[k] = {"bit_identical": a[k] == b[k], "diff": b[k] - a[k]}
    out["all_bit_identical"] = all(v["bit_identical"] for v in out.values())
    return out


def roundtrip_exit(rt: float, R: float, rtol) -> dict:
    """`upoly_verify.roundtrip_production` (target exit) と同じ手順: r_t を問題に入れて生産経路で壁を作り直し、出口半径と目標の差。
    違いは integral_delta_r に同じ腕の rtol を注入することだけ。"""
    p = load_problem(PROB)
    p.spec["r_throat"] = float(rt)
    d = design_chain(p)
    _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"], rtol=rtol)
    W = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    got = rt * float(W.r(np.r_[W.x_e])[0])
    return {"r_t_m": rt, "target_m": R, "achieved_m": got, "diff_m": got - R, "pass": bool(abs(got - R) <= TOL_RT_M),
            "x_throat_rt": float(W.x_throat), "throat_radius_m": rt * float(W.r_throat)}


def floor(arm: str, rtol: float, profiles: dict) -> dict:
    p = load_problem(PROB)
    d = design_chain(p)
    drt = drt_points(RT_C)
    ev = []
    for dd in drt:
        ev.append(evaluate(p, d, RT_C + dd, rtol))
        print(f"  [{arm}] floor drt {dd:+.3e}: throat {ev[-1]['throat_m']!r} exit {ev[-1]['exit_m']!r} nfev {ev[-1]['solve_ivp']['nfev']}", flush=True)
    rts = np.array([e["rt"] for e in ev])
    x = ev[0]["x"]
    if not all(np.array_equal(e["x"], x) for e in ev):
        raise RuntimeError("11 点で δ_r の表の x が違う (r_t に依らないはず)")
    raw = np.array([e["raw"] for e in ev])
    smo = np.array([e["smooth"] for e in ev])
    i0 = int(np.argmin(np.abs(x)))
    out = {"drt_m": drt.tolist(), "r_t_m": rts.tolist(),
           "throat_radius_m": [e["throat_m"] for e in ev], "exit_radius_m": [e["exit_m"] for e in ev],
           "x_throat_rt": [e["x_throat_rt"] for e in ev],
           "spread": {"throat_radius_m": detrend(drt, [e["throat_m"] for e in ev]),
                      "exit_radius_m": detrend(drt, [e["exit_m"] for e in ev])}}
    # 生の δ_r と平滑化後の δ_r: x = 0 (表の補間) と表の終端 (x_F)、実寸 r_t·δ_r [m] と r_t 単位
    pts = {}
    for lab, Y in (("raw", raw), ("smooth", smo)):
        at0 = np.array([np.interp(0.0, x, y) for y in Y])
        atF = Y[:, -1]
        pts[lab] = {"x0_rt_units": detrend(drt, at0), "x0_m": detrend(drt, rts * at0),
                    "xF_rt_units": detrend(drt, atF), "xF_m": detrend(drt, rts * atF),
                    "x0_values_rt": at0.tolist(), "xF_values_rt": atF.tolist()}
        prof = detrend_profile(drt, Y * rts[:, None])
        pts[lab]["profile_m"] = {"max_spread_m": prof["max_spread"], "x_at_max_rt": float(x[prof["i_at_max"]])}
        profiles[f"{arm}_{lab}_spread_per_x_m"] = prof["spread_per_x"]
        profiles[f"{arm}_{lab}"] = Y
    profiles["x_rt"] = x
    profiles[f"{arm}_r_t_m"] = rts
    out["delta_r"] = pts
    out["delta_r_fn_at_x0_and_xe_rt"] = {"x0": [e["delta_r_x0_smooth_fn"] for e in ev], "xe": [e["delta_r_xe_smooth_fn"] for e in ev]}
    out["x_grid"] = {"n": int(len(x)), "x_first": float(x[0]), "x_last": float(x[-1]), "x_near_0": float(x[i0]), "x_e_wall": ev[0]["x_e"]}
    out["solve_ivp_center"] = ev[int(np.flatnonzero(drt == 0.0)[0])]["solve_ivp"]
    out["nfev_range"] = [int(min(e["solve_ivp"]["nfev"] for e in ev)), int(max(e["solve_ivp"]["nfev"] for e in ev))]
    center = ev[int(np.flatnonzero(drt == 0.0)[0])]
    # 同じ入力での再評価: 同じ問題オブジェクトで / 問題を読み直して
    re_same = evaluate(p, d, RT_C, rtol)
    p2 = load_problem(PROB)
    re_fresh = evaluate(p2, design_chain(p2), RT_C, rtol)
    out["reevaluation_same_input"] = {"same_objects": compare_eval(center, re_same), "fresh_load": compare_eval(center, re_fresh)}
    return out, center


def solve(arm: str, rtol: float) -> dict:
    t0 = time.time()
    try:
        r = DL.solve_rt(PROB, R_EXIT, integral_rtol=rtol)
        row = {"converged": True, "r_t_m": r["r_t_m"], "residual_m": r["residual_m"], "n_iter": r["n_iter"],
               "iter_history": r["iter_history"], "delta_r_source": r["delta_r_source"], "physical_throat": r["physical_throat"]}
        rt_ret = r["r_t_m"]
    except DL.SizingNotConverged as e:
        row = {"converged": False, "error": str(e), "iter_history": e.history, "n_iter": len(e.history)}
        rt_ret = None
    row["solve_seconds"] = time.time() - t0
    t1 = time.time()
    if rt_ret is not None:
        row["roundtrip_production"] = roundtrip_exit(rt_ret, R_EXIT, rtol)
    else:
        # 返った寸法が無い (反復の上限で例外) — 記録として最後の反復の寸法で往復を測る (判定には使わない)
        row["roundtrip_production_last_iterate"] = roundtrip_exit(row["iter_history"][-1]["r_t_m"], R_EXIT, rtol)
    row["roundtrip_seconds"] = time.time() - t1
    rp = row.get("roundtrip_production") or row.get("roundtrip_production_last_iterate")
    print(f"  [{arm}] solve_rt 0.775: converged {row['converged']} n_iter {row['n_iter']} r_t {rt_ret!r} roundtrip diff {rp['diff_m']:.3e} m", flush=True)
    return row


def git_info() -> dict:
    def g(*a):
        return subprocess.run(["git", "-C", str(C.parents[1]), *a], capture_output=True, text=True).stdout.strip()
    return {"head": g("rev-parse", "HEAD"), "dirty_design_and_script": g("status", "--porcelain", "--", "design", str(Path(__file__).resolve()))}


def main():
    import scipy
    assert DL.SIZING_MAX_ITER == 30 and DL.SIZING_TOL_M == TOL_RT_M == 1e-9, "登録の反復の上限・許容差と違う"
    u2b = json.loads((UV.UP / "U2b.json").read_text())["rows"]["0.775"]
    # 11 点の定義が U2 の noise_floor と同じであること (U2b の 0.775 の雑音の床は同じ中心で測られている)
    if u2b["iter_history"][-1]["r_t_m"] != RT_C:
        raise RuntimeError("U2b の 0.775 の最後の反復の r_t が RT_C と違う")
    if u2b["noise_floor"]["drt_m"] != drt_points(RT_C).tolist():
        raise RuntimeError("11 点の定義が upoly_verify.noise_floor (U2b の記録) と違う")
    rec = {"item": "U2c 積分の精度だけを変える診断 (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c、CFD 0 step)",
           "registration": {"commit": REG_COMMIT, "plan": "plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §6 U2c",
                            "consult": "notes/reviews/2026-10-07-upoly-sizing-noise-floor-diagnose.md"},
           "criterion": f"B (rtol 1e-10) で物理スロートの半径と出口半径の局所の散らばり (11 点、1 次の傾向を除いた max − min) がともに ≤ {TOL_SPREAD_M:g} m、"
                        f"かつ出口の目標 {R_EXIT} m を固定点反復 (最大 30 回、許容差 1e-9 m) で解いて返った寸法から同じ腕の精度で生産の壁を作り直した"
                        f"出口半径と目標の差 ≤ {TOL_RT_M:g} m → 「今の積分の精度が障害だった」説を支持。どちらかを満たさない → 「精度を上げるだけで直る」案を棄却。"
                        "途中まで改善しただけの結果を PASS にしない",
           "script": str(Path(__file__).resolve()), "script_sha256": sha_file(Path(__file__)),
           "changed_sources_sha256": {s: sha_file(DESIGN / s) for s in CHANGED_SRC}, "git": git_info(),
           "python": sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
           "fixed_inputs": {"problem": str(PROB), "problem_sha256": sha_file(PROB),
                            "deltastar_initializer": load_problem(PROB).raw["deltastar_initializer"],
                            "r_t_center_m": RT_C, "R_exit_target_m": R_EXIT, "drt_m": drt_points(RT_C).tolist(),
                            "sizing_tol_m": DL.SIZING_TOL_M, "sizing_max_iter": DL.SIZING_MAX_ITER,
                            "integral_fixed": "RK45・atol 1e-14・max_step (x1−x0)/400 (integral_bl の既定)、5 次 P-spline (knot 2.0, lam 1.0)"},
           "arms_rtol": ARMS, "arms": {}, "partial": True}
    profiles = {}
    for arm, rtol in ARMS.items():
        print(f"arm {arm} (rtol {rtol:g})", flush=True)
        n0 = len(IVP_CALLS)
        row = {"rtol_injected": rtol}
        t0 = time.time()
        row["floor"], center = floor(arm, rtol, profiles)
        row["floor_seconds"] = time.time() - t0
        n_floor = len(IVP_CALLS) - n0
        if arm == "A":
            # 変更前の記録の再現: A (rtol 1e-6 を明示) の出口半径の 11 点 = U2b (変更前のコード、rtol 未指定) の雑音の床の値
            row["reproduces_U2b_noise_floor"] = {"exit_radius_bit_identical": row["floor"]["exit_radius_m"] == u2b["noise_floor"]["R_m"],
                                                 "spread_A": row["floor"]["spread"]["exit_radius_m"]["detrended_spread"],
                                                 "spread_U2b": u2b["noise_floor"]["detrended_spread_m"]}
            # 注入しない (None) と 1e-6 の明示が同じ
            k0 = len(IVP_CALLS)
            pn = load_problem(PROB)
            e_none = evaluate(pn, design_chain(pn), RT_C, None)
            row["none_vs_explicit_1e-6"] = {**compare_eval(center, e_none), "solve_ivp_rtol_when_none": IVP_CALLS[k0:]}
        t1 = time.time()
        k1 = len(IVP_CALLS)
        row["solve"] = solve(arm, rtol)
        row["solve_and_roundtrip_seconds"] = time.time() - t1
        if arm == "A":
            h_old, h_new = u2b["iter_history"], row["solve"]["iter_history"]
            d_old = (u2b.get("roundtrip_production_last_iterate") or {}).get("diff_m")
            d_new = (row["solve"].get("roundtrip_production_last_iterate") or {}).get("diff_m")
            row["reproduces_U2b_history"] = {"identical": h_old == json.loads(json.dumps(h_new, default=float)),
                                             "n_old": len(h_old), "n_new": len(h_new),
                                             "roundtrip_last_iterate_diff_U2b": d_old, "roundtrip_last_iterate_diff_A": d_new,
                                             "roundtrip_last_iterate_identical": d_old is not None and d_old == d_new}
        calls = IVP_CALLS[n0:]
        if arm == "A":
            calls = IVP_CALLS[n0:n0 + n_floor] + IVP_CALLS[k1:]          # rtol=None の 1 回は別に記録した
        row["solve_ivp_rtol_seen"] = {"n_calls": len(calls), "values": sorted({repr(v) for v in calls}),
                                      "all_equal_injected": all(v == rtol for v in calls),
                                      "solve_rt_record_integral_rtol": (row["solve"].get("delta_r_source") or {}).get("integral_rtol")}
        row["seconds_total"] = time.time() - t0
        rec["arms"][arm] = row
        OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    np.savez(OUT_NPZ, **profiles)
    # 判定 (登録どおり、B だけで決める。A は記録)
    j = {}
    for arm in ARMS:
        r = rec["arms"][arm]
        sp_t = r["floor"]["spread"]["throat_radius_m"]["detrended_spread"]
        sp_e = r["floor"]["spread"]["exit_radius_m"]["detrended_spread"]
        s = r["solve"]
        rt_diff = (s.get("roundtrip_production") or {}).get("diff_m")
        j[arm] = {"throat_radius_spread_m": sp_t, "exit_radius_spread_m": sp_e,
                  "throat_spread_ok": bool(sp_t <= TOL_SPREAD_M), "exit_spread_ok": bool(sp_e <= TOL_SPREAD_M),
                  "solve_converged": s["converged"], "n_iter": s["n_iter"], "r_t_returned_m": s.get("r_t_m"),
                  "roundtrip_diff_m": rt_diff,
                  "roundtrip_last_iterate_diff_m": (s.get("roundtrip_production_last_iterate") or {}).get("diff_m"),
                  "roundtrip_ok": bool(s["converged"] and rt_diff is not None and abs(rt_diff) <= TOL_RT_M),
                  "rtol_effective_ok": r["solve_ivp_rtol_seen"]["all_equal_injected"], "seconds_total": r["seconds_total"]}
        j[arm]["all_conditions"] = bool(j[arm]["throat_spread_ok"] and j[arm]["exit_spread_ok"] and j[arm]["roundtrip_ok"])
    rec["judgement"] = j
    if not (j["A"]["rtol_effective_ok"] and j["B"]["rtol_effective_ok"]):
        rec["verdict"] = "判定不能 (solve_ivp に注入した rtol 以外が渡った呼び出しがある)"
    elif j["B"]["all_conditions"]:
        rec["verdict"] = "支持: B で両半径の局所の散らばり ≤ 1e-10 m かつ出口の往復の誤差 ≤ 1e-9 m (「今の積分の精度が障害だった」説を支持)"
    else:
        miss = [k for k in ("throat_spread_ok", "exit_spread_ok", "roundtrip_ok") if not j["B"][k]]
        rec["verdict"] = f"棄却: B で満たさない条件 {miss} (「精度を上げるだけで直る」案を棄却)"
    rec["partial"] = False
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print("U2c", rec["verdict"], flush=True)


if __name__ == "__main__":
    main()
