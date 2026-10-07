"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d: CFD 前の初期見積もりの許容差 (CFD 0 step、forge は起動しない)。

事前登録 (commit 548778bc、2026-10-07、結果を見る前。諮問 notes/reviews/2026-10-07-upoly-sizing-after-u2c-diagnose.md):
  - 対象は case/45 の今の MOC・k_f・熱条件・`poly`・δ_r の共通経路に限る。rtol 1e-6・平滑化・初期の r_t・反復法・形のゲートは固定。
  - 変えるのは許容差 tol_R_m だけ: A = 1e-9 m、B = スロート 1e-7 m・出口 1e-5 m。CFD 0 step、各最大 30 評価。
  - 未使用の目標: 出口 0.775 × 1.025 m、スロート 0.07676037975901401 × 1.025 m の各 1 点。
  - B の合格条件: 両目標とも 30 評価以内に返り、返った r_t を問題に設定して生産の経路で作り直した壁の絶対残差がそれぞれの許容差以下。
    さらに返った r_t を中心とする U2c と同じ 11 点すべてで、傾向を除く前の絶対の目標残差が同じ許容差以下。全評価で既存の形のゲートを通る。
  - 許容差・反復の上限を実行中に増やさない。最終の残差・実効の許容差・δ_r の出典・評価回数を保存する。
  - 判定: A が失敗し B が全条件を満たす → 「厳密な停止条件が初期見積もりを妨げていた」を支持 / B も失敗 → 「P1 だけで復旧できる」を棄却し、
    P3 (閉包の分岐の発動の計測) へ。追加の緩和はしない / 両腕とも合格 → この目標では厳密な条件が障害だったとは言えない (U2b の FAIL と併記)。

固定入力は U2・U2b・U2c と同じ (`_band_ab/upoly/inputs/u3_A_poly_legacy.problem.yaml`、poly 明示 + legacy、k_f 1.054129…)。
腕 A は tol_R_m = 1e-9 を明示、腕 B は既定 (None) で呼び、関数が記録した実効の許容差が登録値と一致することを確かめる (違えば判定不能)。
rtol は注入しない (`integral_bl` の既定 1e-6)。solve_ivp に実際に渡った rtol を数える。
11 点の定義は `upoly_u2c.drt_points` (= `upoly_verify.noise_floor`) と同じ式で、中心を返った r_t にする。評価は U2c の `evaluate` と同じ
(問題を 1 回読み、`integral_delta_r(scale=r)` → `build_physical_wall`)。往復は `upoly_verify.roundtrip_production` と同じ (spec.r_throat を
書き換えた問題で生産の経路)。11 点の中心と往復の値がビット同一かも記録する。
形のゲート: 作ったすべての物理壁 (逆算の各評価・往復・11 点) で、`poly` の上流のゲート (`upstream_gate.pass`、不合格は構築時に例外) と
`prepare_ns` の物理壁フィルタ (`wall.validate()` が空) を記録し、両方を通ることを条件にする。

usage: python upoly_u2d.py → _band_ab/upoly/U2d.json
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
import upoly_verify as UV  # noqa: E402  (problem_of・UP)
import upoly_u2c as U2C  # noqa: E402  (U2c と同じ 11 点 drt_points・同じ傾向の除き方 detrend)

import forge_design.evaluate.runner_axismach as RA  # noqa: E402
import forge_design.feedback.deltastar_integral as DI  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

drt_points, detrend = U2C.drt_points, U2C.detrend

REG_COMMIT = "548778bc"
R_EXIT_BASE = 0.775
R_THROAT_BASE = 0.07676037975901401          # U2 の基準 (r_t 0.0766539 の生産の壁の物理スロート半径、U2.json targets_m.base)
FAC = 1.025
TARGETS = {"exit": R_EXIT_BASE * FAC, "throat": R_THROAT_BASE * FAC}
ARMS_TOL = {"A": {"exit": 1e-9, "throat": 1e-9}, "B": {"exit": 1e-5, "throat": 1e-7}}
MAX_EVAL = 30
RTOL_EXPECTED = 1e-6
PROB = UV.problem_of("u3_A_poly_legacy")
OUT = UV.UP / "U2d.json"
CHANGED_SRC = ("forge_design/feedback/deltastar_loop.py", "forge_design/feedback/deltastar_integral.py",
               "forge_design/evaluate/runner_axismach.py", "forge_design/geometry/wall_axismach.py")

# --- 記録: solve_ivp に実際に渡った rtol (U2c と同じ方法) -----------------------------------------------------------
# upoly_u2c は import 時に DI.solve_ivp を自分の記録器に差し替えるので、本物 (U2C._REAL_SOLVE_IVP) を包み直す
_REAL_SOLVE_IVP = U2C._REAL_SOLVE_IVP
IVP_CALLS: list = []


def _recording_solve_ivp(*a, **kw):
    IVP_CALLS.append(kw.get("rtol", "<not passed>"))
    return _REAL_SOLVE_IVP(*a, **kw)


DI.solve_ivp = _recording_solve_ivp

# --- 記録: 作ったすべての物理壁の形のゲート ------------------------------------------------------------------------------
# 逆算 (`_fixed_point_rt`) は呼び出しのたびに runner_axismach から build_physical_wall を取り出すので、モジュールの属性を包めば全評価を拾える
_REAL_BUILD = RA.build_physical_wall
WALLS: list = []
PHASE = {"name": None}


def _recording_build(p, d, scale, *a, **kw):
    rec = {"phase": PHASE["name"], "r_t_m": float(scale)}
    try:
        W = _REAL_BUILD(p, d, scale, *a, **kw)
    except Exception as e:  # noqa: BLE001  (ゲート不合格は構築時の例外。記録して投げ直す)
        rec.update(built=False, exception=f"{type(e).__name__}: {e}"[:600], upstream_gate_pass=False, validate_msgs=None)
        WALLS.append(rec)
        raise
    g = getattr(W, "upstream_gate", None)
    msgs = W.validate()
    rec.update(built=True, pw_upstream=W.pw_upstream, upstream_gate_pass=bool(g is not None and g["pass"]),
               max_abs_d2_change_vs_H=(None if g is None else float(g["max_abs_d2_change_vs_H"])),
               max_seam_jump=(None if g is None else float(g["max_seam_jump"])),
               unique_min=(None if g is None else bool(g["unique_min"])), validate_msgs=list(msgs))
    WALLS.append(rec)
    return W


RA.build_physical_wall = _recording_build


def sha_file(f: Path) -> str:
    return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def git_info() -> dict:
    def g(*a):
        return subprocess.run(["git", "-C", str(C.parents[1]), *a], capture_output=True, text=True).stdout.strip()
    return {"head": g("rev-parse", "HEAD"),
            "dirty_design_and_script": g("status", "--porcelain", "--", "design", str(Path(__file__).resolve()))}


def wall_value(W, rt: float, target: str) -> float:
    return rt * (float(W.r_throat) if target == "throat" else float(W.r(np.r_[W.x_e])[0]))


def solve(arm: str, target: str) -> dict:
    """逆算: A は tol_R_m を明示 (1e-9)、B は既定 (None)。反復の上限は既定 (SIZING_MAX_ITER = 30) のまま。"""
    R = TARGETS[target]
    kw = {"tol_R_m": ARMS_TOL[arm][target]} if arm == "A" else {}
    fn = DL.solve_rt if target == "exit" else DL.solve_rt_throat
    PHASE["name"] = f"{arm}/{target}/solve"
    w0 = len(WALLS)
    t0 = time.time()
    try:
        r = fn(PROB, R, **kw)
        row = {"returned": True, "r_t_m": r["r_t_m"], "residual_m": r["residual_m"], "n_eval": r["n_iter"],
               "tol_R_m_effective": r["tol_R_m"], "tol_R_m_source": r["tol_R_m_source"],
               "iter_history": r["iter_history"], "delta_r_source": r["delta_r_source"],
               "physical_throat": r["physical_throat"], "exit_radius_m": r["exit_radius_m"]}
    except DL.SizingNotConverged as e:
        row = {"returned": False, "error": str(e), "iter_history": e.history, "n_eval": len(e.history),
               "tol_R_m_effective": e.tol_R_m, "tol_R_m_source": "(例外のため戻りなし; 実効値は SizingNotConverged.tol_R_m)"}
    except ValueError as e:                      # 形のゲートの不合格 (構築時の例外) など
        row = {"returned": False, "error": f"ValueError: {e}"[:800], "iter_history": None,
               "n_eval": len(WALLS) - w0,
               "tol_R_m_effective": DL._resolve_sizing_tol(kw.get("tol_R_m"), target, cfd_before=True)[0],
               "tol_R_m_source": "(例外のため戻りなし; _resolve_sizing_tol で同じ引数から再現)"}
    row["solve_seconds"] = time.time() - t0
    row["call"] = {"function": fn.__name__, "kwargs": kw, "max_iter": "既定 (SIZING_MAX_ITER)"}
    row["n_walls_built_in_solve"] = len(WALLS) - w0
    return row


def roundtrip(arm: str, target: str, rt: float, tag: str) -> dict:
    """返った r_t を問題に設定し (spec.r_throat を書き換え)、生産の経路 (prepare_ns と同じ: integral_delta_r → build_physical_wall) で
    壁を作り直して目標との絶対残差 (`upoly_verify.roundtrip_production` と同じ手順)。"""
    PHASE["name"] = f"{arm}/{target}/{tag}"
    p = RA.load_problem(PROB)
    p.spec["r_throat"] = float(rt)
    d = RA.design_chain(p)
    _, drx, info = RA.integral_delta_r(p, d, p.raw["deltastar_initializer"])
    W = RA.build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    got = wall_value(W, rt, target)
    tol = ARMS_TOL[arm][target]
    return {"r_t_m": rt, "target_m": TARGETS[target], "achieved_m": got, "abs_residual_m": abs(got - TARGETS[target]),
            "residual_m": got - TARGETS[target], "tol_m": tol, "pass": bool(abs(got - TARGETS[target]) <= tol),
            "x_throat_rt": float(W.x_throat), "throat_radius_m": rt * float(W.r_throat), "exit_radius_m": rt * float(W.r(np.r_[W.x_e])[0]),
            "delta_r_source": {k: info.get(k) for k in ("model", "thermal_bc", "cf_scale", "n_scale", "a_crocco", "closure")}
            | {"route": "integral_delta_r (prepare_ns と同じ経路、rtol 注入なし)", "smooth": info.get("smooth", {}).get("kind")}}


def eleven(arm: str, target: str, rt: float, tag: str) -> dict:
    """返った r_t を中心とする U2c と同じ 11 点で、傾向を除く前の絶対の目標残差 (評価は U2c の evaluate と同じ)。"""
    PHASE["name"] = f"{arm}/{target}/{tag}"
    p = RA.load_problem(PROB)
    d = RA.design_chain(p)
    drt = drt_points(rt)
    vals = []
    for dd in drt:
        r_ = rt + dd
        _, drx, _ = RA.integral_delta_r(p, d, p.raw["deltastar_initializer"], scale=r_)
        W = RA.build_physical_wall(p, d, r_, delta_r_x=drx, offset="radial")
        vals.append(wall_value(W, r_, target))
    vals = np.asarray(vals)
    res = vals - TARGETS[target]
    tol = ARMS_TOL[arm][target]
    return {"center_r_t_m": rt, "drt_m": drt.tolist(), "R_m": vals.tolist(), "residual_m": res.tolist(),
            "max_abs_residual_m": float(np.max(np.abs(res))), "i_at_max": int(np.argmax(np.abs(res))), "tol_m": tol,
            "pass": bool(np.all(np.abs(res) <= tol)), "detrended_spread_info": detrend(drt, vals),
            "center_value_m": float(vals[int(np.flatnonzero(drt == 0.0)[0])])}


def run_row(arm: str, target: str) -> dict:
    t0 = time.time()
    n0 = len(IVP_CALLS)
    w0 = len(WALLS)
    row = {"target": target, "target_m": TARGETS[target], "tol_registered_m": ARMS_TOL[arm][target]}
    row["solve"] = solve(arm, target)
    s = row["solve"]
    if s["returned"]:
        row["roundtrip_production"] = roundtrip(arm, target, s["r_t_m"], "roundtrip")
        row["eleven_points"] = eleven(arm, target, s["r_t_m"], "eleven")
        row["center_vs_roundtrip_bit_identical"] = row["eleven_points"]["center_value_m"] == row["roundtrip_production"]["achieved_m"]
    elif s.get("iter_history"):
        # 返った寸法が無い (反復の上限で例外) — 記録として最後の反復の寸法で往復と 11 点を測る (判定には使わない)
        rt_last = s["iter_history"][-1]["r_t_m"]
        row["roundtrip_production_last_iterate"] = roundtrip(arm, target, rt_last, "roundtrip_last_iterate")
        row["eleven_points_last_iterate"] = eleven(arm, target, rt_last, "eleven_last_iterate")
    walls = WALLS[w0:]
    row["shape_gates"] = {"n_walls": len(walls), "n_built": sum(w["built"] for w in walls),
                          "n_upstream_gate_pass": sum(w["upstream_gate_pass"] for w in walls),
                          "n_validate_empty": sum(1 for w in walls if w["built"] and w["validate_msgs"] == []),
                          "all_pass": bool(walls) and all(w["built"] and w["upstream_gate_pass"] and w["validate_msgs"] == [] for w in walls),
                          "max_abs_d2_change_vs_H": max((w["max_abs_d2_change_vs_H"] or 0.0) for w in walls) if walls else None,
                          "max_seam_jump": max((w["max_seam_jump"] or 0.0) for w in walls) if walls else None,
                          "by_phase": {ph: sum(1 for w in walls if w["phase"] == ph) for ph in sorted({w["phase"] for w in walls})},
                          "failures": [w for w in walls if not (w["built"] and w["upstream_gate_pass"] and w["validate_msgs"] == [])][:10]}
    calls = IVP_CALLS[n0:]
    row["solve_ivp_rtol_seen"] = {"n_calls": len(calls), "values": sorted({repr(v) for v in calls}),
                                  "all_equal_expected": bool(calls) and all(v == RTOL_EXPECTED for v in calls)}
    row["seconds"] = time.time() - t0
    rp, el = row.get("roundtrip_production"), row.get("eleven_points")
    print(f"  [{arm}/{target}] returned {s['returned']} n_eval {s['n_eval']} tol {s['tol_R_m_effective']} r_t {s.get('r_t_m')!r} "
          f"roundtrip {None if rp is None else rp['residual_m']} eleven max {None if el is None else el['max_abs_residual_m']} "
          f"gates {row['shape_gates']['all_pass']} ({row['seconds']:.0f} s)", flush=True)
    return row


def judge_arm(rec: dict, arm: str) -> dict:
    j = {}
    for target in TARGETS:
        r = rec["arms"][arm][target]
        s = r["solve"]
        tol_ok = s["tol_R_m_effective"] == ARMS_TOL[arm][target]
        c = {"returned_within_30": bool(s["returned"] and s["n_eval"] <= MAX_EVAL),
             "roundtrip_ok": bool(s["returned"] and r["roundtrip_production"]["pass"]),
             "eleven_ok": bool(s["returned"] and r["eleven_points"]["pass"]),
             "gates_ok": bool(r["shape_gates"]["all_pass"])}
        j[target] = {**c, "all_conditions": all(c.values()), "tol_effective_matches_registration": bool(tol_ok),
                     "rtol_effective_ok": r["solve_ivp_rtol_seen"]["all_equal_expected"],
                     "n_eval": s["n_eval"], "r_t_returned_m": s.get("r_t_m"), "solver_residual_m": s.get("residual_m"),
                     "roundtrip_residual_m": (r.get("roundtrip_production") or {}).get("residual_m"),
                     "eleven_max_abs_residual_m": (r.get("eleven_points") or {}).get("max_abs_residual_m"),
                     "tol_m": ARMS_TOL[arm][target]}
    j["all_conditions"] = all(j[t]["all_conditions"] for t in TARGETS)
    j["valid"] = all(j[t]["tol_effective_matches_registration"] and j[t]["rtol_effective_ok"] for t in TARGETS)
    return j


def main():
    import scipy
    # 登録の許容差・反復の上限と同じであること (実行中に増やさない)
    assert DL.SIZING_MAX_ITER == MAX_EVAL, "反復の上限が登録 (30) と違う"
    assert DL.SIZING_TOL_M == 1e-9, "共通の SIZING_TOL_M が 1e-9 でない"
    assert (DL.SIZING_TOL_PRE_CFD_EXIT_M, DL.SIZING_TOL_PRE_CFD_THROAT_M) == (ARMS_TOL["B"]["exit"], ARMS_TOL["B"]["throat"]), \
        "CFD 前の既定の許容差が登録 (出口 1e-5・スロート 1e-7) と違う"
    u2 = json.loads((UV.UP / "U2.json").read_text())
    u2c = json.loads((UV.UP / "U2c.json").read_text())
    # 11 点の定義が U2c の記録 (同じ中心) と同じであること
    if drt_points(U2C.RT_C).tolist() != u2c["arms"]["A"]["floor"]["drt_m"]:
        raise RuntimeError("11 点の定義が U2c の記録と違う")
    if PROB != Path(u2c["fixed_inputs"]["problem"]) or sha_file(PROB) != u2c["fixed_inputs"]["problem_sha256"]:
        raise RuntimeError("固定入力の問題が U2c と違う")
    p0 = RA.load_problem(PROB)
    rec = {"item": "U2d CFD 前の初期見積もりの許容差 (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d、CFD 0 step)",
           "registration": {"commit": REG_COMMIT, "plan": "plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §6 U2d",
                            "consult": "notes/reviews/2026-10-07-upoly-sizing-after-u2c-diagnose.md"},
           "criterion": "B (スロート 1e-7 m・出口 1e-5 m): 両目標とも 30 評価以内に返り、返った r_t を問題に設定して生産の経路で作り直した壁の"
                        "絶対残差がそれぞれの許容差以下。さらに返った r_t を中心とする U2c と同じ 11 点すべてで、傾向を除く前の絶対の目標残差が"
                        "同じ許容差以下。全評価で既存の形のゲートを通る。A (1e-9 m) も同じ条件で記録。"
                        "判定: A 失敗・B 全条件 → 「厳密な停止条件が初期見積もりを妨げていた」を支持 / B 失敗 → 「P1 だけで復旧できる」を棄却し P3 へ"
                        " (追加の緩和はしない) / 両腕合格 → この目標では厳密な条件が障害だったとは言えない (U2b の FAIL と併記)",
           "script": str(Path(__file__).resolve()), "script_sha256": sha_file(Path(__file__)),
           "helper_scripts_sha256": {f: sha_file(C / f) for f in ("upoly_u2c.py", "upoly_verify.py")},
           "changed_sources_sha256": {s: sha_file(DESIGN / s) for s in CHANGED_SRC}, "git": git_info(),
           "python": sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
           "fixed_inputs": {"problem": str(PROB), "problem_sha256": sha_file(PROB),
                            "deltastar_initializer": p0.raw["deltastar_initializer"], "pw_upstream": p0.geometry.get("pw_upstream"),
                            "physical_wall_repr": p0.geometry.get("physical_wall_repr"),
                            "r_t_initial_m": float(p0.spec["r_throat"]), "targets_m": TARGETS,
                            "target_bases_m": {"exit": R_EXIT_BASE, "throat": R_THROAT_BASE}, "factor": FAC,
                            "throat_base_equals_U2_base": R_THROAT_BASE == u2["targets_m"]["base"],
                            "problem_same_as_U2c": True, "drt_definition": "upoly_u2c.drt_points (U2c の記録の drt_m と一致を確認済み)、中心 = 返った r_t",
                            "max_eval": MAX_EVAL, "sizing_max_iter": DL.SIZING_MAX_ITER,
                            "integral": "rtol 注入なし (integral_bl の既定 1e-6)・RK45・atol 1e-14・max_step (x1−x0)/400、5 次 P-spline (knot 2.0, lam 1.0)",
                            "iteration": "固定点反復 r_t ← R/r(r_t) (_fixed_point_rt)、初期 r_t = spec.r_throat"},
           "arms_tol_m": ARMS_TOL,
           "arms_call": {"A": "tol_R_m を明示 (1e-9)", "B": "既定 (tol_R_m=None → SIZING_TOL_PRE_CFD_EXIT_M・SIZING_TOL_PRE_CFD_THROAT_M)"},
           "arms": {}, "partial": True}
    t_all = time.time()
    for arm in ARMS_TOL:
        print(f"arm {arm} (tol {ARMS_TOL[arm]})", flush=True)
        rec["arms"][arm] = {}
        for target in TARGETS:
            rec["arms"][arm][target] = run_row(arm, target)
            OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    rec["seconds_total"] = time.time() - t_all
    jA, jB = judge_arm(rec, "A"), judge_arm(rec, "B")
    rec["judgement"] = {"A": jA, "B": jB}
    if not (jA["valid"] and jB["valid"]):
        rec["verdict"] = "判定不能 (実効の許容差が登録値と違う、または solve_ivp に 1e-6 以外の rtol が渡った)"
    elif jB["all_conditions"] and not jA["all_conditions"]:
        rec["verdict"] = "支持: A が失敗し B が全条件を満たす (「厳密な停止条件が初期見積もりを妨げていた」を支持)"
    elif not jB["all_conditions"]:
        miss = {t: [k for k in ("returned_within_30", "roundtrip_ok", "eleven_ok", "gates_ok") if not jB[t][k]] for t in TARGETS}
        rec["verdict"] = f"棄却: B が満たさない条件 {miss} (「P1 だけで復旧できる」を棄却し P3 へ。追加の緩和はしない)"
    else:
        rec["verdict"] = "両腕とも合格: この目標では厳密な条件が障害だったとは言えない (U2b の FAIL と併記)"
    rec["partial"] = False
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print("U2d", rec["verdict"], f"({rec['seconds_total']:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
