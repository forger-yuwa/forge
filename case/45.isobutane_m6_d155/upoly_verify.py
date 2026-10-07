"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6.0 U0〜U3 と、全域 1 本の plan (tooling-nozzle-wall-single-bspline) の
W1・W2・W3・W5・保存 → 復元 → 報告の往復のノット挿入の版でのやり直し (CFD 0 step、forge は起動しない)。判定条件は事前登録どおり。

前提 (_band_ab/upoly/ 以下):
  u0_base_head (wsb_prepare.py を HEAD 1a6f7daa の design/ で、キー無し = ramp) / u0_after_ramp (upoly_prepare.py ramp_explicit) /
  u1_poly_default / u3_A_poly_legacy / u3_B_poly_single_bspline (upoly_prepare.py) /
  u0_solve_rt_before.json (HEAD) / u0_solve_rt_after.json (変更後) (upoly_solve_rt_record.py) /
  tests_before_1a6f7daa.txt・tests_after.txt (design/tests の結果)
U2 は 2 段: `u2solve` (逆算・往復・雑音の床 → u2_solve.json、prepare_ns の往復に要る r_t を表示) → upoly_prepare.py で
u2_prep_* を作る → `u2` (判定)。

usage: [CASE_RUNS=...] FORGE_CONVERTER=... python upoly_verify.py {u0|u1|u2solve|u2|u2b|u3|w1w2|w5|roundtrip|neg}
→ _band_ab/upoly/*.json
"""
import json
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
sys.path.insert(0, str(C))
import wsb_verify as WV  # noqa: E402  (compare_dirs・judge_w3・judge_w2・judge_w5_orig・_tests_summary)

UP = C / "_band_ab" / "upoly"
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
R_T0 = 0.0766539
U2_TOL_M = 1e-9
U3_TOL = {"r": 1e-12, "r1": 1e-10, "r2": 1e-8}
U1_TOL = {"x": 1e-6, "r": 1e-9, "d2": 5e-3, "seam": 1e-8}
W1_TOL = WV.W1_TOL                    # 全域 1 本の plan の W1 (半径 1.3e-7・r′ 1e-7・r″ 1e-5・継ぎ目 1e-8・S6 0.002)
EXPECTED_NEW_INFO_KEYS = {"pw_upstream", "pw_upstream_gate", "sizing"}


def dump(name: str, obj) -> dict:
    (UP / name).write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=float))
    return obj


def problem_of(run: str) -> Path:
    return UP / "inputs" / f"{run}.problem.yaml"


def build(problem: Path, mut=None):
    """問題の写しから prepare_ns と同じ経路で設計壁・δ_r・物理壁 (build_physical_wall) を作る。mut(p) で問題を書き換えてから。"""
    from forge_design.evaluate.runner_axismach import build_physical_wall, design_chain, integral_delta_r, load_problem
    p = load_problem(problem)
    if mut is not None:
        mut(p)
    d = design_chain(p)
    res, drx, info = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    W = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    return p, d, res, drx, W


def piecewise_poly(problem: Path):
    """poly の区分表現 (physical_wall_repr を外した同じ問題) と 1 本の B-spline (同じ壁から作り直す)。"""
    from forge_design.geometry.wall_axismach import SingleBSplinePhysicalWall

    def mut(p):
        p.geometry.pop("physical_wall_repr", None)
    p, d, res, drx, PW = build(problem, mut)
    t0 = time.perf_counter()
    B = SingleBSplinePhysicalWall(PW)
    return p, d, res, drx, PW, B, time.perf_counter() - t0


# ======================================================================================================== U0
def judge_u0(cmp: dict, ia: dict, ib: dict) -> dict:
    """U0: ramp を明示した変更後 (ib) と変更前 (ia) で、物理壁と prepare_info.json の既存の値が完全一致。
    ファイル: prepare_info.json 以外は全部バイト列一致 (nozzle.h5 はデータセット単位と属性)、ファイルの集合も同じ。
    prepare_info.json: 変更前の全キーの値が一致、増えたキーは pw_upstream・pw_upstream_gate・sizing だけ、pw_upstream は ramp の明示。"""
    files = cmp["files"]
    only = [n for n, v in files.items() if "only_in" in v]
    bytes_diff = [n for n, v in files.items() if n not in ("nozzle.h5", "prepare_info.json") and not v.get("identical_bytes", False)]
    h5_ok = cmp["h5"]["all_identical"]
    old_keys_diff = [k for k in ia if ia[k] != ib.get(k, object())]
    new_keys = set(ib) - set(ia)
    pwu = ib.get("pw_upstream")
    ok = (not only) and (not bytes_diff) and h5_ok and (not old_keys_diff) and new_keys == EXPECTED_NEW_INFO_KEYS \
        and pwu == {"value": "ramp", "source": "explicit", "requested": "ramp"} and ib.get("pw_upstream_gate") is None
    return {"pass": bool(ok), "files_only_in_one": only, "files_bytes_diff_except_h5_info": bytes_diff,
            "nozzle_h5_datasets_and_attrs_identical": h5_ok,
            "nozzle_h5_file_bytes_identical": files.get("nozzle.h5", {}).get("identical_bytes"),
            "prepare_info_old_keys_differing": old_keys_diff, "prepare_info_new_keys": sorted(new_keys),
            "prepare_info_pw_upstream": pwu, "prepare_info_pw_upstream_gate_is_none": ib.get("pw_upstream_gate") is None,
            "prepare_info_sizing": ib.get("sizing")}


def judge_solve_rt_ns_after(before: dict, after: dict) -> dict:
    a, b = before["ns_after"], after["ns_after"]
    diff = [k for k in a if a[k] != b.get(k)]
    return {"pass": (not diff) and set(a) <= set(b), "keys_differing": diff, "r_t_before": a["r_t_m"], "r_t_after": b["r_t_m"]}


def run_u0():
    from forge_design.evaluate.runner_axismach import _gam_or_gas, design_chain, integral_delta_r, load_problem
    from forge_design.feedback.deltastar_integral import integral_bl
    base, after = UP / "u0_base_head", UP / "u0_after_ramp"
    cmp = WV.compare_dirs(base, after)
    ia, ib = (json.loads((r / "prepare_info.json").read_text()) for r in (base, after))
    j = judge_u0(cmp, ia, ib)
    sb, sa = (json.loads((UP / f"u0_solve_rt_{w}.json").read_text()) for w in ("before", "after"))
    jn = judge_solve_rt_ns_after(sb, sa)
    tests = WV._tests_summary(UP / "tests_after.txt", UP / "tests_after_logs")
    tests_before = WV._tests_summary(UP / "tests_before_1a6f7daa.txt", UP / "tests_before_logs")
    new_test = tests.get("rows", {}).get("run_pw_upstream_poly_tests.py")
    tests_ok = bool(tests.get("pass")) and new_test is not None and new_test["rc"] == 0
    # キー無しが poly になる (prepare_ns の端から端まで)
    idf = json.loads((UP / "u1_poly_default" / "prepare_info.json").read_text())
    default_ok = idf.get("pw_upstream") == {"value": "poly", "source": "default", "requested": None} \
        and (idf.get("pw_upstream_gate") or {}).get("pass") is True and idf.get("pw_ramp_gate") is None
    # CFD 前の solve_rt の経路の違い (一致は求めない; 同じ r_t での δ_r を旧経路・新経路で記録 — レビューの独立計算の再現)
    p = load_problem(problem_of("u0_after_ramp"))
    d = design_chain(p)
    old = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), R_T0)
    _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"], scale=R_T0)
    x_e = float(d["wall"].x_e)
    rec_path = {"r_t_m": R_T0, "old_integral_bl_direct": {"delta_r_x0_rt": float(np.interp(0.0, old["x"], old["delta_r"])),
                                                         "delta_r_exit_rt": float(old["delta_r"][-1])},
                "new_prepare_ns_path": {"delta_r_x0_rt": float(drx(np.r_[0.0])[0]), "delta_r_exit_rt": float(drx(np.r_[x_e])[0])},
                "review_values": {"x0_old": 0.0010135, "x0_new": 0.0013897, "exit_old": 0.70149, "exit_new": 0.73276, "exit_radius_diff_mm": 2.40}}
    rec_path["exit_radius_diff_mm"] = (rec_path["new_prepare_ns_path"]["delta_r_exit_rt"] - rec_path["old_integral_bl_direct"]["delta_r_exit_rt"]) * R_T0 * 1e3
    out = {"item": "U0 旧経路のビット同一 (plan §6.0)",
           "criterion": "pw_upstream: ramp を明示すると、今の物理壁・prepare_info.json の既存の値が変更前 (HEAD 1a6f7daa) と完全一致 / "
                        "solve_rt の NS 後の経路の結果が変更前と一致 / design/tests の既存テストは FAIL 0 (変更前でも失敗する 3 本を除く) / "
                        "キー無しが poly になる。CFD 前の solve_rt は変えるので一致は求めない (記録だけ)",
           "verdict": "PASS" if (j["pass"] and jn["pass"] and tests_ok and default_ok) else "FAIL",
           "artifacts": j, "solve_rt_ns_after": jn, "tests_after": tests, "tests_before_1a6f7daa": tests_before,
           "new_test_row": new_test, "key_absent_is_poly_prepare_ns": {"pass": default_ok, "pw_upstream": idf.get("pw_upstream"),
                                                                       "gate_pass": (idf.get("pw_upstream_gate") or {}).get("pass")},
           "solve_rt_cfd_before_record": {"before": {k: sb["cfd_before"].get(k) for k in ("r_t_m", "delta_exit_rt", "source")},
                                          "after": {k: sa.get("cfd_before", {}).get(k) for k in ("r_t_m", "delta_exit_rt", "source", "residual_m", "error")}
                                          if isinstance(sa.get("cfd_before"), dict) else sa.get("cfd_before")},
           "delta_r_path_old_vs_new_same_rt": rec_path,
           "compare_files": cmp["files"], "prepare_info_diff": cmp["prepare_info_diff"]}
    return dump("U0.json", out)


# ======================================================================================================== U1
def brute_min(W, lo=-1.0, hi=1.0):
    xx = np.linspace(lo, hi, 400001)
    i = int(np.argmin(W.r(xx)))
    for _ in range(4):
        xx = np.linspace(xx[max(i - 2, 0)], xx[min(i + 2, len(xx) - 1)], 4001)
        i = int(np.argmin(W.r(xx)))
    return float(xx[i])


def judge_u1_identity(PP, PR, extra_x=()) -> dict:
    """[0, x_e] の物理壁が ramp とビット同一 (密な点 + ノット + 継ぎ目、r・r′・r″・r‴)。"""
    xs = np.unique(np.r_[np.linspace(0.0, PP.x_e, 400001), PP.design._spl.t, PP._dr.spline.t, list(extra_x)])
    xs = xs[(xs >= 0.0) & (xs <= PP.x_e)]
    same = {f"d{n}": bool(np.array_equal(PP.r(xs, n), PR.r(xs, n))) for n in range(4)}
    return {"pass": all(same.values()), "n_points": int(len(xs)), "identical": same}


def judge_u1_gate(PP) -> dict:
    """判定側でゲートを測り直す (壁の属性 upstream_gate を信用しない): |Q″ − H″| (係数から厳密 + 密な点)、継ぎ目の跳び (左右の極限)、
    大域最小の位置 (総当たり)、単調性 (密な点の r′ の符号、2,000,001 点)。"""
    from forge_design.geometry.wall import _poly_eval
    w, L = PP.design, PP.L_U
    P = np.polynomial.Polynomial
    d2 = P(np.asarray(PP._q_c) - np.asarray(w.up._c)).deriv(2) / L ** 2
    crit = np.r_[0.0, 1.0, [z.real for z in d2.deriv().roots() if abs(z.imag) <= 1e-12 and 0 < z.real < 1]]
    d2max = float(np.max(np.abs(d2(crit))))
    xs = np.r_[np.linspace(-L, 0.0, 200001)[:-1], -1e-13]
    d2dense = float(np.max(np.abs(PP.r(xs, 2) - w.r(xs, 2))))
    qL = [float(_poly_eval(PP._q_c, -L, 0.0, np.r_[-L], n)[0]) for n in range(3)]
    q0 = [float(_poly_eval(PP._q_c, -L, 0.0, np.r_[0.0], n)[0]) for n in range(3)]
    e0 = [float(w.r(np.r_[0.0], n)[0] + PP._dr(np.r_[0.0], n)[0]) for n in range(3)]
    seams = {"x=-L_U": [abs(qL[0] - PP.r_U), abs(qL[1]), abs(qL[2])], "x=0": [abs(q0[n] - e0[n]) for n in range(3)]}
    seam_max = max(max(v) for v in seams.values())
    xb = brute_min(PP)
    xd = np.linspace(PP.x_in, PP.x_e, 2000001)
    r1 = PP.r(xd, 1)
    before, after = xd < PP.x_throat, xd > PP.x_throat
    mono = {"max_r1_before": float(r1[before].max()), "min_r1_after": float(r1[after].min()),
            "x_min_r1_after": float(xd[after][int(np.argmin(r1[after]))])}
    ok = (d2max <= U1_TOL["d2"] and seam_max <= U1_TOL["seam"] and mono["max_r1_before"] <= 1e-12 and mono["min_r1_after"] >= -1e-12
          and abs(xb - PP.x_throat) <= 2e-6 and PP.r_throat > 0)
    return {"pass": bool(ok), "abs_d2_change_vs_H_max_exact": d2max, "abs_d2_change_vs_H_max_dense": d2dense, "limit_d2": U1_TOL["d2"],
            "seam_jumps": seams, "seam_max": seam_max, "limit_seam": U1_TOL["seam"],
            "throat_brute_force_x": xb, "throat_x": PP.x_throat, "monotone_dense_2M": mono,
            "wall_attr_gate_pass": PP.upstream_gate["pass"], "wall_attr_gate_d2": PP.upstream_gate["max_abs_d2_change_vs_H"]}


def run_u1():
    from forge_design.evaluate.runner_axismach import _gam_or_gas, delta_r_from_table
    from forge_design.geometry.wall_axismach import PhysicalNozzleWall, wall_global_min
    out = {"item": "U1 上流の多項式の形 (plan §6.0、case/45 の単調壁の生産問題、CFD 0 step)", "checks": {}}
    ck = out["checks"]
    p, d, res, drx, PP = build(problem_of("u1_poly_default"))
    pr, dr_, _, drx_r, PR = build(problem_of("u0_after_ramp"))
    S = float(p.spec["r_throat"])
    um = S * 1e6
    same_inputs = bool(np.array_equal(drx.spline.t, drx_r.spline.t) and np.array_equal(drx.spline.c, drx_r.spline.c)
                       and np.array_equal(d["wall"]._spl.t, dr_["wall"]._spl.t) and np.array_equal(d["wall"]._spl.c, dr_["wall"]._spl.c))
    ck["same_design_and_delta_r"] = {"pass": same_inputs, "note": "poly と ramp の問題の写しで設計壁 S と δ_r が同一 (違うのは上流の作り方だけ)"}
    ck["downstream_bit_identical"] = judge_u1_identity(PP, PR)
    dx, drr = PP.x_throat - PR.x_throat, PP.r_throat - PR.r_throat
    ck["throat_vs_ramp"] = {"pass": bool(abs(dx) <= U1_TOL["x"] and abs(drr) <= U1_TOL["r"]), "dx_rt": dx, "dr_rt": drr,
                            "tol": {"x": U1_TOL["x"], "r": U1_TOL["r"]}, "poly": {"x": PP.x_throat, "r": PP.r_throat, "kappa": PP.kappa_throat},
                            "ramp": {"x": PR.x_throat, "r": PR.r_throat, "kappa": PR.kappa_throat},
                            "poly_m": {"x_mm": PP.x_throat * S * 1e3, "r_m": PP.r_throat * S}}
    ck["gate"] = judge_u1_gate(PP)
    # 今の物理壁 (ramp [−11, −6]) との差 (記録): 半径の差と断面積差 100·(Q²/W²−1)
    xs = np.linspace(PP.x_in, 0.0, 240001)[:-1]
    dq = PP.r(xs) - PR.r(xs)
    i = int(np.argmin(dq))
    area = 100.0 * (PP.r(xs) ** 2 / PR.r(xs) ** 2 - 1.0)
    ia = int(np.argmin(area))
    a7 = float(100.0 * (PP.r(np.r_[-7.085])[0] ** 2 / PR.r(np.r_[-7.085])[0] ** 2 - 1.0))
    ck["diff_vs_current_wall"] = {"max_abs_um": float(abs(dq[i]) * um), "signed_um": float(dq[i] * um), "x_rt": float(xs[i]),
                                  "max_positive_um": float(dq.max() * um), "area_diff_pct_at_x_minus_7_085": a7,
                                  "area_diff_pct_min": float(area[ia]), "x_area_min": float(xs[ia]),
                                  "within_1rt_max_abs_um": float(np.abs(dq[xs >= -1.0]).max() * um),
                                  "within_0p2rt_max_abs_um": float(np.abs(dq[xs >= -0.2]).max() * um),
                                  "plan_prior_estimate": "−0.52 mm (x = −7.1)、面積 −0.247 % (x = −7.085)・最小 −0.250 % (x = −6.84)、1 r_t 以内 ≤ 6 µm、0.2 r_t 以内 ≤ 0.06 µm"}
    ck["diff_vs_current_wall"]["pass"] = None
    # prepare_ns の記録との一致 (u1_poly_default)
    idf = json.loads((UP / "u1_poly_default" / "prepare_info.json").read_text())
    ck["prepare_info_matches_rebuilt"] = {"pass": bool(idf["throat_physical"]["x"] == PP.x_throat and idf["throat_physical"]["r"] == PP.r_throat
                                                       and idf["pw_upstream_gate"]["max_abs_d2_change_vs_H"] == PP.upstream_gate["max_abs_d2_change_vs_H"]),
                                          "prepare_info_throat": idf["throat_physical"], "pw_upstream": idf["pw_upstream"]}
    # 一般性: 標準の L_U 3.5・r_U 2.5 (端条件も作り直す: δ_r を同じ経路で作り直す)。通らなければそう記録し、自動で ramp にしない
    def mut35(q):
        # CFD ピンの初期線は凍結源 (L_U 12) の形に縛られるので、縮流部を変える一般性の試験は Hall の初期線で行う
        q.geometry["L_U"], q.geometry["r_inlet"] = 3.5, 2.5
        q.geometry["initial_line"] = "hall"
        q.geometry.pop("initial_line_run", None)
        q.geometry.pop("initial_line_res", None)
    try:
        p35, d35, _, dr35, W35 = build(problem_of("u1_poly_default"), mut35)
        g = W35.upstream_gate
        ck["generality_LU3p5_rU2p5"] = {"pass": None, "gate_pass": g["pass"], "pw_upstream": W35.pw_upstream,
                                        "abs_d2_change_vs_H": g["max_abs_d2_change_vs_H"], "throat": g["throat"],
                                        "monotone": [g["monotone_before"], g["monotone_after"]], "seam_max": g["max_seam_jump"],
                                        "x_in": W35.x_in, "L_U": W35.L_U, "r_U": W35.r_U,
                                        "note": "記録 (判定なし)。δ_r はこの形状で積分法を回し直したもの"}
    except ValueError as e:
        ck["generality_LU3p5_rU2p5"] = {"pass": None, "gate_pass": False, "error": str(e), "note": "poly のゲートで止まった (ramp に切り替えない)"}
    # 大域最小の探索器の負例 (生産の設計壁、δ_r の表の x は生産の積分法の 1500 点)
    xt = res["x"]
    a1 = float(drx(np.r_[0.0], 1)[0])
    grow = lambda x: 1e-6 * np.maximum(x - 10.0, 0.0) ** 3  # noqa: E731
    cases = {
        "delta_r_prime0_positive (生産の δ_r)": lambda x: drx(x),
        "delta_r_prime0_zero (生産の δ_r − δ_r′(0)·x·exp(−x²))": lambda x: drx(x) - a1 * x * np.exp(-x ** 2),
        "delta_r_prime0_negative (生産の δ_r − 2δ_r′(0)·x·exp(−x²))": lambda x: drx(x) - 2 * a1 * x * np.exp(-x ** 2),
        "review_slope_wide (0.0014 − 0.000776·x·exp(−x²) + 下流の伸び; 表の間隔で鈍らない幅)": lambda x: 0.0014 - 0.000776 * x * np.exp(-x ** 2) + grow(x),
        "review_counterexample (0.0014 − 0.000776·x·exp(−(x/0.1)²) + 下流の伸び)": lambda x: 0.0014 - 0.000776 * x * np.exp(-(x / 0.1) ** 2) + grow(x),
        "review_counterexample_exact (0.0014 − 0.000776·x·exp(−(x/0.1)²)、伸びなし)": lambda x: 0.0014 - 0.000776 * x * np.exp(-(x / 0.1) ** 2),
        "extra_downstream_extremum (生産の δ_r − 0.4·exp(−((x−3)/0.3)²))": lambda x: drx(x) - 0.4 * np.exp(-((x - 3.0) / 0.3) ** 2),
    }
    args = (d["wall"], d["wall_inv"], S, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
    neg = {}
    for lab, fun in cases.items():
        f = delta_r_from_table(xt, fun(xt))
        row = {"delta_r_prime0": float(f(np.r_[0.0], 1)[0])}
        try:
            W = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f)
            xb = brute_min(W)
            row.update(built=True, x_t=W.x_throat, brute_x=xb, unique=W.upstream_gate["unique_min"],
                       monotone=[W.upstream_gate["monotone_before"], W.upstream_gate["monotone_after"]], gate_pass=W.upstream_gate["pass"],
                       position_ok=bool(abs(W.x_throat - xb) <= 2e-6))
        except ValueError as e:
            # ゲートで止まった壁は、同じ式の ramp の壁 (x ≥ 0 は同一) と Q を探索器に直接かけて位置を記録する
            W0 = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, ramp=(-11.0, -6.0), upstream="ramp")
            br = np.unique(np.r_[0.0, d["wall"]._spl.t, f.spline.t, f.x_range, W0.x_e])
            br = br[(br >= 0.0) & (br <= W0.x_e)]
            s0 = wall_global_min(W0.r, br)
            row.update(built=False, error=str(e)[:300], downstream_search={k: s0[k] for k in ("x_t", "monotone_after", "min_r1_after", "x_min_r1_after", "unique")})
        neg[lab] = row
    sign_ok = (neg["delta_r_prime0_positive (生産の δ_r)"].get("x_t", 1) < 0
               and abs(neg["delta_r_prime0_zero (生産の δ_r − δ_r′(0)·x·exp(−x²))"].get("x_t", 1)) <= 1e-6
               and neg["delta_r_prime0_negative (生産の δ_r − 2δ_r′(0)·x·exp(−x²))"].get("x_t", -1) > 0)
    pos_ok = all(r.get("position_ok", True) for r in neg.values() if r.get("built"))
    dip = neg["extra_downstream_extremum (生産の δ_r − 0.4·exp(−((x−3)/0.3)²))"]
    dip_ok = (not dip.get("built")) and (not dip["downstream_search"]["monotone_after"]) and 2.0 < dip["downstream_search"]["x_min_r1_after"] < 3.2
    ck["global_min_negatives"] = {"pass": bool(sign_ok and pos_ok and dip_ok), "cases": neg, "sign_ok": sign_ok, "position_ok": pos_ok,
                                  "extra_extremum_detected": dip_ok,
                                  "review_value_x_t": 0.00155091}
    judged = [k for k, v in ck.items() if v.get("pass") is not None]
    out["verdict"] = "PASS" if all(ck[k]["pass"] for k in judged) else "FAIL"
    out["judged"] = judged
    return dump("U1.json", out)


# ======================================================================================================== U2
def _supplier_A(p):
    """腕 A: 旧 solve_rt 型の δ_r の供給 (`integral_bl` を直接、既定の引数 = 未較正・熱条件既定・未平滑化) を壁に渡す関数にする。"""
    from forge_design.evaluate.runner_axismach import _gam_or_gas, delta_r_from_table
    from forge_design.feedback.deltastar_integral import integral_bl

    def sup(d, rt):
        r = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt)
        return delta_r_from_table(r["x"], r["delta_r"]), {"kind": "arm A: integral_bl 直接 (未較正・未平滑化、旧 solve_rt 型)"}
    return sup


def roundtrip_production(problem: Path, rt: float, R: float, target: str = "throat", delta_csv=None) -> dict:
    """往復: 解いた r_t を問題に入れ (spec.r_throat を書き換えた新しい問題)、生産経路 (prepare_ns と同じ: integral_delta_r か
    delta_r_csv の表 → build_physical_wall) で壁を作り直して、連続な壁の最小半径 (または出口半径) と目標の差を測る。"""
    from forge_design.evaluate.runner_axismach import build_physical_wall, delta_r_from_table, design_chain, integral_delta_r, load_problem
    p = load_problem(problem)
    p.spec["r_throat"] = float(rt)
    d = design_chain(p)
    if delta_csv is None:
        _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    else:
        tbl = np.loadtxt(delta_csv, delimiter=",", skiprows=1)
        drx = delta_r_from_table(tbl[:, 0], tbl[:, 1])
    W = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    got = rt * (W.r_throat if target == "throat" else float(W.r(np.r_[W.x_e])[0]))
    return {"r_t_m": rt, "target_m": R, "achieved_m": got, "diff_m": got - R, "pass": bool(abs(got - R) <= U2_TOL_M),
            "x_throat_rt": W.x_throat, "exit_radius_m": rt * float(W.r(np.r_[W.x_e])[0])}


def noise_floor(problem: Path, rt: float, target: str, supplier=None) -> dict:
    """残差の雑音の床: r_t を ±1e-12・1e-11・1e-10 m と ulp の数倍ずらして R(r_t) を評価し、1 次の傾向を除いた散らばり。"""
    from forge_design.evaluate.runner_axismach import build_physical_wall, design_chain, integral_delta_r, load_problem
    p = load_problem(problem)
    d = design_chain(p)
    if supplier is None:
        def supplier(d_, r_):
            _, f, _ = integral_delta_r(p, d_, p.raw["deltastar_initializer"], scale=r_)
            return f, None
    drt = np.r_[-1e-10, -1e-11, -3e-12, -1e-12, -3 * np.spacing(rt), 0.0, 3 * np.spacing(rt), 1e-12, 3e-12, 1e-11, 1e-10]
    vals = []
    for dd in drt:
        f, _ = supplier(d, rt + dd)
        W = build_physical_wall(p, d, rt + dd, delta_r_x=f, offset="radial")
        vals.append((rt + dd) * (W.r_throat if target == "throat" else float(W.r(np.r_[W.x_e])[0])))
    vals = np.asarray(vals)
    slope = np.polyfit(drt, vals, 1)
    res = vals - np.polyval(slope, drt)
    return {"drt_m": drt.tolist(), "R_m": vals.tolist(), "slope": float(slope[0]), "detrended_spread_m": float(res.max() - res.min()),
            "detrended_rms_m": float(np.sqrt(np.mean(res ** 2)))}


def run_u2solve():
    from forge_design.feedback.deltastar_loop import SizingNotConverged, solve_rt_throat
    from forge_design.evaluate.runner_axismach import load_problem
    prob = problem_of("u3_A_poly_legacy")          # poly (明示) + legacy。寸法の逆算は壁の表現に依らない (区分表現で解く)
    p = load_problem(prob)
    _, _, _, _, W0 = build(prob)
    R_b = R_T0 * W0.r_throat
    targets = {"base": R_b, "plus5": 1.05 * R_b, "minus5": 0.95 * R_b}
    out = {"problem": str(prob), "r_t0_m": R_T0, "base_throat_radius_m": R_b, "targets_m": targets, "tol_m": U2_TOL_M, "arms": {}}
    for arm in ("B", "A"):
        rows = {}
        for tag, R in targets.items():
            t0 = time.time()
            try:
                kw = {} if arm == "B" else {"_delta_supplier": _supplier_A(p)}
                r = solve_rt_throat(prob, R, **kw)
                row = {"converged": True, "r_t_m": r["r_t_m"], "residual_m": r["residual_m"], "n_iter": r["n_iter"],
                       "iter_history": r["iter_history"], "physical_throat": r["physical_throat"], "exit_radius_m": r["exit_radius_m"],
                       "delta_r_source": r["delta_r_source"]}
                row["roundtrip_production"] = roundtrip_production(prob, r["r_t_m"], R)
                rt_last = r["r_t_m"]
            except SizingNotConverged as e:
                row = {"converged": False, "error": str(e), "iter_history": e.history}
                rt_last = e.history[-1]["r_t_m"]
                row["roundtrip_production_last_iterate"] = roundtrip_production(prob, rt_last, R)
            except ValueError as e:
                row = {"converged": False, "error": f"ValueError: {e}"[:400]}
                rt_last = None
            row["seconds"] = time.time() - t0
            if rt_last is not None and arm == "B":
                row["noise_floor"] = noise_floor(prob, rt_last, "throat")
            rows[tag] = row
            print(arm, tag, row.get("converged"), row.get("r_t_m"), row.get("residual_m"),
                  (row.get("roundtrip_production") or row.get("roundtrip_production_last_iterate") or {}).get("diff_m"), flush=True)
        out["arms"][arm] = rows
    # NS 後: run_0092 の帯 E の抽出 (_extract_edge/delta_r_next.csv) の δ_E の全分布
    prev = RUNS / "run_0092_ns_c2pin_pass2"
    dn = prev / "_extract_edge" / "delta_r_next.csv"
    ns = {}
    for tag, R in targets.items():
        csv = UP / f"u2_ns_{tag}_delta_r.csv"
        try:
            r = solve_rt_throat(prob, R, prev_run=prev, delta_next=dn, delta_r_out=csv)
            row = {"converged": True, "r_t_m": r["r_t_m"], "residual_m": r["residual_m"], "n_iter": r["n_iter"],
                   "iter_history": r["iter_history"], "delta_r_source": r["delta_r_source"], "exit_radius_m": r["exit_radius_m"],
                   "physical_throat": r["physical_throat"], "delta_r_out": str(csv)}
            row["roundtrip_production"] = roundtrip_production(prob, r["r_t_m"], R, delta_csv=csv)
        except (SizingNotConverged, ValueError) as e:
            row = {"converged": False, "error": str(e)[:400]}
        ns[tag] = row
        print("NS", tag, row.get("converged"), row.get("r_t_m"), (row.get("roundtrip_production") or {}).get("diff_m"), flush=True)
    out["ns_after"] = {"prev_run": str(prev), "delta_next": str(dn), "rows": ns}
    try:
        solve_rt_throat(prob, 0.9 * R_b, max_iter=1)
        out["iteration_cap"] = {"raised": False}
    except SizingNotConverged as e:
        out["iteration_cap"] = {"raised": True, "error": str(e)}
    dump("u2_solve.json", out)
    print("prepare_ns の往復に使う r_t:", json.dumps({t: out["arms"]["B"][t].get("r_t_m") for t in targets}),
          json.dumps({t: ns[t].get("r_t_m") for t in targets}))
    return out


def run_u2():
    s = json.loads((UP / "u2_solve.json").read_text())
    rows = {}
    for arm in ("B", "A"):
        for tag, R in s["targets_m"].items():
            r = s["arms"][arm][tag]
            rt = r.get("roundtrip_production") or {}
            rows[f"{arm}/{tag}"] = {"converged": r["converged"], "r_t_m": r.get("r_t_m"), "solver_residual_m": r.get("residual_m"),
                                    "n_iter": r.get("n_iter"), "roundtrip_diff_m": rt.get("diff_m"),
                                    "pass": bool(r["converged"] and rt.get("pass")),
                                    "last_iterate_roundtrip_diff_m": (r.get("roundtrip_production_last_iterate") or {}).get("diff_m"),
                                    "noise_floor_spread_m": (r.get("noise_floor") or {}).get("detrended_spread_m"),
                                    "error": r.get("error")}
    # prepare_ns の端から端まで (B の 3 目標 + NS 後の基準): spec.sizing の記録の残差
    e2e = {}
    for tag in s["targets_m"]:
        for kind in ("B", "NS"):
            run = UP / f"u2_prep_{kind}_{tag}"
            if (run / "prepare_info.json").exists():
                sz = json.loads((run / "prepare_info.json").read_text())["sizing"]
                e2e[f"{kind}/{tag}"] = {"run": str(run), "sizing": sz, "pass": bool(sz["residual_m"] is not None and abs(sz["residual_m"]) <= U2_TOL_M)}
    ns = {tag: {"converged": r["converged"], "r_t_m": r.get("r_t_m"), "solver_residual_m": r.get("residual_m"), "n_iter": r.get("n_iter"),
                "roundtrip_diff_m": (r.get("roundtrip_production") or {}).get("diff_m"),
                "pass": bool(r["converged"] and (r.get("roundtrip_production") or {}).get("pass")), "error": r.get("error")}
          for tag, r in s["ns_after"]["rows"].items()}
    B_pass = all(rows[f"B/{t}"]["pass"] for t in s["targets_m"])
    A_pass = all(rows[f"A/{t}"]["pass"] for t in s["targets_m"])
    if B_pass and not A_pass:
        interp = "A だけ外れ B が通る → 共通経路の採用を支持 (事前登録)"
    elif B_pass and A_pass:
        interp = "両方通る → 「経路の違いで外れる」仮説を棄却 (事前登録)"
    else:
        interp = "B も外れる → 事前登録どおり、反復の停止条件・寸法の反映・壁の生成の不一致を調べて止める"
    out = {"item": "U2 スロート径から決める: 往復検証 (plan §6.0、判別 A/B、CFD 0 step)",
           "criterion": f"B (prepare_ns と共通の δ_r の経路) が全目標 (基準と ±5 %) で、逆算した r_t から生産経路で作り直した壁の最小半径と目標の差 ≤ {U2_TOL_M} m。"
                        "A (integral_bl 直接) も記録。NS 後の経路も同じ往復で ≤ 1e-9 m。反復の上限で不合格になること",
           "arms": {"A": "solve_rt 型の integral_bl 直接 (未較正・未平滑化)", "B": "prepare_ns と共通 (integral_delta_r)"},
           "rows": rows, "B_pass": B_pass, "A_pass": A_pass, "interpretation_preregistered": interp,
           "prepare_ns_end_to_end": e2e, "ns_after": ns, "ns_after_pass": all(v["pass"] for v in ns.values()),
           "iteration_cap": s.get("iteration_cap"), "targets_m": s["targets_m"]}
    out["verdict"] = ("B: " + ("PASS" if B_pass else "FAIL") + " / A: " + ("PASS" if A_pass else "FAIL")
                      + " / NS 後: " + ("PASS" if out["ns_after_pass"] else "FAIL")
                      + " / 反復の上限: " + ("PASS" if (s.get("iteration_cap") or {}).get("raised") else "FAIL"))
    return dump("U2.json", out)


def run_u2b():
    from forge_design.feedback.deltastar_loop import SizingNotConverged, solve_rt
    prob = problem_of("u3_A_poly_legacy")
    rows = {}
    for tag, R in (("0.775", 0.775), ("plus5", 0.775 * 1.05), ("minus5", 0.775 * 0.95)):
        t0 = time.time()
        try:
            r = solve_rt(prob, R)
            row = {"converged": True, "r_t_m": r["r_t_m"], "residual_m": r["residual_m"], "n_iter": r["n_iter"], "iter_history": r["iter_history"],
                   "delta_r_source": r["delta_r_source"], "physical_throat": r["physical_throat"]}
            row["roundtrip_production"] = roundtrip_production(prob, r["r_t_m"], R, target="exit")
            rt_last = r["r_t_m"]
        except SizingNotConverged as e:
            row = {"converged": False, "error": str(e), "iter_history": e.history}
            rt_last = e.history[-1]["r_t_m"]
            row["roundtrip_production_last_iterate"] = roundtrip_production(prob, rt_last, R, target="exit")
        row["noise_floor"] = noise_floor(prob, rt_last, "exit")
        row["seconds"] = time.time() - t0
        row["pass"] = bool(row["converged"] and row["roundtrip_production"]["pass"])
        rows[tag] = row
        print(tag, row["converged"], rt_last, (row.get("roundtrip_production") or row.get("roundtrip_production_last_iterate"))["diff_m"],
              row["noise_floor"]["detrended_spread_m"], flush=True)
    out = {"item": "U2b 出口径から決める CFD 前の経路の往復検証 (plan §6.0)",
           "criterion": "変えた solve_rt で出口半径の目標 (0.775 m と ±5 %) から r_t を解き、生産経路で壁を作り直して、出口半径と目標の差 ≤ 1e-9 m",
           "rows": rows, "verdict": "PASS" if all(r["pass"] for r in rows.values()) else "FAIL"}
    return dump("U2b.json", out)


# ======================================================================================================== U3 (+ W1・W2)
def judge_u3(B, PW) -> dict:
    """U3: 1 本の B-spline と poly の区分表現の差 (各ノット区間の内部の密な点 + 区間多項式の極値、継ぎ目の左右の極限)。
    判定側で測り直す (B.fit_diag を信用しない)。"""
    from forge_design.geometry.wall_axismach import bspline_piece_limits, bspline_wall_errors
    spl = B.spline
    dist = np.unique(np.asarray(spl.t))
    iv = np.c_[dist[:-1], dist[1:]]
    err = bspline_wall_errors(spl, PW.r, iv, n_dense=61)
    emax = {n: float(err[n][0].max()) for n in range(4)}
    xat = {n: float(err[n][1][int(np.argmax(err[n][0]))]) for n in range(4)}
    L = PW.L_U
    jl, jr = bspline_piece_limits(spl, [-L, 0.0], nmax=2)
    src = {-L: ([PW.r_U, 0.0, 0.0], [float(PW._q_eval(PW._q_c, -L, 0.0, np.r_[-L], n)[0]) for n in range(3)]),
           0.0: ([float(PW._q_eval(PW._q_c, -L, 0.0, np.r_[0.0], n)[0]) for n in range(3)],
                 [float(PW.design.r(np.r_[0.0], n)[0] + PW._dr(np.r_[0.0], n)[0]) for n in range(3)])}
    lim = []
    for q, xj in enumerate((-L, 0.0)):
        sl, sr = src[xj]
        lim.append({"x": xj, **{f"left_d{n}": float(abs(jl[n, q] - sl[n])) for n in range(3)},
                    **{f"right_d{n}": float(abs(jr[n, q] - sr[n])) for n in range(3)}})
    lmax = {n: max(max(r_[f"left_d{n}"], r_[f"right_d{n}"]) for r_ in lim) for n in range(3)}
    tot = {n: max(emax[n], lmax[n]) for n in range(3)}
    ok = tot[0] <= U3_TOL["r"] and tot[1] <= U3_TOL["r1"] and tot[2] <= U3_TOL["r2"]
    return {"pass": bool(ok), "tol": U3_TOL, "max_err_interior": {f"d{n}": emax[n] for n in range(4)},
            "x_at_max_err_interior": {f"d{n}": xat[n] for n in range(4)}, "joint_limits": lim,
            "max_err_joint_limits": {f"d{n}": lmax[n] for n in range(3)}, "max_err_total": {f"d{n}": tot[n] for n in range(3)},
            "n_intervals": int(len(iv))}


def judge_w1_poly(B, PW, mono) -> dict:
    """全域 1 本の plan の W1 (ノット挿入の版): 許容誤差 (半径 1.3e-7・r′ 1e-7・r″ 1e-5)・継ぎ目の跳び ≤ 1e-8・区間内 C⁴ (重複度)・
    形のゲート S6 (r″ の最大増加 ≤ 0.002・validate 空) と、ランプのゲートの代わりに上流のゲート (判定側で測り直す: |r″ − H″|・
    単調性を密な点で)。"""
    from forge_design.geometry.wall_axismach import bspline_piece_limits, bspline_wall_errors, r3_piecewise_exact
    spl = B.spline
    t = np.asarray(spl.t)
    dist, mult = np.unique(t, return_counts=True)
    iv = np.c_[dist[:-1], dist[1:]]
    err = bspline_wall_errors(spl, PW.r, iv)
    emax = {n: float(err[n][0].max()) for n in range(4)}
    tol_ok = emax[0] <= W1_TOL["r"] and emax[1] <= W1_TOL["r1"] and emax[2] <= W1_TOL["r2"]
    jl, jr = bspline_piece_limits(spl, B.joints, nmax=3)
    jumps = [{"x": float(x), **{f"d{n}": float(abs(jr[n, q] - jl[n, q])) for n in range(4)}} for q, x in enumerate(B.joints)]
    jump_ok = all(j["d0"] <= W1_TOL["joint_jump_d1_d2"] and j["d1"] <= W1_TOL["joint_jump_d1_d2"] and j["d2"] <= W1_TOL["joint_jump_d1_d2"] for j in jumps)
    inner = (dist > B.x_in) & (dist < B.x_e)
    is_joint = np.isin(dist, np.asarray(B.joints))
    simple = inner & ~is_joint
    struct_ok = bool(np.all(mult[simple] == 1) and np.all(mult[is_joint] == 3) and mult[0] == 6 and mult[-1] == 6)
    sl, sr = bspline_piece_limits(spl, dist[simple], nmax=4)
    rel = np.abs(sr - sl) / np.maximum(1.0, np.maximum(np.abs(sl), np.abs(sr)))
    m6 = r3_piecewise_exact(spl, None, *mono)
    v = B.validate()
    s6_ok = m6["r2_max_increase"] <= W1_TOL["S6_r2_max_increase"] and v == []
    L = B.L_U
    xs = np.r_[np.linspace(-L, 0.0, 200001)[:-1], -1e-13]
    d2 = float(np.max(np.abs(spl(xs, 2) - B.design.r(xs, 2))))
    xd = np.linspace(B.x_in, B.x_e, 2000001)
    r1 = spl(xd, 1)
    up_ok = bool(d2 <= 5e-3 and r1[xd < B.x_throat].max() <= 1e-12 and r1[xd > B.x_throat].min() >= -1e-12)
    ok = tol_ok and jump_ok and struct_ok and s6_ok and up_ok
    return {"pass": bool(ok), "tol": dict(W1_TOL), "max_err": {f"d{n}": emax[n] for n in range(4)}, "tol_pass": bool(tol_ok),
            "joint_jumps": jumps, "joint_jump_pass": bool(jump_ok), "c4_structure_pass": struct_ok,
            "n_simple_interior_knots": int(simple.sum()), "joint_mults": mult[is_joint].tolist(),
            "c4_numeric_record": {f"max_rel_jump_d{n}": float(np.nanmax(rel[n])) for n in range(5)},
            "S6": {"pass": bool(s6_ok), "interval": list(mono), "r2_max_increase": m6["r2_max_increase"], "r2_max": m6["r2_max"],
                   "limit": W1_TOL["S6_r2_max_increase"], "validate": v},
            "upstream_gate_rejudged": {"pass": up_ok, "abs_d2_change_vs_H_dense": d2, "max_r1_before_dense": float(r1[xd < B.x_throat].max()),
                                       "min_r1_after_dense": float(r1[xd > B.x_throat].min()), "wall_attr": B.upstream_gate["pass"]},
            "counts": {"n_coef": int(len(spl.c)), "n_knots": int(len(t)), "n_distinct_knots": int(len(dist)),
                       "min_knot_gap_rt": float(np.diff(dist).min()), "joints": list(B.joints)},
            "construction": {k: B.fit_diag[k] for k in ("construction", "insertion", "knot_removal", "build_seconds")}}


def run_u3():
    from forge_design.geometry.wall_axismach import load_wall_file
    p, d, res, drx, PW, B, sec = piecewise_poly(problem_of("u3_B_poly_single_bspline"))
    j = judge_u3(B, PW)
    Wf = load_wall_file(UP / "u3_B_poly_single_bspline")
    same = bool(np.array_equal(Wf["physical"].spline.t, B.spline.t) and np.array_equal(Wf["physical"].spline.c, B.spline.c))
    cmp = WV.compare_dirs(UP / "u3_A_poly_legacy", UP / "u3_B_poly_single_bspline")
    j3 = WV.judge_w3(cmp)
    ha, hb = (json.loads((UP / f"{n}.hashes.json").read_text()) for n in ("u3_A_poly_legacy", "u3_B_poly_single_bspline"))
    env_same = {k: ha.get(k) == hb.get(k) for k in ("design_dir", "python", "numpy", "scipy", "converter_sha256", "FORGE_ALLOW_UNVERIFIED_SPECIES")}
    # 参考: キー無し (既定 poly、壁ファイルなし) と腕 A のソルバ入力
    extra = WV.compare_dirs(UP / "u1_poly_default", UP / "u3_A_poly_legacy")
    out = {"item": "U3 全域 1 本の B-spline (ノット挿入、plan §6.0) と W3 のやり直し (全域 1 本の plan §6.0)",
           "criterion_u3": f"poly の物理壁と 1 本の B-spline の差 (各ノット区間の内部の密な点と継ぎ目の左右の極限): 半径 ≤ {U3_TOL['r']}、"
                           f"1 階 ≤ {U3_TOL['r1']}、2 階 ≤ {U3_TOL['r2']}。所要時間を記録",
           "u3": {**j, "build_seconds_in_verify": sec, "build_seconds_in_class": B.fit_diag["build_seconds"],
                  "same_as_run_wall_file": same, "run_wall_file": Wf["path"]},
           "criterion_w3": "同じ poly の壁を両方の表現 (A = legacy、B = single_bspline) で prepare_ns まで作り、ソルバ入力 (nozzle.h5 の全データセットと"
                           "属性・solverConfig.yaml・bcondConfig.yaml) が全項目ビット同一なら支持、1 つでも違えば棄却して差の経路を特定して止める",
           "w3": {**j3, "verdict": "PASS (全項目ビット同一)" if j3["pass"] else "FAIL (入力が違う — 差の経路を特定して止める)",
                  "same_environment": env_same, "problem_copy_diff": WV._problem_diff(ha["problem_copy"], hb["problem_copy"]),
                  "nozzle_msh_text": cmp["nozzle_msh"], "files": cmp["files"], "prepare_info_diff": cmp["prepare_info_diff"][:40],
                  "h5_datasets_differing": {k: v for k, v in cmp["h5"]["datasets"].items() if not v["identical"]}},
           "extra_default_vs_A": {"note": "参考: キー無し (既定 poly、壁ファイルなし) と腕 A (poly 明示 + legacy) のソルバ入力",
                                  "h5_all_identical": extra["h5"]["all_identical"],
                                  "solver_config_identical": all(v["identical_bytes"] for v in extra["solver_config"].values()),
                                  "files": extra["files"]}}
    out["verdict"] = ("U3: " + ("PASS" if j["pass"] else "FAIL") + " / W3: " + ("PASS" if j3["pass"] else "FAIL"))
    return dump("U3.json", out)


def run_w1w2():
    p, d, res, drx, PW, B, _ = piecewise_poly(problem_of("u3_B_poly_single_bspline"))
    mono = tuple(float(v) for v in p.geometry["wall_fit_mono_r2"])
    j1 = judge_w1_poly(B, PW, mono)
    out1 = {"item": "W1 作り直しの精度 (全域 1 本の plan §6.0、ノット挿入の版でやり直し)", "verdict": "PASS" if j1["pass"] else "FAIL", **j1}
    dump("W1.json", out1)
    j2 = WV.judge_w2(B, PW)
    out2 = {"item": "W2 スロート量 (全域 1 本の plan §6.0、ノット挿入の版)", "verdict": "PASS" if j2["pass"] else "FAIL", **j2,
            "note": "judge_w2 は旧来の囲い込み (−0.3, 0.2) の Brent 法で判定側が求め直す (壁は大域最小の探索器)"}
    dump("W2.json", out2)
    return out1, out2


def run_w5():
    from forge_design.export.wall_step import export_run, read_step, sample_params, step_curve_data
    from forge_design.geometry.wall_axismach import load_wall_file
    RB = UP / "u3_B_poly_single_bspline"
    outd = UP / "step"
    if outd.exists():
        shutil.rmtree(outd)
    sc = export_run(RB, outd, n_per_interval=3)
    W = load_wall_file(RB)
    data = step_curve_data(W["record"])
    pts = sample_params(data, 3)
    rd = read_step(outd / "wall_physical.step", [u for u, _, _ in pts], revolve=False)
    p, d, res, drx, PW, B, _ = piecewise_poly(problem_of("u3_B_poly_single_bspline"))
    same = bool(np.array_equal(B.spline.c, W["physical"].spline.c) and np.array_equal(B.spline.t, W["physical"].spline.t))
    orig = WV.judge_w5_orig(data, rd, pts, B, PW, data["scale_mm_per_rt"])
    tr = sc["readback"]["transfer"]
    rv = sc["readback"]["revolve"] or {}
    out = {"item": "W5 STEP (全域 1 本の plan §6.0、ノット挿入の版でやり直し)", "step": str(outd / "wall_physical.step"),
           "sidecar": str(outd / "wall_physical_step.json"), "same_as_rebuilt_wall": same,
           "structure": tr["structure"], "structure_pass": all(tr["structure"].values()),
           "transfer": {k: tr[k] for k in tr if k != "structure"}, "orig": orig,
           "revolve": {**rv, "pass": bool(rv.get("is_valid") is True and rv.get("area_mm2", 0) > 0)},
           "freecad_version": sc["readback"].get("freecad_version"), "chord_minus_curve_um": sc["cad_vs_cfd"].get("chord_minus_curve"),
           "receiving_cad": "未確認 (受け取り側の CAD の種類・版での読み込みはユーザ側で確認)"}
    ok = out["structure_pass"] and tr["pass"] and orig["pass"] and orig["derivative_units_pass"] and out["revolve"]["pass"] and same
    out["verdict"] = "PASS" if ok else "FAIL"
    return dump("W5.json", out)


def run_roundtrip():
    """保存 → 復元 → 報告の往復 (poly の区分表現 legacy と single_bspline の両方): 壁ファイルから復元した物理壁が作り直した壁と
    ビット一致、prepare_info.json の physical_wall と壁ファイルの一致、報告 (fig_wall_shape) が保存した係数から評価する。"""
    from forge_design.geometry.wall_axismach import WALL_FILE, load_wall_file
    from forge_design.report.nozzle_report import fig_wall_shape
    res = {}
    for run, kind in (("u3_A_poly_legacy", "legacy"), ("u3_B_poly_single_bspline", "single_bspline")):
        R = UP / run
        p, d, _, drx, PW, B, _ = piecewise_poly(problem_of(run))
        live = PW if kind == "legacy" else B
        W = load_wall_file(R)
        info = json.loads((R / "prepare_info.json").read_text())
        S = float(info["scale_m"])
        xq = np.unique(np.r_[np.linspace(W["domain"][0], W["domain"][1], 200001), np.asarray(B.spline.t)])
        nder = 4 if kind == "legacy" else 3
        bit = {f"d{n}": bool(np.array_equal(W["physical"].r(xq, n), live.r(xq, n))) for n in range(nder)}
        with tempfile.TemporaryDirectory() as td:
            o = fig_wall_shape(R, {"S": S}, Path(td) / "fig.png")
            shutil.copy(Path(td) / "fig.png", UP / f"roundtrip_fig_wall_shape_{kind}.png")
        rec_ok = (W["record"]["version"] == 2 and W["pw_upstream"] == "poly" and W["kind"] == kind
                  and info["physical_wall"]["physical_wall"] == W["record"]["physical_wall"] and info["physical_wall"]["sha256"] == WV.sha_file(R / WALL_FILE)
                  and W["record"]["physical_wall"].get("upstream_poly", {}).get("coef") == [float(v) for v in PW._q_c])
        rep_ok = (o.get("wall_source") == "saved_coefficients" and o.get("wall_repr") == kind
                  and o["exit_radius_m"] == float(W["physical"].r(np.r_[W["domain"][1]])[0] * S))
        res[kind] = {"run": str(R), "restored_bit_identical": bit, "record_ok": bool(rec_ok), "report_ok": bool(rep_ok),
                     "report": {k: o.get(k) for k in ("wall_source", "wall_repr", "exit_radius_m", "r2_highfreq_max_x_gt2")},
                     "throat_record": W["throat"], "pass": bool(all(bit.values()) and rec_ok and rep_ok)}
    out = {"item": "保存 → 復元 → 報告の往復 (plan §4.1 M1、legacy (poly) と single_bspline の両方)", "rows": res,
           "verdict": "PASS" if all(r["pass"] for r in res.values()) else "FAIL"}
    return dump("roundtrip.json", out)


# ======================================================================================================== 負例
def run_neg():
    import copy
    import h5py
    from scipy.interpolate import BSpline
    from forge_design.evaluate.runner_axismach import delta_r_from_table
    from forge_design.export.wall_step import read_step, sample_params, step_curve_data, transfer_check, write_step
    from forge_design.geometry.wall_axismach import WALL_FILE, PhysicalNozzleWall, load_wall_file
    res = {}
    with tempfile.TemporaryDirectory(dir=str(Path(os.environ.get("TMPDIR", "/tmp")))) as td:
        td = Path(td)
        base, after = UP / "u0_base_head", UP / "u0_after_ramp"
        ia = json.loads((base / "prepare_info.json").read_text())
        # U0: 既存の値を 1 ulp 変える / キーを 1 つ増やす / nozzle.h5 の 1 要素を 1 ulp
        a1 = td / "after_val"; shutil.copytree(after, a1)
        ib = json.loads((a1 / "prepare_info.json").read_text()); ib["throat_physical"]["x"] = float(np.nextafter(ib["throat_physical"]["x"], 1.0))
        (a1 / "prepare_info.json").write_text(json.dumps(ib, indent=1))
        j = judge_u0(WV.compare_dirs(base, a1), ia, ib)
        res["U0_existing_value_one_ulp"] = {"detected": not j["pass"], "why": j["prepare_info_old_keys_differing"]}
        ib2 = json.loads((after / "prepare_info.json").read_text()); ib2["extra"] = 1
        j = judge_u0(WV.compare_dirs(base, after), ia, ib2)
        res["U0_extra_key"] = {"detected": not j["pass"], "why": j["prepare_info_new_keys"]}
        a3 = td / "after_h5"; shutil.copytree(after, a3)
        with h5py.File(a3 / "nozzle.h5", "r+") as f:
            v = f["/VALUE/ro"][:]; v[777] = np.nextafter(v[777], np.float32(np.inf)); f["/VALUE/ro"][:] = v
        j = judge_u0(WV.compare_dirs(base, a3), ia, json.loads((after / "prepare_info.json").read_text()))
        res["U0_nozzle_h5_one_ulp"] = {"detected": not j["pass"], "why": j["nozzle_h5_datasets_and_attrs_identical"]}
        ib4 = json.loads((after / "prepare_info.json").read_text()); ib4["pw_upstream"]["source"] = "default"
        j = judge_u0(WV.compare_dirs(base, after), ia, ib4)
        res["U0_pw_upstream_not_explicit"] = {"detected": not j["pass"]}
        # U1: [0, x_e] のビット同一の判定が δ_r の 1 点の 1 ulp を拾う / ゲートの判定が r′ < 0 を拾う
        p, d, rr, drx, PP = build(problem_of("u1_poly_default"))
        y = drx(rr["x"]).copy(); y[800] = np.nextafter(y[800], 1.0)
        f2 = delta_r_from_table(rr["x"], y)
        from forge_design.evaluate.runner_axismach import _gam_or_gas
        args = (d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
        PR2 = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f2, ramp=(-11.0, -6.0), upstream="ramp")
        j = judge_u1_identity(PP, PR2)
        res["U1_identity_delta_one_ulp"] = {"detected": not j["pass"], "why": j["identical"]}
        Pm = copy.copy(PP); Pm.x_throat = PP.x_throat + 1e-5
        j = judge_u1_gate(Pm)
        res["U1_gate_throat_attr_shifted"] = {"detected": not j["pass"], "why": {"brute": j["throat_brute_force_x"], "attr": j["throat_x"]}}
        # U2: 往復の判定が r_t の 1e-8 m のずれを拾う
        s = json.loads((UP / "u2_solve.json").read_text())
        rB = s["arms"]["B"]["base"]
        if rB.get("converged"):
            j = roundtrip_production(problem_of("u3_A_poly_legacy"), rB["r_t_m"] + 1e-8, s["targets_m"]["base"])
            res["U2_roundtrip_rt_plus_1e-8m"] = {"detected": not j["pass"], "why": j["diff_m"]}
        # U3: 係数を 1e-11 ずらす → 半径の誤差 > 1e-12 / W1: 係数を 1e-6 / 継ぎ目のノットを単純ノットに戻せない (重複度) の検出
        p3, d3, _, _, PW, B, _ = piecewise_poly(problem_of("u3_B_poly_single_bspline"))
        c = np.array(B.spline.c); c[700] += 1e-11
        Bx = copy.copy(B); Bx.spline = BSpline(B.spline.t, c, 5)
        j = judge_u3(Bx, PW)
        res["U3_coef_plus_1e-11"] = {"detected": not j["pass"], "why": j["max_err_total"]}
        c = np.array(B.spline.c); c[700] += 1e-6
        Bx = copy.copy(B); Bx.spline = BSpline(B.spline.t, c, 5)
        mono = tuple(float(v) for v in p3.geometry["wall_fit_mono_r2"])
        j = judge_w1_poly(Bx, PW, mono)
        res["W1_coef_plus_1e-6"] = {"detected": (not j["pass"]) and (not j["tol_pass"]), "why": j["max_err"]}
        # W3: 腕 B の写しの VALUE/wall_dist を 1 ulp
        b1 = td / "B_wd"; shutil.copytree(UP / "u3_B_poly_single_bspline", b1)
        with h5py.File(b1 / "nozzle.h5", "r+") as f:
            v = f["/VALUE/wall_dist"][:]; v[-1] = np.nextafter(v[-1], np.float32(np.inf)); f["/VALUE/wall_dist"][:] = v
        j = WV.judge_w3(WV.compare_dirs(UP / "u3_A_poly_legacy", b1))
        res["W3_wall_dist_one_ulp"] = {"detected": not j["pass"], "why": [k for k, v in j["items"].items() if not v["identical"]]}
        # 往復: 壁ファイルの Q の係数を 1e-6 → 読み込みが例外 / 1 本の係数を 1e-12 → 復元が作り直しとビット一致しない
        rA = json.loads((UP / "u3_A_poly_legacy" / WALL_FILE).read_text())
        rA["physical_wall"]["upstream_poly"]["coef"][2] += 1e-6
        q = td / "rtA"; q.mkdir(); (q / WALL_FILE).write_text(json.dumps(rA))
        try:
            load_wall_file(q); res["roundtrip_Q_coef_broken_raises"] = False
        except ValueError:
            res["roundtrip_Q_coef_broken_raises"] = True
        rBf = json.loads((UP / "u3_B_poly_single_bspline" / WALL_FILE).read_text())
        rBf["physical_wall"]["c"][900] += 1e-12
        q = td / "rtB"; q.mkdir(); (q / WALL_FILE).write_text(json.dumps(rBf))
        Wq = load_wall_file(q)
        xq = np.linspace(Wq["domain"][0], Wq["domain"][1], 200001)
        res["roundtrip_sb_coef_plus_1e-12_not_identical"] = {"detected": not np.array_equal(Wq["physical"].r(xq), B.r(xq))}
        # W5: 制御点 1 個を 1e-5 mm ずらした STEP
        W = load_wall_file(UP / "u3_B_poly_single_bspline")
        data = step_curve_data(W["record"])
        pts = sample_params(data, 1)
        dd = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v) for k, v in data.items()}
        dd["poles_mm"][800, 1] += 1e-5
        stp = td / "neg_pole.step"
        write_step(dd, stp)
        j = transfer_check(data, read_step(stp, [u for u, _, _ in pts]), pts)
        res["W5_pole_shift_1e-5mm"] = {"detected": not j["pass"], "why": j.get("pos_max_mm")}
    det = {k: (v["detected"] if isinstance(v, dict) else bool(v)) for k, v in res.items()}
    out = {"item": "判定スクリプトの負例 (FAIL を FAIL と出すか; detected = 意図した検査で FAIL・例外になった)", "cases": res, "detected": det,
           "verdict": "PASS (全負例を FAIL・例外として検出)" if all(det.values()) else "FAIL (検出できない負例がある)"}
    return dump("negatives.json", out)


if __name__ == "__main__":
    cmd = sys.argv[1]
    r = globals()[f"run_{cmd}"]()
    for x in (r if isinstance(r, tuple) else (r,)):
        print(cmd, x.get("verdict"))
