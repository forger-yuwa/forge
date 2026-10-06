"""plan discretization-moc-axis-limit-and-corrector §6 V4 (case/45 の形状の試算) と V6 (所要時間と収束)。CFD 0 step。

対象: 単調壁の生産問題 problem_d155_ns_finemesh_recal_final_mono.yaml (n_start 41・n_axis_inv 2400・axis_dx0 0.03)。
比較の基準 = 現行の単調壁 (キー無し = legacy + fixed2)。新 = geometry.moc_axis_limit: analytic + moc_corrector: converge。
参考として legacy + converge (修正子だけ収束) も同じ量を記録する (軸端の源項と修正子の寄与を分けるため。判定には使わない)。

記録 (§6 V4): 第 1 点の角度差 Δθ₁ = θ₁ − atan(x₁/R)、単調拘束なしの r″ の山 (joint 当てはめ、[0, 0.3] の max r″ − 1/R)、
AXIS_LIMIT_FRAC の発火の数と位置、未収束の対の数 / 設計壁の変化 (最大・位置・出口)、試験部 (x 40〜94 r_t) の壁角の変化、出口半径の変化。
候補にする条件 (結果を見る前に plan §6 V4 に登録されたもの。ここで変えない):
  未収束の対 0、かつ単調壁の形状ゲート S1・S2・S6・S7・S8 が PASS (`throat_mono_shape_gate.py` と同じ式・同じ閾値、
  S8 の比較相手は現行の単調壁)、かつ置き換えた S3 (x ≥ 0.3 での設計壁の変化 ≤ 10 µm、出口での変化 ≤ 0.1 µm) と
  S4 (新しい点群から求めた解析的下限 θ₁ − atan(x₁/R) と第 1 点の当てはめのずれの差 ≤ 0.002°)。外れたら諮問。
  S7 は本変更の「既定のビット同一」= キー無しの設計 (MOC 点群・当てはめ後の係数とノット) が実装前のコミットと完全一致。
V6: design_chain 1 回の所要時間 (各 3 回の中央値、現行・新) と新の修正子の反復回数の分布。目安は現行の 5 倍以下。

usage: [CASE_RUNS=<run_0062 のある case dir>] python3 moc_axis_limit_v4_case45.py [OUT_JSON]
       既定の CASE_RUNS はローカル主ツリー /home/sano/work/forge/case/45.isobutane_m6_d155。出力 _band_ab/moc_axis_limit_v4_case45.json
"""
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate.runner_axismach import (_gam_or_gas, design_chain, integral_delta_r,  # noqa: E402
                                                   load_problem)
from forge_design.geometry.wall_axismach import PhysicalNozzleWall, joint_fit_wall, r3_piecewise_exact  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/moc_axis_limit_v4_case45.json"
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
PROB = C / "problem_d155_ns_finemesh_recal_final_mono.yaml"
BASE = "170b0d75924fa695134be0840e32fb37e5754917"      # 実装前 (S7 = 既定のビット同一の基準)
MONO = [0.0, 1.5]
TS = (40.0, 94.0)                                       # 試験部の壁角の比較区間 [r_t] (throat_moc_fix_probe.py と同じ)
ARMS = {"current": None, "new": ("analytic", "converge"), "ref_legacy_converge": ("legacy", "converge")}
N_TIME = 3


def load(keys):
    p = load_problem(PROB)
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    if keys is not None:
        p.geometry["moc_axis_limit"], p.geometry["moc_corrector"] = keys
    return p


def run(keys):
    p = load(keys)
    t0 = time.time()
    d = design_chain(p)
    return p, d, time.time() - t0


p0 = load(None)
assert (int(p0.geometry["n_start"]), int(p0.geometry["n_axis_inv"]), float(p0.geometry["axis_dx0"])) == (41, 2400, 0.03), \
    "問題 YAML の分解能が §6 V4 の条件 (n_start 41・n_axis 2400・dx0 0.03) と違う"
assert list(p0.geometry["wall_fit_mono_r2"]) == MONO
R = float(p0.geometry["R"])
RT_UM = float(p0.spec["r_throat"]) * 1e6
res = {}
for arm, keys in ARMS.items():
    p, d, dt = run(keys)
    res[arm] = {"p": p, "d": d, "seconds_first": dt}
    m = d["moc"]
    print(f"{arm}: {dt:.1f} s  対 {m['pairs']}  反復 max {m['iters']['max']} mean {m['iters']['mean']:.3f}  "
          f"gate {m['gate']}", flush=True)


def cloud(d):
    tb = d["wall_inv"]
    return tb, tb[:, 0], tb[:, 1], tb[:, 2]


rec = {}
for arm in ARMS:
    d = res[arm]["d"]
    tb, x, r, th = cloud(d)
    s0, _ = joint_fit_wall(tb, R)                                     # 単調拘束なし
    bump = r3_piecewise_exact(s0, None, 0.0, 0.3)
    W = d["wall"]
    m = d["moc"]
    rec[arm] = {"x1": float(x[1]), "dtheta1_deg": float(np.degrees(th[1] - np.arctan(x[1] / R))),
                "bump_unconstrained_r2max_minus_1overR": float(bump["r2_max"] - 1.0 / R),
                "x_bump_unconstrained": float(bump["x_r2_max"]),
                "mono_first_point_dev_deg": float(np.degrees(np.arctan(W._spl(np.r_[x[1]], 1)[0]) - th[1])),
                "moc_exit_x_rt": float(x[-1]), "moc_exit_r_rt": float(r[-1]), "wall_x_e": float(W.x_e),
                "wall_r_e_rt": float(W.r(np.r_[W.x_e])[0]),
                "moc_pairs": m["pairs"], "moc_iters": m["iters"], "moc_resid": m["resid"], "moc_branch": m["branch"],
                "axis_limit_frac": m["axis_limit_frac"], "moc_gate": m["gate"], "theta_r": m["theta_r"],
                "unconverged_pairs": int(m["pairs"]["iter_nonfinite"] + m["pairs"]["iter_maxiter"]),
                "seconds_first": res[arm]["seconds_first"]}

# --- 設計壁の変化 (新・参考 − 現行) -----------------------------------------------------------------
W0 = res["current"]["d"]["wall"]
for arm in ("new", "ref_legacy_converge"):
    W = res[arm]["d"]["wall"]
    xe = min(W.x_e, W0.x_e)
    xs = np.unique(np.r_[np.linspace(-12.0, xe, 400001), W._spl.t[W._spl.t <= xe], W0._spl.t[W0._spl.t <= xe]])
    dr = (W.r(xs) - W0.r(xs)) * RT_UM
    dth = np.degrees(np.arctan(W.r(xs, 1)) - np.arctan(W0.r(xs, 1)))
    g3 = xs >= 0.3
    ts = (xs >= TS[0]) & (xs <= TS[1])
    rec[arm]["wall_change"] = {
        "max_abs_dr_um": float(np.abs(dr).max()), "x_max_abs_dr": float(xs[np.argmax(np.abs(dr))]),
        "max_abs_dr_um_x_ge_0p3": float(np.abs(dr[g3]).max()), "x_max_abs_dr_x_ge_0p3": float(xs[g3][np.argmax(np.abs(dr[g3]))]),
        "dr_at_common_exit_um": float(dr[-1]), "x_common_exit": float(xe),
        "exit_radius_change_um": float((W.r(np.r_[W.x_e])[0] - W0.r(np.r_[W0.x_e])[0]) * RT_UM),
        "exit_x_shift_rt": float(W.x_e - W0.x_e), "exit_x_shift_um": float((W.x_e - W0.x_e) * RT_UM),
        "max_abs_dtheta_test_section_deg": float(np.abs(dth[ts]).max()), "test_section_rt": list(TS),
        "max_abs_dtheta_all_deg": float(np.abs(dth).max()), "x_max_abs_dtheta_all": float(xs[np.argmax(np.abs(dth))])}
    print(f"{arm}: 設計壁の変化 {rec[arm]['wall_change']}", flush=True)

# --- 形状ゲート (新しい単調壁) ---------------------------------------------------------------------
dN, dC = res["new"]["d"], res["current"]["d"]
Wn, Wc = dN["wall"], dC["wall"]
tbn, xn, rn, thn = cloud(dN)
gate = {}
# S1 単調性 ([0, 1.5]、区間多項式で厳密)
m1 = r3_piecewise_exact(Wn._spl, None, *MONO)
gate["S1"] = {"status": "PASS" if (m1["r3_max"] <= 1e-6 and m1["r2_max_increase"] <= 1e-7 and m1["r2_max"] <= 1.0 / R + 1e-9)
              else "FAIL", "r3_max": m1["r3_max"], "r2_max_increase": m1["r2_max_increase"], "r2_max": m1["r2_max"],
              "limits": {"r3_max": 1e-6, "r2_max_increase": 1e-7, "r2_max": "1/R + 1e-9"}}
# S2 接続 (x = 0)
z = np.array([0.0])
jump = [float(abs(Wn.up.r(z, k)[0] - Wn._spl(z, k)[0])) for k in (0, 1, 2)]
gate["S2"] = {"status": "PASS" if max(jump) <= 1e-8 else "FAIL", "jump_r": jump[0], "jump_r1": jump[1], "jump_r2": jump[2],
              "limit": 1e-8}
# S3' (置き換え): x ≥ 0.3 の設計壁の変化 ≤ 10 µm、出口での変化 ≤ 0.1 µm
wc = rec["new"]["wall_change"]
ex = max(abs(wc["exit_radius_change_um"]), abs(wc["dr_at_common_exit_um"]))
gate["S3p"] = {"status": "PASS" if (wc["max_abs_dr_um_x_ge_0p3"] <= 10.0 and ex <= 0.1) else "FAIL",
               "max_abs_dr_um_x_ge_0p3": wc["max_abs_dr_um_x_ge_0p3"], "x": wc["x_max_abs_dr_x_ge_0p3"],
               "exit_change_um": ex, "exit_radius_change_um": wc["exit_radius_change_um"],
               "dr_at_common_exit_um": wc["dr_at_common_exit_um"], "limits_um": {"x_ge_0p3": 10.0, "exit": 0.1}}
# S4' (置き換え): 新しい点群の解析的下限 θ₁ − atan(x₁/R) と第 1 点の当てはめのずれ
bound = float(np.degrees(thn[1] - np.arctan(xn[1] / R)))
dth1 = float(np.degrees(np.arctan(Wn._spl(np.r_[xn[1]], 1)[0]) - thn[1]))
gap = abs(abs(dth1) - bound)
gate["S4p"] = {"status": "PASS" if gap <= 0.002 else "FAIL", "x1": float(xn[1]), "bound_deg": bound,
               "dtheta_first_deg": dth1, "gap_deg": gap, "limit_deg": 0.002}
# S6 物理壁 (生産経路: integral_delta_r → PhysicalNozzleWall、pw_ramp)。積分法の設定は throat_mono_shape_gate と同じ
kf = float(json.loads((C / "c2pin_solve_recal.json").read_text())["k_f"])
init_cfg = {"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}


def physical(p, d):
    _, drx, info = integral_delta_r(p, d, init_cfg)
    pw = p.geometry.get("pw_ramp")
    PW = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]),
                            _gam_or_gas(p), p.cp, dstar_x=None, offset="radial", delta_r_x=drx,
                            ramp=(None if pw is None else tuple(float(v) for v in pw)))
    return PW, drx, info


PWn, drn, infn = physical(res["new"]["p"], dN)
PWc, drc, infc = physical(res["current"]["p"], dC)
assert PWn.analytic and PWc.analytic, "物理壁が解析経路でない"
m6 = r3_piecewise_exact(lambda xq, k: PWn.r(xq, k), np.r_[Wn._spl.t, drn.spline.t], *MONO)
m6c = r3_piecewise_exact(lambda xq, k: PWc.r(xq, k), np.r_[Wc._spl.t, drc.spline.t], *MONO)
vn = PWn.validate()
xg6 = np.linspace(max(PWn.x_in, PWc.x_in), min(PWn.x_e, PWc.x_e), 800001)
d6 = (PWn.r(xg6) - PWc.r(xg6)) * RT_UM
gate["S6"] = {"status": "PASS" if (m6["r2_max_increase"] <= 0.002 and vn == []) else "FAIL",
              "r2_max_increase_phys": m6["r2_max_increase"], "limit": 0.002, "r2_max_phys": m6["r2_max"],
              "x_r2_max_phys": m6["x_r2_max"], "validate": vn,
              "current_physical": {"r2_max_increase": m6c["r2_max_increase"], "r2_max": m6c["r2_max"]},
              "physical_wall_change_um": {"max_abs": float(np.abs(d6).max()), "x": float(xg6[np.argmax(np.abs(d6))]),
                                          "at_exit": float(d6[-1])},
              "delta_r_exit": {"new": infn["delta_r_exit"], "current": infc["delta_r_exit"]},
              "throat": {"new": {"x": PWn.x_throat, "r": PWn.r_throat}, "current": {"x": PWc.x_throat, "r": PWc.r_throat}},
              "init_cfg": init_cfg}
# S8 滑らかさ ([0, 0.3]、厳密): 新しい単調壁 ≤ 現行の単調壁
m8, m8c = r3_piecewise_exact(Wn._spl, None, 0.0, 0.3), r3_piecewise_exact(Wc._spl, None, 0.0, 0.3)
gate["S8"] = {"status": "PASS" if (m8["int_r3sq"] <= m8c["int_r3sq"] and m8["r4_absmax"] <= m8c["r4_absmax"]) else "FAIL",
              "int_r3sq_new": m8["int_r3sq"], "int_r3sq_current": m8c["int_r3sq"],
              "r4_absmax_new": m8["r4_absmax"], "r4_absmax_current": m8c["r4_absmax"]}
# S7 既定のビット同一 (キー無しの設計が実装前のコミットと完全一致)
_SNIP = r'''
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from forge_design.evaluate.runner_axismach import design_chain, load_problem
import forge_design
assert forge_design.__file__.startswith(sys.argv[1]), forge_design.__file__
p = load_problem(sys.argv[3]); p.geometry["initial_line_run"] = sys.argv[4]
d = design_chain(p)
np.savez(sys.argv[2], wall_inv=d["wall_inv"], c=d["wall"]._spl.c, t=d["wall"]._spl.t)
'''
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    subprocess.run(f"git -C {ROOT} archive {BASE} design/forge_design solver_density_cuda/data | tar -x -C {td}",
                   shell=True, check=True, capture_output=True)
    (td / "snip.py").write_text(_SNIP)
    r7 = subprocess.run([sys.executable, str(td / "snip.py"), str(td / "design"), str(td / "base.npz"), str(PROB),
                         str((RUNS / p0.geometry["initial_line_run"]).resolve())], capture_output=True, text=True)
    if r7.returncode != 0:
        raise RuntimeError("S7: 実装前のコミットでの設計に失敗\n" + r7.stderr[-2000:])
    b = np.load(td / "base.npz")
    same = {"wall_inv": bool(np.array_equal(b["wall_inv"], dC["wall_inv"])),
            "coef": bool(np.array_equal(b["c"], Wc._spl.c)), "knots": bool(np.array_equal(b["t"], Wc._spl.t))}
gate["S7"] = {"status": "PASS" if all(same.values()) else "FAIL", "identical": same, "base_commit": BASE}
unconv = rec["new"]["unconverged_pairs"]
gate["unconverged_pairs"] = {"status": "PASS" if unconv == 0 else "FAIL", "n": unconv}
for k, v in gate.items():
    print(f"{k:5s} {v['status']}  " + "  ".join(f"{a}={b:.6g}" if isinstance(b, float) else f"{a}={b}" for a, b in v.items()
                                              if a != "status" and not isinstance(b, (dict, list))), flush=True)
cand = all(v["status"] == "PASS" for v in gate.values())

# --- V6: 所要時間 (各 N_TIME 回の中央値) と反復回数の分布 ----------------------------------------------
tim = {}
for arm in ("current", "new"):
    ts_ = []
    for _ in range(N_TIME):
        _, _, dt = run(ARMS[arm])
        ts_.append(dt)
    tim[arm] = {"seconds": ts_, "median": float(np.median(ts_))}
ratio = tim["new"]["median"] / tim["current"]["median"]
v6 = {"timing": tim, "ratio_new_over_current": ratio, "within_5x": bool(ratio <= 5.0),
      "iters_new": dN["moc"]["iters"], "iters_ref_legacy_converge": res["ref_legacy_converge"]["d"]["moc"]["iters"],
      "unconverged_pairs_new": unconv,
      "note": "design_chain 1 回 (CFD ピン初期線の読込・逆 MOC・壁 QA・joint 当てはめ・単調拘束) の壁時計。同じ機械で連続して計測"}
print(f"V6: 現行 {tim['current']['median']:.2f} s / 新 {tim['new']['median']:.2f} s (比 {ratio:.2f})", flush=True)

commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                       capture_output=True, text=True).stdout.strip()
out = {"plan": "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V4・V6", "problem": PROB.name,
       "commit": commit, "design_tree_dirty": dirty, "case_runs": str(RUNS), "R": R, "r_t_um": RT_UM,
       "arms": {k: (None if v is None else {"moc_axis_limit": v[0], "moc_corrector": v[1]}) for k, v in ARMS.items()},
       "record": rec, "gates": gate,
       "candidate_conditions": "未収束の対 0、S1・S2・S6・S7・S8 PASS、S3′ (x≥0.3 ≤10 µm・出口 ≤0.1 µm)、S4′ (≤0.002°) — plan §6 V4",
       "VERDICT_candidate": "候補の条件を満たす" if cand else "候補の条件を満たさない (諮問)",
       "fail_items": [k for k, v in gate.items() if v["status"] != "PASS"], "V6": v6}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
print(f"VERDICT (V4 候補): {out['VERDICT_candidate']}" + (f" (外れ: {', '.join(out['fail_items'])})" if out["fail_items"] else "")
      + f"  -> {OUT}")
