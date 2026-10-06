"""plan tooling-nozzle-throat-monotone-r2 の result 段諮問 (notes/reviews/2026-10-07-throat-mono-result-interpretation-diagnose.md) の判別 A/B (0 step)。
変えるのは境界層初期化 (initializer) の供給経路だけ。単調壁の生産問題 YAML・環境は同じ。
  A: 検証時の経路 = prep_c2pin.py と同じ initializer {"model": "contur", "a_crocco": 1, "cf_scale": k_f, "n_scale": 1} (k_f は c2pin_solve_recal.json)
  B: 生産の入口 `deltastar_loop --init-integral` が渡す initializer (resolve_integral_initializer で YAML から解決。
     2026-10-07 の修正前は常に {"model": "contur"} で、出口半径が 2.4 mm ずれた)
  参考: run_0147 が実際に使った物理壁 (wall_physical.csv、AWS で作成。ローカルとは numpy の版で約 1e-5 r_t の差が既知)
比較: 設計 spline (A・B で同一のはず)、δ_r(x)、物理壁 r(x) の差 (最大・位置)、出口半径。
usage: [CASE_RUNS=<run_0062 のある case dir>] python3 throat_mono_entry_ab.py [RUN_0147_DIR] → _band_ab/throat_mono_entry_ab.json
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, integral_delta_r, _gam_or_gas  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
R147 = Path(sys.argv[1]) if len(sys.argv) > 1 else RUNS / "run_0147_ns_mono_final"
KF = json.loads((C / "c2pin_solve_recal.json").read_text())["k_f"]
p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
d = design_chain(p)
S = float(p.spec["r_throat"])
INIT = {"A_prep_c2pin": {"model": "contur", "a_crocco": 1.0, "cf_scale": KF, "n_scale": 1.0},
        "B_deltastar_loop": None}   # 生産の入口 deltastar_loop.resolve_integral_initializer で解決する (下)
from forge_design.feedback.deltastar_loop import resolve_integral_initializer  # noqa: E402
INIT["B_deltastar_loop"] = resolve_integral_initializer(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
W = {}
for k, init in INIT.items():
    _, drx, info = integral_delta_r(p, d, init)
    W[k] = (PhysicalNozzleWall(d["wall"], d["wall_inv"], S, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp,
                               offset="radial", delta_r_x=drx, ramp=tuple(float(v) for v in p.geometry["pw_ramp"])), drx, info)
x = np.linspace(-12.0, float(d["wall_inv"][-1, 0]), 200001)
rA, rB = W["A_prep_c2pin"][0].r(x), W["B_deltastar_loop"][0].r(x)
out = {"k_f": KF, "initializers": INIT,
       "A_vs_B": {"max_abs_dr_m": float(np.abs(rB - rA).max() * S), "x_of_max_rt": float(x[np.argmax(np.abs(rB - rA))]),
                  "exit_radius_A_m": float(rA[-1] * S), "exit_radius_B_m": float(rB[-1] * S),
                  "delta_r_exit_A_rt": float(W["A_prep_c2pin"][1](np.r_[x[-1]])[0]), "delta_r_exit_B_rt": float(W["B_deltastar_loop"][1](np.r_[x[-1]])[0])}}
f = R147 / "wall_physical.csv"
if f.is_file():
    hdr = f.read_text().splitlines()[0].strip()
    if hdr != "x_m,r_m":
        raise SystemExit(f"{f}: 列名 {hdr!r} が想定 (x_m,r_m) と違う")
    w = np.loadtxt(f, delimiter=",", skiprows=1)
    xw, rw = w[:, 0] / S, w[:, 1] / S          # m → r_t
    m = (xw >= x[0]) & (xw <= x[-1])
    for k in INIT:
        rr = W[k][0].r(xw[m])
        out[f"{k}_vs_run0147"] = {"max_abs_dr_m": float(np.abs(rr - rw[m]).max() * S), "x_of_max_rt": float(xw[m][np.argmax(np.abs(rr - rw[m]))]),
                                  "n_points": int(m.sum())}
    out["run0147_wall_physical_header"] = f.read_text().splitlines()[0]
(C / "_band_ab").mkdir(exist_ok=True)
(C / "_band_ab/throat_mono_entry_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
