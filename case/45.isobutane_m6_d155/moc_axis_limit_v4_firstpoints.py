"""plan discretization-moc-axis-limit-and-corrector §6 V4 の補足 (2026-10-07): S4′ の余裕 (0.00004°) が小さい理由を見る。
始点付近の MOC 点の角度・位置と、当てはめ (拘束なし・単調拘束) を、円弧 r = 1 + x²/2R からの差で並べる。
usage: [CASE_RUNS=<run_0062 のある case dir>] python3 moc_axis_limit_v4_firstpoints.py > _band_ab/moc_axis_limit_v4_firstpoints.txt"""
import os
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.geometry.wall_axismach import joint_fit_wall  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
R = 2.0
for arm, keys in (("current", {}), ("new", dict(moc_axis_limit="analytic", moc_corrector="converge"))):
    p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    p.geometry.update(keys)
    d = design_chain(p)
    tb = d["wall_inv"]; x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    s0, _ = joint_fit_wall(tb, R)
    Wm = d["wall"]
    circ_r = 1 + x ** 2 / (2 * R)          # 円弧 (r″ = 1/R) の近似 (x ≪ R)
    circ_th = np.arctan(x / R)
    print(f"--- {arm}")
    print(" i      x        MOC dθ[deg]  MOC dr[um]   free dθ   free dr   mono dθ   mono dr   (d = 円弧からの差; free/mono は当てはめ − 円弧)")
    for i in range(0, 12):
        fr, fth = s0(x[i]), np.degrees(np.arctan(s0(x[i], 1)))
        mr, mth = Wm._spl(x[i]), np.degrees(np.arctan(Wm._spl(x[i], 1)))
        print(f"{i:2d} {x[i]:9.5f} {np.degrees(th[i] - circ_th[i]):11.5f} {(r[i] - circ_r[i]) * 76653.9:11.4f}"
              f" {fth - np.degrees(circ_th[i]):9.5f} {(fr - circ_r[i]) * 76653.9:9.4f}"
              f" {mth - np.degrees(circ_th[i]):9.5f} {(mr - circ_r[i]) * 76653.9:9.4f}")
    xg = np.linspace(0, 0.3, 3001)
    print("mono r″ at x = 0, 0.0125, 0.025, 0.05, 0.1:", [round(float(Wm._spl(v, 2)), 6) for v in (0, 0.0125, 0.025, 0.05, 0.1)])
    print("free r″ at same:", [round(float(s0(v, 2)), 6) for v in (0, 0.0125, 0.025, 0.05, 0.1)])
