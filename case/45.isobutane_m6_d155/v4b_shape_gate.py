"""V4b の形状ゲート (事前登録: plan verification-m6-axis-wave-mesh-su2 §5.1 #16、diagnostician 2026-10-05 2 回目)。CFD 0 step。
G1 x ≤ −L_b で V0 と |Δr| ≤ 1e-12 / G2 接合 x=−L_b の r, r′, r″ の跳び ≤ 1e-8 / G3 blend 内で r″ に内部極値なし・r‴ ≤ 1.0 / G4 x<x_min で r′ ≤ 0、x_min ∈ [−0.01, 0]、1−r_min ≤ 1e-5
G5 x ≥ 0 のスプライン係数が V4 と一致 / G6 (併記) blend 内の max|Δr| vs V0 と Δr(−0.25/−0.5/−1)。
usage: design/.venv-opt/bin/python v4b_shape_gate.py → _band_ab/v4b_shape_gate.json
"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from wall_v4 import build_v4, build_v4b  # noqa: E402

d = design_chain(load_problem(C / "problem_d155_euler_c2final_n2400.yaml")); V0 = d["wall"]; V4 = build_v4(d); W = build_v4b(d); Lb = W.up.L_b
xl = np.linspace(-V0.up.L_U, -Lb, 20001); e = 1e-9
G = {}
G["G1"] = dict(maxdr=float(np.abs(W.r(xl) - V0.r(xl)).max())); G["G1"]["ok"] = G["G1"]["maxdr"] <= 1e-12
jumps = [float(W.r(np.r_[-Lb + e], k)[0] - W.r(np.r_[-Lb - e], k)[0]) for k in range(3)]
G["G2"] = dict(jumps=jumps, ok=bool(max(abs(v) for v in jumps) <= 1e-8))   # x=−L_b ± 1e-9 で評価
xb = np.linspace(-Lb, 0, 30001); r2 = W.r(xb, 2); r3 = W.r(xb, 3)
G["G3"] = dict(r2_interior_extrema=int(np.count_nonzero(np.diff(np.sign(np.diff(r2))) != 0)), r3_max=float(np.abs(r3).max()), r2_range=[float(r2.min()), float(r2.max())])
G["G3"]["ok"] = G["G3"]["r2_interior_extrema"] == 0 and G["G3"]["r3_max"] <= 1.0
rb = W.r(xb); i = int(np.argmin(rb))
G["G4"] = dict(x_min=float(xb[i]), one_minus_rmin=float(1 - rb[i]), r1_max_before=float(W.r(xb[:i], 1).max()))
G["G4"]["ok"] = bool(-0.01 <= xb[i] <= 0 and 1 - rb[i] <= 1e-5 and W.r(xb[:i], 1).max() <= 0)
G["G5"] = dict(same_spline=bool(np.array_equal(W._spl.c, V4._spl.c))); G["G5"]["ok"] = G["G5"]["same_spline"]
dr = W.r(xb) - V0.r(xb); j = int(np.argmax(np.abs(dr)))
G["G6"] = dict(blend_maxdr_vs_V0=float(dr[j]), x_at=float(xb[j]), dr_at={str(x): float(W.r(np.r_[x])[0] - V0.r(np.r_[x])[0]) for x in (-0.25, -0.5, -1.0)},
               V4_maxdr_vs_V0_upstream=float(np.abs(V4.r(xl) - V0.r(xl)).max()))
G["all_ok"] = all(G[k]["ok"] for k in ("G1", "G2", "G3", "G4", "G5"))
(C / "_band_ab/v4b_shape_gate.json").write_text(json.dumps(G, indent=1, ensure_ascii=False))
print(json.dumps(G, indent=1, ensure_ascii=False))
