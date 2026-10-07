"""plan tooling-nozzle-wall-single-bspline §4.1 の事前試算 (2026-10-07、実装前、CFD 0 step)。
run_0147 の物理壁を全域 1 本の x の 5 次 B-spline に作り直したときの誤差と係数の数を測る。
厳密に書ける区間 (直管・H・H+δ_r・S+δ_r) はノットの和集合、ランプ [−11,−6] (10 次) は同じ空間への最良近似。
usage: [DESIGN_DIR=<design/ の写し>] [CASE_RUNS=<run_0147・run_0062 のある case dir>] python3 wall_single_bspline_probe.py
  → _band_ab/wall_single_bspline_probe.json
2026-10-07 の実測は、同じツリーで MOC の実装が進行中だったため HEAD (934d3086) の design/ を git archive した写しで回した。"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import lstsq

C = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("DESIGN_DIR", str(C.parents[1] / "design")))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, delta_r_from_table, _gam_or_gas  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
d = design_chain(p)
scale = float(p.spec["r_throat"])
tbl = np.loadtxt(RUNS / "run_0147_ns_mono_final/delta_r_initial.csv", delimiter=",", skiprows=1)
f = delta_r_from_table(tbl[:, 0], tbl[:, 1])
W = PhysicalNozzleWall(d["wall"], d["wall_inv"], scale, float(p.spec["Pt"]), float(p.spec["Tt"]),
                       _gam_or_gas(p), p.cp, offset="radial", delta_r_x=f,
                       ramp=tuple(float(v) for v in p.geometry["pw_ramp"]), upstream="ramp")   # 2026-10-07 に既定が poly

# 生産壁の再現確認 (run_0147 の wall_physical.csv、m 単位)
wp = np.loadtxt(RUNS / "run_0147_ns_mono_final/wall_physical.csv", delimiter=",", skiprows=1)
rep = float(np.abs(W.r(wp[:, 0] / scale) * scale - wp[:, 1]).max())
x_in, x_e = float(W.x_in), float(W.x_e)
dw = d["wall"]
print(f"x_in {x_in} x_e {x_e} L_U {dw.up.L_U} r_U {dw.up.r_U}  reproduce wall_physical.csv max|dr| = {rep:.3e} m")

k = 5
ramp = tuple(float(v) for v in p.geometry["pw_ramp"])
junctions = [-float(dw.up.L_U), ramp[0], ramp[1], 0.0]
tS = np.asarray(dw._spl.t)
tD = np.asarray(f.spline.t)
inner = lambda t, a, b: np.unique(t[(t > a) & (t < b)])   # noqa: E731
segs = [(x_in, junctions[0], np.array([])),                          # 直管
        (junctions[0], junctions[1], np.array([])),                  # H
        (junctions[1], junctions[2], inner(tD, *junctions[1:3])),    # ランプ (当てはめ)
        (junctions[2], junctions[3], inner(tD, junctions[2], 0.0)),  # H + δ_r
        (0.0, x_e, np.union1d(inner(tS, 0.0, x_e), inner(tD, 0.0, x_e)))]  # S + δ_r
knots = [x_in] * (k + 1)
for i, (a, b, ins) in enumerate(segs):
    knots += list(ins)
    knots += [b] * (3 if i < len(segs) - 1 else k + 1)
t = np.array(knots)
nc = len(t) - k - 1
dist = np.unique(t)
print(f"coefficients {nc}, distinct knots {len(dist)}, min knot gap {np.diff(dist).min():.3e}")

# 各ノット区間の Gauss 点 (8 点) で最小二乗
gx, _ = np.polynomial.legendre.leggauss(8)
xs = np.concatenate([(a + b) / 2 + (b - a) / 2 * gx for a, b in zip(dist[:-1], dist[1:])])
A = BSpline.design_matrix(xs, t, k).toarray()
c, *_ = lstsq(A, W.r(xs))
B = BSpline(t, c, k)

out = {"reproduce_wall_physical_csv_m": rep, "n_coef": int(nc), "n_distinct_knots": int(len(dist)),
       "min_knot_gap_rt": float(np.diff(dist).min()), "segments": []}
names = ["pipe", "H", "ramp (fit)", "H+dr", "S+dr"]
for (a, b, _), nm in zip(segs, names):
    m = (dist >= a) & (dist <= b)
    dd = dist[m]
    # 区間ごとに 41 点 (端を少し内側に) + Gauss 点
    q = np.concatenate([np.linspace(lo, hi, 41)[1:-1] for lo, hi in zip(dd[:-1], dd[1:])])
    row = {"seg": nm, "range": [a, b]}
    for n in range(4):
        e = np.abs(B(q, n) - W.r(q, n)) if n <= 3 else None
        row[f"max_d{n}"] = float(e.max())
        row[f"x_at_max_d{n}"] = float(q[np.argmax(e)])
    out["segments"].append(row)
    print(nm, {kk: (f"{v:.3e}" if isinstance(v, float) else v) for kk, v in row.items() if kk.startswith("max")})
# 継ぎ目の跳び (左右極限)
eps = 1e-12
for xj in junctions:
    jumps = [float(abs(B(xj + eps, n) - B(xj - eps, n))) for n in range(4)]
    print(f"junction {xj}: jump d0..d3 = " + " ".join(f"{v:.2e}" for v in jumps))
out["rt_um"] = scale * 1e6
(C / "_band_ab").mkdir(exist_ok=True)
(C / "_band_ab/wall_single_bspline_probe.json").write_text(json.dumps(out, indent=1))
