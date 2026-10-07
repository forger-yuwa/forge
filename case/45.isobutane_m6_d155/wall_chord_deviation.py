"""plan tooling-nozzle-wall-single-bspline §4.3 の「CAD の形と CFD の形の関係」の実測 (2026-10-07)。
run_0147 の物理壁と、同じメッシュ設定の壁節点を直線でつないだ多角形 (CFD が解く形) の半径方向の差
(弦 − 曲線、符号付き。正 = 弦が曲線より外 = 流路が広い側) を、全壁辺について測る。
あわせて、δ_r の表の範囲がランプ開始から出口までを覆うか (範囲外は端値延長で導関数 0) を確かめる。
usage: [CASE_RUNS=<run_0147・run_0062 のある case dir>] python3 wall_chord_deviation.py → _band_ab/wall_chord_deviation.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, delta_r_from_table, _gam_or_gas, mesh_params  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402
from forge_design.meshing.mesh2d import _x_stations  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
RUN = RUNS / "run_0147_ns_mono_final"
p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
d = design_chain(p)
scale = float(p.spec["r_throat"])
tbl = np.loadtxt(RUN / "delta_r_initial.csv", delimiter=",", skiprows=1)
f = delta_r_from_table(tbl[:, 0], tbl[:, 1])
ramp = tuple(float(v) for v in p.geometry["pw_ramp"])
W = PhysicalNozzleWall(d["wall"], d["wall_inv"], scale, float(p.spec["Pt"]), float(p.spec["Tt"]),
                       _gam_or_gas(p), p.cp, offset="radial", delta_r_x=f, ramp=ramp, upstream="ramp")   # 2026-10-07 に既定が poly
wp = np.loadtxt(RUN / "wall_physical.csv", delimiter=",", skiprows=1)
rep = float(np.abs(W.r(wp[:, 0] / scale) * scale - wp[:, 1]).max())

info = json.loads((RUN / "prepare_info.json").read_text())["mesh"]
mp = mesh_params(p, scale, ni=int(info["ni"]), nj=int(info["nj"]), wall_first_frac=float(info["wall_first_frac"]))
xs = _x_stations(W.x_in, W.x_e, mp.ni, mp.throat_refine, mp.throat_width, mp.local_center, mp.local_refine, mp.local_width)
rw = W.r(xs)
um = scale * 1e6
dev_min, dev_max, rows = [], [], []
for i in range(len(xs) - 1):
    xa, xb = xs[i], xs[i + 1]
    q = np.linspace(xa, xb, 201)
    chord = rw[i] + (rw[i + 1] - rw[i]) * (q - xa) / (xb - xa)
    e = (chord - W.r(q)) * um
    rows.append((xa, xb, float(e.min()), float(q[np.argmin(e)]), float(e.max()), float(q[np.argmax(e)])))
rows = np.array(rows)
i_min, i_max = int(np.argmin(rows[:, 2])), int(np.argmax(rows[:, 4]))
h = np.diff(xs)
out = {
    "reproduce_wall_physical_csv_m": rep,
    "mesh": {"ni": mp.ni, "nj": mp.nj, "n_wall_edges": int(len(xs) - 1),
             "dx_min_rt": float(h.min()), "dx_max_rt": float(h.max()),
             "dx_at_throat_rt": float(h[np.argmin(np.abs(xs[:-1]))])},
    "chord_minus_curve_um": {
        "min": float(rows[i_min, 2]), "x_of_min_rt": float(rows[i_min, 3]), "edge_dx_at_min_rt": float(rows[i_min, 1] - rows[i_min, 0]),
        "max": float(rows[i_max, 4]), "x_of_max_rt": float(rows[i_max, 5]), "edge_dx_at_max_rt": float(rows[i_max, 1] - rows[i_max, 0])},
    "by_region_um": {},
    "delta_r_table": {"x_lo": float(tbl[0, 0]), "x_hi": float(tbl[-1, 0]), "ramp_lo": ramp[0], "x_e": float(W.x_e),
                      "covers_ramp_to_exit": bool(tbl[0, 0] <= ramp[0] and tbl[-1, 0] >= W.x_e),
                      "x_hi_minus_x_e": float(tbl[-1, 0] - W.x_e)},
}
for name, lo, hi in (("pipe+contraction [-12.5,0)", -12.5, 0.0), ("throat [0,1)", 0.0, 1.0), ("expansion [1,x_e]", 1.0, 1e9)):
    m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
    out["by_region_um"][name] = {"min": float(rows[m, 2].min()), "max": float(rows[m, 4].max()),
                                 "dx_max_rt": float((rows[m, 1] - rows[m, 0]).max())}
(C / "_band_ab").mkdir(exist_ok=True)
(C / "_band_ab/wall_chord_deviation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
