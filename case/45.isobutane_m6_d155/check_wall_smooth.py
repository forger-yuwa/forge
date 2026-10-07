"""壁の入力を滑らかにした修正 (plan verification-m6-axis-wave-mesh-su2 §5.1 #8a) の事前登録基準 ①〜④ を確かめる。
B1 の表 (_band_ab/edge_tests_v2/T5_F3_c1.25/delta_r_next.csv) で、旧 (np.interp) と新 (delta_r_from_table) と δ=0 の壁を作って比べる。
usage: design/.venv-opt/bin/python check_wall_smooth.py → _band_ab/wall_smooth_check.json, _band_ab/points_B1_smoothwall*.csv
"""
import json, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas, delta_r_from_table
from forge_design.geometry.wall_axismach import PhysicalNozzleWall

C = Path(__file__).resolve().parent
p = load_problem(C / "problem_d155_ns_rt77p02.yaml"); d = design_chain(p); S = float(p.spec["r_throat"])
t = np.loadtxt(C / "_band_ab/edge_tests_v2/T5_F3_c1.25/delta_r_next.csv", delimiter=",", skiprows=1)
mk = lambda f: PhysicalNozzleWall(d["wall"], d["wall_inv"], S, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp,
                                  offset="radial", delta_r_x=f)
W = {"old_interp": mk(lambda x: np.interp(x, t[:, 0], t[:, 1])), "new_quintic": mk(delta_r_from_table(t[:, 0], t[:, 1])),
     "delta0": mk(lambda x: np.zeros_like(np.asarray(x, dtype=float)))}
xw = np.arange(0.05, 93.0, 0.001)
ri2 = d["wall"].r(xw, 2)
k = int(round(0.5 / 0.001))


def hf(w):
    dd = w.r(xw, 2) - ri2
    return dd - np.convolve(dd, np.ones(k) / k, "same")


H = {n: hf(w) for n, w in W.items()}
out = {"hf_max": {}}
for n, h in H.items():
    out["hf_max"][n] = {f"[{a},{b})": float(np.abs(h[(xw >= a) & (xw < b) & (xw > xw[0] + 0.25) & (xw < xw[-1] - 0.25)]).max())
                        for a, b in ((0.3, 1), (1, 2), (2, 6), (6, 20), (20, 60), (60, 93))}
new, d0 = out["hf_max"]["new_quintic"], out["hf_max"]["delta0"]
c1 = (max(new["[2,6)"], new["[6,20)"], new["[20,60)"], new["[60,93)"]) <= 1e-5 and new["[1,2)"] <= 1e-4 and new["[0.3,1)"] <= 1.5 * d0["[0.3,1)"])
# ② 節点位置差 (B1 の格子の x 断面)
with h5py.File(C / "run_0046_ns_band_edge/nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
X = (nc[:, 0] / S).reshape(1250, 97)[:, 0]
dpos = np.abs(W["new_quintic"].r(X) - W["old_interp"].r(X))
c2 = float(dpos.max()) <= 5e-6
out["node_pos_diff_max_rt"] = float(dpos.max()); out["node_pos_diff_x"] = float(X[np.argmax(dpos)])


# ③ 格子節点の離散曲率の高周波 (参考)
def disc_hf(w):
    R = w.r(X); Ri = d["wall"].r(X)
    c = lambda R_: np.r_[0, 2 * ((R_[2:] - R_[1:-1]) / (X[2:] - X[1:-1]) - (R_[1:-1] - R_[:-2]) / (X[1:-1] - X[:-2])) / (X[2:] - X[:-2]), 0]
    dc = c(R) - c(Ri); sm = np.convolve(dc, np.ones(11) / 11, "same"); m = (X >= 6) & (X < 90)
    return float(np.abs((dc - sm)[m]).max())


out["mesh_disc_curv_hf_6_90"] = {n: disc_hf(w) for n, w in W.items()}
# ④ 点列 (格子断面と 0.01 r_t の密な点) を新壁から作り、密な点列で ① を再検査 (5 次補間で戻して r″)
from scipy.interpolate import make_interp_spline
xd = np.arange(float(W["new_quintic"].x_throat), float(W["new_quintic"].x_e), 0.01)
rd = W["new_quintic"].r(xd)
np.savetxt(C / "_band_ab/points_B1_smoothwall_dense.csv", np.c_[xd * S, rd * S, xd, rd], delimiter=",",
           header="x_m,r_m,x_over_rt,r_over_rt  (B1 の表を 5 次補間で載せた物理壁; スロート〜x_F を 0.01 r_t 間隔; r_t 0.07702 m)", comments="")
np.savetxt(C / "_band_ab/points_B1_smoothwall_mesh.csv", np.c_[X * S, W["new_quintic"].r(X) * S, X, W["new_quintic"].r(X)], delimiter=",",
           header="x_m,r_m,x_over_rt,r_over_rt  (同じ壁を格子 1250 断面で)", comments="")
sp = make_interp_spline(xd, rd, k=5)
m4 = (xw >= 2) & (xw < xd[-1] - 0.5)
dd4 = sp(xw[m4], 2) - ri2[m4]; h4 = dd4 - np.convolve(dd4, np.ones(k) / k, "same")
out["points_dense_hf_2_end"] = float(np.abs(h4[k:-k]).max())
c4 = out["points_dense_hf_2_end"] <= 1e-5
out["criteria"] = dict(c1_hf=c1, c2_node_pos=c2, c3_reference_only=True, c4_points=c4, all=bool(c1 and c2 and c4))
(C / "_band_ab/wall_smooth_check.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
