"""G1x の x の粗さが Q_w の壁に沿った積分 (後処理) に与える誤差: G0 の解の q_w(x) を G1x の壁の x に線形補間して積分し直す。"""
import sys, numpy as np, h5py
from pathlib import Path
sys.argv = ["core_grid_mesh.py", "select"]
sys.path.insert(0, "/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155")
import core_grid_mesh as CG
import cold_xcheck as XC
I = Path("/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155")
XC.OUTD = I / "_band_ab/cold_pair"; yb_x, yb = XC.common_yb()
wall = CG.Wall(I / "_band_ab/prod_confirm/prep"); S = wall.scale
with h5py.File(I / "run_0353_m9_L5/res_115000.h5") as h:
    xy = np.asarray(h["MESH/COORD"][:], dtype=float).reshape(-1, 3)[:, :2]
    v = [np.asarray(h["VALUE/" + k][:], dtype=float) for k in ("ro", "Ux", "Uy", "T", "k")]
o = XC.reduce_fields(xy, *v, False, 4719, 121, S, yb_x, yb)
xw = o["x_w"]; qw = o["q_w"]; rw = xy.reshape(4719, 121, 2)[:, -1, 1] / S
def Qint(x, q, r):
    P = np.c_[x, r] * S; ds = np.linalg.norm(np.diff(P, axis=0), axis=1); f = 2 * np.pi * q * r * S
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * ds))
Q0 = Qint(xw, qw, rw)
P, _ = CG.generate(wall, CG.mesh_block(75, 0.03, xcoarse=True))
x1 = P[:, -1, 0]; r1 = P[:, -1, 1]
Q1 = Qint(x1, np.interp(x1, xw, qw), r1)
print(f"Q_w: G0 の station {Q0:.6e} W、G1x の station に補間 {Q1:.6e} W、差 {100 * (Q1 / Q0 - 1):+.4f} %")
for a, b in ((-12.5, -2), (-2, 2), (2, 12), (12, 96)):
    m0 = (xw >= a) & (xw < b); m1 = (x1 >= a) & (x1 < b)
    print(f"  x∈[{a},{b}): {100 * (Qint(x1[m1], np.interp(x1[m1], xw, qw), r1[m1]) / Qint(xw[m0], qw[m0], rw[m0]) - 1):+.4f} %")
