#!/usr/bin/env python3
"""翼列の入口から前縁までの自由流乱れ Tu(x) と粘性比 mu_t/mu(x) を出す (ピッチ方向平均)。遷移モデルの入口条件の感度を読むための診断。

usage: tu_upstream.py RUN_DIR [RUN_DIR ...]
"""
import sys, glob, os
import numpy as np, h5py
for run in sys.argv[1:]:
    fs = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
    with h5py.File(fs[-1], "r") as h:
        c = h["MESH/COORD"][:].reshape(-1, 3); V = {k: h["VALUE"][k][:].astype(float) for k in ("k", "Ux", "Uy", "vis_lam", "vis_turb", "wall_dist")}
    xw = c[V["wall_dist"] <= 0.0, 0].min(); x0 = c[:, 0].min()
    print(f"{run}: inlet x={x0:.4f}  leading edge x={xw:.4f}  ({fs[-1]})")
    for fr in (0.0, 0.25, 0.5, 0.75, 0.9, 0.97):
        xs = x0 + fr * (xw - x0); m = np.abs(c[:, 0] - xs) < 0.004 * (xw - x0) + 5e-4
        U = np.hypot(V["Ux"][m], V["Uy"][m]); tu = 100 * np.sqrt(2 * V["k"][m] / 3) / U
        print(f"   x-x_in = {xs-x0:7.4f} m ({fr:4.0%})  Tu = {tu.mean():5.2f} %   mu_t/mu = {(V['vis_turb'][m]/V['vis_lam'][m]).mean():7.2f}   |U| = {U.mean():6.1f}  n={m.sum()}")
