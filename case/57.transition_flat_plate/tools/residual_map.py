#!/usr/bin/env python3
"""残差の下げ止まりが**どこ**にあるかを出す (output.level 2 の run。res_* 場を使う)。

各保存量について、残差二乗和への寄与を領域別 (前縁近傍 / 入口・出口の列 / 壁第一層 / それ以外) に分け、最大の節点の位置を出す。
usage: residual_map.py RUN_DIR [--le 0.01]
"""
import argparse, glob, os
import h5py, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--le", type=float, default=0.01, help="前縁近傍とみなす半径 [m]")
a = ap.parse_args()
f = sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda p: int(os.path.basename(p)[4:-3]))[-1]
with h5py.File(f, "r") as h:
    c = h["MESH/COORD"][:].reshape(-1, 3); V = {k: h["VALUE"][k][:].astype(float) for k in h["VALUE"] if k.startswith("res_")}; wd = h["VALUE/wall_dist"][:]
x, y = c[:, 0], c[:, 1]
reg = {"前縁 r<%g m" % a.le: np.hypot(x, y) < a.le, "入口列": x <= x.min() + 1e-9, "出口列": x >= x.max() - 1e-9, "上面列": y >= y.max() - 1e-9}
print(f"{f}")
print(f"{'残差':<13}{'rms':>11}" + "".join(f"{k:>14}" for k in reg) + f"{'上位 10 節点':>14}   最大の位置 (x, y) [m]")
for k, r in V.items():
    r2 = r ** 2; tot = r2.sum()
    if tot == 0: continue
    i = np.argsort(r2)[::-1]
    print(f"{k:<13}{np.sqrt(tot/len(r)):11.3e}" + "".join(f"{r2[m].sum()/tot:14.1%}" for m in reg.values()) + f"{r2[i[:10]].sum()/tot:14.1%}   ({x[i[0]]:.4f}, {y[i[0]]:.2e})")
