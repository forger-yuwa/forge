#!/usr/bin/env python3
"""3D 半幅 run の res_*.h5 系列で、対称面輪郭壁の p/p0 (x=16/46/85/95 mm) と角線 T−Tt を表にする。
usage: wall_series.py RUN_DIR [RUN_DIR2 ...]  (複数は step を通しで並べる)"""
import sys, glob, h5py, numpy as np
Tt, p0 = 286.65, 59070.0
rows = []; base = 0
for d in sys.argv[1:]:
    files = sorted(glob.glob(d.rstrip('/') + "/res_[0-9]*.h5"), key=lambda s: int(s.split('_')[-1][:-3]))
    f0 = h5py.File(files[-1], "r"); c = np.array(f0["MESH/COORD"]).reshape(-1, 3); wd = np.array(f0["VALUE/wall_dist"])
    zmax = c[:, 2].max(); xs = np.round(c[:, 0], 7); ux = np.unique(xs); wall = wd <= 0
    sel = []
    for xv in [0.0164, 0.0456, 0.0848, 0.0946]:
        xv = ux[np.argmin(np.abs(ux - xv))]; col = np.where(xs == xv)[0]; cw = col[wall[col]]
        symc = cw[np.abs(c[cw, 2] - zmax) < 1e-7]; sel.append(("w%.0f" % (xv * 1e3), symc[np.argmax(c[symc, 1])]))
        sel.append(("c%.0f" % (xv * 1e3), cw[np.argmax(c[cw, 1] * 1e3 - c[cw, 2] * 1e4)]))
    for fn in files:
        st = int(fn.split('_')[-1][:-3])
        if st == 0 and base > 0: continue
        f = h5py.File(fn, "r"); T = np.array(f["VALUE/T"]); P = np.array(f["VALUE/P"])
        rows.append((base + st, d.split('/')[-1][:8], [P[j] / p0 for k, j in sel if k[0] == 'w'], [T[j] - Tt for k, j in sel if k[0] == 'c']))
    base += int(files[-1].split('_')[-1][:-3])
print("  step   run      | wall p/p0 @x=16   46     85     95   | corner T-Tt @16   46    85    95")
for st, r, pw, tc in rows:
    print(f"{st:6d} {r:9s} | " + " ".join(f"{v:.4f}" for v in pw) + " | " + " ".join(f"{v:+5.1f}" for v in tc))
