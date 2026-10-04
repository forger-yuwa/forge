#!/usr/bin/env python3
"""(a) ドリフト方向と run 間オフセットのコサイン。閉包試験は v7_closure.py。"""
import sys, h5py, numpy as np
sys.path.insert(0,".")
from v7_dist import load, wall_mask
MESH="run_0447_cfl02_ext/sern.h5"
S02="run_0439_3d_junction_outflow_steady/res_36000.h5"; E02="run_0447_cfl02_ext/res_10000.h5"
S04="run_0442_3d_junction_cfl04_long/res_27000.h5";     E04="run_0448_cfl04_ext/res_10000.h5"
a,b,e,f = load(S02), load(S04), load(E02), load(E04)
wm = wall_mask(MESH, len(a["ro"]))

print("=== (a) cos(δ, D)  δ = 両 run のドリフトの差, D = 種同士のオフセット ===")
for k in ("ro","roe","P"):
    D  = (b[k]-a[k])[wm]                       # 種同士のオフセット
    d  = ((f[k]-b[k]) - (e[k]-a[k]))[wm]       # ドリフトの差
    c  = float(D@d/(np.linalg.norm(D)*np.linalg.norm(d)))
    print(f"  {k:4s} cos {c:+.4f}   |D| {np.linalg.norm(D):.4e}  |δ| {np.linalg.norm(d):.4e}  |δ|/|D| {np.linalg.norm(d)/np.linalg.norm(D):.4f}")
print("  判定: |cos| > 0.7 ならオフセットはドリフトモードそのもの -> 未決")

# 閉包試験は v7_closure.py に分離 (ここに書いた版はプレースホルダ行が残って落ちた)
