#!/usr/bin/env python3
"""G0 (並進周期) の入力: 線形場を node の VALUE に焼き込む。

plan boundary-node-periodic-gradient-fix §6 G0。Ux = 10 + y、k = 1 + 0.1 y、omega = 100 + 2 y (密度・圧力一様)。
y は周期方向なので y 継ぎ目では線形場が跳ぶ。評価は y・z の継ぎ目から離れた節点の x 継ぎ目だけで行う (aggregate.py)。

使い方: python3 bake_linear.py <変換済み Taylor-Green.h5>   (VALUE を上書き)
"""
import sys
import h5py
import numpy as np

GAMMA, CP, P0 = 1.4, 1038.8, 101325.0

with h5py.File(sys.argv[1], "r+") as f:
    xyz = f["MESH/COORD"][()].reshape(-1, 3).astype(np.float64)
    y = xyz[:, 1]
    v = f["VALUE"]
    n = v["ro"].shape[0]
    assert n == xyz.shape[0], (n, xyz.shape)
    ro = float(v["ro"][0])
    ux = 10.0 + y
    k = 1.0 + 0.1 * y
    om = 100.0 + 2.0 * y
    e = P0 / (GAMMA - 1.0) + 0.5 * ro * ux**2     # roe は k を含めない (sstEnergyIncludesK 既定 0)
    for name, a in (("roUx", ro * ux), ("roUy", 0 * y), ("roUz", 0 * y),
                    ("roK", ro * k), ("roOmega", ro * om), ("roe", e)):
        v[name][...] = a.astype(np.float32)
print("baked linear field:", sys.argv[1])
