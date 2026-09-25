#!/usr/bin/env python3
"""G0 集計: x 継ぎ目の節点と内部節点の勾配を比べる (plan boundary-node-periodic-gradient-fix §6 G0)。

真値: dUx/dy = 1、dK/dy = 0.1、dOmega/dy = 2、他成分 0。y・z の継ぎ目から 2.5 格子以上離れた節点だけを使う
(y 方向は周期なので y 継ぎ目で線形場が跳ぶ)。

使い方: python3 aggregate.py <res_*.h5>   (勾配は output.level 2 で出る)
"""
import sys
import h5py
import numpy as np

L = 2.0 * np.pi
with h5py.File(sys.argv[1], "r") as f:
    xyz = f["MESH/COORD"][()].reshape(-1, 3).astype(np.float64)
    V = {k: f["VALUE/" + k][()].astype(np.float64) for k in f["VALUE"] if k.startswith("d")}
n1 = round(len(xyz) ** (1 / 3))
h = L / (n1 - 1)
x, y, z = xyz.T
tol = 1e-3 * h
far = lambda c: (c > 2.5 * h) & (c < L - 2.5 * h)
ok = far(y) & far(z)
seam = ok & ((np.abs(x) < tol) | (np.abs(x - L) < tol))
inner = ok & far(x)
truth = {"dUxdy": 1.0, "dKdy": 0.1, "dOmegady": 2.0}
# 場の最大値 (ゼロ成分の閾値 1e-6·|φ|/h 用、bake_linear.py の式)
phimax = {"Ux": 10.0 + L, "K": 1.0 + 0.1 * L, "Omega": 100.0 + 2.0 * L}
print(f"G0/G1 (並進周期、線形場、h={h:.5f}): x 継ぎ目 {seam.sum()} 点 / 内部 {inner.sum()} 点")
verdict = {}
for var, scheme in (("Ux", "LSQ (G0)"), ("K", "GG (G1 参考)"), ("Omega", "GG (G1 参考)")):
    ok_all = True
    for c in "xyz":
        name = f"d{var}d{c}"
        a = V[name]
        t = truth.get(name, 0.0)
        s, i = a[seam], a[inner]
        err = max(np.abs(s - t).max(), np.abs(i - t).max())
        # plan §6 G0: 非退化方向は相対 1e-5、ゼロ成分は絶対 1e-6·|φ|/h
        tol = 1e-5 * abs(t) if t != 0.0 else 1e-6 * phimax[var] / h
        ok = err <= tol
        ok_all &= ok
        print(f"  {name:9s}: 継ぎ目 平均 {s.mean():.6f} [{s.min():.6f}, {s.max():.6f}]  内部 平均 {i.mean():.6f} "
              f"[{i.min():.6f}, {i.max():.6f}]  真値 {t}  誤差 {err:.2e} / 閾値 {tol:.2e} {'ok' if ok else 'NG'}")
    verdict[var] = ok_all
    print(f"  -> {var} {scheme}: {'PASS' if ok_all else 'FAIL'}")
print("VERDICT (G0, LSQ):", "PASS" if verdict["Ux"] else "FAIL")
print("G1 参考 (GG、正式判定は CPU double 参照で行う):",
      "PASS" if verdict["K"] and verdict["Omega"] else "FAIL")
