#!/usr/bin/env python3
"""第一内部ノード基準の局所 y1+ (plan §4.5 / G10)。

ソルバの `ypls` は node で壁ノードが壁面に乗るため退化する ([tooling-convergence-and-wall-resolution-gates])。
`check_wall_resolution.py` は viscMethod 2 (kinetic theory) を解決できず「判定不能」を返すので、
場から直接計算する: y1+ = y1 sqrt(ρ_w τ_w) / μ_w,  τ_w = μ_w |du/dn|_w (粘性底層)。
"""
import argparse, sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]


def y1p_line(c, U, V, ro, mu, fixed_axis, fixed_val, span_axis, lo, hi, inward, tol=1e-9):
    axis_vals = np.unique(np.round(c[:, fixed_axis], 9))
    cand = (axis_vals[axis_vals > fixed_val + tol][:1] if inward > 0
            else axis_vals[axis_vals < fixed_val - tol][-1:])
    sel = (np.abs(c[:, fixed_axis] - fixed_val) < tol) & \
          (c[:, span_axis] > lo - tol) & (c[:, span_axis] < hi + tol)
    idx = np.where(sel)[0]
    out = []
    for i in idx:
        sv = c[i, span_axis]
        m = np.where((np.abs(c[:, fixed_axis] - cand[0]) < tol) & (np.abs(c[:, span_axis] - sv) < tol))[0]
        if not len(m):
            continue
        j = int(m[0])
        y1 = abs(c[j, fixed_axis] - fixed_val)
        # 壁接線方向速度 (壁面は no-slip なので第一ノードの接線成分で勾配を作る)
        ut = abs(V[j]) if fixed_axis == 0 else abs(U[j])
        tau = mu[i] * ut / y1
        out.append((c[i, span_axis], y1 * np.sqrt(max(ro[i], 0) * tau) / max(mu[i], 1e-30)))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--res", default=None)
    a = ap.parse_args()
    rd = CASE / a.run
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else files[-1]
    with h5py.File(res, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        U = f["/VALUE/Ux"][:].astype(float); V = f["/VALUE/Uy"][:].astype(float)
        ro = f["/VALUE/ro"][:].astype(float); mu = f["/VALUE/vis_lam"][:].astype(float)
    import json
    g = json.loads((CASE / "geometry.json").read_text())
    w = g["cavity"]["widths"]["0.05"] * 1e-3; d = g["cavity"]["depth"] * 1e-3
    xr = g["cavity"]["x_rear_wall_from_le"] * 1e-3; xf = xr - w
    print(f"run {a.run}  res {res.name}")
    walls = [("平板 (y=0, x>0)", lambda: y1p_line(c, U, V, ro, mu, 1, 0.0, 0, 1e-9, 0.45, +1))]
    if (c[:, 1] < -1e-9).any():
        walls += [("すきま 後壁 (x=xr)", lambda: y1p_line(c, U, V, ro, mu, 0, xr, 1, -d, -1e-9, -1)),
                  ("すきま 前壁 (x=xf)", lambda: y1p_line(c, U, V, ro, mu, 0, xf, 1, -d, -1e-9, +1)),
                  ("すきま 床 (y=-d)",   lambda: y1p_line(c, U, V, ro, mu, 1, -d, 0, xf, xr, +1))]
    ok = True
    for name, fn in walls:
        A = fn()
        if not len(A):
            print(f"  {name:22} 点なし"); ok = False; continue
        yp = A[:, 1]
        over = 100.0 * np.count_nonzero(yp > 1.0) / len(yp)
        ok &= over == 0.0
        print(f"  {name:22} {len(yp):5d} 点  y1+ 平均 {yp.mean():7.3f} / p99 {np.percentile(yp,99):7.3f} "
              f"/ 最大 {yp.max():7.3f}   >1 の点 {over:5.1f} %")
    print(f"\nVERDICT: {'PASS' if ok else 'FAIL'} (局所 y1+ <= 1 を全点で要求)")


if __name__ == "__main__":
    main()
