#!/usr/bin/env python3
"""すきま中心線の速度・温度超過の深さ方向減衰率を測り、Stokes 極限 (Papkovich–Fadle) と比べる。

2D 平行壁の間の Stokes 流は、端部の駆動から exp(-Re(λ1) z / a) で減衰する
(λ1 = 4.2124 ± 2.2507i, a = W/2) → **exp(-8.42 z/W)** = 3.66 decade / すきま幅。
実際は慣性 (Re_W = ρ u W/μ) があるので減衰は緩むが、桁のオーダーはここで決まる。
"""
import sys, json, glob
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]


def centerline(run, res=None):
    rd = CASE / run
    if list(rd.glob("res_nan_*.h5")):
        raise SystemExit(f"REFUSED: {run} は発散している")
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    f = rd / res if res else files[-1]
    geom = json.loads((CASE / "geometry.json").read_text())
    setup = json.loads((rd / "case_setup.json").read_text())
    w = geom["cavity"]["widths"][str(setup["series"]["w_over_d"])] * 1e-3
    d = geom["cavity"]["depth"] * 1e-3
    xr = geom["cavity"]["x_rear_wall_from_le"] * 1e-3
    xc = xr - 0.5 * w
    with h5py.File(f, "r") as h:
        c = h["/MESH/COORD"][:].reshape(-1, 3)
        U = h["/VALUE/Ux"][:].astype(float); V = h["/VALUE/Uy"][:].astype(float)
        T = h["/VALUE/T"][:].astype(float)
    m = (np.abs(c[:, 0] - xc) < 0.55 * w / 60) & (c[:, 1] < -1e-6)
    idx = np.where(m)[0]
    idx = idx[np.argsort(-c[idx, 1])]
    z = -c[idx, 1]
    speed = np.hypot(U[idx], V[idx])
    dT = T[idx] - setup["Tw"]
    return z / w, speed, dT, f.name, w


def decay(zw, v, lo=0.3, hi=3.0):
    m = (zw > lo) & (zw < hi) & (v > 0)
    if m.sum() < 5:
        return float("nan")
    p = np.polyfit(zw[m], np.log(v[m]), 1)
    return -p[0]


if __name__ == "__main__":
    print("Stokes 極限: 速度は exp(-8.42 z/W) = 3.66 decade / すきま幅\n")
    print(f"{'run':>34} {'res':>12} {'速度 減衰率':>12} {'dec/W':>7} {'ΔT 減衰率':>11} {'dec/W':>7}")
    for run in sys.argv[1:]:
        try:
            zw, sp, dT, name, w = centerline(run)
        except SystemExit as e:
            print(f"{run:>34} {e}"); continue
        a_u, a_T = decay(zw, sp), decay(zw, np.abs(dT))
        print(f"{run:>34} {name:>12} {a_u:12.3f} {a_u/np.log(10):7.2f} {a_T:11.3f} {a_T/np.log(10):7.2f}")
