#!/usr/bin/env python3
"""V5 (#10d-2): ΔP を設計で与えた人工状態を h5 に書く。

plan convection-slau-wall-normal-chi §6 V5。収束場は使わない — 収束した ZPG 平板は
壁隣接面の ΔP が最小のケースで、χ_n の信号が出ない (V3 で Cf 差 0.000 %)。

  Ux = U_t (M_t≈0.6, chi が 0/1 にクリップされない)   Uy = 0.2 U_t r1   Uz = 0
  P  = P_inf (1 + 0.05 r2)    T = T_inf 一様    ro = P/(R T)    壁ノードは u=0
  r1, r2 ∈ [-1,1] は**節点 index からの決定的疑似乱数** (run 間で同一)。
  x が BAND_X の帯では r2 = 0 → 「対象面だが ΔP=0」の集合を意図的に作る (P2 用)。
"""
import argparse
from pathlib import Path
import h5py, numpy as np

GAM, CP = 1.4, 1004.5
R = CP * (GAM - 1) / GAM
T_INF, P_INF = 283.0, 5037.4
M_T = 0.6
BAND_X = (0.40, 0.60)     # この x 帯は P 一様 (ΔP=0 の対象面を作る)

def det_rand(idx, salt):
    """節点 index からの決定的疑似乱数 [-1,1) (splitmix64)。"""
    x = (idx.astype(np.uint64) + np.uint64(salt)) * np.uint64(0x9E3779B97F4A7C15)
    x ^= x >> np.uint64(30); x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27); x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return (x >> np.uint64(11)).astype(np.float64) / float(1 << 53) * 2.0 - 1.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    a = ap.parse_args()
    with h5py.File(a.h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        idx = np.arange(n, dtype=np.int64)
        r1, r2 = det_rand(idx, 1), det_rand(idx, 2)
        band = (c[:, 0] >= BAND_X[0]) & (c[:, 0] <= BAND_X[1])
        r2 = np.where(band, 0.0, r2)

        u_t = M_T * np.sqrt(GAM * R * T_INF)
        P = P_INF * (1.0 + 0.05 * r2)
        T = np.full(n, T_INF)
        ro = P / (R * T)
        Ux = np.full(n, u_t); Uy = 0.2 * u_t * r1; Uz = np.zeros(n)
        wall = (wd <= 0.0) & (c[:, 0] >= -1e-9)
        Ux[wall] = 0.0; Uy[wall] = 0.0

        roe = P / (GAM - 1) + 0.5 * ro * (Ux**2 + Uy**2 + Uz**2)
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * Ux).astype(np.float32)
        f["/VALUE/roUy"][:] = (ro * Uy).astype(np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name in ("roK", "roOmega"):
            if f"/VALUE/{name}" in f: del f[f"/VALUE/{name}"]
        print(f"V5 state: n={n} wall={int(wall.sum())} band(ΔP=0)={int(band.sum())} "
              f"U_t={u_t:.1f} m/s  P {P.min():.1f}..{P.max():.1f} Pa  ro {ro.min():.5f}..{ro.max():.5f}")
if __name__ == "__main__":
    main()
