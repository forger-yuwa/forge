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
BAND_FRAC = (0.40, 0.60)  # x 範囲のこの割合の帯は P 一様 (ΔP=0 の対象面を作る)

def det_rand(key, salt):
    """**座標**からの決定的疑似乱数 [-1,1) (splitmix64)。

    節点 index から作ると、周期 partner (別 index・同一 DOF) に別の値が乗り
    「同一 DOF に 2 つの状態」になる (plan §6 V6 の前提条件)。周期方向は
    折り返した座標を量子化してから混ぜるので、partner はビット一致する。"""
    x = (key.astype(np.uint64) + np.uint64(salt)) * np.uint64(0x9E3779B97F4A7C15)
    x ^= x >> np.uint64(30); x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27); x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return (x >> np.uint64(11)).astype(np.float64) / float(1 << 53) * 2.0 - 1.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--period", default="", help="周期長を 'dx,dy,dz' で (0 = 非周期)。指定すると座標を折り返してからシードする")
    ap.add_argument("--axisymmetric", action="store_true", help="軸 (y=0) のノードは壁扱いにしない (slip/axis)")
    a = ap.parse_args()
    per = [float(v) for v in a.period.split(",")] if a.period else [0.0, 0.0, 0.0]
    with h5py.File(a.h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        # 座標キー: 周期方向は折り返し、1 nm 刻みで量子化 (float32 の丸め違いを吸収)
        q = np.empty((n, 3), dtype=np.int64)
        for d in range(3):
            v = c[:, d].astype(np.float64)
            if per[d] > 0.0: v = np.mod(v, per[d])
            q[:, d] = np.rint(v / 1e-9).astype(np.int64)
        key = (q[:, 0] * np.int64(73856093)) ^ (q[:, 1] * np.int64(19349663)) ^ (q[:, 2] * np.int64(83492791))
        r1, r2 = det_rand(key, 1), det_rand(key, 2)
        x0, x1 = float(c[:, 0].min()), float(c[:, 0].max())
        bl, bh = x0 + BAND_FRAC[0] * (x1 - x0), x0 + BAND_FRAC[1] * (x1 - x0)
        band = (c[:, 0] >= bl) & (c[:, 0] <= bh)
        r2 = np.where(band, 0.0, r2)

        u_t = M_T * np.sqrt(GAM * R * T_INF)
        P = P_INF * (1.0 + 0.05 * r2)
        T = np.full(n, T_INF)
        ro = P / (R * T)
        Ux = np.full(n, u_t); Uy = 0.2 * u_t * r1; Uz = np.zeros(n)
        wall = (wd <= 0.0)
        if a.axisymmetric: wall &= (np.abs(c[:, 1]) > 1e-9)   # 軸ノードは壁でない
        Ux[wall] = 0.0; Uy[wall] = 0.0

        roe = P / (GAM - 1) + 0.5 * ro * (Ux**2 + Uy**2 + Uz**2)
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * Ux).astype(np.float32)
        f["/VALUE/roUy"][:] = (ro * Uy).astype(np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name in ("roK", "roOmega"):
            if f"/VALUE/{name}" in f: del f[f"/VALUE/{name}"]
        print(f"V5 state: n={n} wall={int(wall.sum())} band(ΔP=0)={int(band.sum())} [x {bl:.4g}..{bh:.4g}] "
              f"U_t={u_t:.1f} m/s  P {P.min():.1f}..{P.max():.1f} Pa  ro {ro.min():.5f}..{ro.max():.5f}")
if __name__ == "__main__":
    main()
