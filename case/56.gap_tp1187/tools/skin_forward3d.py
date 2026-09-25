#!/usr/bin/env python3
"""3D の前向き壁 q_conv を**薄肉過渡法と同じ土俵**に乗せてから実測と比べる。

TP-1187 は薄板 0.08 cm・304SS の温度上昇率から q = ρcτ dT/dt を出している。
薄板内の横方向伝導は無補正なので、**CFD の生の q_conv を実測と直接比べてはいけない**
(`acceptance.json` の T4-3D observation_model)。ここでは case/50 と同じ順方向解法で

    dT/dt = q_conv(s)/(ρcτ) + α_s ∂²T/∂s²,   α_s = λ_s/(ρc)

を解き、実測と同じ演算 ρcτ dT/dt をかけて**見かけ熱流束**にする。逆問題は解かない。

**帯はリップの熱的接続の両極で作る** (事前登録):
  adiabatic  : 両端断熱。外へ熱を逃がさない = 深部を**高めに**出す側
  sink       : 浅い端 (リップ) を T=0 固定。母材タイルへ逃がす = 深部を**低めに**出す側

    python3 tools/skin_forward3d.py --run run_0052_ab3d_open_long [--step 60000] [--t 0.25]
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
LAM_S = 15.0                 # 304SS 熱伝導率 [W/(m K)]
RHO_C = 7900 * 500.0         # ρc [J/(m³ K)]
ALPHA_S = LAM_S / RHO_C
TAU_TP1187 = 0.08e-2         # 薄板厚 [m] (conditions.json: skin_thickness_cm 0.08)
DEPTH_CM = [0.25, 0.51, 0.76, 1.52, 2.54, 3.81]
TC = [92, 91, 90, 89, 88, 87]


def wall_profile(run, step, W, r, D, zband, tol=1e-6):
    """前向き壁 (x=+W/2) の q_conv を**リップからの深さ**の関数として返す (z 帯平均)。"""
    f = CASE / run / f"res_gap_6_{step}.h5"
    with h5py.File(f) as h:
        c = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
        q = np.asarray(h["VALUE/qwall"], dtype=float)
    sel = ((np.abs(c[:, 0] - 0.5 * W) < tol) & (c[:, 1] < -r + tol)
           & (c[:, 1] > -D - tol) & (c[:, 2] <= zband))
    c, q = c[sel], q[sel]
    zs = np.unique(np.round(c[:, 2], 9))
    ys = np.unique(np.round(c[:, 1], 9))[::-1]          # 浅い -> 深い
    prof = np.full((len(zs), len(ys)), np.nan)
    for i, z0 in enumerate(zs):
        m = np.abs(c[:, 2] - z0) < 1e-9
        yy, qq = c[m, 1], q[m]
        o = np.argsort(yy)
        prof[i] = np.interp(ys, yy[o], qq[o])
    return (-ys - r), np.nanmean(prof, axis=0)          # リップからの深さ [m], q [W/m²]


def solve(s, qc, t_end, tau, end="adiabatic"):
    """薄板過渡を陽解法で解き、実測と同じ ρcτ dT/dt を返す。"""
    ds = s[1] - s[0]
    nst = max(1, int(np.ceil(t_end / (0.4 * ds ** 2 / ALPHA_S))))
    dt = t_end / nst
    T = np.zeros_like(qc)
    src = qc / (RHO_C * tau)

    def lap_of(T):
        lap = np.zeros_like(T)
        lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
        if end == "adiabatic":
            lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2
        else:                                            # sink: リップ端 T=0
            lap[0] = (T[1] - 2.0 * T[0]) / ds ** 2
        lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2        # 深い端は常に断熱 (袋小路)
        return lap

    for _ in range(nst):
        T += dt * (src + ALPHA_S * lap_of(T))
        if end != "adiabatic":
            T[0] = 0.0
    return RHO_C * tau * (src + ALPHA_S * lap_of(T))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--step", type=int, default=None, help="既定は最新")
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--zband", type=float, default=6.0e-3)
    ap.add_argument("--tau", type=float, default=TAU_TP1187, help="薄板厚 [m]")
    ap.add_argument("--t", type=float, default=0.25, help="解析時刻 [s] (平均化窓)")
    ap.add_argument("--times", default=None, help="時刻感度を見る (カンマ区切り)")
    ap.add_argument("--n", type=int, default=1601, help="薄板方向の格子点")
    ap.add_argument("--qfp-forge", type=float, default=77.34e3,
                    help="forge の 2D 平板 q (位置 II)。既定は run_0002_fp_t8_long")
    ap.add_argument("--qfp-meas", type=float, default=63.44e3,
                    help="実測 q_FP (Table III run 8 位置 II)")
    a = ap.parse_args()

    if a.step is None:
        fs = [int(re.search(r"res_gap_6_(\d+)\.h5$", p).group(1))
              for p in glob.glob(str(CASE / a.run / "res_gap_6_*.h5"))]
        if not fs:
            raise SystemExit(f"{a.run}: res_gap_6_*.h5 が無い")
        a.step = max(fs)

    d, q = wall_profile(a.run, a.step, a.w, a.r, a.depth, a.zband)
    s = np.linspace(d.min(), d.max(), a.n)
    qc = np.interp(s, d, q)
    print(f"[{a.run}] step {a.step}  前向き壁 {len(d)} 点  深さ {d.min()*1e3:.2f}–{d.max()*1e3:.1f} mm")
    print(f"  薄板 τ={a.tau*1e3:.2f} mm  α_s={ALPHA_S:.3e} m²/s  "
          f"横方向拡散長 √(α_s t)={np.sqrt(ALPHA_S*a.t)*1e3:.2f} mm @ t={a.t:.2f} s")

    res = {e: solve(s, qc, a.t, a.tau, e) for e in ("adiabatic", "sink")}
    print(f"\n  TC  深さ[cm]   生 q_conv     観測モデル後 (見かけ q)        q/q_FP (forge 分母)"
          f"      q/q_FP (実測分母)")
    print(f"                  [W/m²]      断熱端      リップ吸熱        断熱端   吸熱端"
          f"        断熱端   吸熱端")
    for tc, dc in zip(TC, DEPTH_CM):
        # 深さは上面 y=0 基準なので、リップ (y=-r) からの深さに直す
        sd = dc * 1e-2 - a.r
        if sd < s.min() - 1e-9:
            print(f"  {tc:3d} {dc:7.2f}   (リップより浅い — 半径上)"); continue
        raw = np.interp(sd, d, q)
        va = np.interp(sd, s, res["adiabatic"])
        vs = np.interp(sd, s, res["sink"])
        print(f"  {tc:3d} {dc:7.2f} {raw:11.1f} {va:11.1f} {vs:13.1f}"
              f"   {abs(va)/a.qfp_forge:9.3f} {abs(vs)/a.qfp_forge:8.3f}"
              f"   {abs(va)/a.qfp_meas:9.3f} {abs(vs)/a.qfp_meas:8.3f}")

    if a.times:
        print("\n  時刻感度 (断熱端、q/q_FP forge 分母)")
        for t in [float(x) for x in a.times.split(",")]:
            r2 = solve(s, qc, t, a.tau, "adiabatic")
            vals = [abs(np.interp(dc * 1e-2 - a.r, s, r2)) / a.qfp_forge
                    for dc in DEPTH_CM if dc * 1e-2 - a.r >= s.min() - 1e-9]
            print(f"    t={t:.2f} s: " + " ".join(f"{v:.3f}" for v in vals))


if __name__ == "__main__":
    main()
