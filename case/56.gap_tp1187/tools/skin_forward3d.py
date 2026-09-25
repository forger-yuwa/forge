#!/usr/bin/env python3
"""3D の壁熱流束を**薄肉過渡法と同じ土俵**に乗せてから実測と比べる。

TP-1187 は薄板 0.08 cm・304SS の温度上昇率から q = ρcτ dT/dt を出しており、薄板内の
横方向伝導を補正していない。**CFD の生の q_conv を実測と直接比べてはいけない**
(`acceptance.json` の T4-3D observation_model)。ここでは case/50 と同じ順方向解法で

    dT/dt = q_conv(ℓ)/(ρcτ) + α_s ∂²T/∂ℓ²,   α_s = λ_s/(ρc)

を解き、実測と同じ演算 ρcτ dT/dt をかけて**見かけ熱流束**にする。逆問題は解かない。

**経路 ℓ は濡れ面をひとつなぎに取る** (2026-09-26 に原報 Fig 5 で確定):

    ℓ=0: すきま底 (y=-D) → 前向き壁を上る → エッジ円弧 → 上面を下流へ

熱電対は Fig 5 の寸法どおり **上面 y=0 からの鉛直深さ**で置く:
  前向き壁 0.25/0.51/0.76/1.52/2.54/3.81 cm = TC 92/91/90/89/88/87
  上面側は角から下流へ 0.25/0.51 cm = TC 93/94 (r=0.25 cm なので 93 は円弧の上端接点)

端条件は両極で挟む: `adiabatic` (両端断熱 = 深部を高めに) と `sink` (両端 T=0 = 低めに)。
"""
import argparse, glob, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
LAM_S, RHO_C = 15.0, 7900 * 500.0        # 304SS
ALPHA_S = LAM_S / RHO_C
TAU_TP1187 = 0.08e-2                     # 薄板厚 [m]
WALL_TC = [(92, 0.25), (91, 0.51), (90, 0.76), (89, 1.52), (88, 2.54), (87, 3.81)]
TOP_TC = [(93, 0.25), (94, 0.51)]        # 角から下流へ [cm]


def path_profile(run, step, W, r, D, zband, tol=1e-6):
    """濡れ面に沿った (ℓ, q_conv) を返す。ℓ=0 がすきま底、上面下流が正の端。"""
    xd = 0.5 * W
    xvd = xd + r
    with h5py.File(CASE / run / f"res_gap_6_{step}.h5") as h:
        cg = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
        qg = np.asarray(h["VALUE/qwall"], dtype=float)
    pf = CASE / run / f"res_plate_4_{step}.h5"
    if pf.exists():
        with h5py.File(pf) as h:
            cp = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
            qp = np.asarray(h["VALUE/qwall"], dtype=float)
    else:
        cp, qp = np.zeros((0, 3)), np.zeros(0)

    def band(c, q, m):
        m = m & (c[:, 2] <= zband)
        return c[m], q[m]

    segs = []
    # (1) 前向き壁: x=xd, -D<=y<=-r。ℓ = y + D
    c1, q1 = band(cg, qg, (np.abs(cg[:, 0] - xd) < tol) & (cg[:, 1] < -r + tol))
    if len(c1):
        segs.append((c1[:, 1] + D, c1, q1, "wall"))
    # (2) 円弧 arc_d: 中心 (xvd,-r) 半径 r。ℓ = (D-r) + r*θ、θ は壁側 0 → 上面側 π/2
    onarc = (np.abs(np.hypot(cg[:, 0] - xvd, cg[:, 1] + r) - r) < 5e-5) & (cg[:, 1] > -r - tol)
    c2, q2 = band(cg, qg, onarc)
    if len(c2):
        th = np.arctan2(c2[:, 1] + r, -(c2[:, 0] - xvd))     # 壁側で 0、上面側で π/2
        segs.append(((D - r) + r * np.clip(th, 0.0, np.pi / 2), c2, q2, "arc"))
    # (3) 上面: y=0, x>=xvd。ℓ = (D-r) + rπ/2 + (x - xvd)
    if len(cp):
        c3, q3 = band(cp, qp, (np.abs(cp[:, 1]) < tol) & (cp[:, 0] > xvd - tol))
        if len(c3):
            segs.append(((D - r) + r * np.pi / 2 + (c3[:, 0] - xvd), c3, q3, "top"))

    ls, qs = [], []
    for l, c, q, _ in segs:
        # z 帯の平均を ℓ ごとに取る
        for lv in np.unique(np.round(l, 9)):
            m = np.abs(l - lv) < 1e-9
            ls.append(lv); qs.append(float(np.nanmean(q[m])))
    o = np.argsort(ls)
    return np.array(ls)[o], np.array(qs)[o], {s[3]: int(len(np.unique(np.round(s[0], 9)))) for s in segs}


def solve(l, qc, t_end, tau, end="adiabatic"):
    dl = l[1] - l[0]
    nst = max(1, int(np.ceil(t_end / (0.4 * dl ** 2 / ALPHA_S))))
    dt = t_end / nst
    T = np.zeros_like(qc)
    src = qc / (RHO_C * tau)

    def lap_of(T):
        lp = np.zeros_like(T)
        lp[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / dl ** 2
        if end == "adiabatic":
            lp[0] = 2.0 * (T[1] - T[0]) / dl ** 2
            lp[-1] = 2.0 * (T[-2] - T[-1]) / dl ** 2
        else:
            lp[0] = (T[1] - 2.0 * T[0]) / dl ** 2
            lp[-1] = (T[-2] - 2.0 * T[-1]) / dl ** 2
        return lp

    for _ in range(nst):
        T += dt * (src + ALPHA_S * lap_of(T))
        if end != "adiabatic":
            T[0] = T[-1] = 0.0
    return RHO_C * tau * (src + ALPHA_S * lap_of(T))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--zband", type=float, default=6.0e-3)
    ap.add_argument("--tau", type=float, default=TAU_TP1187)
    ap.add_argument("--t", type=float, default=0.25, help="解析時刻 [s] (平均化窓 0.25 s)")
    ap.add_argument("--n", type=int, default=2001)
    ap.add_argument("--qfp-forge", type=float, default=77.34e3)
    ap.add_argument("--qfp-meas", type=float, default=63.44e3)
    a = ap.parse_args()

    if a.step is None:
        fs = [int(re.search(r"res_gap_6_(\d+)\.h5$", p).group(1))
              for p in glob.glob(str(CASE / a.run / "res_gap_6_*.h5"))]
        if not fs:
            raise SystemExit(f"{a.run}: res_gap_6_*.h5 が無い")
        a.step = max(fs)

    l, q, segn = path_profile(a.run, a.step, a.w, a.r, a.depth, a.zband)
    ll = np.linspace(l.min(), l.max(), a.n)
    qq = np.interp(ll, l, q)
    print(f"[{a.run}] step {a.step}  経路 ℓ {l.min()*1e3:.2f}–{l.max()*1e3:.1f} mm  区間点数 {segn}")
    print(f"  薄板 τ={a.tau*1e3:.2f} mm  √(α_s t)={np.sqrt(ALPHA_S*a.t)*1e3:.2f} mm @ t={a.t:.2f} s")
    res = {e: solve(ll, qq, a.t, a.tau, e) for e in ("adiabatic", "sink")}

    print("\n  TC  位置              生 q_conv    観測後(断熱) 観測後(吸熱)   q/q_FP forge   q/q_FP 実測")
    rows = ([(tc, (a.depth - d * 1e-2), f"壁 深さ{d:.2f}cm") for tc, d in WALL_TC]
            + [(tc, (a.depth - a.r) + a.r * np.pi / 2 + max(dx * 1e-2 - a.r, 0.0),
                f"上面 角から{dx:.2f}cm") for tc, dx in TOP_TC])
    for tc, lv, lbl in sorted(rows, key=lambda x: -x[1]):
        if not (ll.min() - 1e-9 <= lv <= ll.max() + 1e-9):
            print(f"  {tc:3d}  {lbl:16s} (経路外)"); continue
        raw = np.interp(lv, l, q)
        va, vs = (np.interp(lv, ll, res[e]) for e in ("adiabatic", "sink"))
        print(f"  {tc:3d}  {lbl:16s} {raw:11.1f} {va:12.1f} {vs:12.1f}"
              f"   {abs(va)/a.qfp_forge:5.3f}-{abs(vs)/a.qfp_forge:5.3f}"
              f"   {abs(va)/a.qfp_meas:5.3f}-{abs(vs)/a.qfp_meas:5.3f}")


if __name__ == "__main__":
    main()
