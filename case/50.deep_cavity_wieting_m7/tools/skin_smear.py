#!/usr/bin/env python3
"""薄肉過渡法の横方向伝導を **順方向に解いて**、CFD と実測を同じ土俵に乗せる。

W70 (TN D-5908) は後壁・床に 0.305 mm の 304SS 薄板を貼り、挿入後の温度上昇率から
    q_meas(s) = rho c tau (dT/dt)_s
で局所熱流束を出す。本文は「表面伝導・放射の補正をしていない」と明記している。
薄板 (背面断熱) の実際の式は

    dT/dt = q_conv(s)/(rho c tau) + alpha_s d2T/ds2,      alpha_s = lam_s/(rho c)

で、alpha_s の横方向拡散長は 304SS・0.8 s で **1.74 mm = すきま幅の 1.4 倍**。
これは q_conv(s) の減衰長 (0.5 幅程度) より長いので、**局所補正 (q + alpha t q'') では
足りず** (深部では q''≈0 なので何も起きない)、過渡方程式をそのまま解く必要がある。

ここでは CFD の q_conv(s) を熱源として T(s,t) を解き、実測と同じ演算
(rho c tau dT/dt) を適用して **見かけ熱流束**にしてから比較する。逆問題 (後退拡散) は
解かない。両端は断熱とする — 外から熱を入れない**保守側**の仮定で、
横方向伝導は上端の熱を深部へ再配分するだけになる。
"""
import argparse, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent

LAM_S = 15.0            # 304SS 熱伝導率 [W/(m K)]
RHO_C = 7900 * 500.0    # rho c [J/(m^3 K)]
TAU = 0.305e-3          # 薄板厚 [m] (0.012 in)
ALPHA_S = LAM_S / RHO_C


def read_ref(path):
    rows = []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        try:
            rows.append([float(v) for v in ln.split(",")[:2]])
        except ValueError:
            continue
    return np.array(rows)


def solve(s, qc, t_end, n_sub=None):
    """dT/dt = qc/(rho c tau) + alpha_s T''  を陽解法で t_end まで。両端断熱。"""
    ds = s[1] - s[0]
    dt_max = 0.4 * ds ** 2 / ALPHA_S
    nst = int(np.ceil(t_end / dt_max)) if n_sub is None else n_sub
    dt = t_end / nst
    T = np.zeros_like(qc)
    src = qc / (RHO_C * TAU)
    for _ in range(nst):
        lap = np.zeros_like(T)
        lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
        lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2          # 断熱端
        lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2
        T += dt * (src + ALPHA_S * lap)
    # 実測と同じ演算: q_meas = rho c tau dT/dt (終端時刻の瞬時値)
    lap = np.zeros_like(T)
    lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
    lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2
    lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2
    return T, RHO_C * TAU * (src + ALPHA_S * lap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--ref", default="ref/w70_fig6a_wd0063_rear.csv")
    ap.add_argument("--qfp", type=float, default=None, help="既定は case_setup.json の qfp_kW")
    ap.add_argument("--depth", type=float, default=None, help="既定は geometry.json")
    ap.add_argument("--width", type=float, default=None, help="既定は geometry.json + 系列の w/d")
    ap.add_argument("--times", default="0.4,0.8,1.6")
    ap.add_argument("--n", type=int, default=1201)
    a = ap.parse_args()

    import json
    setup = json.loads((CASE / a.run / "case_setup.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    if a.qfp is None:
        a.qfp = setup["series"]["qfp_kW"] * 1e3
    if a.depth is None:
        a.depth = geom["cavity"]["depth"] * 1e-3
    if a.width is None:
        a.width = geom["cavity"]["widths"][str(setup["series"]["w_over_d"])] * 1e-3

    arr = np.loadtxt(CASE / a.run / "cavity_rear.csv", delimiter=",", skiprows=1)
    xd, q = arr[:, 0], arr[:, 1]
    o = np.argsort(xd); xd, q = xd[o], q[o]
    s = np.linspace(0.0, xd.max() * a.depth, a.n)
    qc = np.interp(s, xd * a.depth, q)

    ref = read_ref(CASE / a.ref)
    ts = [float(t) for t in a.times.split(",")]
    out = {t: solve(s, qc, t) for t in ts}

    print(f"run {a.run}  w/d={setup['series']['w_over_d']}  W={a.width*1e3:.3f} mm  "
          f"q_fp={a.qfp*1e-3:.1f} kW/m²   alpha_s={ALPHA_S:.3e} m²/s")
    print("横方向拡散長: " + ", ".join(
        f"t={t}s → {np.sqrt(ALPHA_S*t)*1e3:.2f} mm ({np.sqrt(ALPHA_S*t)/a.width:.1f} W)" for t in ts))
    print(f"\n{'z/W':>6} {'x/d':>6} {'実測':>9} {'CFD 生':>9} " +
          "".join(f"{'見かけ t='+str(t):>14}" for t in ts))
    for zt in (0.5, 1, 2, 3, 4, 6, 8, 12):
        sz = zt * a.width
        if sz > s[-1]:
            continue
        i = int(np.argmin(np.abs(s - sz)))
        xdi = sz / a.depth
        rq = np.interp(xdi, ref[:, 0], ref[:, 1]) if ref[:, 0].min() <= xdi <= ref[:, 0].max() else np.nan
        print(f"{zt:6.1f} {xdi:6.3f} {rq:9.4f} {qc[i]/a.qfp:9.4f} " +
              "".join(f"{out[t][1][i]/a.qfp:14.4f}" for t in ts))
    t0 = ts[len(ts)//2]
    print(f"\n温度上昇 (t={t0}s): 上端 {out[t0][0][0]:.2f} K, z/W=2 {out[t0][0][int(np.argmin(np.abs(s-2*a.width)))]:.2f} K, "
          f"最深 {out[t0][0][-1]:.3f} K   (本文: 計測区間の最大温度上昇 < 22 K)")
    print("値はすべて q/q_fp。『見かけ』= CFD の q_conv を熱源に薄板を解き、実測と同じ rho c tau dT/dt を取ったもの。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
