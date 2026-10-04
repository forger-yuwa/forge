#!/usr/bin/env python3
"""TP-1187 の深さ分布が薄板内伝導に支配されるかを見積もる (3D に進むかの判断材料)。

case/50 で、W70 の深部 $q$ が**薄板内の横方向伝導**でほぼ説明できることが分かった
(`case/50/tools/skin_smear.py`)。TP-1187 も同じ薄肉過渡法なので、同じ検査を先にかける。
汚染されているなら、3D 計算の照合先そのものが対流熱流束ではないことになる。

TP-1187 の計測側 (Fig 2 / 本文):
    薄板 0.08 cm 304SS (W70 の 0.0305 cm より 2.6 倍厚い)
    平均化 0.25 s、挿入走査 2 秒強、解析時刻 = 模型が中心線に達した時点
    すきま幅 W = 0.10 / 0.18 / 0.30 / 0.41 cm

真の対流分布 $q_{conv}(z)$ は 3D 計算をしないと出せないので、ここでは **case/55 の
forge 乱流すきま ($M$5, $W$ 2.5 mm, $D$ 50 mm) の深さ減衰形をそのまま借り**、
$z/W$ を保って TP-1187 の $W$ に縮尺し、$q_{FP}$ を掛けて熱源にする。
(形だけ借りる粗い見積り。桁の議論にのみ使う。)
"""
import argparse, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
C50 = CASE.parents[0] / "50.deep_cavity_wieting_m7" / "tools"
C55 = CASE.parents[0] / "55.gap_turbulent_m5"
sys.path.insert(0, str(C50))
import importlib.util as _ilu
_sp = _ilu.spec_from_file_location("skin_smear", C50 / "skin_smear.py")
ss = _ilu.module_from_spec(_sp); _sp.loader.exec_module(ss)

TAU_TP = 0.08e-2        # TP-1187 薄板厚 [m]
Q_FP = 63.44e3          # run 8 位置 II の分母 [W/m2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.18e-2, help="すきま幅 [m]")
    ap.add_argument("--depth", type=float, default=6.35e-2, help="タイル厚 = すきま深さ [m]")
    ap.add_argument("--times", default="0.25,1.0,2.0")
    ap.add_argument("--src", default=str(C55 / "run_0010_T1_cav_yp1" / "cavity_rear.csv"))
    ap.add_argument("--src-w", type=float, default=2.5e-3)
    ap.add_argument("--src-d", type=float, default=50e-3)
    ap.add_argument("--src-qfp", type=float, default=83.698e3)
    ap.add_argument("--n", type=int, default=1501)
    a = ap.parse_args()

    arr = np.loadtxt(a.src, delimiter=",", skiprows=1)
    xd, q = arr[:, 0], arr[:, 1]
    o = np.argsort(xd); xd, q = xd[o], q[o]
    zw_src = xd * a.src_d / a.src_w          # z/W (case/55)
    qn = q / a.src_qfp                        # q/q_fp

    s = np.linspace(0.0, a.depth, a.n)
    zw = s / a.w
    qc = np.interp(zw, zw_src, qn, left=qn[0], right=0.0) * Q_FP

    ss.TAU = TAU_TP                           # 薄板厚を差し替える
    ts = [float(t) for t in a.times.split(",")]
    out = {t: ss.solve(s, qc, t) for t in ts}

    print(f"TP-1187 薄板 {TAU_TP*1e3:.2f} mm (W70 は {0.305:.3f} mm), W = {a.w*1e3:.2f} mm, "
          f"D = {a.depth*1e2:.2f} cm, q_FP = {Q_FP*1e-3:.2f} kW/m²")
    print("横方向拡散長: " + ", ".join(
        f"t={t}s → {np.sqrt(ss.ALPHA_S*t)*1e3:.2f} mm ({np.sqrt(ss.ALPHA_S*t)/a.w:.2f} W)" for t in ts))
    print(f"\n{'z/W':>6} {'z [cm]':>7} {'q_conv':>9} " + "".join(f"{'見かけ t='+str(t):>14}" for t in ts)
          + "   ← 見かけ/真")
    for zt in (1, 1.4, 2.8, 4.2, 8.4, 14.1, 21.2):    # 90° 配列の熱電対深さ (W=0.18cm)
        sz = zt * a.w
        if sz > s[-1]:
            continue
        i = int(np.argmin(np.abs(s - sz)))
        row = f"{zt:6.1f} {sz*1e2:7.2f} {qc[i]/Q_FP:9.4f} "
        row += "".join(f"{out[t][1][i]/Q_FP:14.4f}" for t in ts)
        rat = [out[t][1][i] / qc[i] if qc[i] > 1e-9 else np.inf for t in ts]
        row += "   " + " / ".join(f"{r:.0f}" if np.isfinite(r) and r < 1e4 else "∞" for r in rat)
        print(row)
    t0 = ts[1]
    print(f"\n温度上昇 (t={t0}s): 上端 {out[t0][0][0]:.2f} K, 最深 {out[t0][0][-1]:.3f} K")
    print("値は q/q_FP。『見かけ』= q_conv を熱源に TP-1187 の薄板を解き rho c tau dT/dt を取ったもの。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
