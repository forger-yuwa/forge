#!/usr/bin/env python3
"""W70 の観測モデル — **薄板全体を 1 本の経路として解き、全点・全幅に同一適用する**。

`skin_smear.py` (2026-09-20 初版) は後壁だけを両端断熱で解いていた。codex result M2 の指摘:

  * 計測時刻が固定されておらず、深さ 4 幅で 0.4-1.6 s の振りが **9.6 倍**の幅を生む
  * 両端断熱・後壁単独・一定熱源は挿入中の熱履歴や床との接続を再現していない
  * $w/d$=0.211 の $z/W$=0.5 だけ補正前の CFD を支持材料に使っていた (不整合)

ここでは観測モデルを 1 つに固定し、**全比較点・全幅に同じものを掛ける**:

  経路   : 後壁 (リップ s=0 → 床隅) + 床 (隅 → 前壁側)。W70 は**後壁と床**に薄板を貼る。
  方程式 : dT/dt = q_conv(s)/(rho c tau) + alpha_s d2T/ds2
  リップ端: 模型表面に繋がる。表面は q_fp を受けるので、同じ時間だけ温まった温度を
           Dirichlet で与える (断熱より現実的で、伝導寄与を**過小評価しない**)
  遠端    : 前壁側は断熱 (前壁は薄板でない)
  計測    : q_meas = rho c tau dT/dt を **t_eval で評価**。t_eval は原典の
           「挿入 0.8 s・中心線到達時に評価」に固定し、感度は 0.4/1.6 s で別途示す。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
LAM_S, RHO_C, TAU = 15.0, 7900 * 500.0, 0.305e-3
# [W70] の報告された精度限界 (標準偏差ではなく区間)。plan §4.10。
Q_PRECISION_W_M2 = 0.68e3
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


def build_path(run, depth, width, n):
    """後壁 (s: 0→D) + 床 (s: D→D+W) の 1 本の経路に q_conv を並べる。"""
    rd = CASE / run
    rear = np.loadtxt(rd / "cavity_rear.csv", delimiter=",", skiprows=1)
    flo = np.loadtxt(rd / "cavity_floor.csv", delimiter=",", skiprows=1)
    o = np.argsort(rear[:, 0]); rear = rear[o]
    o = np.argsort(flo[:, 0]); flo = flo[o]
    s = np.linspace(0.0, depth + width, n)
    q = np.where(s <= depth,
                 np.interp(np.clip(s, 0, depth) / depth, rear[:, 0], rear[:, 1]),
                 # 床 csv の座標は **y/d** (幅 W ではない)。W で割ると床の大半が端値の外挿になる
                 # (2026-09-20 codex result-2 M8)。
                 np.interp(np.clip(s - depth, 0, width) / depth, flo[:, 0], flo[:, 1]))
    return s, q


def solve(s, qc, t_end, T_lip, lip="dirichlet"):
    """リップ端の扱いを選べる。lip="adiabatic" なら断熱 (伝導寄与の下限側)。"""
    ds = s[1] - s[0]
    nst = int(np.ceil(t_end / (0.4 * ds ** 2 / ALPHA_S)))
    dt = t_end / nst
    T = np.zeros_like(qc)
    src = qc / (RHO_C * TAU)
    for k in range(nst):
        lap = np.zeros_like(T)
        lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
        lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2          # 遠端断熱
        if lip == "adiabatic":
            lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2
        T += dt * (src + ALPHA_S * lap)
        if lip == "dirichlet":
            T[0] = T_lip * (k + 1) * dt / t_end
    # **時間更新と同じ境界演算子**を使う。lap[0]=lap[1] で代用していたため
    # 全経路積分が 1.017 倍ずれ、それを「保存」と誤って報告していた
    # (2026-09-20 codex result-2 M8)。同じ演算子なら厳密に保存する。
    lap = np.zeros_like(T)
    lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
    lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2
    if lip == "adiabatic":
        lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2
    else:
        lap[0] = lap[1]          # Dirichlet 端では壁側の値を使わない
    return T, RHO_C * TAU * (src + ALPHA_S * lap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-eval", type=float, default=0.8, help="評価時刻 [s] (原典: 挿入 0.8 s)")
    ap.add_argument("--t-sens", default="0.4,1.6", help="感度用の時刻")
    ap.add_argument("--n", type=int, default=1201)
    a = ap.parse_args()
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    d = geom["cavity"]["depth"] * 1e-3
    cases = [("run_0006_T1_wd0063_long", 0.063, "ref/w70_fig6a_wd0063_rear.csv"),
             ("run_0008_T1_wd0211_long", 0.211, "ref/w70_fig6b_wd0211_rear.csv")]
    ts = [a.t_eval] + [float(x) for x in a.t_sens.split(",")]
    for run, wd, ref in cases:
        setup = json.loads((CASE / run / "case_setup.json").read_text(encoding="utf-8"))
        W = geom["cavity"]["widths"][str(wd)] * 1e-3
        qfp = setup["series"]["qfp_kW"] * 1e3
        s, qc = build_path(run, d, W, a.n)
        # リップの温度上昇: 表面は q_fp を t だけ受ける
        out = {t: solve(s, qc, t, qfp * t / (RHO_C * TAU)) for t in ts}
        adb = solve(s, qc, a.t_eval, 0.0, lip="adiabatic")
        R = read_ref(CASE / ref)
        print(f"\n=== {run}  w/d={wd}  W={W*1e3:.3f} mm  q_fp={qfp*1e-3:.1f} kW/m²  "
              f"拡散長({a.t_eval}s)={np.sqrt(ALPHA_S*a.t_eval)/W:.2f} W ===")
        eps = Q_PRECISION_W_M2 / qfp
        print(f"  実測の精度限界 ±{Q_PRECISION_W_M2*1e-3:.2f} kW/m² = ±{eps:.4f} (q/q_fp)")
        print(f"{'z/W':>6} {'x/d':>6} {'実測±限界':>17} {'CFD 生':>9} | "
              f"{'リップ断熱':>10} {'リップ Dir':>10} | {'帯が重なるか':>13}")
        for zt in (0.5, 1, 2, 3, 4, 6, 8):
            sz = zt * W
            if sz > d:
                continue
            i = int(np.argmin(np.abs(s - sz)))
            xdi = sz / d
            rq = (np.interp(xdi, R[:, 0], R[:, 1])
                  if R[:, 0].min() <= xdi <= R[:, 0].max() else np.nan)
            lo = adb[1][i] / qfp; hi = out[a.t_eval][1][i] / qfp
            lo, hi = min(lo, hi), max(lo, hi)
            # **実測にも精度限界の区間がある**。両方の区間が重なるかで判定する
            # (2026-09-20 codex result-2 M9: 点値で「届かない」と判定していた)。
            mlo, mhi = rq - eps, rq + eps
            ov = "YES" if (lo <= mhi and mlo <= hi) else ("低すぎ" if mhi < lo else "高すぎ")
            print(f"{zt:6.1f} {xdi:6.3f} {rq:8.4f}±{eps:.4f} {qc[i]/qfp:9.4f} | "
                  f"{lo:10.4f} {hi:10.4f} | {ov:>13}")
        print(f"  リップ温度上昇 ({a.t_eval}s) = {qfp*a.t_eval/(RHO_C*TAU):.2f} K"
              f"  (本文: 計測区間の最大温度上昇 < 22 K)")
    print("\n値は q/q_fp。**全点・全幅に同一の観測モデル**を掛け、リップ端の両極 "
          "(断熱 / 表面温度 Dirichlet) で挟んでいる。端条件は原典から決まらないので、"
          "この帯が観測モデルの不定性そのもの。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
