#!/usr/bin/env python3
r"""case/65 B3 (リミッタの有無 `space.limiter` 2 / 0) の判定 (plan `boundary-cht-conjugate-flat-plate.md` §4.4、事前登録どおり)。

    python3 ab_judge_lim.py <run_A (limiter 2)> <run_B (limiter 0)> <対照の元 (run_0016)>

- 区間: A・B とも 4000 step × 3。
- 能動な残差: `rms_ro/rms_roUx/rms_roUy/rms_roe` の区間算術平均 (`residual_history.csv` の outer_end 行)。
- 下流 slip の P の振幅: `ab_series.npz` の下流 slip 群の P (q∞ で規格化) の区間内標準偏差の節点最大。
- 前提 (対照の再現): A の最終区間の残差 (全列) と下流 slip の P の振幅が、元 run の最後の 4000 step の値の ±10 % 以内。
- 支持: B の最後の 2 区間で全残差と下流 P の振幅が A の最終区間の 1/10 以下、かつ B の 3 区間で厳密減少
  → リミッタを含む再構成への依存を支持。
- 棄却: B の全区間で全指標が A の最終区間の 1/2 以上、かつ B の 3 区間の max/min ≤ 1.1 → 「リミッタ単独が原因」を棄却。
- それ以外は中間 (判定不能)。
- 併記 (判定外): 残差二乗和の領域別割合 (前縁帯・後縁帯・その他、時間積算)、全域 RMS と CSV の照合、
  最終 1/3 区間のリミッタ係数の時間標準偏差が大きい節点の位置、res² の大きい節点の位置。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plate_common as pc  # noqa: E402
from ab_judge_cfl import COLS, Q_INF, fmt  # noqa: E402
from ab_series import refuse  # noqa: E402

NI, LI = 3, 4000


def resid(run, lo, hi, li):
    import csv
    v = {c: [] for c in COLS}
    with open(Path(run) / "residual_history.csv") as f:
        for row in csv.DictReader(f):
            if row["phase"] != "outer_end":
                continue
            s = int(row["step"])
            if lo <= s < hi:
                for c in COLS:
                    x = float(row[c])
                    if not np.isfinite(x):
                        refuse(f"{run}: {c} に非有限値 (step {s})")
                    v[c].append(x)
    out = {}
    for c in COLS:
        a = np.array(v[c])
        if len(a) != hi - lo:
            refuse(f"{run}: {c} の行数 {len(a)} (期待 {hi - lo})")
        out[c] = a.reshape(-1, li).mean(axis=1)
    return out, {c: np.array(v[c]) for c in COLS}


def down_amp(run, lo, hi, li):
    d = np.load(Path(run) / "ab_series.npz", allow_pickle=True)
    st = d["step"]; g = d["group"]
    m = (st > lo) & (st <= hi)
    if m.sum() != hi - lo:
        refuse(f"{run}: 系列の点数 {m.sum()} (期待 {hi - lo})")
    P = d["P"][m][:, g == "down"] / Q_INF
    return np.array([np.std(P[k * li:(k + 1) * li], axis=0).max() for k in range(len(P) // li)]), d


def main():
    rA, rB, r0 = (Path(p) for p in sys.argv[1:4])
    print(f"=== B3 判定 (limiter 2 vs 0、cfl_pseudo 0.5)  A={rA.name}  B={rB.name}  対照の元={r0.name} (最後の 4000 step)")
    n0 = int(np.load(r0 / "ab_series.npz", allow_pickle=True)["step"][-1])
    R0, _ = resid(r0, n0 - LI, n0, LI)
    D0, _ = down_amp(r0, n0 - LI, n0, LI)
    RA, rawA = resid(rA, 0, NI * LI, LI); RB, rawB = resid(rB, 0, NI * LI, LI)
    DA, dA = down_amp(rA, 0, NI * LI, LI); DB, dB = down_amp(rB, 0, NI * LI, LI)
    print("  残差の区間平均:")
    for c in COLS:
        print(f"    {c:9s} 元 {R0[c][0]:.3e} | A {fmt(RA[c])} | B {fmt(RB[c])}")
    print(f"    下流 slip の P 振幅  元 {D0[0]:.3e} | A {fmt(DA)} | B {fmt(DB)}")
    rep = all(abs(RA[c][-1] / R0[c][0] - 1) <= 0.10 for c in COLS) and abs(DA[-1] / D0[0] - 1) <= 0.10
    print(f"  対照の再現 (A 最終区間 / 元, ±10 %): {'OK' if rep else 'NG'} ("
          + ", ".join(f"{c} {RA[c][-1]/R0[c][0]:.3f}" for c in COLS) + f", 下流 P {DA[-1]/D0[0]:.3f})")
    allA = {**RA, "down_P": DA}; allB = {**RB, "down_P": DB}
    print("  B/A (A 最終区間基準): " + ", ".join(f"{k} " + "/".join(f"{v:.3f}" for v in allB[k] / allA[k][-1]) for k in allB))
    low = all((allB[k][1:] <= allA[k][-1] / 10).all() and allB[k][0] > allB[k][1] > allB[k][2] for k in allB)
    high = all((allB[k] >= allA[k][-1] / 2).all() and allB[k].max() / allB[k].min() <= 1.1 for k in allB)
    print("  併記 (判定外):")
    for tag, d, raw in (("A", dA, rawA), ("B", dB, rawB)):
        if "region_res2" not in d:
            print(f"    {tag}: 領域統計なし"); continue
        rr = d["region_res2"]  # [step, var, region]
        tot = rr.sum(axis=0)
        frac = tot / tot.sum(axis=1, keepdims=True)
        names = [str(x) for x in d["region_names"]]; cnt = d["region_count"]
        for j, nm in enumerate(d["res_names"]):
            print(f"    {tag} {str(nm):9s} 残差二乗和の割合: " + ", ".join(f"{names[i]} {frac[j, i]*100:.1f} % ({cnt[i]} 節点)" for i in range(len(names))))
        # 全域 RMS と CSV の照合 (出力 step n は残差 step n-1 と対応)
        nn = int(d["n_nodes"])
        rms = np.sqrt(rr.sum(axis=2) / nn)  # [step, var]
        csvv = np.array([raw[c] for c in COLS]).T  # [step, var]
        m = min(len(rms), len(csvv))
        for off in (0, 1):
            a_, b_ = rms[:m - 1], csvv[off:off + m - 1]
            rel = np.abs(a_ / b_ - 1).max(axis=0)
            print(f"    {tag} 全域 RMS / CSV (ずれ {off}) の最大相対差: " + ", ".join(f"{c} {v:.2e}" for c, v in zip(COLS, rel)))
        c = d["coord"]; cnt_t = int(d["tail_count"])
        mu = d["tail_lim_sum"] / cnt_t; sd = np.sqrt(np.maximum(d["tail_lim_sumsq"] / cnt_t - mu ** 2, 0))
        for j, nm in enumerate(d["lim_names"]):
            o = np.argsort(sd[j])[::-1][:4]
            print(f"    {tag} {str(nm):10s} 時間標準偏差の上位: " + ", ".join(f"({c[i,0]/pc.L:.3f},{c[i,1]/pc.L:.4f}) {sd[j,i]:.2f}" for i in o))
        r2 = d["tail_res2_sum"]
        for j, nm in enumerate(d["res_names"]):
            o = np.argsort(r2[j])[::-1][:4]
            print(f"    {tag} {str(nm):9s} res² 上位節点 (x/L,y/L,割合): " + ", ".join(f"({c[i,0]/pc.L:.3f},{c[i,1]/pc.L:.4f},{r2[j,i]/r2[j].sum()*100:.0f}%)" for i in o))
    if not rep:
        print("VERDICT: 判定不能 (対照 A が run_0016 の状態を再現していない)")
        return 3
    if low:
        print("VERDICT: リミッタを含む再構成への依存を支持")
    elif high:
        print("VERDICT: 「リミッタ単独が原因」を棄却")
    else:
        print("VERDICT: 中間 (判定不能)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
