#!/usr/bin/env python3
r"""case/65 B2′ (擬似時間刻み `cfl_pseudo` 2 ↔ 0.5 の A/B) の判定 (plan `boundary-cht-conjugate-flat-plate.md` §4.3、事前登録どおり)。

    python3 ab_judge_cfl.py <run_A (cfl 2, 12000 step)> <run_B (cfl 0.5, 48000 step)> <対照の元 (run_0012)>

- 区間: A は 4000 step × 3、B は 16000 step × 3 (CFL × step を揃える)。
- 残差: `residual_history.csv` の outer_end 行の能動列 `rms_ro/rms_roUx/rms_roUy/rms_roe` の区間算術平均。
- 主要な変動: `ab_series.npz` (毎 step・FP64) から、前縁帯の界面 q (評価窓平均 |q| で規格化)・下流 slip の P・上流 slip の P
  (q∞ で規格化) の、区間内の標準偏差の節点最大。
- 前提 (対照の再現): A の残差の区間平均が、元 (run_0012) の同じ区間 (0–12000 step、4000 × 3) の 1/2〜2 倍、かつ A が横ばい (max/min ≤ 1.1)。
- 支持: B の最後の 2 区間で残差の全列と主要な変動がすべて A の最終区間の 1/10 以下、かつ B の 3 区間で厳密減少
  → 「空間的な欠陥だけで同程度の振動が不可避」を棄却 (有限刻み依存を支持)。
- 棄却: B の全区間で残差の全列と主要な変動が A の最終区間の 1/2 以上、かつ B の区間平均・振幅の max/min ≤ 1.1
  → 「CFL 2 → 0.5 の縮小で主要な停滞を解消できる」を棄却。
- それ以外は中間 (判定不能)。発散 (非有限) は REFUSED。
- 併記 (判定外): 観測節点の累積 `dt_local` (中央値)、slip 上 `Ux` の空間の交互成分、リミッタ値の時間変動。
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plate_common as pc  # noqa: E402
from ab_series import q_scale, refuse  # noqa: E402

COLS = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe")
NI = 3
Q_INF = 709.275


def resid_means(run, li):
    vals = {c: [[] for _ in range(NI)] for c in COLS}
    with open(Path(run) / "residual_history.csv") as f:
        for row in csv.DictReader(f):
            if row["phase"] != "outer_end":
                continue
            s = int(row["step"])
            if s >= NI * li:
                break
            for c in COLS:
                v = float(row[c])
                if not np.isfinite(v):
                    refuse(f"{run}: {c} に非有限値 (step {s})")
                vals[c][s // li].append(v)
    for c in COLS:
        for k in range(NI):
            if len(vals[c][k]) != li:
                refuse(f"{run}: {c} 区間 {k} の行数 {len(vals[c][k])} (期待 {li})")
    return {c: np.array([np.mean(v) for v in vals[c]]) for c in COLS}


def series(run, li, qs):
    d = np.load(Path(run) / "ab_series.npz", allow_pickle=True)
    st = d["step"]
    if len(st) < NI * li or (np.diff(st[:NI * li]) != 1).any():
        refuse(f"{run}: 系列が {len(st)} 点・不連続 (期待 ≥ {NI * li})")
    g = d["group"]; wx = d["wall_x"] / pc.L
    def iv(a, k):
        return a[k * li:(k + 1) * li]
    main = {
        "前縁帯の界面 q": np.array([np.std(iv(d["q"][:, wx <= 0.02 + 1e-9] / qs, k), axis=0).max() for k in range(NI)]),
        "下流 slip の P": np.array([np.std(iv(d["P"][:, g == "down"] / Q_INF, k), axis=0).max() for k in range(NI)]),
        "上流 slip の P": np.array([np.std(iv(d["P"][:, g == "up"] / Q_INF, k), axis=0).max() for k in range(NI)]),
    }
    extra = {}
    if "dt_local" in d:
        extra["累積 dt_local 中央値 [s]"] = float(np.median(d["dt_local"][:NI * li].sum(axis=0)))
    if "Ux" in d:
        for grp in ("up", "down"):
            u = d["Ux"][:NI * li][:, g == grp] / pc.U_INF
            alt = np.abs(u[:, 1:-1] - 0.5 * (u[:, :-2] + u[:, 2:])).mean(axis=0).max()
            extra[f"{grp} slip の Ux 交互成分 (時間平均の節点最大、/U∞)"] = float(alt)
    for nm in ("limiter_ro", "limiter_P", "limiter_T", "limiter_Ux", "limiter_Uy"):
        if nm in d:
            a = d[nm][NI * li - li:NI * li]
            extra[f"{nm} 最終区間の時間標準偏差 (節点最大)"] = float(np.std(a, axis=0).max())
            extra[f"{nm} 最終区間の最小値"] = float(a.min())
    for nm in d.files:
        if d[nm].dtype.kind == "f" and not np.isfinite(d[nm]).all():
            refuse(f"{run}: 系列 {nm} に非有限値")
    return main, extra


def fmt(a):
    return " ".join(f"{v:.3e}" for v in a)


def main():
    rA, rB, r0 = (Path(p) for p in sys.argv[1:4])
    qs = q_scale()
    LA, LB = 4000, 16000
    print(f"=== B2′ 判定 (cfl_pseudo 2 vs 0.5)  A={rA.name} (4000×3)  B={rB.name} (16000×3)  対照の元={r0.name}")
    RA, RB, R0 = resid_means(rA, LA), resid_means(rB, LB), resid_means(r0, LA)
    print("  残差の区間平均:")
    for c in COLS:
        print(f"    {c:9s} 元 {fmt(R0[c])} | A {fmt(RA[c])} | B {fmt(RB[c])}")
    rep = all(0.5 <= RA[c][k] / R0[c][k] <= 2.0 for c in COLS for k in range(NI))
    flatA = all(RA[c].max() / RA[c].min() <= 1.1 for c in COLS)
    print(f"  対照の再現: A/元 ∈ [0.5,2] {'OK' if rep else 'NG'}、A 横ばい {'OK' if flatA else 'NG'} ("
          + ", ".join(f"{c} {RA[c].max()/RA[c].min():.3f}" for c in COLS) + ")")
    MA, XA = series(rA, LA, qs); MB, XB = series(rB, LB, qs)
    print("  主要な変動 (規格化標準偏差の節点最大、区間 1/2/3):")
    for k in MA:
        print(f"    {k:14s} A {fmt(MA[k])} | B {fmt(MB[k])}  B最終/A最終 {MB[k][-1]/MA[k][-1]:.3f}")
    print("  B/A (A 最終区間基準) 残差: " + ", ".join(f"{c} " + "/".join(f"{v:.3f}" for v in RB[c] / RA[c][-1]) for c in COLS))
    allB = {**{c: RB[c] for c in COLS}, **MB}; allA = {**{c: RA[c] for c in COLS}, **MA}
    low = all((allB[k][1:] <= allA[k][-1] / 10).all() and allB[k][0] > allB[k][1] > allB[k][2] for k in allB)
    high = all((allB[k] >= allA[k][-1] / 2).all() and allB[k].max() / allB[k].min() <= 1.1 for k in allB)
    print("  併記 (判定外):")
    for k in XA:
        print(f"    {k}: A {XA[k]:.4g} | B {XB.get(k, float('nan')):.4g}")
    if not (rep and flatA):
        print("VERDICT: 判定不能 (対照 A が run_0012 の停滞を再現していない)")
        return 3
    if low:
        print("VERDICT: 「空間的な欠陥だけで同程度の振動が不可避」を棄却 (有限刻み依存を支持)")
    elif high:
        print("VERDICT: 「CFL 2 → 0.5 の縮小で主要な停滞を解消できる」を棄却")
    else:
        print("VERDICT: 中間 (判定不能)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
