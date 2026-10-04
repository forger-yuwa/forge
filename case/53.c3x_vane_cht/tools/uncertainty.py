#!/usr/bin/env python3
r"""NASA CR-168015 が公表している**実験不確かさ** (表 V / VI / VII) — ページ画像から転記。

出典 (PDF ページ番号は 0 起点の索引、括弧内は報告の印字ページ):

- **表 V** (PDF idx 32 / 印字 p.26): Mark II カスケードの熱伝達係数の不確かさ。弧長 % の帯ごと。
- **表 VI** (PDF idx 33 / 印字 p.27): C3X カスケードの同じ表。
- **表 VII** (PDF idx 34 / 印字 p.28): 試験パラメータ (Re, MN, Tw/Tg, Tu) の不確かさ。
- **図 17** (印字 p.27): Mark II run 46 の $h$ 分布に各点の不確かさ幅を描いた図。**本ファイルの
  帯データはその図の根拠と同じもの**。

報告本文 (p.28) の重要な但し書き:

> The uncertainties presented in this subsection are intended to provide the analyst with an
> indication of **the uncertainty in absolute level** in utilizing the data for verification
> purposes. In comparing data from runs for a given cascade (that is, looking for Re trends,
> etc.), the uncertainty in the comparisons is considerably less than the values in Tables V
> and VI.

CFD を実測の**絶対値**と比べる用途はまさに前者なので、表 V / VI をそのまま使ってよい。
run 間のトレンド比較に使うと過大評価になる。

**原典の誤植**: 表 V の負圧面 7 行目は "73-32" と印字されているが、前後が 63-73 / 82-91 なので
**73-82** の誤り。ここでは 73-82 として扱う (`MARKII_SS` のコメント参照)。

usage: python3 case/53.c3x_vane_cht/tools/uncertainty.py   (自己検査)
"""
from __future__ import annotations

import sys

# ---- 表 V: Mark II (弧長 % の下限, 上限, ±%) ----
MARKII_PS = [(0, 20, 8.4), (20, 29, 6.9), (29, 42, 8.4), (42, 55, 10.0),
             (55, 67, 16.7), (67, 78, 14.4), (78, 88, 18.8), (88, 100, 18.2)]
MARKII_SS = [(0, 18, 9.0), (18, 32, 8.1), (32, 42, 7.1), (42, 52, 7.7),
             (52, 63, 10.0), (63, 73, 12.6), (73, 82, 10.4),   # 原典 "73-32" は誤植
             (82, 91, 15.8), (91, 100, 15.4)]

# ---- 表 VI: C3X ----
C3X_PS = [(0, 16, 6.8), (16, 23, 6.3), (23, 34, 6.6), (34, 45, 7.3), (45, 55, 8.9),
          (55, 66, 13.3), (66, 78, 11.6), (78, 89, 20.1), (89, 100, 23.5)]
C3X_SS = [(0, 8, 6.8), (8, 31, 6.7), (31, 39, 6.2), (39, 49, 6.5), (49, 58, 7.1),
          (58, 67, 8.6), (67, 76, 11.9), (76, 85, 10.9), (85, 94, 15.8), (94, 100, 23.5)]

TABLES = {("markii", "PS"): MARKII_PS, ("markii", "SS"): MARKII_SS,
          ("c3x", "PS"): C3X_PS, ("c3x", "SS"): C3X_SS}

# ---- 表 VII: 試験パラメータの不確かさ (±%) ----
PARAM = {"Re": 3.1, "M": 0.9, "Tw_over_Tg": 2.0, "Tu": 10.0}
PRESSURE_ABS_KPA = 0.7    # Scanivalve + Mensor 校正。報告 p.27 "+0.7 kPA (0.1 psi)"


def h_uncertainty_pct(vane: str, is_suction, s_norm):
    """弧長 $s/S$ (0–1) に対する $h$ の不確かさ [±%]。配列でもスカラでも受ける。"""
    import numpy as np
    s = np.atleast_1d(np.asarray(s_norm, dtype=float))
    ss = np.atleast_1d(np.asarray(is_suction))
    if ss.size == 1 and s.size > 1:
        ss = np.repeat(ss, s.size)
    out = np.full(s.shape, np.nan)
    for side in ("PS", "SS"):
        tbl = TABLES[(vane, side)]
        m = ss if side == "SS" else ~ss
        pct = s * 100.0
        for lo, hi, u in tbl:
            sel = m & (pct >= lo) & (pct <= hi)
            out[sel] = u
    # 100 % を超える節点 (後縁の丸め) は最後の帯に寄せる
    for side in ("PS", "SS"):
        tbl = TABLES[(vane, side)]
        m = (ss if side == "SS" else ~ss) & np.isnan(out)
        out[m] = tbl[-1][2]
    return out if np.ndim(s_norm) else float(out[0])


def validate() -> int:
    bad = 0
    for key, tbl in TABLES.items():
        # 帯が 0 から 100 まで隙間なく連続しているか
        if tbl[0][0] != 0 or tbl[-1][1] != 100:
            print(f"[uncertainty] {key}: 帯が 0–100 を覆っていない"); bad += 1
        for (l0, h0, _), (l1, h1, _) in zip(tbl, tbl[1:]):
            if h0 != l1:
                print(f"[uncertainty] {key}: 帯が不連続 {h0} -> {l1}"); bad += 1
            if h0 <= l0:
                print(f"[uncertainty] {key}: 帯が逆順 {l0}-{h0}"); bad += 1
        for _, _, u in tbl:
            if not (5.0 <= u <= 25.0):
                print(f"[uncertainty] {key}: ±{u} % は想定範囲外"); bad += 1
    # 両翼とも後縁側で不確かさが増える (報告本文: 翼厚が薄くなるため)
    for vane in ("markii", "c3x"):
        for side in ("PS", "SS"):
            tbl = TABLES[(vane, side)]
            if tbl[-1][2] <= tbl[0][2]:
                print(f"[uncertainty] {vane}/{side}: 後縁側が前縁側より小さい"); bad += 1
    print("[uncertainty] 表 V/VI/VII 検査: " + ("PASS" if bad == 0 else f"FAIL ({bad})"))
    return bad


if __name__ == "__main__":
    sys.exit(1 if validate() else 0)
