#!/usr/bin/env python3
r"""収束の判定 (案) — 「残差が準定常・流量の収支が取れる・目的の量が準定常」の AND (plan tooling-nozzle-card-local-mesh §4.4・§5.1 #4b、
ユーザ決定 2026-10-15「残さが準定常になり、流量バランスが取れ、目的とする物理量が準定常担ったら収束と判定しないか」)。

この版は**残差の準定常**の部分だけ (流量の収支と目的の量は forge の出力の追加 #4a の後に足す)。低下の桁数は見ない。

残差の基準 (codex diagnose 2026-10-15 m55-settle-and-monitor M1: log10 の残差に `check_quasisteady` (系列の平均で割る) を使わない):
- 列ごとに y = log10(rms) を、本段の outer_end の行から取る (step の間隔は履歴のまま)。全履歴に非有限があれば DIVERGED。
- 末尾の窓 W [step] (既定 20,000) に点が `min_points` (既定 11) 以上あること (足りなければ判定不能)。
- 傾き: 窓の y の step に対する 1 次の最小二乗の傾き s で |s|·W ≤ `slope_dec` (既定 log10(1.05) ≈ 0.0212 桁)。
- 振幅: 窓の y の p95 − p5 ≤ `amp_dec` (既定 0.1 桁)。
- 持続的な上昇: 直前の窓 (W) の平均から末尾の窓の平均への上昇 ≤ `rise_dec` (既定 log10(1.10) ≈ 0.0414 桁; 4 万 step で 10 % を超える上昇を通さない)。
- 恒等ゼロの列は非活動として外す。途中からゼロ・一部ゼロは異常として判定不能。
閾値は**検証済みの既定ではない** (試して決める; plan §4.4 の修正 M1)。

使い方:
  python3 check_settled.py <run_dir|residual_history.csv[.gz]> [--window 20000] [--scan 2000]
  --scan N: N step ごとに「その時点までの履歴」で判定し、最初に残差の準定常を満たした step を報告する (到達の時点の比較用)。
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

import numpy as np

SLOPE_DEC = math.log10(1.05)
RISE_DEC = math.log10(1.10)
AMP_DEC = 0.1


def read_residuals(path) -> tuple[np.ndarray, dict]:
    """residual_history.csv (gz 可) の outer_end の行 → (steps, {列: 値})。同じ step が重なれば最後の行。"""
    p = Path(path)
    if p.is_dir():
        p = p / "residual_history.csv"
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        r = csv.reader(f)
        h = next(r)
        ic = [i for i, k in enumerate(h) if k.startswith("rms_") and not k.startswith("rms_dq_")]
        rows = {}
        for row in r:
            if len(row) < 3 or row[2] != "outer_end":
                continue
            try:
                rows[int(row[0])] = [float(row[i]) for i in ic]
            except ValueError:
                rows[int(row[0])] = [float("nan")] * len(ic)
    st = np.array(sorted(rows))
    v = np.array([rows[s] for s in st]) if len(st) else np.zeros((0, len(ic)))
    return st, {h[i]: v[:, k] for k, i in enumerate(ic)}


def residual_settled(steps, cols: dict, window: int = 20000, min_points: int = 11, slope_dec: float = SLOPE_DEC,
                     amp_dec: float = AMP_DEC, rise_dec: float = RISE_DEC) -> dict:
    """残差の準定常の判定。戻り: {verdict: SETTLED|NOT_SETTLED|DIVERGED|INDETERMINATE, columns: {...}, reasons: [...]}。"""
    steps = np.asarray(steps, dtype=float)
    out = {"window": window, "thresholds": {"slope_dec": slope_dec, "amp_dec": amp_dec, "rise_dec": rise_dec, "min_points": min_points},
           "columns": {}, "reasons": []}
    if len(steps) == 0:
        out.update(verdict="INDETERMINATE", reasons=["残差の行が無い"]); return out
    last = steps[-1]
    w1 = steps > last - window
    w0 = (steps > last - 2 * window) & ~w1
    verdicts = []
    for k, v in cols.items():
        v = np.asarray(v, dtype=float)
        if not np.all(np.isfinite(v)):
            out["columns"][k] = {"verdict": "DIVERGED"}; verdicts.append("DIVERGED"); continue
        if np.all(v == 0.0):
            out["columns"][k] = {"verdict": "INACTIVE"}; continue
        if np.any(v <= 0.0):
            out["columns"][k] = {"verdict": "INDETERMINATE", "why": "途中・一部がゼロ"}; verdicts.append("INDETERMINATE"); continue
        y = np.log10(v)
        if w1.sum() < min_points:
            out["columns"][k] = {"verdict": "INDETERMINATE", "why": f"窓の点 {int(w1.sum())} < {min_points}"}; verdicts.append("INDETERMINATE"); continue
        s = float(np.polyfit(steps[w1], y[w1], 1)[0])
        amp = float(np.percentile(y[w1], 95) - np.percentile(y[w1], 5))
        rise = float(np.mean(y[w1]) - np.mean(y[w0])) if w0.sum() >= min_points else float("nan")
        ok_s = abs(s) * window <= slope_dec
        ok_a = amp <= amp_dec
        ok_r = (not np.isfinite(rise)) or rise <= rise_dec
        why = [w for w, ok in (("傾き", ok_s), ("振幅", ok_a), ("上昇", ok_r)) if not ok]
        vd = "SETTLED" if not why else "NOT_SETTLED"
        if not np.isfinite(rise):
            why_r = "直前の窓の点が足りず上昇は見ていない"
        else:
            why_r = None
        out["columns"][k] = {"verdict": vd, "slope_dec_per_window": s * window, "amp_dec": amp, "rise_dec": rise,
                             "fail": why, **({"note": why_r} if why_r else {}), "last": float(v[-1])}
        verdicts.append(vd)
    if not verdicts:
        out.update(verdict="INDETERMINATE", reasons=["活動している残差の列が無い"]); return out
    if "DIVERGED" in verdicts:
        out["verdict"] = "DIVERGED"
    elif "INDETERMINATE" in verdicts:
        out["verdict"] = "INDETERMINATE"
    elif all(v == "SETTLED" for v in verdicts):
        out["verdict"] = "SETTLED"
    else:
        out["verdict"] = "NOT_SETTLED"
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src")
    ap.add_argument("--window", type=int, default=20000)
    ap.add_argument("--min-points", type=int, default=11)
    ap.add_argument("--scan", type=int, default=0)
    a = ap.parse_args(argv)
    st, cols = read_residuals(a.src)
    res = residual_settled(st, cols, a.window, a.min_points)
    print(f"=== {a.src}  [last step {int(st[-1]) if len(st) else '-'}]  残差の準定常 (窓 {a.window}): {res['verdict']} ===")
    for k, c in res["columns"].items():
        if c["verdict"] in ("INACTIVE",):
            print(f"  {k:12s}: 非活動 (恒等ゼロ)"); continue
        if "slope_dec_per_window" not in c:
            print(f"  {k:12s}: {c['verdict']} {c.get('why', '')}"); continue
        print(f"  {k:12s}: 末尾 {c['last']:.2e}  傾き {c['slope_dec_per_window']:+.4f} 桁/窓  振幅 {c['amp_dec']:.4f} 桁  上昇 {c['rise_dec']:+.4f} 桁  {c['verdict']}"
              + (f"  (不成立: {', '.join(c['fail'])})" if c['fail'] else ""))
    if a.scan:
        first = None
        for s in range(a.scan, int(st[-1]) + 1, a.scan):
            m = st <= s
            r = residual_settled(st[m], {k: v[m] for k, v in cols.items()}, a.window, a.min_points)
            if r["verdict"] == "SETTLED":
                first = s; break
        print(f"最初に残差の準定常を満たした step (--scan {a.scan}): {first}")
    print(f"VERDICT: {res['verdict']}")
    return 0 if res["verdict"] == "SETTLED" else 1


if __name__ == "__main__":
    sys.exit(main())
