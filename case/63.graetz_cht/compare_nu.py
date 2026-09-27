#!/usr/bin/env python3
r"""case/63 Graetz の run 間比較 (plan `boundary-cht-axisymmetric-graetz.md` §6 V-g3 / V-g4 / V-g5 / 固体層感度)。

入力は `eval_graetz.py snap` が書く `graetz_nu_<step>.csv` (各 run の差し引き Nu と、その格子の区間平均基準 Nu_ref)。

    pair A.csv B.csv --tol T        : 同じ格子の 2 run。窓内の全壁節点で max |Nu_A/Nu_B − 1| ≤ T (V-g4 0.3 %、V-g5 0.5 %、固体層 0.05 %)
    grid C16.csv C32.csv C64.csv    : 入れ子の 3 格子。正規化誤差 e_h = Nu_h/Nu_ref,h − 1 の**格子間差**
                                      |e16 − e32| > |e32 − e64| を 4 観測点 (x⁺ 3e-3, 1e-2, 3e-2, 0.1、対数補間) と
                                      共通節点 (N_r=16 の窓内壁節点) の最大で要求する (V-g3)。観測次数は参考

e_h の差を取るのは、基準解のモデル差 (有限 Pe) が共通に消え、各格子で双対面の区間幅が違う (= 基準の区間平均が違う)
ことも吸収するため (codex plan レビュー M2)。
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

PTS = (3e-3, 1e-2, 3e-2, 0.1)


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d


def window(d, lo, hi):
    return (d["xplus"] >= lo) & (d["xplus"] <= hi)


def at(d, p, key):
    m = d["xplus"] > 0
    return float(np.interp(np.log(p), np.log(d["xplus"][m]), d[key][m]))


def cmd_pair(a):
    A, B = load(a.a), load(a.b)
    if len(A) != len(B) or not np.allclose(A["x"], B["x"], rtol=0, atol=1e-9):
        print("REFUSED: 2 run の壁節点が一致しない (同じ格子でない)"); print("VERDICT: REFUSED"); return 2
    W = window(A, a.lo, a.hi) & window(B, a.lo, a.hi)
    r = A["Nu"][W] / B["Nu"][W] - 1
    if not np.isfinite(r).all():
        print("VERDICT: FAIL (非有限)"); return 1
    i = int(np.argmax(np.abs(r)))
    v = float(np.abs(r).max())
    print(f"pair: {a.a}\n      {a.b}\n  窓内 {W.sum()} 節点、max |Nu_A/Nu_B − 1| = {v:.4e} (x⁺ {A['xplus'][W][i]:.4g})、許容 {a.tol:g}")
    for p in PTS:
        print(f"  参考 x⁺={p:g}: {at(A, p, 'Nu') / at(B, p, 'Nu') - 1:+.4e}")
    ok = v <= a.tol
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def cmd_grid(a):
    C = [load(p) for p in (a.c16, a.c32, a.c64)]
    e = [d["err"] for d in C]
    bad = []
    print(f"grid: {a.c16} / {a.c32} / {a.c64}")
    print("  x⁺        e16        e32        e64       |e16−e32|  |e32−e64|  観測次数(参考)")
    for p in PTS:
        v = [at(d, p, "err") for d in C]
        d1, d2 = abs(v[0] - v[1]), abs(v[1] - v[2])
        order = np.log2(d1 / d2) if d2 > 0 and d1 > 0 else np.nan
        ok = d1 > d2
        bad += [] if ok else [f"x⁺={p:g}"]
        print(f"  {p:8.1e} {v[0]:+.3e} {v[1]:+.3e} {v[2]:+.3e}  {d1:.3e}  {d2:.3e}  {order:6.2f}  {'ok' if ok else 'NG'}")
    # 共通節点: N_r=16 の窓内壁節点 (入れ子なので 32・64 にも同じ x の節点がある)
    W = window(C[0], a.lo, a.hi)
    xs = C[0]["x"][W]
    def pick(d):
        idx = [int(np.argmin(np.abs(d["x"] - x))) for x in xs]
        if max(abs(d["x"][i] - x) for i, x in zip(idx, xs)) > 1e-7:
            raise SystemExit("REFUSED: 共通節点が入れ子になっていない")
        return d["err"][idx]
    e16, e32, e64 = pick(C[0]), pick(C[1]), pick(C[2])
    d1, d2 = np.abs(e16 - e32).max(), np.abs(e32 - e64).max()
    ok = d1 > d2
    bad += [] if ok else ["共通節点の最大"]
    print(f"  共通 {len(xs)} 節点の最大: |e16−e32| {d1:.3e} > |e32−e64| {d2:.3e} ? {'ok' if ok else 'NG'}  (観測次数 参考 {np.log2(d1/d2):.2f})")
    print("  (観測次数は 3 格子から 1 つしか出ないので漸近域の確認ではない。合否は格子間差の減少のみ)")
    print(f"VERDICT: {'PASS' if not bad else 'FAIL'}" + (f"  ({', '.join(bad)})" if bad else ""))
    return 0 if not bad else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("pair"); p.add_argument("a"); p.add_argument("b"); p.add_argument("--tol", type=float, required=True)
    g = sub.add_parser("grid"); g.add_argument("c16"); g.add_argument("c32"); g.add_argument("c64")
    for q in (p, g):
        q.add_argument("--lo", type=float, default=3e-3); q.add_argument("--hi", type=float, default=0.1)
    a = ap.parse_args()
    return cmd_pair(a) if a.mode == "pair" else cmd_grid(a)


if __name__ == "__main__":
    sys.exit(main())
