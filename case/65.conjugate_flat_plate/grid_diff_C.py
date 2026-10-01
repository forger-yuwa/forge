#!/usr/bin/env python3
r"""case/65 C の forge 同士の格子間差 (後継 plan §4.6.2 m3)。隣接格子の run の `eval_conj_150000_L2.csv` から、
共通節点 (x の一致 1e-9) で窓 x/L∈[0.2,0.9] の max|Δq|/(細格子の窓内平均 |q|) と max|ΔT| を出す。参照解との差とは別の量。

    python3 grid_diff_C.py
"""
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
L = 0.01


def load(r):
    lines = [l for l in open(HERE / r / "eval_conj_150000_L2.csv") if l.strip()]
    hdr = lines[0].split("#")[0].strip().split(",")
    a = np.array([[float(v) for v in l.split(",")] for l in lines[1:]])
    d = {h: a[:, i] for i, h in enumerate(hdr)}
    o = np.argsort(d["x"])
    return d["x"][o], d["q_forge"][o], d["T_forge"][o]


for case, runs in (("C1", ["run_0024_c1_n16_lim0", "run_0026_c1_n32_lim0", "run_0022_c1_n64_lim0"]),
                   ("C2", ["run_0025_c2_n16_lim0", "run_0027_c2_n32_lim0", "run_0023_c2_n64_lim0"])):
    out = []
    for a, b in zip(runs[:-1], runs[1:]):
        xa, qa, Ta = load(a); xb, qb, Tb = load(b)
        idx = np.array([int(np.argmin(np.abs(xb - x))) for x in xa]); ok = np.abs(xb[idx] - xa) < 1e-9
        w = ok & (xa >= 0.2 * L - 1e-12) & (xa <= 0.9 * L + 1e-12)
        qm = np.mean(np.abs(qb[idx][w]))
        out.append(f"{a} → {b}: 共通節点 {w.sum()}、max|Δq|/q̄ {np.abs(qa[w] - qb[idx][w]).max() / qm * 100:.3f} %、max|ΔT| {np.abs(Ta[w] - Tb[idx][w]).max():.3e} K")
    print(case); print("  " + "\n  ".join(out))
