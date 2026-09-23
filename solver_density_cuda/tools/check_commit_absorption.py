#!/usr/bin/env python3
"""commit の丸め吸収を場で測る — その run は「足しても消える」領域にいるか。

    python3 solver_density_cuda/tools/check_commit_absorption.py <res_with_dq.h5> [--split y|x|none]

**背景**: 定常陰解法の commit は `Q = Q_N + dq` で `Q` は float32。
$\\lvert dq\\rvert < \\tfrac12\\mathrm{ULP}(Q)$ になった加算は**丸めで消える**ので、
残差が残っているのに場が動かなくなる (plans/active/time_integration-fp64-accumulator.md)。

**使い方**: 収束済みの場から **1 step だけ**回し、`output.extraFields` に
`[dq_block_old_0, dq_block_old_1, dq_block_old_2, dq_block_old_3, dq_block_old_4]`
を指定して `res_1.h5` を出す。本ツールはそれを読んで $\\lvert dq\\rvert/\\mathrm{ULP}(Q)$ を出す。

**読み方**: 0.5 未満の割合が大きいほど「その領域は丸めで止まっている」。
`case/56.gap_tp1187` の深いすきまでは **99.8 %** が 0.5 未満だった (中央値 0.146)。
一方 `case/48` の平板のように流れが活きている領域では 1 を大きく超える。

⚠ **`rms_dq_*` 列は使えない**: `residual_history.csv` に列はあるが `main.cpp` が常に 0 を書く。
"""
import argparse, sys
import h5py
import numpy as np

CONS = ["ro", "roUx", "roUy", "roUz", "roe"]

def ulp32(a):
    a = np.abs(np.asarray(a, dtype=np.float32))
    out = np.full(a.shape, np.inf, dtype=np.float64)
    m = a > 0
    out[m] = (np.nextafter(a[m], np.float32(np.inf)) - a[m]).astype(np.float64)
    return out

ap = argparse.ArgumentParser()
ap.add_argument("res", help="dq_block_old_* を含む res_*.h5 (1 step 回したもの)")
ap.add_argument("--split", default="y", choices=["y", "x", "none"], help="深さ方向の分割軸")
ap.add_argument("--bins", type=int, default=4, help="分割数")
a = ap.parse_args()

with h5py.File(a.res) as h:
    if "VALUE/dq_block_old_0" not in h:
        sys.exit("dq_block_old_* が無い。output.extraFields に 5 本を指定して 1 step 回すこと")
    c = h["MESH/COORD"][:].reshape(-1, 3)
    Q = {v: np.asarray(h[f"VALUE/{v}"][:]) for v in CONS if f"VALUE/{v}" in h}
    D = {v: np.asarray(h[f"VALUE/dq_block_old_{i}"][:]).astype(np.float64)
         for i, v in enumerate(CONS) if f"VALUE/dq_block_old_{i}" in h}

order = [v for v in CONS if v in Q and v in D]
n = len(next(iter(Q.values())))
axis = {"y": 1, "x": 0}.get(a.split)
if axis is None:
    groups = [("全域", np.ones(n, bool))]
else:
    co = c[:n, axis]
    edges = np.quantile(co, np.linspace(0, 1, a.bins + 1))
    groups = [("全域", np.ones(n, bool))]
    for k in range(a.bins):
        m = (co >= edges[k]) & (co <= edges[k + 1] if k == a.bins - 1 else co < edges[k + 1])
        groups.append((f"{'xyz'[axis]} {edges[k]:+.3e}..{edges[k+1]:+.3e}", m))

print(f"=== {a.res}  ({n} CV) ===")
print("  |dq| / ULP(Q) — 0.5 未満は**その加算が丸めで消える**\n")
hdr = "%-30s %8s" % ("領域", "CV 数")
for v in order: hdr += " %10s" % v
print(hdr + "   " + " ".join("%6s" % (v + "<.5") for v in order))
worst = 0.0
for lbl, m in groups:
    if m.sum() < 5: continue
    line = "%-30s %8d" % (lbl[:30], m.sum())
    fr = []
    for v in order:
        r = np.abs(D[v][m]) / ulp32(Q[v][m])
        r = r[np.isfinite(r)]
        if len(r) == 0: line += " %10s" % "-"; fr.append("-"); continue
        line += " %10.3f" % np.median(r)
        f = (r < 0.5).mean()
        fr.append("%5.1f%%" % (f * 100))
        if lbl != "全域": worst = max(worst, f)
    print(line + "   " + " ".join("%6s" % x for x in fr))
print(f"\n  最悪の帯で 0.5 未満だった割合: **{worst*100:.1f} %**")
print("  VERDICT:", "吸収している (丸めで止まっている領域がある)" if worst > 0.5
      else "吸収は目立たない" if worst > 0.1 else "吸収していない")
