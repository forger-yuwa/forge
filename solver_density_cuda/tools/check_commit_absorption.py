#!/usr/bin/env python3
"""commit の丸め吸収を場で測る — その run は「足しても消える」領域にいるか。

    python3 solver_density_cuda/tools/check_commit_absorption.py <res_with_dq.h5> [--split y|x|none]

**背景**: 定常陰解法の commit は `Q = Q_N + dq` で `Q` は float32。
$\\lvert dq\\rvert < \\tfrac12\\mathrm{ULP}(Q)$ になった加算は**丸めで消える**ので、
残差が残っているのに場が動かなくなる (plans/active/time_integration-fp64-accumulator.md)。

**使い方**: 収束済みの場から **1 step だけ**回し、`output.extraFields` に
`[dq_block_old_0, dq_block_old_1, dq_block_old_2, dq_block_old_3, dq_block_old_4]`
を指定して `res_1.h5` を出す。本ツールはそれを読んで $\\lvert dq\\rvert/\\mathrm{ULP}(Q)$ を出す。

**読み方 — これは「壊れているか」の判定ではない** (2026-09-24 に訂正):
6 case で測ったところ、**健全な `case/48` 冷却平板でも最悪の帯で 92.2 %** が 0.5 未満だった。
**吸収は例外でなく、float32 で「収束した」run の normal な姿**である
(収束すれば定義上 $\lvert dq\rvert$ は小さくなり、いずれ $\tfrac12$ULP を下回る)。

実測 (最悪の帯で $\lvert dq_\rho\rvert<\tfrac12$ULP だった CV の割合 / 中央値):
`case/56` すきま **99.8 % / 0.146** ・ `case/53` 翼 CHT **93.5 % / 0.177** ・
`case/48` 平板 **92.2 % / 0.189** ・ `case/55` 乱流すきま 52.4 % / 0.473 ・
`case/50` 深キャビティ 22.6 % / 4.537。

**正しい問いは「吸収する前に物理の答えに着いていたか」**で、本ツールだけでは決まらない。
**判別するには `qAccumulatorFP64` の ON/OFF で答えが動くかを見る** (動けばその run は床律速だった)。
本ツールはその**入口**として「どの領域が床に当たっているか」を示すもの。

⚠ **`rms_dq_*` 列は使えない**: `residual_history.csv` に列はあるが `main.cpp` が常に 0 を書く。
"""
import argparse, sys
import h5py
import numpy as np

CONS = ["ro", "roUx", "roUy", "roUz", "roe"]

def ulp32(a):
    """**補助指標**。符号を見ない片側 ULP なので、判定には使わない (下の absorbed を使う)。"""
    a = np.abs(np.asarray(a, dtype=np.float32))
    out = np.full(a.shape, np.inf, dtype=np.float64)
    m = a > 0
    out[m] = (np.nextafter(a[m], np.float32(np.inf)) - a[m]).astype(np.float64)
    return out

def absorbed(q, dq):
    """**直接検査**: float32 で `q + dq` を実際に計算し、値が変わらなければ吸収。

    ⚠ ULP 比で判定してはいけない (2026-09-24, codex rollout plan M7):
      - 2 のべき乗の境界では**下向きの ULP が半分**なので、正方向 ULP で割ると誤判定する
        (`Q=1, dq=-4e-8` は比 0.336 だが **値は 0.99999994 に変わる**)。
      - `Q=0` の ULP を無限大にすると、非ゼロ増分でも比が 0 になって吸収と誤判定する。
    戻り値は (吸収したか, 増分がゼロか)。**増分ゼロは「吸収」と別に数える**
      (更新する気が無いだけで、丸めで消されたのではない)。
    """
    q = np.asarray(q, dtype=np.float32)
    d = np.asarray(dq, dtype=np.float32)
    zero = (d == np.float32(0.0))
    same = (np.float32(q) + d) == q
    return (same & ~zero), zero

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
# **全 CV で dq がちょうど 0 の成分は判定から外す** (2026-09-24)。
# 疑似 2D の roUz などは 0/ULP = 0 になり「0.5 未満」に数えられて VERDICT を汚す。
zero = [v for v in order if not np.any(np.asarray(D[v]) != 0.0)]
order = [v for v in order if v not in zero]
if zero:
    print(f"  ⚠ 全 CV で dq=0 のため判定から除外: {zero}")
if not order:
    sys.exit("判定できる成分が無い (全部 dq=0)")
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
print("  **直接検査**: float32 で q+dq を計算し値が変わらなければ吸収 (ULP 比は補助)\n")
hdr = "%-30s %8s" % ("領域", "CV 数")
for v in order: hdr += " %9s" % (v + " 吸収")
for v in order: hdr += " %8s" % (v + " dq=0")
print(hdr)
worst = 0.0
for lbl, m in groups:
    if m.sum() < 5: continue
    line = "%-30s %8d" % (lbl[:30], m.sum())
    zs = []
    for v in order:
        ab, z = absorbed(Q[v][m], D[v][m])
        nz = (~z).sum()
        f = ab.sum() / nz if nz > 0 else float("nan")
        line += " %8.1f%%" % (f * 100 if f == f else float("nan"))
        zs.append(z.mean() * 100)
        if lbl != "全域" and f == f: worst = max(worst, f)
    for x in zs: line += " %7.2f%%" % x
    print(line)
print(f"\n  最悪の帯の吸収率 (非ゼロ増分のうち丸めで消えた割合): **{worst*100:.1f} %**")
print("  (吸収率は **増分が非ゼロの CV だけ**を分母にした。増分ゼロは別列)")
print("  VERDICT:", "床に当たっている領域が広い" if worst > 0.5
      else "一部が床に当たっている" if worst > 0.1 else "床には当たっていない")
print("  ⚠ これは「結果が壊れている」判定ではない。**収束した run では正常にこうなる**")
print("     (健全な case/48 冷却平板でも最悪の帯で 92.2 %)。")
print("     答えが床律速だったかは **qAccumulatorFP64 の ON/OFF で場が動くか**で判別する。")
