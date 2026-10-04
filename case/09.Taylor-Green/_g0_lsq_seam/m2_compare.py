#!/usr/bin/env python3
"""plan gradient-scalar-lsq-unification §5.1 #5g (codex result-1 M2): 実装直前版 `36d8ba03` と現行 gg の case/48 比較。

判定 (投入前に固定):
  (1) 延長段の全活動残差列の末尾 20 % の |値| 平均 f で、各列 2/3 ≤ f_old/f_new ≤ 1.5。
      列の活動状態が両者で違う・NaN/Inf は判定不能。
  (2) Cf・q_w・δ*・θ (x 0.3/0.6/0.9) の**延長段 step 12000–24000 のスナップショット平均**で、旧版を分母に相対差 ≤ 0.1 %。
      同じ区間での各量の (max−min)/|mean| が 0.1 % を超える量があれば、0.1 % の判別ができないので判定不能。
  (3) 系列 STEADY・既存物理ゲート・NaN/RISING は `cooled_plate_eval.py --series`・`check_convergence.py` の出力を別に貼る。

  python3 m2_compare.py --old RUN_old_ext --new RUN_new_ext [--out M2_case48.txt]
"""
import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "case", "48.flat_plate_cooled_m4", "tools"))
sys.path.insert(0, os.path.join(REPO, "design"))
import cooled_plate_eval as CPE  # noqa: E402

STATIONS = (0.3, 0.6, 0.9)
KEYS = ("Cf", "qw", "dstar", "theta")
STEP_LO, STEP_HI = 12000, 24000


def resid_tail(run, frac=0.2):
    rows = [r for r in csv.DictReader(open(os.path.join(run, "residual_history.csv"))) if r.get("phase") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq_")]
    out = {}
    for c in cols:
        ser = [float(r[c]) for r in rows]
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            out[c] = float("nan"); continue
        if not any(v != 0.0 for v in ser):
            out[c] = 0.0; continue
        a = ser[int(len(ser) * (1 - frac)):]
        out[c] = sum(abs(v) for v in a) / len(a)
    return out


def phys_tail(run):
    fs = [f for f in glob.glob(os.path.join(run, "res_[0-9]*.h5"))
          if STEP_LO <= int(os.path.basename(f)[4:-3]) <= STEP_HI]
    fs.sort(key=lambda f: int(os.path.basename(f)[4:-3]))
    vals = {(x, k): [] for x in STATIONS for k in KEYS}
    for f in fs:
        D = CPE.load_forge(f)
        for x in STATIONS:
            st = CPE.station(D, x)
            for k in KEYS:
                vals[(x, k)].append(st[k])
    return {kk: np.array(v) for kk, v in vals.items()}, [os.path.basename(f) for f in fs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
    ap.add_argument("--out", default="M2_case48.txt")
    a = ap.parse_args()
    L = [f"# #5g 実装直前版 vs 現行 gg (case/48)\n- old: {a.old}\n- new: {a.new}"]
    ok, undecided = True, False
    fo, fn = resid_tail(a.old), resid_tail(a.new)
    L.append("\n## (1) 残差床比 f_old/f_new (延長段 末尾 20 % の |値| 平均、許容 [2/3, 1.5])")
    L.append("| 列 | f_old | f_new | 比 | 判定 |"); L.append("| --- | --- | --- | --- | --- |")
    for c in sorted(set(fo) | set(fn)):
        o, n = fo.get(c), fn.get(c)
        if o is None or n is None or (o == 0.0) != (n == 0.0) or math.isnan(o or 0) or math.isnan(n or 0):
            L.append(f"| {c} | {o} | {n} | - | 判定不能 |"); undecided = True; continue
        if o == 0.0 and n == 0.0:
            L.append(f"| {c} | 0 | 0 | - | 両方非活動 |"); continue
        r = o / n; good = 2.0 / 3.0 <= r <= 1.5; ok &= good
        L.append(f"| {c} | {o:.3e} | {n:.3e} | {r:.3f} | {'ok' if good else 'NG'} |")
    po, so = phys_tail(a.old); pn, sn = phys_tail(a.new)
    L.append(f"\n## (2) 物理量 (延長段 step {STEP_LO}–{STEP_HI} の平均、旧版を分母、許容 0.1 %)")
    L.append(f"old snapshots {so}\nnew snapshots {sn}")
    if len(so) != len(sn) or not so:
        L.append("スナップショットが揃わない → 判定不能"); undecided = True
    L.append("| x | 量 | old 平均 | new 平均 | 相対差 % | old 変動 % | new 変動 % | 判定 |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for x in STATIONS:
        for k in KEYS:
            o, n = po[(x, k)], pn[(x, k)]
            mo, mn = o.mean(), n.mean()
            d = (mn - mo) / mo * 100
            vo = (o.max() - o.min()) / abs(mo) * 100; vn = (n.max() - n.min()) / abs(mn) * 100
            if vo > 0.1 or vn > 0.1:
                verdict = "判定不能 (変動 > 0.1 %)"; undecided = True
            else:
                good = abs(d) <= 0.1; ok &= good; verdict = "ok" if good else "NG"
            L.append(f"| {x} | {k} | {mo:.6e} | {mn:.6e} | {d:+.4f} | {vo:.4f} | {vn:.4f} | {verdict} |")
    v = "判定不能" if undecided else ("条件 (1)(2) 成立" if ok else "B (床比または物理差が上限超過)")
    L.append(f"\nVERDICT #5g (1)(2): {v}  — (3) の STEADY・物理ゲート・NaN/RISING と、旧版が起点床を外れるかは別出力で確認して A/B を決める")
    txt = "\n".join(L) + "\n"
    open(a.out, "w").write(txt); print(txt)


if __name__ == "__main__":
    main()
