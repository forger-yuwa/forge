#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B1 (ii): SERN 2D の力係数 5 列を、chi 明示 0 と省略 (auto 1) で差区間比較する。

規則は usage-rule plan (accepted `convection-slau-wall-normal-chi-usage-rule.md`) §4.2・§4.3 をそのまま使う (測る前に固定済み):
- 末尾 40 % の窓。各列を `check_quasisteady.classify` で判定 (drift/osc は列ごと: C_T 系 0.0002/0.0005、C_L 系・C_M 0.0014/0.0014)。
  DRIFTING / TRANSIENT-UNSETTLED / NONFINITE → 「判定保留」。STEADY / OSCILLATING の列だけ差区間で判定。
- 差区間 I_Δ = [min C1 − max C0, max C1 − min C0] (C0 = 明示 0、C1 = 省略)、許容帯 C_T 系・C_L 系 ±0.002、C_M ±0.05。
  I_Δ ⊂ 帯 → 帯内 / I_Δ ∩ 帯 = ∅ → 差が残る / 他 → 判定不能。全体は最悪列。

  python3 chi_default_band.py FLAG0.csv FLAG1.csv [--label m4_off] [--out FILE]
CSV は `v3sern_ct.py` (= 本番 runner と同じ `sern_forces.force_history`) の出力。
"""
import argparse
import csv
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solver_density_cuda", "tools"))
import check_quasisteady as QS  # noqa: E402

COLS = [("C_T", 0.0002, 0.0005, 0.002), ("C_T_with_shear", 0.0002, 0.0005, 0.002),
        ("C_L", 0.0014, 0.0014, 0.002), ("C_L_with_shear", 0.0014, 0.0014, 0.002), ("C_M", 0.0014, 0.0014, 0.05)]
TAIL = 0.4


def load(p):
    rows = list(csv.DictReader(open(p)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("flag0"); ap.add_argument("flag1"); ap.add_argument("--label", default=""); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    R0, R1 = load(a.flag0), load(a.flag1)
    L = [f"# chi 既定化 B1 (ii) 帯判定 {a.label}: 明示 0 = {a.flag0} / 省略 = {a.flag1}",
         f"末尾 {int(TAIL * 100)} % の窓、差区間 I_Δ = [min C1 − max C0, max C1 − min C0]"]
    worst = "帯内"
    rank = {"帯内": 0, "判定不能": 1, "判定保留": 1, "差が残る": 2}
    for c, dr, osc, tol in COLS:
        s0 = [float(r["step"]) for r in R0]; v0 = [float(r[c]) for r in R0]
        s1 = [float(r["step"]) for r in R1]; v1 = [float(r[c]) for r in R1]
        l0, d0, _ = QS.classify(s0, v0, TAIL, dr, osc, 4)
        l1, d1, _ = QS.classify(s1, v1, TAIL, dr, osc, 4)
        k0 = max(3, int(np.ceil(TAIL * len(v0)))); k1 = max(3, int(np.ceil(TAIL * len(v1))))
        t0, t1 = np.array(v0[-k0:]), np.array(v1[-k1:])
        lo, hi = t1.min() - t0.max(), t1.max() - t0.min()
        if any(l in ("DRIFTING", "TRANSIENT-UNSETTLED", "NONFINITE") for l in (l0, l1)):
            v = "判定保留"
        elif -tol <= lo and hi <= tol:
            v = "帯内"
        elif hi < -tol or lo > tol:
            v = "差が残る"
        else:
            v = "判定不能"
        if rank[v] > rank[worst]:
            worst = v
        L.append(f"{c:15s} flag0 {l0:12s} flag1 {l1:12s} 平均差 {t1.mean() - t0.mean():+.6f}  I_Δ [{lo:+.6f}, {hi:+.6f}]  許容 ±{tol}  → {v}")
    L.append(f"全体 (最悪列): {worst}")
    txt = "\n".join(L) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
