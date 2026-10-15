#!/usr/bin/env python3
"""plan tooling-nozzle-card-local-mesh §6 V5c (事前登録 2026-10-15、ユーザ了承「よい」): 試験部の一様性の量を全場の出力ごとに測る。

窓 x 4.0〜7.95 m、線 r/r_w = 0.1 (判定) と r = 0 (軸の節点、併記)。各線で
  maxabs  = max |100 (M/5.5 − 1)|           (目標からの最大絶対偏差 [%]、許容 ≤ 0.1)
  halfamp = 100 (max M − min M) / (2·5.5)   (半振幅 [%]、許容 ≤ 0.05)
を測り、CSV (step ごと) に書く。線の値は `nozzle_report.eta_line` (列の中で r/r_w に補間し列の間を x で補間) と同じ。
両量の時系列は `check_quasisteady.py --series-csv` に渡す (停止規則の系列は波高を含まないので別に要る; codex diagnose 2026-10-15 ③)。
使い方: PYTHONPATH=<design> python v5c_series.py <run_dir> <out.csv>
"""
import re
import sys
from pathlib import Path

import numpy as np

from forge_design.report.nozzle_report import eta_line, load_field

MD, X0_M, X1_M = 5.5, 4.0, 7.95
run, out = Path(sys.argv[1]), Path(sys.argv[2])
files = sorted((f for f in run.iterdir() if re.fullmatch(r"res_\d+\.h5", f.name) and f.name != "res_0.h5"),
               key=lambda f: int(re.findall(r"\d+", f.name)[0]))
rows = []
for f in files:
    F = load_field(run, f.name)
    S = float(F["info"]["scale_m"])
    xq = np.linspace(X0_M / S, X1_M / S, 2001)
    rec = [int(re.findall(r"\d+", f.name)[0])]
    for eta in (0.1, 0.0):
        M = eta_line(F, "M", eta, xq)
        rec += [float(np.max(np.abs(100.0 * (M / MD - 1.0)))), float(100.0 * (M.max() - M.min()) / (2.0 * MD))]
    rows.append(rec)
    print(f"{f.name}: eta0.1 maxabs {rec[1]:.4f} % halfamp {rec[2]:.4f} % | axis maxabs {rec[3]:.4f} % halfamp {rec[4]:.4f} %", flush=True)
np.savetxt(out, np.array(rows), delimiter=",", comments="", fmt=["%d", "%.8g", "%.8g", "%.8g", "%.8g"],
           header="step,maxabs_eta0.1,halfamp_eta0.1,maxabs_axis,halfamp_axis")
last = rows[-1]
v = {"maxabs_eta0.1": last[1] <= 0.1, "halfamp_eta0.1": last[2] <= 0.05}
print(f"V5c (最後の出力 step {last[0]}): 最大絶対偏差 {last[1]:.4f} % (≤ 0.1: {'PASS' if v['maxabs_eta0.1'] else 'FAIL'})、"
      f"半振幅 {last[2]:.4f} % (≤ 0.05: {'PASS' if v['halfamp_eta0.1'] else 'FAIL'})。定常は check_quasisteady の VERDICT で別に判定")
