"""case/70: run のスナップショットごとの出口コア M と無次元の流量の系列を作り、check_quasisteady --series-csv で判定する
(plan tooling-nozzle-deltastar-band-root-zone §6 の参照場の受入れゲート、codex 諮問 2026-10-13 ns-rebuild-route)。

- 出口コア M: 最後の列の η = r/r_w ∈ [0.05, 0.7] の M の平均 (nozzle_report の exit_core と同じ定義)。
- 流量 (無次元、r_t 単位): 2π∫ρU_x r dr を x ∈ (−2.5, x_E) の列で計算した中央値 (nozzle_report の mdot_ratio_vs_euler の分子・分母と同じ定義)。
判定は末尾 5 枚 (--tail 4.5/n) で check_quasisteady の既定の許容 (drift 5 %・osc 10 %)。**出口 M の合否 (plan tooling-nozzle-deltastar-band-root-zone §6 V5) には緩すぎる**
(codex result レビュー 2026-10-13 M2) ので、その判定は `--tail 0.4 --drift 5e-5 --osc 5e-5` (M ≈ 5 で絶対 2.5e-4) を渡す。出力 <run>/quantities_series.csv・<run>/quasisteady_series.txt。
usage (case dir、design/.venv-opt):  python qs_series.py RUN [--every K] [--tail F --drift F --osc F]   (K 枚ごとに間引く。既定 1)
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--every", type=int, default=1)
ap.add_argument("--tail", type=float, default=None); ap.add_argument("--drift", type=float, default=None); ap.add_argument("--osc", type=float, default=None)
a = ap.parse_args()
run = (C / a.run).resolve()
files = _res_files(run)[::a.every]
rows = []
for r in files:
    F = load_field(run, r)
    eta = F["R"][-1] / F["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    xs = F["X"][:, 0]; xE = float(F["info"].get("x_E", 40.0))
    q = F["V"]["ro"] * F["V"]["Ux"]
    m = 2 * np.pi * np.trapezoid(q * F["R"], F["R"], axis=1)
    sel = (xs > -2.5) & (xs < xE)
    step = int(''.join(ch for ch in Path(str(r)).stem.split("_")[-1] if ch.isdigit()))
    rows.append({"step": step, "exit_core_M": float(F["V"]["M"][-1][core].mean()), "mdot_rt": float(np.median(m[sel]))})
rows.sort(key=lambda d: d["step"])
out = run / "quantities_series.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["step", "exit_core_M", "mdot_rt"]); w.writeheader(); w.writerows(rows)
n = len(rows)
if n < 6:
    raise SystemExit(f"{out}: スナップショットが {n} 枚 (末尾 5 枚の判定には 6 枚以上)")
q = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(out),
                    "--series-cols", "exit_core_M,mdot_rt", "--tail", f"{(a.tail if a.tail is not None else 4.5 / n):.6f}"]
                   + (["--drift", str(a.drift)] if a.drift is not None else []) + (["--osc", str(a.osc)] if a.osc is not None else []),
                   capture_output=True, text=True)
txt = q.stdout + q.stderr
(run / "quasisteady_series.txt").write_text(txt)
tail = rows[-5:]
print(txt.strip().splitlines()[-6:] if txt.strip() else txt)
print({"n": n, "tail5_exit_core_M": [round(t["exit_core_M"], 6) for t in tail], "tail5_mdot_rt": [round(t["mdot_rt"], 6) for t in tail]})
