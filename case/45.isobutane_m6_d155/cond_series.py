"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11 ④: 凝縮 ON run の凝縮 4 量の時系列 (準定常判定用)。
usage: python3 cond_series.py RUN → RUN/cond_series.csv (+ check_quasisteady --series-csv の VERDICT を表示)
量の定義は nozzle_report.metrics と同じ: 軸の凝縮開始位置 (g_0 > 1e-4 になる最初の x)、過飽和度 S の最大、出口コアの g_0 平均、出口コア M (η∈[0.05,0.7] 平均)。"""
import csv, subprocess, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field, eta_line  # noqa: E402
run = C / sys.argv[1]; rows = []
for res in _res_files(run):
    G = load_field(run, res); xF = float(G["X"][-1, 0]); xq = np.linspace(float(G["X"][0, 0]), xF, 2401)
    g0 = eta_line(G, "g_0", 0.0, xq); on = np.where(g0 > 1e-4)[0]
    eta = G["R"][-1] / G["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    rows.append(dict(step=int(res.split("_")[1].split(".")[0]), onset_x_axis=float(xq[on[0]]) if len(on) else float("nan"),
                     S_max=float(np.nanmax(G["V"]["condS_0"])) if "condS_0" in G["V"] else float("nan"),
                     exit_g_core=float(G["V"]["g_0"][-1][core].mean()), exit_core_M=float(G["V"]["M"][-1][core].mean())))
with open(run / "cond_series.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows[-5:]: print(r)
q = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/check_quasisteady.py"), str(run), "--series-csv", str(run / "cond_series.csv"),
                    "--series-cols", "onset_x_axis,S_max,exit_g_core,exit_core_M"], capture_output=True, text=True)
print("\n".join(l for l in q.stdout.splitlines() if "csv" in l or any(k in l for k in ("onset", "S_max", "exit_g", "exit_core"))))
