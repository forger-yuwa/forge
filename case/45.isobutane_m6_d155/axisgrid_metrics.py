"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11h (事前登録): 軸側上限 A/B の量の時系列。
usage: python3 axisgrid_metrics.py RUN → RUN/axisgrid_series.csv
各スナップショットで η=0・0.05・0.1 の Mach 波 (P-spline knot 10 の残差の最大絶対値、試験窓 [x_E+2, x_F−1]、全長 2401 点)・オーバーシュート (x ≥ x_E−15 の最大)・傾き
(試験窓の 1 次当てはめ × 窓幅) と最大位置。値は 100(M/6−1) [%]。nozzle_report.metrics と同じ定義。"""
import csv, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field, eta_line, pspline  # noqa: E402
run = C / sys.argv[1]; rows = []
for res in _res_files(run):
    G = load_field(run, res); info = G["info"]; Md = float(info.get("Md", 6.0)); xE = float(info.get("x_E", 40.0)); xF = float(G["X"][-1, 0])
    xq = np.linspace(float(G["X"][0, 0]), xF, 2401); wt = (xq >= xE + 2) & (xq <= xF - 1); wo = (xq >= xE - 15) & (xq <= xF)
    r = dict(step=int(res.split("_")[1].split(".")[0]))
    for eta in (0.0, 0.05, 0.1):
        d = 100 * (eta_line(G, "M", eta, xq) / Md - 1); xx, v = xq[wt], d[wt]; res_ = np.abs(v - pspline(xx, v)); i = int(np.argmax(res_))
        tag = f"e{int(round(eta * 100)):02d}"
        r.update({f"wave_{tag}": float(res_[i]), f"xwave_{tag}": float(xx[i]), f"os_{tag}": float(d[wo].max()),
                  f"xos_{tag}": float(xq[wo][np.argmax(d[wo])]), f"slope_{tag}": float(np.polyfit(xx, v, 1)[0] * (xx[-1] - xx[0]))})
    rows.append(r)
with open(run / "axisgrid_series.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(run.name, len(rows), "rows")
