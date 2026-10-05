"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11g (codex diagnose 2026-10-06 wave-drifting, 事前登録): 波 η0.1 の最大値探索の標本間隔 A/B (CFD 0 step)。
usage (AWS, case dir): python3 wave_sampling_ab.py RUN [RUN ...] → <RUN>/wave_sampling_series.csv、_band_ab/wave_sampling_ab.json
腕 A: nozzle_report.metrics と同じ (全長 2401 点、試験窓 [x_E+2, x_F−1]、P-spline knot 10・lam 1 を A の標本で当てはめ、残差の最大絶対値)。
腕 B: 同じ時刻・同じ P-spline 曲線 (A の標本で当てはめた係数を再利用、B では再フィットしない) を、試験窓の標本間隔を 4 分割した点で評価して最大値を探す。"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field, eta_line  # noqa: E402


def pspline_obj(x, v, knot=10.0, lam=1.0, k=3, order=2):   # metrics.deltastar.pspline_uniform と同じ当てはめ、曲線オブジェクトを返す
    lo, hi = float(x.min()), float(x.max()); n_int = max(int(np.ceil((hi - lo) / knot)) - 1, 1)
    t = np.concatenate([[lo] * (k + 1), np.linspace(lo, hi, n_int + 2)[1:-1], [hi] * (k + 1)]); n = len(t) - k - 1
    B = BSpline.design_matrix(x, t, k).toarray(); D = np.diff(np.eye(n), n=order, axis=0)
    return BSpline(t, np.linalg.solve(B.T @ B + lam * (D.T @ D), B.T @ v), k, extrapolate=True)


out = {}
for name in sys.argv[1:]:
    run = C / name; rows = []
    for res in _res_files(run):
        G = load_field(run, res); info = G["info"]; Md = float(info.get("Md", 6.0)); xE = float(info.get("x_E", 40.0)); xF = float(G["X"][-1, 0])
        xq = np.linspace(float(G["X"][0, 0]), xF, 2401); w = (xq >= xE + 2) & (xq <= xF - 1)
        xa = xq[w]; va = 100 * (eta_line(G, "M", 0.1, xa) / Md - 1); sp = pspline_obj(xa, va)
        ra = np.abs(va - sp(xa)); ia = int(np.argmax(ra))
        xb = np.linspace(xa[0], xa[-1], 4 * (len(xa) - 1) + 1); vb = 100 * (eta_line(G, "M", 0.1, xb) / Md - 1)
        rb = np.abs(vb - sp(xb)); ib = int(np.argmax(rb))
        rows.append(dict(step=int(res.split("_")[1].split(".")[0]), wave_A=float(ra[ia]), x_A=float(xa[ia]), wave_B=float(rb[ib]), x_B=float(xb[ib]),
                         B_minus_A=float(rb[ib] - ra[ia])))
    with open(run / "wave_sampling_series.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
    t = rows[-5:]; st = np.array([r["step"] for r in t], float)
    out[name] = dict(tail5=t, max_abs_BmA_all=float(max(abs(r["B_minus_A"]) for r in rows)),
                     trend_A=float(np.polyfit(st, [r["wave_A"] for r in t], 1)[0] * (st[-1] - st[0])),
                     trend_B=float(np.polyfit(st, [r["wave_B"] for r in t], 1)[0] * (st[-1] - st[0])))
    print(name, json.dumps({k: v for k, v in out[name].items() if k != "tail5"}))
    for r in t: print("  ", r)
(C / "_band_ab" / "wave_sampling_ab.json").write_text(json.dumps(out, indent=1))
