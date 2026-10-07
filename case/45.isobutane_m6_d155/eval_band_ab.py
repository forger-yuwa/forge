"""帯修正の CFD 確認 A/B (plan verification-m6-axis-wave-mesh-su2 §4.4・§6) の判定量。
P(η) = max − min_{x∈[42,94]} ΔM(x,η)、ΔM = M/M_Euler(run_0037) − 1、b70(η)、各 res の時系列、B1/B0 の直接差。
usage: python3 eval_band_ab.py → _band_ab/cfd_ab.json
"""
import json, re
from pathlib import Path
import numpy as np, h5py
from axis_wave_compare import load, eta_profile, res_files, bump, XQ

C = Path(__file__).resolve().parent
EU = "run_0037_euler_rt77p02"
RUNS = {"B0": "run_0045_ns_band_adaptive", "B1": "run_0046_ns_band_edge", "A0'(旧壁)": "run_0042_ns_restart_ctrl"}
W = (XQ >= 42) & (XQ <= 94)
E = {eta: eta_profile(*load(EU), XQ, eta) for eta in (0.0, 0.1)}


def metrics(run, res=None):
    F = load(run, res); out = {}
    for eta in (0.0, 0.1):
        d = 100 * (eta_profile(*F, XQ, eta) / E[eta] - 1)
        out[f"P{eta}"] = float(d[W].max() - d[W].min())
        out[f"b70_{eta}"] = bump(d)[0]
        out[f"argmax{eta}"] = float(XQ[W][np.argmax(d[W])]); out[f"argmin{eta}"] = float(XQ[W][np.argmin(d[W])])
    return out


R = {}
for k, run in RUNS.items():
    ser = []
    for r in res_files(C / run):
        ser.append((int(re.findall(r"\d+", r)[0]), metrics(run, r)))
    fin = ser[-1][1]
    tail = ser[-5:]
    fin["tail5_range"] = {q: float(max(m[q] for _, m in tail) - min(m[q] for _, m in tail)) for q in ("P0.0", "P0.1", "b70_0.0")}
    with h5py.File(C / run / res_files(C / run)[-1]) as f:
        fin["nonfinite"] = int(sum((~np.isfinite(f["/VALUE/" + v][:])).sum() for v in ("ro", "P", "T", "Ux", "Uy", "roe")))
    fin["series"] = [(s, {q: round(m[q], 5) for q in ("P0.0", "P0.1", "b70_0.0")}) for s, m in ser]
    R[k] = fin
# 直接差 B1/B0
FA, FB = load(RUNS["B0"]), load(RUNS["B1"])
for eta in (0.0, 0.05, 0.1):
    D = 100 * (eta_profile(*FB, XQ, eta) / eta_profile(*FA, XQ, eta) - 1)
    R[f"D_B1_B0_eta{eta}"] = dict(median_42_94=float(np.median(D[W])), maxabs_dev_42_94=float(np.abs(D[W] - np.median(D[W])).max()),
                                  x_at=float(XQ[W][np.argmax(np.abs(D[W] - np.median(D[W])))]))
P1, P0 = R["B1"]["P0.1"], R["B0"]["P0.1"]
R["verdict"] = dict(effective=bool(P1 <= 0.06 and P1 <= 0.5 * P0 and R["B1"]["b70_0.0"] <= 0.08),
                    not_effective=bool(P1 >= 0.8 * P0))
(C / "_band_ab" / "cfd_ab.json").write_text(json.dumps(R, indent=1))
print(json.dumps({k: (v if not isinstance(v, dict) or "series" not in v else {q: v[q] for q in v if q != "series"}) for k, v in R.items()}, indent=1))
