"""CONTUR (積分法) の閉包を CFD 場の抽出 δ に較正する (plan verification-m6-axis-wave-mesh-su2 §5.1 #8b、オフライン)。
係数: k_f (C_f 倍率), k_N (N 倍率), a (Crocco 係数 0〜1)。目的: x∈[8,90] で ln(δ_E/δ_CONTUR) の最小二乗。
対象: B1 (run_0046) の場から方式 E で測った δ (_band_ab/B1_reextract_edge/delta_r_equiv.csv)。
usage: design/.venv-opt/bin/python calibrate_contur.py → _band_ab/contur_calibration.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
from forge_design.feedback.deltastar_integral import integral_bl
from forge_design.metrics.deltastar import pspline_uniform

C = Path(__file__).resolve().parent
p = load_problem(C / "problem_d155_ns_rt77p02.yaml"); d = design_chain(p); S = float(p.spec["r_throat"])
E = np.loadtxt(C / "_band_ab/B1_reextract_edge/delta_r_equiv.csv", delimiter=",", skiprows=1)
xe, de = E[:, 0], E[:, 1]
m_fit = (xe >= 8) & (xe <= 90) & np.isfinite(de) & (de > 0)
m_lo = (xe >= 1) & (xe < 8) & np.isfinite(de) & (de > 0)


def contur(kf, kN, a):
    r = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), S,
                    thermal_bc=p.wall_thermal_bc_integral, a_crocco=a, cf_scale=kf, n_scale=kN, x_out=xe[m_fit | m_lo])
    return np.interp(xe, r["x"], r["delta_r"])


t0 = time.time(); base = contur(1.0, 1.0, 1.0); t1 = time.time()
res_fun = lambda q: np.log(de[m_fit] / contur(*q)[m_fit])
sol = least_squares(res_fun, x0=[1.0, 1.0, 1.0], bounds=([0.5, 0.5, 0.0], [2.0, 2.0, 1.0]), diff_step=1e-3, xtol=1e-6)
cal = contur(*sol.x)


def report(model):
    lr = np.log(de / model)
    f = lr[m_fit]; xx = xe[m_fit]
    loc = f - pspline_uniform(xx, f, knot=10.0, lam=1.0)
    return dict(rms_fit_pct=float(100 * np.sqrt(np.mean(f ** 2))), mean_fit_pct=float(100 * f.mean()),
                local_detrended_max_pct=float(100 * np.abs(loc).max()), x_local=float(xx[np.argmax(np.abs(loc))]),
                lt8_max_pct=float(100 * np.abs(lr[m_lo]).max()), x_lt8=float(xe[m_lo][np.argmax(np.abs(lr[m_lo]))]),
                profile={f"{x}": float(100 * lr[np.argmin(abs(xe - x))]) for x in (2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 90)})


out = dict(coeffs=dict(k_f=float(sol.x[0]), k_N=float(sol.x[1]), a=float(sol.x[2])), nfev=int(sol.nfev), sec_per_eval=t1 - t0,
           baseline=report(base), calibrated=report(cal))
out["criteria"] = dict(rms_le_1pct=out["calibrated"]["rms_fit_pct"] <= 1.0, local_le_0p3pct=out["calibrated"]["local_detrended_max_pct"] <= 0.3,
                       lt8_le_5pct_report=out["calibrated"]["lt8_max_pct"] <= 5.0)
(C / "_band_ab/contur_calibration.json").write_text(json.dumps(out, indent=1))
np.savetxt(C / "_band_ab/contur_calibration_profile.csv", np.c_[xe, de, base, cal], delimiter=",", header="x_rt,delta_E,delta_contur,delta_contur_cal", comments="")
print(json.dumps(out, indent=1))
