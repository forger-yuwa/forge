"""CFD ピン §6 V2 (凍結の自己整合): 凍結元 (V0 場) と pass 1 の Euler 場から抽出した線・アンカー・m* の差 D₁ を、D₀ (V0 場の線 − Hall の線) と比べる。
登録: D₁ max|ΔM| ≤ 5e-4・max|Δθ| ≤ 0.01°・|Δm*/m*| ≤ 2e-4 → 採用、超過は諮問、D₁ > D₀ (各許容差で無次元化した最大) は棄却。
usage: python3 cfdpin_v2.py FROZEN_RUN FROZEN_RES PASS1_RUN PASS1_RES → _band_ab/cfdpin_v2.json"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from cfd_initial_line import load_run_field, CFDPinnedThroat  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.evaluate.runner_axismach import load_problem, _gam_or_gas  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu  # noqa: E402
from forge_design.geometry.moc_inverse import _flux_along, _Pt  # noqa: E402
G = 1.2735422978978068
gas = _gam_or_gas(load_problem(C / "problem_d155_euler_c2final_n2400.yaml")); gg = gas if hasattr(gas, "nu") else G


def line(t):
    x, r, M, th = t.throat_characteristic(n=161)
    ms = _flux_along([_Pt(float(a), float(b), float(d), float(pm_nu(float(c), gg)), gg) for a, b, c, d in zip(x, r, M, th)], gg)[-1]
    return dict(x=x, M=M, th=th, mstar=float(ms))


L0 = line(CFDPinnedThroat(2.0, G, *load_run_field(C / sys.argv[1], [sys.argv[2]])[:4]))
L1 = line(CFDPinnedThroat(2.0, G, *load_run_field(C / sys.argv[3], [sys.argv[4]])[:4]))
LH = line(HallThroat(R=2.0, gamma=G))
tol = dict(M=5e-4, th=0.01, mstar=2e-4)


def D(a, b):
    return dict(M=float(np.abs(a["M"] - b["M"]).max()), th=float(np.degrees(np.abs(a["th"] - b["th"])).max()),
                mstar=float(abs(a["mstar"] / b["mstar"] - 1)), x=float(np.abs(a["x"] - b["x"]).max()))


D1, D0 = D(L1, L0), D(L0, LH)
n1 = max(D1[k] / tol[k] for k in tol); n0 = max(D0[k] / tol[k] for k in tol)
out = dict(D1=D1, D0=D0, D1_norm=n1, D0_norm=n0, adopt=bool(n1 <= 1.0), reject=bool(n1 > n0))
(C / "_band_ab/cfdpin_v2.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
