"""CFD ピンの凍結入力の時間安定性 (codex plan レビュー M2 採用): 同じ run の複数スナップショットから線・アンカー・m* を抽出し、ばらつきを測る。
usage: python3 cfdpin_line_stability.py RUN res_A.h5 res_B.h5 ... → _band_ab/cfdpin_line_stability_<run>.json"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from cfd_initial_line import load_run_field, CFDPinnedThroat  # noqa: E402
from forge_design.evaluate.runner_axismach import load_problem, _gam_or_gas  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu  # noqa: E402
from forge_design.geometry.moc_inverse import _flux_along, _Pt  # noqa: E402

run = C / sys.argv[1]; resl = sys.argv[2:]
p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); gas = _gam_or_gas(p); gg = gas if hasattr(gas, "nu") else 1.2735422978978068
L = {}
for r in resl:
    X, R, M, T, _ = load_run_field(run, [r]); t = CFDPinnedThroat(2.0, 1.2735422978978068, X, R, M, T)
    x, rr, Mm, th = t.throat_characteristic(n=161)
    ms = _flux_along([_Pt(float(a), float(b), float(d), float(pm_nu(float(c), gg)), gg) for a, b, c, d in zip(x, rr, Mm, th)], gg)[-1]
    L[r] = dict(x=x, M=Mm, th=th, x0=t.x0_cfd, anchor=t.axis_anchor(t.x0_cfd), mstar=float(ms))
ref = L[resl[-1]]
out = {"run": sys.argv[1], "res": resl,
       "line_dM_max": float(max(np.abs(v["M"] - ref["M"]).max() for v in L.values())),
       "line_dtheta_max_deg": float(max(np.degrees(np.abs(v["th"] - ref["th"])).max() for v in L.values())),
       "line_dx_max": float(max(np.abs(v["x"] - ref["x"]).max() for v in L.values())),
       "x0_range": float(np.ptp([v["x0"] for v in L.values()])),
       "M_A_range": float(np.ptp([v["anchor"][0] for v in L.values()])), "Mp_A_range": float(np.ptp([v["anchor"][1] for v in L.values()])),
       "mstar_rel_range": float(np.ptp([v["mstar"] for v in L.values()]) / ref["mstar"])}
(C / f"_band_ab/cfdpin_line_stability_{sys.argv[1][:8]}.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
