"""CFD ピン抽出器の単体試験 = plan tooling-nozzle-cfd-pinned-initial-line §6 V0 (事前登録):
合成 Hall 場 (case/45 の Euler 格子に Hall を載せる) から 線 max|ΔM| ≤ 2e-4・max|Δθ| ≤ 0.004°・|Δx₀| ≤ 2e-4、アンカー |ΔM| ≤ 2e-5・|ΔM′| ≤ 1e-4・|ΔM″| ≤ 1e-3、m* の再現 ≤ 1e-4 (相対)、ν↔M 往復 ≤ 2e-4。
usage: design/.venv-opt/bin/python test_cfd_initial_line.py [格子の run] → _band_ab/cfdpin_V0_test.json"""
import json, sys
from pathlib import Path
import numpy as np, h5py
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu, pm_mach  # noqa: E402
from forge_design.geometry.moc_inverse import _flux_along, _Pt  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas  # noqa: E402
from cfd_initial_line import CFDPinnedThroat  # noqa: E402

run = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k")
info = json.loads((run / "prepare_info.json").read_text()); S = float(info["scale_m"]); ni = int(info["mesh"]["ni"])
with h5py.File(run / "nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
X = (nc[:, 0] / S).reshape(ni, -1); R = (nc[:, 1] / S).reshape(ni, -1)
g = 1.2735422978978068; ht = HallThroat(R=2.0, gamma=g)
m = (X[:, 0] > -2.2) & (X[:, 0] < 2.7)
Mh = np.ones(X.shape); Th = np.zeros(X.shape); Mh[m] = ht.mach(X[m], R[m]); Th[m] = ht.theta(X[m], R[m])
cfd = CFDPinnedThroat(2.0, g, X, R, Mh, Th)
xc, rc, Mc, tc = cfd.throat_characteristic(n=41); xh, rh, Mhl, thl = ht.throat_characteristic(n=41)
aC = cfd.axis_anchor(cfd.x0_cfd); aH = ht.axis_anchor(cfd.x0_cfd)
p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); gas = _gam_or_gas(p)
gg = gas if hasattr(gas, "nu") else g
flux = lambda x, r, M, t: _flux_along([_Pt(float(a), float(b), float(d), float(pm_nu(float(c), gg)), gg) for a, b, c, d in zip(x, r, M, t)], gg)[-1]
mC, mH = flux(xc, rc, Mc, tc), flux(xh, rh, Mhl, thl)
Ms = np.linspace(1.0005, 1.3, 200); rt = np.array([float(pm_mach(float(pm_nu(float(v), gg)), gg)) for v in Ms])
out = dict(line_dM=float(np.abs(Mc - Mhl).max()), line_dtheta_deg=float(np.degrees(np.abs(tc - thl)).max()), dx0=float(xc[0] - xh[0]),
           anchor_dM=aC[0] - aH[0], anchor_dMp=aC[1] - aH[1], anchor_dMpp=aC[2] - aH[2], mstar_rel=float(mC / mH - 1), nu_roundtrip=float(np.abs(rt - Ms).max()))
gates = dict(line_dM=2e-4, line_dtheta_deg=0.004, dx0=2e-4, anchor_dM=2e-5, anchor_dMp=1e-4, anchor_dMpp=1e-3, mstar_rel=1e-4, nu_roundtrip=2e-4)
out["pass"] = {k: bool(abs(out[k]) <= v) for k, v in gates.items()}; out["all_pass"] = all(out["pass"].values())
(C / "_band_ab/cfdpin_V0_test.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
