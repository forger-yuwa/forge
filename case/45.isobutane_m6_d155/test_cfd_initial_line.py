"""CFD ピン抽出器の単体試験 (plan tooling-nozzle-cfd-pinned-initial-line §6 V0 案): Hall 場を case/45 の Euler 格子に載せた合成場から抽出し、
HallThroat.throat_characteristic を再現できるか (M・θ の最大差、軸着地 x)。usage: design/.venv-opt/bin/python test_cfd_initial_line.py [格子の run]"""
import json, sys
from pathlib import Path
import numpy as np, h5py
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from cfd_initial_line import CFDThroatCharacteristic  # noqa: E402

run = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k")
info = json.loads((run / "prepare_info.json").read_text()); S = float(info["scale_m"]); ni = int(info["mesh"]["ni"])
with h5py.File(run / "nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
X = (nc[:, 0] / S).reshape(ni, -1); R = (nc[:, 1] / S).reshape(ni, -1)
ht = HallThroat(R=2.0, gamma=1.2735422978978068)
# 合成場: 壁は格子の壁 (≈ 設計壁)。Hall 場は R=2 の円弧前提なので壁近傍で厳密には壁に沿わないが、抽出器の検査にはそのまま使う
m = (X[:, 0] > -2.2) & (X[:, 0] < 2.7)
Mh = np.full(X.shape, 1.0); Th = np.zeros(X.shape)
Mh[m] = ht.mach(X[m], R[m]); Th[m] = ht.theta(X[m], R[m])
Mh[~m] = 1.0
cfd = CFDThroatCharacteristic(X, R, Mh, Th)
out = {}
for n in (41, 161):
    xc, rc, Mc, tc = cfd.throat_characteristic(n=n)
    xh, rh, Mhl, thl = ht.throat_characteristic(n=n)
    out[n] = dict(axis_x=[float(xc[0]), float(xh[0])], dM_max=float(np.abs(Mc - Mhl).max()), dtheta_max_deg=float(np.degrees(np.abs(tc - thl)[1:-1]).max()),
                  dx_max=float(np.abs(xc - xh).max()), wall_r=[float(rc[-1]), float(rh[-1])])
ok = all(v["dM_max"] <= 1e-4 and v["dtheta_max_deg"] <= 0.002 for v in out.values())
out["pass"] = bool(ok)
print(json.dumps(out, indent=1))
