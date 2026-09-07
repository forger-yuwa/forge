#!/usr/bin/env python3
"""0-D 着火検証 (有限速度化学) の IC: 量論 H2-air, T0, p0, u=0 を node 周期箱 (cube.h5) に書く。
roe は forge と同じ NASA-9 + sensible datum (thermoHrefTemp=298.15) で組む。
  .venv-chem/bin/python make_h2_ignition_ic.py run_dir [--T0 1200] [--p 101325]"""
import argparse, h5py, numpy as np, yaml, pathlib
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--T0", type=float, default=1200.0)
ap.add_argument("--p", type=float, default=101325.0); ap.add_argument("--Tref", type=float, default=298.15); a = ap.parse_args()
d = pathlib.Path(a.run_dir); Ru = 8.314462618
names = ["H2", "O2", "H", "O", "OH", "H2O", "HO2", "H2O2", "N2"]   # solverConfig.yaml の species 順
db = yaml.safe_load(open(d / "species_db.yaml"))
def h_molar(sp, T):
    c = sp["nasa9_low"] if T < sp["Tmid"] else sp["nasa9_high"]
    return Ru*T*(-c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T)
X = np.zeros(9); X[0] = 2; X[1] = 1; X[8] = 3.76; X /= X.sum()
W = np.array([db[n]["MW"] for n in names]); Wmix = (X*W).sum(); Y = X*W/Wmix
ro = a.p*Wmix/(Ru*a.T0); e = sum(Y[i]*(h_molar(db[n], a.T0) - h_molar(db[n], a.Tref))/W[i] for i, n in enumerate(names)) - (Ru/Wmix)*a.T0
with h5py.File(d / "cube.h5", "r+") as f:
    v = f["VALUE"]; N = v["ro"].shape[0]
    v["ro"][:] = ro; v["roUx"][:] = 0; v["roUy"][:] = 0; v["roUz"][:] = 0; v["roe"][:] = ro*e
    for i in range(9):
        k = f"roY{i}"
        if k in v: v[k][:] = ro*Y[i]
        else: v.create_dataset(k, data=np.full(N, ro*Y[i], dtype=v["ro"].dtype))
print(f"IC: T0={a.T0} p={a.p} ro={ro:.5f} e={e:.2f} J/kg Y={np.round(Y, 5)} N={N}")
