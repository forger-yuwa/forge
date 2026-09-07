#!/usr/bin/env python3
"""forge 同梱の Jachimowski 機構 YAML を Cantera で読み、量論 H2-air の定積着火遅れを GRI 系 h2o2.yaml と比較する。
実行: .venv-chem/bin/python notes/investigations/scripts/chem_ignition_check.py  (repo root で)"""
import cantera as ct, numpy as np, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[3]
MECH = [str(ROOT / "solver_density_cuda/tools/mechanisms/h2air_jachimowski1988_13sp33r.yaml"),
        str(ROOT / "solver_density_cuda/tools/mechanisms/h2o2_jachimowski1988_9sp20r.yaml"),
        "h2o2.yaml"]

def tau_ign(mech, T0, p_atm=1.0, phi=1.0):
    g = ct.Solution(mech); g.TP = T0, p_atm * ct.one_atm
    g.set_equivalence_ratio(phi, "H2", "O2:1, N2:3.76")
    r = ct.IdealGasReactor(g, clone=False); net = ct.ReactorNet([r]); t = 0; Ts = []; ts = []
    while t < 5e-3:
        t = net.step(); ts.append(t); Ts.append(r.T)
        if r.T > T0 + 800: break
    ts = np.array(ts); T = np.array(Ts)
    return ts[np.argmax(np.gradient(T, ts))]

print("loaded reactions:", [ct.Solution(m).n_reactions for m in MECH])
print("ignition delay [us], stoich H2-air, const-V (max dT/dt)")
print(" T0[K]  p[atm]  Jach13   Jach9   Cantera-h2o2")
for T0, p in [(1000, 1), (1100, 1), (1200, 1), (1500, 1), (1200, 2), (2000, 1)]:
    r = [tau_ign(m, T0, p) * 1e6 for m in MECH]
    print(f" {T0:5d}  {p:4.1f}   {r[0]:7.1f} {r[1]:7.1f} {r[2]:7.1f}")
