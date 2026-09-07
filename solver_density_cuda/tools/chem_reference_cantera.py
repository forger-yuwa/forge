#!/usr/bin/env python3
"""tools/test_chemistry.cpp の参照解: 同じ機構 YAML (NASA-9 熱力学込み) で Cantera の定積断熱反応器を積分し
time,T,Y_* の CSV を出力する。
  .venv-chem/bin/python solver_density_cuda/tools/chem_reference_cantera.py MECH.yaml OUT.csv [T0=1200] [p_atm=1] [t_end=1e-3]"""
import sys, cantera as ct, numpy as np
mech, out = sys.argv[1], sys.argv[2]
T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1200.0
p0 = (float(sys.argv[4]) if len(sys.argv) > 4 else 1.0) * 101325.0
tEnd = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-3
g = ct.Solution(mech); g.TP = T0, p0
g.set_equivalence_ratio(1.0, "H2", "O2:1, N2:3.76")
r = ct.IdealGasReactor(g, clone=False); net = ct.ReactorNet([r])
net.rtol, net.atol = 1e-10, 1e-18
rows = [(0.0, r.T, *r.thermo.Y)]
t = 0.0
while t < tEnd:
    t = min(net.step(), tEnd) if net.time < tEnd else tEnd
    net.advance(t) if net.time < t else None
    rows.append((net.time, r.T, *r.thermo.Y))
    if net.time >= tEnd: break
with open(out, "w") as f:
    f.write("time,T," + ",".join("Y_" + s for s in g.species_names) + "\n")
    for row in rows: f.write(",".join(f"{v:.10e}" for v in row) + "\n")
T = np.array([r[1] for r in rows]); ts = np.array([r[0] for r in rows])
i = np.argmax(np.gradient(T, ts)); print(f"cantera: rho={g.density:.4f} tau_ign={ts[i]*1e6:.3f} us T_end={T[-1]:.2f} K rows={len(rows)}")
