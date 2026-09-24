#!/usr/bin/env python3
"""#10b 補助: 各 run **単体**のランプ量の時系列 (差ではない)。揺れが差にあるか run にあるかの切り分け。"""
import argparse, sys
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v3sern_series import dumps, ramp_geom

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--out", required=True)
ap.add_argument("--ramp-id", type=int, default=4); ap.add_argument("--mesh", default=None)
a = ap.parse_args()
mesh = a.mesh or str(Path(a.run)/"sern.h5")
nodes, x, ipl, S, own = ramp_geom(mesh, a.ramp_id)
d = dumps(a.run); steps = sorted(s for s in d if s > 0)
rows = []
for s in steps:
    with h5py.File(d[s]) as g: P = np.asarray(g["/VALUE/P"], np.float64)
    p = P[nodes]; j = int(np.argmax(np.abs(np.gradient(p, x))))
    F = (P[own[ipl]][:,None] * S[ipl]).sum(axis=0)
    rows.append((s, F[0], F[1], x[j], float(np.sqrt(np.mean(p**2)))))
with open(a.out,"w") as f:
    f.write("step,Fx,Fy,x_foot,p_rms\n")
    for r in rows: f.write(f"{r[0]},{r[1]:.9e},{r[2]:.9e},{r[3]:.9e},{r[4]:.9e}\n")
print(f"# {a.out}: {len(rows)} 行  終端 Fx {rows[-1][1]:.6e}  Fy {rows[-1][2]:.6e}  x_foot {rows[-1][3]:.6f}")
