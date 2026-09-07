#!/usr/bin/env python3
"""run_0049 系: res_*.h5 の T(t) (内部/面/角ノード平均) を Cantera 参照 CSV と重ねて図にする。
  .venv-chem/bin/python plot_h2_ignition.py run_dir ref.csv [--dt 5e-9]"""
import argparse, glob, re, h5py, numpy as np, pathlib, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("ref"); ap.add_argument("--dt", type=float, default=5e-9); a = ap.parse_args()
d = pathlib.Path(a.run_dir)
files = sorted(glob.glob(str(d / "res_*.h5")), key=lambda f: int(re.findall(r"\d+", pathlib.Path(f).name)[0]))
with h5py.File(d / "cube.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
onface = (np.isclose(xyz, 0) | np.isclose(xyz, 1)).sum(axis=1)
sel = {"interior": onface == 0, "face": onface == 1, "corner": onface == 3}
t, T, Yh2o = [], {k: [] for k in sel}, []
for f in files:
    n = int(re.findall(r"\d+", pathlib.Path(f).name)[0])
    with h5py.File(f) as h:
        Tf = h["VALUE/T"][:]; t.append(n*a.dt*1e6)
        for k, m in sel.items(): T[k].append(Tf[m].mean())
        Yh2o.append(h["VALUE/Y5"][:][sel["interior"]].mean())
ref = np.loadtxt(a.ref, delimiter=",", skiprows=1)
t = np.array(t); Ti = np.array(T["interior"]); i = np.argmax(np.gradient(Ti, t)); ir = np.argmax(np.gradient(ref[:, 1], ref[:, 0]))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(ref[:, 0]*1e6, ref[:, 1], "k-", lw=2, label="Cantera (same YAML)")
for k, ls in (("interior", "o"), ("face", "s"), ("corner", "^")): ax[0].plot(t, T[k], ls, ms=4, mfc="none", label=f"forge node {k}")
ax[0].set_xlabel("t [us]"); ax[0].set_ylabel("T [K]"); ax[0].legend(); ax[0].grid(alpha=.3)
ax[0].set_title(f"stoich H2-air 1200 K 1 atm: tau_ign forge {t[i]:.1f} / Cantera {ref[ir,0]*1e6:.1f} us, T_end {Ti[-1]:.0f}/{ref[-1,1]:.0f} K")
ax[1].plot(ref[:, 0]*1e6, ref[:, [c for c in range(2, ref.shape[1])][5]], "k-", lw=2, label="Cantera Y_H2O")
ax[1].plot(t, Yh2o, "o", ms=4, mfc="none", label="forge Y_H2O (interior)"); ax[1].set_xlabel("t [us]"); ax[1].legend(); ax[1].grid(alpha=.3)
fig.tight_layout(); out = d / "h2_ignition_vs_cantera.png"; fig.savefig(out, dpi=130); print("wrote", out, f"tau forge {t[i]:.1f} us ref {ref[ir,0]*1e6:.2f} us T_end {Ti[-1]:.1f} ref {ref[-1,1]:.1f}")
