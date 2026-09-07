#!/usr/bin/env python3
# KEEP 陰解法CFL掃引の KE/エントロピー保存解析。各 run の res_*.h5 から K(t)=Σ½ρ|u|²V,
# S(t)=Σρ cv ln(P/ρ^γ)V を計算し、explicit ref と比較。dt は run ごとに solverConfig から読む。
import glob, os, re
import numpy as np, h5py, yaml
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _cell_volume(f, path):
    """セル体積: res には既定 (output.level 1) で volume が無いので、mesh h5 (solverConfig の meshFileName) の CELLS/volume を読む。"""
    import os, yaml
    if "VALUE/volume" in f:
        return f["VALUE/volume"][:].astype(np.float64)
    rd = os.path.dirname(os.path.abspath(path))
    mesh = yaml.safe_load(open(os.path.join(rd, "solverConfig.yaml")))["mesh"]["meshFileName"]
    with h5py.File(os.path.join(rd, mesh), "r") as m:
        return m["CELLS/volume"][:].astype(np.float64)


gamma = 1.4; cp = 0.4; cv = cp / gamma
def step_of(p):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(p)); return int(m.group(1)) if m else -1

def totals(path):
    with h5py.File(path, "r") as f:
        V = _cell_volume(f, path)
        ro = f["VALUE/ro"][:].astype(np.float64)
        P  = f["VALUE/P"][:].astype(np.float64)
        ux = f["VALUE/Ux"][:].astype(np.float64); uy=f["VALUE/Uy"][:].astype(np.float64); uz=f["VALUE/Uz"][:].astype(np.float64)
    K = np.sum(0.5*ro*(ux*ux+uy*uy+uz*uz)*V)
    S = np.sum(ro*(cv*np.log(P/np.power(ro, gamma)))*V)
    Pmin = float(P.min()); romin = float(ro.min())
    return K, S, Pmin, romin

def dt_of(run_dir):
    with open(os.path.join(run_dir, "solverConfig.yaml")) as fh:
        c = yaml.safe_load(fh)
    return float(c["time"]["deltaT"]["dt"])

def series(run_dir):
    dt = dt_of(run_dir)
    fs = sorted([p for p in glob.glob(os.path.join(run_dir,"res_*.h5")) if step_of(p)>=0 and "nan" not in p], key=step_of)
    t=[];K=[];S=[];Pm=[];rm=[]
    for p in fs:
        k,s,pmin,romin = totals(p); t.append(step_of(p)*dt); K.append(k); S.append(s); Pm.append(pmin); rm.append(romin)
    return np.array(t),np.array(K),np.array(S),np.array(Pm),np.array(rm)

RUNS = [
    ("explicit RK4 (CFL.05)",        "run_0011_cell_keep_expl_ref",   "k",  "-"),
    ("implicit dual-time (CFL.05)",  "run_0012_cell_keep_impl_cfl005","C0", "--"),
    ("implicit dual-time (CFL.2)",   "run_0013_cell_keep_impl_cfl02", "C1", "--"),
    ("implicit dual-time (CFL.5)",   "run_0014_cell_keep_impl_cfl05", "C2", "--"),
    ("implicit dual-time (CFL1)",    "run_0015_cell_keep_impl_cfl1",  "C3", "--"),
    ("implicit dual-time (CFL2)",    "run_0016_cell_keep_impl_cfl2",  "C4", "--"),
    ("implicit dual-time (CFL4)",    "run_0017_cell_keep_impl_cfl4",  "C5", "-."),
    ("implicit dual-time (CFL8)",    "run_0018_cell_keep_impl_cfl8",  "C6", "-."),
    ("implicit dual-time (CFL16)",   "run_0019_cell_keep_impl_cfl16", "C7", "-."),
]
base = os.path.dirname(os.path.abspath(__file__))
fig, ax = plt.subplots(1,2, figsize=(13,5))
print(f"{'run':30s}{'dt':>7}{'t_end':>7}{'K/K0_end':>11}{'(S-S0)/|S0|':>14}{'Pmin_end':>11}")
for nm,d,c,ls in RUNS:
    t,K,S,Pm,rm = series(os.path.join(base,d))
    if len(t)==0: print(f"{nm:30s} (no data)"); continue
    ax[0].plot(t, K/K[0], ls, color=c, lw=1.6, marker="o", ms=3, label=nm)
    ax[1].plot(t, (S-S[0])/abs(S[0]), ls, color=c, lw=1.6, marker="o", ms=3, label=nm)
    print(f"{nm:30s}{dt_of(os.path.join(base,d)):7.3f}{t[-1]:7.2f}{K[-1]/K[0]:11.4f}{(S[-1]-S[0])/abs(S[0]):14.3e}{Pm[-1]:11.4f}")
ax[0].axhline(1.0, ls=":", c="gray", lw=0.8); ax[0].set_xlabel("time"); ax[0].set_ylabel("K(t)/K(0)")
ax[0].set_title("Total kinetic energy (inviscid: ideal=const)"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
ax[1].axhline(0.0, ls=":", c="gray", lw=0.8); ax[1].set_xlabel("time"); ax[1].set_ylabel("(S(t)-S(0))/|S(0)|")
ax[1].set_title("Total entropy change"); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)
fig.suptitle("Taylor-Green 32^3 KEEP: explicit vs implicit dual-time, physical-CFL sweep")
fig.tight_layout(); out=os.path.join(base,"implicit_sweep_ke_entropy.png"); fig.savefig(out, dpi=120)
print("saved", out)
