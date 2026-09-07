#!/usr/bin/env python3
# TGV の保存量 (運動エネルギー K, エントロピー S, 全エネルギー E, 質量 M) を
# res_*.h5 スナップショット時系列から計算し、初期値からの相対変化を出力/作図する。
import sys, glob, os, re
import numpy as np
import h5py

def _cell_volume(f, path):
    """セル体積: res には既定 (output.level 1) で volume が無いので、mesh h5 (solverConfig の meshFileName) の CELLS/volume を読む。"""
    import os, yaml
    if "VALUE/volume" in f:
        return f["VALUE/volume"][:].astype(np.float64)
    rd = os.path.dirname(os.path.abspath(path))
    mesh = yaml.safe_load(open(os.path.join(rd, "solverConfig.yaml")))["mesh"]["meshFileName"]
    with h5py.File(os.path.join(rd, mesh), "r") as m:
        return m["CELLS/volume"][:].astype(np.float64)


gamma = 1.4
cp = 0.4
cv = cp / gamma   # calorically perfect
dt = 0.007

def step_of(p):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(p))
    return int(m.group(1)) if m else -1

def totals(path):
    with h5py.File(path, "r") as f:
        V = _cell_volume(f, path)
        ro = f["VALUE/ro"][:].astype(np.float64)
        ux = f["VALUE/Ux"][:].astype(np.float64)
        uy = f["VALUE/Uy"][:].astype(np.float64)
        uz = f["VALUE/Uz"][:].astype(np.float64)
        P  = f["VALUE/P"][:].astype(np.float64)
        roe = f["VALUE/roe"][:].astype(np.float64)
    K = np.sum(0.5 * ro * (ux*ux + uy*uy + uz*uz) * V)
    # 比エントロピー s = cv ln(P/ro^gamma) (+const)。総エントロピー S = sum ro*s*V
    s = cv * np.log(P / np.power(ro, gamma))
    S = np.sum(ro * s * V)
    E = np.sum(roe * V)
    M = np.sum(ro * V)
    return K, S, E, M

def series(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "res_*.h5")), key=step_of)
    files = [p for p in files if step_of(p) >= 0 and "nan" not in p]
    rows = []
    for p in files:
        st = step_of(p)
        K, S, E, M = totals(p)
        rows.append((st, st*dt, K, S, E, M))
    return np.array(rows)

def report(name, run_dir):
    a = series(run_dir)
    if len(a) == 0:
        print(f"\n## {name}: no res files"); return None
    K0, S0, E0, M0 = a[0,2], a[0,3], a[0,4], a[0,5]
    print(f"\n## {name}  ({run_dir})")
    print(f"  {'step':>6} {'time':>8} {'K/K0':>12} {'(S-S0)/|S0|':>14} {'(E-E0)/E0':>13} {'(M-M0)/M0':>13}")
    for st, t, K, S, E, M in a:
        print(f"  {int(st):6d} {t:8.3f} {K/K0:12.6f} {(S-S0)/abs(S0):14.3e} {(E-E0)/E0:13.3e} {(M-M0)/M0:13.3e}")
    return a

if __name__ == "__main__":
    base = os.path.dirname(os.path.abspath(__file__))
    runs = [
        ("node  pure-KEEP", "run_0002_node_keep_explicit"),
        ("cell  pure-KEEP", "run_0001_cell_keep_explicit"),
        ("cell  KEEP+diss", "run_0003_cell_keep_diss"),
    ]
    data = {}
    for nm, d in runs:
        data[nm] = report(nm, os.path.join(base, d))

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13,5))
        for nm, a in data.items():
            if a is None or len(a)==0: continue
            t = a[:,1]; K=a[:,2]; S=a[:,3]
            ax[0].plot(t, K/K[0], "-o", ms=3, label=nm)
            ax[1].plot(t, (S-S[0])/abs(S[0]), "-o", ms=3, label=nm)
        ax[0].set_xlabel("time [s]"); ax[0].set_ylabel("K / K0"); ax[0].set_title("Total kinetic energy")
        ax[0].grid(alpha=0.3); ax[0].legend()
        ax[1].set_xlabel("time [s]"); ax[1].set_ylabel("(S - S0)/|S0|"); ax[1].set_title("Total entropy change")
        ax[1].grid(alpha=0.3); ax[1].legend()
        fig.tight_layout()
        out = os.path.join(base, "conservation_compare.png")
        fig.savefig(out, dpi=110)
        print(f"\nsaved {out}")
    except Exception as e:
        print("plot skipped:", e)
