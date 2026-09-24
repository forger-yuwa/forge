#!/usr/bin/env python3
"""#10b: SERN 2D の flag0/flag1 差の時系列を作る。plan convection-slau-wall-normal-chi §6 V3 / §6.2。

ノルムの定義は §6 V3 で 2026-09-23 に固定済み (壁圧を並べる前に決めたもの):
  x_foot  = flag0 の run で max|dp_w/dx| を与える位置
  窓      = ランプ上 x ∈ [x_foot − 2t, x_foot + 5t]   (t = カウル板厚 2 mm)
  foot_L2_pct = 窓内の p_w/p_ref の相対 L2 差 [%]
  dx_foot     = x_foot(flag1) − x_foot(flag0)
  dFx_pct/dFy_pct = ランプ面の圧力積分 F = Σ P_f S_f の差 [%]

usage: v3sern_series.py RUN0 RUN1 --out-prefix PFX [--ramp-id 4] [--t-cowl 0.002]
"""
import argparse, glob, re, sys
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / "cad"))
from diag_wall_cv_budget import parse_struct

def dumps(d):
    out = {}
    for f in glob.glob(str(Path(d)/"res_[0-9]*.h5")):
        out[int(re.findall(r"res_(\d+)\.h5", Path(f).name)[0])] = f
    return out

def ramp_geom(mesh, rid):
    with h5py.File(mesh) as f:
        ic = np.asarray(f["BCONDS"][str(rid)]["iCells"]).ravel()
        ip = np.asarray(f["BCONDS"][str(rid)]["iPlanes"]).ravel()
        cc = np.asarray(f["CELLS/centCoords"], np.float64).reshape(-1,3)
        S  = np.asarray(f["PLANES/surfVect"], np.float64).reshape(-1,3)
        st = np.asarray(f["PLANES/STRUCT"])
    own, nei, _ = parse_struct(st, len(S))
    ic = ic[ic >= 0]
    order = np.argsort(cc[ic,0])
    return ic[order], cc[ic[order],0], ip[ip>=0], S, own.astype(int)

def foot_and_force(P, nodes, x, ipl, S, own):
    p = P[nodes]
    dpdx = np.gradient(p, x)
    j = int(np.argmax(np.abs(dpdx)))
    F = (P[own[ipl]][:,None] * S[ipl]).sum(axis=0)     # Σ P_f S_f (1 次: owner ノードの P)
    return x[j], p, F

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run0"); ap.add_argument("run1")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ramp-id", type=int, default=4)
    ap.add_argument("--t-cowl", type=float, default=0.002)
    ap.add_argument("--mesh", default=None)
    a = ap.parse_args()
    mesh = a.mesh or str(Path(a.run0)/"sern.h5")
    nodes, x, ipl, S, own = ramp_geom(mesh, a.ramp_id)
    d0, d1 = dumps(a.run0), dumps(a.run1)
    steps = sorted(set(d0) & set(d1)); steps = [s for s in steps if s > 0]
    print(f"# ramp ノード {len(nodes)}  面 {len(ipl)}  共通 dump {len(steps)} (step {steps[0]}..{steps[-1]})")
    rows_f, rows_F = [], []
    for s in steps:
        with h5py.File(d0[s]) as g: P0 = np.asarray(g["/VALUE/P"], np.float64)
        with h5py.File(d1[s]) as g: P1 = np.asarray(g["/VALUE/P"], np.float64)
        xf0, p0, F0 = foot_and_force(P0, nodes, x, ipl, S, own)
        xf1, p1, F1 = foot_and_force(P1, nodes, x, ipl, S, own)
        w = (x >= xf0 - 2*a.t_cowl) & (x <= xf0 + 5*a.t_cowl)
        l2 = 100.0*np.sqrt(np.mean((p1[w]-p0[w])**2))/np.sqrt(np.mean(p0[w]**2))
        rows_f.append((s, l2, xf1-xf0))
        rows_F.append((s, 100*(F1[0]-F0[0])/abs(F0[0]), 100*(F1[1]-F0[1])/abs(F0[1])))
    with open(a.out_prefix+"_foot.csv","w") as f:
        f.write("step,foot_L2_pct,dx_foot\n")
        for r in rows_f: f.write(f"{r[0]},{r[1]:.6f},{r[2]:.9f}\n")
    with open(a.out_prefix+"_force.csv","w") as f:
        f.write("step,dFx_pct,dFy_pct\n")
        for r in rows_F: f.write(f"{r[0]},{r[1]:.6f},{r[2]:.6f}\n")
    print(f"# 書き出し: {a.out_prefix}_foot.csv / {a.out_prefix}_force.csv")
    print(f"# 終端: foot_L2 {rows_f[-1][1]:.4f} %  dx_foot {rows_f[-1][2]:.3e}  dFx {rows_F[-1][1]:+.4f} %  dFy {rows_F[-1][2]:+.4f} %")
