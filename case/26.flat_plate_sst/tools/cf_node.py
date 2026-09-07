#!/usr/bin/env python3
"""node (median-dual) 平板 SST 結果の Cf を第一内点の線形勾配から出す (postprocess_wall_law.py は cell 用)。
τ_w = μ_w u_1/y_1 (y1+ ≲ 1 前提), Cf = τ_w/(½ρ∞U∞²), Schlichting 0.0592 Re_x^-0.2 と比較。
usage: cf_node.py RUN_DIR [x1 x2 ...]"""
import sys, glob, h5py, numpy as np
MU, GAM, CP = 1.8e-5, 1.4, 1004.5; R = CP*(GAM-1)/GAM; PT, TT, MACH = 100000.0, 288.15, 0.2
PS = PT*(1+0.5*(GAM-1)*MACH**2)**(-GAM/(GAM-1)); TS = TT/(1+0.5*(GAM-1)*MACH**2); RHO = PS/(R*TS); U = MACH*np.sqrt(GAM*R*TS)
d = sys.argv[1]; xs = [float(v) for v in sys.argv[2:]] or [0.3, 0.6, 0.9]
fn = sorted(glob.glob(d + "/res_[0-9]*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))[-1]
f = h5py.File(fn, "r"); c = np.array(f["MESH/COORD"]).reshape(-1, 3); Ux = np.array(f["VALUE/Ux"]); wd = np.array(f["VALUE/wall_dist"]); vl = np.array(f["VALUE/vis_lam"])
xr = np.round(c[:, 0], 7); ux = np.unique(xr)
out = []
for x in xs:
    xv = ux[np.argmin(np.abs(ux - x))]; col = np.where(xr == xv)[0]
    wall = col[wd[col] <= 0]
    if len(wall) == 0: continue
    yw = c[wall, 1].min(); inner = col[wd[col] > 0]; j = inner[np.argmin(wd[inner])]
    tau = vl[j] * Ux[j] / wd[j]; cf = tau / (0.5 * RHO * U * U); rex = RHO * U * xv / MU; cfs = 0.0592 * rex ** -0.2
    out.append((xv, wd[j], cf, cfs)); print(f"x={xv:.3f} y1={wd[j]*1e6:.2f}um Re_x={rex:.3e} Cf={cf:.5f} Schlichting={cfs:.5f} ratio={cf/cfs:.4f}")
print("file", fn.split("/")[-1])
