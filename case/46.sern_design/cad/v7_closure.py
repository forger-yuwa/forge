#!/usr/bin/env python3
"""等温ピン拘束の閉包試験: r-1 == -(ΔR/R)/(ΔP/P) か。

壁ノードは T=Tw, P=ρ R(Y) Tw にピンされる (nodeWallDirichlet_d.cu:84-86)。
よって ΔP/P = Δρ/ρ + ΔR/R、すなわち r := (Δρ/ρ)/(ΔP/P) は
r - 1 = -(ΔR/R)/(ΔP/P) を厳密に満たすはず。float32 丸めで閉じるかを見る。
"""
import sys, h5py, numpy as np, yaml
sys.path.insert(0,".")
from v7_dist import load, wall_mask
MESH="run_0447_cfl02_ext/sern.h5"
E02="run_0447_cfl02_ext/res_10000.h5"; E04="run_0448_cfl04_ext/res_10000.h5"
e,f = load(E02), load(E04)
wm = wall_mask(MESH, len(e["ro"]))
db = yaml.safe_load(open("run_0447_cfl02_ext/species_db.yaml"))
names = [k for k in db.keys()]
MW = [float(db[k]["MW"]) for k in names]
print(f"species {names}  MW {MW}")
Ru = 8.314462618
R = [Ru/m for m in MW]
def y(fn, ro):
    with h5py.File(fn) as g:
        return [np.asarray(g[f"/VALUE/roY{i}"],np.float64)[wm]/ro for i in range(len(MW))]
ye, yf = y(E02, e["ro"][wm]), y(E04, f["ro"][wm])
Re = sum(yi*Ri for yi,Ri in zip(ye,R)); Rf = sum(yi*Ri for yi,Ri in zip(yf,R))
dR = (Rf-Re)/Re
dP = (f["P"][wm]-e["P"][wm])/e["P"][wm]
dr = (f["ro"][wm]-e["ro"][wm])/e["ro"][wm]
print(f"R: {R[0]:.2f} / {R[1]:.2f} J/kgK   ΔR/R の rms {np.sqrt(np.mean(dR**2)):.3e}")
sel = np.abs(dP) > np.percentile(np.abs(dP), 90)
lhs = dr[sel]/dP[sel] - 1.0
rhs = -dR[sel]/dP[sel]
res = np.abs(lhs - rhs)
print(f"|ΔP/P| 上位 10% ({int(sel.sum())} 節点):")
print(f"  |(r-1) + (ΔR/R)/(ΔP/P)| 中央値 {np.median(res):.3e}  95%点 {np.percentile(res,95):.3e}  最大 {res.max():.3e}")
print(f"  参考: |r-1| 中央値 {np.median(np.abs(lhs)):.3e}  -> 拘束式で説明できた割合 "
      f"{100*(1-np.median(res)/np.median(np.abs(lhs))):.2f} %")
