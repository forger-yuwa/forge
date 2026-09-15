"""凝縮固有量の時系列 CSV (check_quasisteady.py --series-csv 用): 全 res_*.h5 から step, g_max, exit_g_avg (質量流束平均, x_max−2 r_t),
onset_x_wall (壁流線 g>1e-3 Y_w), x30_j4 (壁から 4 番目の線で g=30 % Y_w 到達), wall_M_exit, exit_M_avg。
usage: python3 cond_series_csv.py RUN_DIR [--Yw 0.0377] → RUN_DIR/cond_series.csv"""
import sys, glob, os, argparse, h5py, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--Yw", type=float, default=0.0377); a=ap.parse_args()
mesh=h5py.File(os.path.join(a.run,"nozzle.h5")); nc=mesh["MESH/COORD"][:].reshape(-1,3); S=0.205711
x=nc[:,0]/S; r=nc[:,1]/S
res=sorted(glob.glob(os.path.join(a.run,"res_[0-9]*.h5")), key=lambda f:int(os.path.basename(f)[4:-3]))
f0=h5py.File(res[-1]); n=len(f0["VALUE/ro"]); node=(n==len(x))
if node:
    cols=np.unique(nc[:,0]); xc=cols[np.argmin(abs(cols-(nc[:,0].max()-2*S)))]; ex=np.abs(nc[:,0]-xc)<1e-9
    lines={}
    for xv in cols:
        idx=np.where(np.abs(nc[:,0]-xv)<1e-9)[0]; idx=idx[np.argsort(r[idx])]
        for jj in (1,4): lines.setdefault(jj,[]).append(idx[-jj])
    wall=np.array(lines[1]); j4=np.array(lines[4]); xl=cols/S; rr=r
else:
    cc=mesh["CELLS/centCoords"][:].reshape(-1,3); xcell=cc[:,0]/S; rcell=cc[:,1]/S
    # cell: 出口断面 = x_max−2 r_t に最も近い列 (x で丸め), 壁線 = 各 x 列で r 最大, j4 = 4 番目
    xr=np.round(xcell,4); cols=np.unique(xr); xc=cols[np.argmin(abs(cols-(xcell.max()-2)))]; ex=(xr==xc)
    lines={}
    for xv in cols:
        idx=np.where(xr==xv)[0]; idx=idx[np.argsort(rcell[idx])]
        if len(idx)>=4:
            for jj in (1,4): lines.setdefault(jj,[]).append(idx[-jj])
    wall=np.array(lines[1]); j4=np.array(lines[4]); xl=np.array([xcell[i] for i in wall]); rr=rcell
rows=[]
for p in res:
    f=h5py.File(p)["VALUE"]; step=int(os.path.basename(p)[4:-3])
    g=f["g_0"][:]; ro=f["ro"][:]; ux=f["Ux"][:]; M=np.hypot(ux,f["Uy"][:])/f["sonic"][:]
    w=(ro*ux*rr)[ex]; o=np.argsort(rr[ex]); ga=np.trapz((w*g[ex])[o],rr[ex][o])/np.trapz(w[o],rr[ex][o]); Ma=np.trapz((w*M[ex])[o],rr[ex][o])/np.trapz(w[o],rr[ex][o])
    i=np.where(g[wall]>1e-3*a.Yw)[0]; xo=xl[i[0]] if len(i) else np.nan
    k=np.where(g[j4]>0.3*a.Yw)[0]; x30=xl[k[0]] if len(k) else np.nan
    rows.append((step, g.max(), ga, xo, x30, M[wall][-1], Ma))
arr=np.array(rows)
# onset 前のスナップショットは onset 列が NaN になるので、全列が有限になった最初の行以降だけを出す (末尾窓の判定には影響しない)
fin=np.isfinite(arr).all(axis=1); first=int(np.argmax(fin)) if fin.any() else 0
arr=arr[first:] if fin.any() else arr
out=os.path.join(a.run,"cond_series.csv")
np.savetxt(out, arr, delimiter=",", header="step,g_max,exit_g_avg,onset_x_wall,x30_j4,wall_M_exit,exit_M_avg", comments="", fmt="%.8g")
print(out); print(open(out).read())
