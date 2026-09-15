"""Wysłouzil 凝縮 run の報告量の時系列 CSV (check_quasisteady.py --series-csv 用): step, onset_c_g1e3 [mm] (中心線 g>1e-3), onset_c_g1e4, g_exit, T_exit, M_exit, g_max。
usage: python3 wys_cond_series.py RUN_DIR → RUN_DIR/cond_series.csv"""
import sys, glob, os, h5py, numpy as np
run=sys.argv[1]
res=sorted(glob.glob(os.path.join(run,"res_[0-9]*.h5")), key=lambda f:int(os.path.basename(f)[4:-3]))
rows=[]
for p in res:
    f=h5py.File(p); c=np.array(f["MESH/COORD"]).reshape(-1,3); V=f["VALUE"]
    x=c[:,0]*1000.0; y=c[:,1]
    g=V["g_0"][:]; T=V["T"][:]; M=np.hypot(V["Ux"][:],V["Uy"][:])/V["sonic"][:]
    cl=np.abs(y)<1e-9 if (np.abs(y)<1e-9).sum()>10 else (np.abs(y)<=np.sort(np.abs(y))[len(y)//200])
    xc=x[cl]; o=np.argsort(xc); xc=xc[o]; gc=g[cl][o]; Tc=T[cl][o]; Mc=M[cl][o]
    def onset(thr):
        i=np.where(gc>thr)[0]; return xc[i[0]] if len(i) else np.nan
    rows.append((int(os.path.basename(p)[4:-3]), onset(1e-3), onset(1e-4), gc[-1], Tc[-1], Mc[-1], g.max()))
arr=np.array(rows); fin=np.isfinite(arr).all(axis=1); first=int(np.argmax(fin)) if fin.any() else 0; arr=arr[first:]
out=os.path.join(run,"cond_series.csv"); np.savetxt(out, arr, delimiter=",", header="step,onset_c_g1e3_mm,onset_c_g1e4_mm,g_exit,T_exit,M_exit,g_max", comments="", fmt="%.8g"); print(out)
