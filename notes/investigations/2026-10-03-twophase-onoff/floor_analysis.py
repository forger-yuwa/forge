# condensation-two-phase-transport #4l: residual floor ratio L1/L0 for wet & dry, and where the L1 residual sits (limiter_P)
import pandas as pd, numpy as np, h5py, glob, os
C="/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
R={"wet_L1":"run_0522_floor_wet_L1","wet_L0":"run_0523_floor_wet_L0","dry_L1":"run_0524_floor_dry_L1","dry_L0":"run_0525_floor_dry_L0"}
cols=["rms_ro","rms_roUx","rms_roUy","rms_roe","rms_roK","rms_roOmega","rms_roY0","rms_roY1"]
fl={}
for k,r in R.items():
    d=pd.read_csv(C+r+"/residual_history.csv")
    d=d.groupby("step").last().reset_index()
    n=len(d); t=d.iloc[int(0.9*n):]
    fl[k]={c:float(np.median(t[c])) for c in cols if c in d}
    fl[k]["nstep"]=int(d.step.max())
F=pd.DataFrame(fl).T; print(F.to_string(float_format=lambda v:f"{v:.3e}"))
for w in ("wet","dry"):
    print(w, "R=L1/L0:", " ".join(f"{c[4:]} {fl[w+'_L1'][c]/fl[w+'_L0'][c]:.1f}" for c in cols if c in fl[w+"_L0"]))
for k in ("wet_L1","dry_L1","wet_L0","dry_L0"):
    fs=sorted(glob.glob(C+R[k]+"/res_*.h5"),key=lambda s:int(s.split("_")[-1][:-3]))
    a=h5py.File(fs[-1])["VALUE"]; b=h5py.File(fs[-2])["VALUE"]
    res=np.abs(np.array(a["res_ro"][...]).ravel()); lp=np.array(a["limiter_P"][...]).ravel(); lp0=np.array(b["limiter_P"][...]).ravel()
    g=np.array(a["g_0"][...]).ravel() if "g_0" in a else np.zeros_like(res)
    top=res>=np.quantile(res,0.99); act=lp<0.99; mov=np.abs(lp-lp0)>=0.05
    print(f"{k} last {os.path.basename(fs[-1])}: limiter<0.99 all {act.mean():.3f} top1% {act[top].mean():.3f} | moved>=0.05 all {mov.mean():.4f} top1% {mov[top].mean():.3f} | g>1e-6 all {(g>1e-6).mean():.3f} top1% {(g>1e-6)[top].mean():.3f} | top1% share of sum res^2 {np.sum(res[top]**2)/np.sum(res**2):.3f}")
