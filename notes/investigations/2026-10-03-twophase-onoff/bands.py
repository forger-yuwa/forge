# condensation-two-phase-transport #4j: ON(run_0521)-OFF(run_0520) differences by wall-distance band and x region (for the diagnose brief)
import h5py, numpy as np
C="/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
def ld(r):
    f=h5py.File(C+r+"/res_16000.h5"); V=f["VALUE"]
    d={k:np.array(V[k][...]).ravel().astype(float) for k in ["T","P","g_0","Y0","Y1","wall_dist","volume","vis_turb","vis_lam","condS_0","Ux","ro"]}
    xyz=np.array(f["MESH"]["COORD"][...]).reshape(-1,3)
    return d,xyz,list(f["MESH"].keys())
a,xa,mk=ld("run_0520_twophase_off_16k"); b,xb,_=ld("run_0521_twophase_dplur_nn0_cont")
print("MESH keys",mk, "n",a["T"].size)
print("vis_turb max/ratio", a["vis_turb"].max(), (a["vis_turb"]/a["vis_lam"]).max(), b["vis_turb"].max())
x=xa[:,0]*1e3 if xa is not None and xa.shape[0]==a["T"].size else None
y=xa[:,1]*1e3 if x is not None else None
print("y range", y.min(), y.max(), "x range", x.min(), x.max())
wd=a["wall_dist"]*1e3; V=a["volume"]
for nm in ["T","g_0","Y0","Y1","P"]:
    a[nm+"_d"]=b[nm]-a[nm]
edges=[0,0.02,0.05,0.1,0.2,0.4,0.8,1.6,10]
print("band[mm] nodes  vol%  mean dT  mean dg  mean dY0 mean dY1  mean gOFF mean S_OFF  mean T_OFF  onlyB")
onlyB=(b["g_0"]>1e-6)&(a["g_0"]<=1e-6)
for lo,hi in zip(edges[:-1],edges[1:]):
  for xr in [(10,35),(35,70)]:
    m=(wd>=lo)&(wd<hi)&(x>=xr[0])&(x<xr[1])
    if m.sum()==0: continue
    w=V[m]/V[m].sum()
    print(f"{lo:5.2f}-{hi:5.2f} x{xr} {m.sum():6d} {100*V[m].sum()/V.sum():5.2f} {np.sum(w*a['T_d'][m]):+8.3f} {np.sum(w*a['g_0_d'][m]):+.3e} {np.sum(w*a['Y0_d'][m]):+.3e} {np.sum(w*a['Y1_d'][m]):+.3e} {np.sum(w*a['g_0'][m]):.3e} {np.sum(w*a['condS_0'][m]):.3e} {np.sum(w*a['T'][m]):7.2f} {onlyB[m].sum()}")
print("onlyB: wall_dist mm pct 5/50/95", np.percentile(wd[onlyB],[5,50,95]), "x pct", np.percentile(x[onlyB],[5,50,95]), "g_B pct", np.percentile(b["g_0"][onlyB],[50,95]), "S_OFF", np.percentile(a["condS_0"][onlyB],[5,50,95]), "T_OFF",np.percentile(a["T"][onlyB],[5,50,95]))
# total water (vapor+liquid) per species: which Y is water? print means in core
core=(wd>1.6)&(x>40)
print("core means Y0,Y1,g OFF:", a["Y0"][core].mean(), a["Y1"][core].mean(), a["g_0"][core].mean())
# max dT location
i=np.argmax(np.abs(a["T_d"])); print("max|dT|", a["T_d"][i], "x,y,wd", x[i], y[i], wd[i], "g A/B", a["g_0"][i], b["g_0"][i])
# exit section profile
ex=x>x.max()-0.3
o=np.argsort(wd[ex])
for j in o[::max(1,len(o)//15)]:
    k=np.where(ex)[0][j]; print(f"exit wd {wd[k]:7.3f} T {a['T'][k]:7.2f} dT {a['T_d'][k]:+7.3f} g {a['g_0'][k]:.3e} dg {a['g_0_d'][k]:+.3e} Ux {a['Ux'][k]:7.1f}")
