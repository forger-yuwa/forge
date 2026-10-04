# limiter-inlet-oscillation: compare scaled Venkatakrishnan eps-hat with the normalised gradient increment and float32 rounding at the residual-dominant nodes
import h5py, numpy as np, sys
C="/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
r=sys.argv[1] if len(sys.argv)>1 else "run_0524_floor_dry_L1"
f=h5py.File(C+r+"/res_48000.h5"); V=f["VALUE"]; xyz=np.array(f["MESH/COORD"][...]).reshape(-1,3)
g=lambda k: np.array(V[k][...]).ravel().astype(np.float64)
K,L=0.05,0.1570673715; ref={"ro":0.5279511755,"P":43128.33253,"Ux":330.8270384,"Uy":330.8270384}
vol=g("volume"); zs=np.unique(np.round(xyz[:,2],9)); dz=(zs.max()-zs.min()) if zs.size>1 else 1.0
print("z levels", zs.size, "dz", dz)
for name,area in (("sqrt(volume)",vol),("sqrt(volume/dz)",vol/dz)):
    h=np.sqrt(area); epsh=(K*h/L)**1.5
    print(name, "h pct 1/50/99", np.percentile(h,[1,50,99]), "eps_hat pct 1/50/99", np.percentile(epsh,[1,50,99]))
h=np.sqrt(vol/dz); epsh=(K*h/L)**1.5
res=np.abs(g("res_ro")); top=res>=np.quantile(res,0.99)
for q,dk in (("P","dPd"),("ro","drod"),("Ux","dUxd"),("Uy","dUyd")):
    gr=np.hypot(g(dk+"x"),g(dk+"y"))*h/ref[q]; val=np.abs(g(q))/ref[q]; ulp=val*2**-23
    lim=g("limiter_"+q)
    for nm,m in (("top1%",top),("all",np.ones_like(top))):
        print(f"{q:3s} {nm:6s} |grad|h/qref p50 {np.median(gr[m]):.2e} | eps_hat p50 {np.median(epsh[m]):.2e} | ulp(q)/qref p50 {np.median(ulp[m]):.2e} | frac(gr<eps) {np.mean(gr[m]<epsh[m]):.2f} frac(gr<100ulp) {np.mean(gr[m]<100*ulp[m]):.2f} | limiter p50 {np.median(lim[m]):.3f} frac<0.99 {np.mean(lim[m]<0.99):.2f}")
o=np.argsort(res)[::-1][:8]
for i in o:
    print(f"node x {xyz[i,0]*1e3:7.2f} y {xyz[i,1]*1e3:7.2f} h {h[i]:.2e} eps_hat {epsh[i]:.2e} | P: grh {np.hypot(g('dPdx')[i],g('dPdy')[i])*h[i]/ref['P']:.2e} lim {g('limiter_P')[i]:.3f} | ro: grh {np.hypot(g('drodx')[i],g('drody')[i])*h[i]/ref['ro']:.2e} lim {g('limiter_ro')[i]:.3f} | Ux lim {g('limiter_Ux')[i]:.3f} Uy lim {g('limiter_Uy')[i]:.3f} T lim {g('limiter_T')[i]:.3f}")
