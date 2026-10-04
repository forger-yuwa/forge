# total-water / liquid / mass flux through x-planes, ON(run_0521) vs OFF(run_0520), trapezoid over node columns (planar 2D)
import h5py, numpy as np
C="/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
def ld(r):
    f=h5py.File(C+r+"/res_16000.h5"); V=f["VALUE"]
    d={k:np.array(V[k][...]).ravel().astype(float) for k in ["ro","Ux","Y1","g_0","T","h0"]}
    d["xyz"]=np.array(f["MESH/COORD"][...]).reshape(-1,3); return d
a=ld("run_0520_twophase_off_16k"); b=ld("run_0521_twophase_dplur_nn0_cont")
x=a["xyz"][:,0]*1e3; y=a["xyz"][:,1]; z=a["xyz"][:,2]
z0=z==z.min()
xs=np.unique(np.round(x[z0],4))
for xt in [-50,0,20,40,60,80,94.9]:
    xc=xs[np.argmin(abs(xs-xt))]; m=z0&(abs(x-xc)<1e-4); o=np.argsort(y[m]); yy=y[m][o]
    out=[]
    for d in (a,b):
        q=lambda v: np.trapezoid(v[m][o],yy)
        rU=d["ro"]*d["Ux"]; out.append((q(rU),q(rU*d["Y1"]),q(rU*d["g_0"]),q(rU*d["h0"])))
    (m0,w0,l0,h0),(m1,w1,l1,h1)=out
    print(f"x={xc:7.2f} n={m.sum():3d} mdot {m0:.6e} dRel {(m1-m0)/m0:+.2e} | water {w0:.6e} dRel {(w1-w0)/w0:+.2e} | liq {l0:.4e} dRel {(l1-l0)/max(l0,1e-30):+.2e} | H {h0:.6e} dRel {(h1-h0)/h0:+.2e} | Y1w/m OFF {w0/m0:.6e} ON {w1/m1:.6e}")
