import sys, numpy as np, h5py
mesh_h5, res_h5 = sys.argv[1:3]; targets=[int(a) for a in sys.argv[3:]]
with h5py.File(mesh_h5,"r") as f:
    xn=f["MESH/COORD"][...].astype(np.float64).reshape(-1,3); st=f["PLANES/STRUCT"][...].tolist(); nP=f["PLANES/surfArea"].shape[0]
N=xn.shape[0]
with h5py.File(res_h5,"r") as f:
    V=f["VALUE"]; ro=V["ro"][...].astype(np.float64); ux=V["roUx"][...]/ro; T=V["T"][...] if "T" in V else None
nb={}; bd={}
pos=0
for ip in range(nP):
    nn=st[pos]; pos+=1+nn; nc=st[pos]; cells=st[pos+1:pos+1+nc]; pos+=1+nc
    if nc==2 and cells[1]<N and cells[0]<N:
        a,b=cells; nb.setdefault(a,[]).append(b); nb.setdefault(b,[]).append(a)
    else:
        bd[cells[0]]=bd.get(cells[0],0)+1
for i in targets:
    print(f"node {i} x={xn[i]} Ux={ux[i]:.1f} T={T[i] if T is not None else float('nan'):.1f} bnd_halffaces={bd.get(i,0)}")
    for j in nb.get(i,[]):
        d=xn[j]-xn[i]
        print(f"   -> {j:8d} d=({d[0]*1e3:+8.4f},{d[1]*1e3:+8.4f},{d[2]*1e3:+8.4f}) mm  Ux={ux[j]:8.1f}  T={T[j] if T is not None else float('nan'):7.1f} bnd={bd.get(j,0)}")
