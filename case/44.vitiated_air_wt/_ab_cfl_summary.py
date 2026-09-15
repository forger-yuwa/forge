import h5py, numpy as np, sys, glob
n=h5py.File('run_0127_va3_M4.19_Lc8_noneq_inletTt_merged/nozzle.h5'); nc=n['MESH/COORD'][:].reshape(-1,3); S=0.205711
x=nc[:,0]/S; r=nc[:,1]/S
cols=np.unique(nc[:,0]); lines={}
for xc in cols:
    idx=np.where(np.abs(nc[:,0]-xc)<1e-9)[0]; idx=idx[np.argsort(r[idx])]
    for jj in range(1,9): lines.setdefault(jj,[]).append(idx[-jj])
xl=cols/S; xc2=cols[np.argmin(abs(cols-(nc[:,0].max()-2*S)))]; ex=np.abs(nc[:,0]-xc2)<1e-9
wall=np.array(lines[1]); j4=np.array(lines[4])
for run in sys.argv[1:]:
    res=sorted(glob.glob(run+'/res_[0-9]*.h5'), key=lambda f:int(''.join(c for c in f.split('/')[-1] if c.isdigit())))[-1]
    f=h5py.File(res); g=f['VALUE/g_0'][:]; L=f['VALUE/condLim_0'][:]; dt=f['VALUE/dt_local'][:]; ro=f['VALUE/ro'][:]; ux=f['VALUE/Ux'][:]; M=np.hypot(ux,f['VALUE/Uy'][:])/f['VALUE/sonic'][:]; T=f['VALUE/T'][:]; Sv=f['VALUE/condS_0'][:]
    w=(ro*ux*r)[ex]; ga=np.trapz(w*g[ex],r[ex])/np.trapz(w,r[ex]); Ma=np.trapz(w*M[ex],r[ex])/np.trapz(w,r[ex])
    def onset(line): 
        i=np.where(g[line]>3.77e-5)[0]; return xl[i[0]] if len(i) else np.nan
    def x_g(line,frac):
        i=np.where(g[line]>frac*0.0377)[0]; return xl[i[0]] if len(i) else np.nan
    cond=g>1e-4
    print(f'{run[:9]} {res.split("/")[-1]}: dt_local(j4,x20) {dt[j4][np.argmin(abs(xl-20))]:.2e}  condLim median(cond) {np.median(L[cond]):.2f} min {L[cond].min():.2f} | onset wall {onset(wall):.2f} j4 {onset(j4):.2f} | x(g=10%Yv) wall {x_g(wall,0.1):.2f} j4 {x_g(j4,0.1):.2f} | x(g=30%Yv) wall {x_g(wall,0.3):.2f} j4 {x_g(j4,0.3):.2f} | g_max {g.max():.5f} exit g_avg {ga:.5f} exit M {Ma:.4f} T_min {T.min():.1f}')
