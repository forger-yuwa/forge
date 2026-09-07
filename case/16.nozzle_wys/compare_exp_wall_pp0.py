import numpy as np, glob, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager
for fp in glob.glob('/home/sano/.fonts/*Noto*CJK*'): font_manager.fontManager.addfont(fp)
plt.rcParams['font.family']='Noto Sans CJK JP'
S='./wall_pp0_csv/'
exp=np.array([[float(v) for v in l.split(',')[:3]] for l in open('wyslouzil_fig3_pp0.csv').read().strip().split('\n')[1:]])
xe, iso, cond = exp[:,0]*10, exp[:,1], exp[:,2]   # mm
runs={ # label: (csv, dry/cond)
 '2D Euler node (run_0193)':('wall_run_0193_user_node_euler_dry.csv','dry'),
 '2D SST node (run_0213)':('wall_run_0213_user_node_sst_dry_outflow_cont.csv','dry'),
 '2D SST cell (run_0196)':('wall_run_0196_user_cell_sst_dry_cont.csv','dry'),
 '3D laminar node (run_0221, 対称面)':('wall_run_0221_user_node3d_lam_dry_half.csv','dry'),
 '3D SST node (run_0228, 対称面, 出口バッファ+Ps)':('wall_run_0228_user_node3d_sst_dry_half_ext_ps_relguard.csv','dry'),
 '2D Euler cond node (run_0194)':('wall_run_0194_user_node_euler_cond.csv','cond'),
 '2D SST cond cell (run_0197)':('wall_run_0197_user_cell_sst_cond_cont.csv','cond'),
 '3D laminar cond node (run_0224, 対称面)':('wall_run_0224_user_node3d_lam_cond_half.csv','cond'),
}
def load(fn):
    a=np.genfromtxt(S+fn, delimiter=',', names=True); return a['x_mm'], a['contour_wall']
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
ax[0].plot(xe, iso, 'ko', ms=6, label='実験 dry (isentrope, Wyslouzil Fig.3)'); ax[1].plot(xe, cond, 'ks', ms=6, label='実験 凝縮 p_v0=1.00 kPa (Fig.3)')
rows=[]
for lab,(fn,kind) in runs.items():
    x,p=load(fn); m=(x>=-1)&(x<=95); a=ax[0] if kind=='dry' else ax[1]
    a.plot(x[m],p[m],'-',lw=1.4,label=lab)
    ref = iso if kind=='dry' else cond
    pi=np.interp(xe, x, p); dev=(pi-ref)/ref*100
    rows.append((lab, kind, pi, dev))
for a,t in zip(ax,['dry: 壁 p/p0 vs isentrope','凝縮 (H2O 1.0 kPa): 壁 p/p0 vs 実験']):
    a.set_xlim(-2,95); a.set_ylim(0.15,0.55); a.set_xlabel('x from throat [mm]'); a.set_ylabel('p_wall / p0'); a.grid(alpha=.3); a.legend(fontsize=8); a.set_title(t)
fig.tight_layout(); fig.savefig(S+'compare_exp_wall_pp0.png', dpi=120)
print('x_exp [mm]:', ' '.join(f'{v:6.1f}' for v in xe))
print('exp iso    :', ' '.join(f'{v:6.3f}' for v in iso)); print('exp cond   :', ' '.join(f'{v:6.3f}' for v in cond))
for lab,kind,pi,dev in rows:
    print(f'{lab:44s} [{kind}] p/p0: '+' '.join(f'{v:6.3f}' for v in pi))
    print(f'{"":44s}        dev%: '+' '.join(f'{v:+6.1f}' for v in dev)+f'   | mean|dev| x>=10mm: {np.mean(np.abs(dev[xe>=10])):.2f} %')
