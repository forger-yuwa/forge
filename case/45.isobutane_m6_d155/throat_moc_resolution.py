"""plan tooling-nozzle-throat-monotone-r2 §3 仮説 H 候補 (3) の切り分け: MOC の分解能 (n_start, n_axis_inv) を振り、第 1 区間の曲率と x=0.025 の流れ角超過を見る (CFD 0 step)。
usage: python3 throat_moc_resolution.py  (凍結元 run はローカル主ツリー /home/sano/work/forge/case/45.isobutane_m6_d155 を参照)
"""
import sys, numpy as np, time
sys.path.insert(0,'/home/sano/work/forge-integ-1005/design')
from forge_design.evaluate.runner_axismach import design_chain, load_problem
M="/home/sano/work/forge/case/45.isobutane_m6_d155/"
for prob in ("problem_d155_euler_c2final_n2400.yaml","problem_d155_ns_finemesh_recal_final.yaml"):
    for ns, na in ((41,2400),(81,2400),(161,2400),(41,4800),(81,4800)):
        p=load_problem(prob); p.geometry["n_start"]=ns; p.geometry["n_axis_inv"]=na
        if p.geometry.get("initial_line")=="cfd": p.geometry["initial_line_run"]=M+p.geometry["initial_line_run"]
        t=time.time()
        try:
            d=design_chain(p)
        except Exception as e:
            print(prob[:24],ns,na,"ERR",str(e)[:100]); continue
        tb=d['wall_inv']; x,th=tb[:,0],tb[:,2]; k=np.diff(np.tan(th))/np.diff(x)
        # 区間幅に依らない比較: x=0.025 での流れ角 (点間を tan θ 線形補間)
        th25=np.arctan(np.interp(0.025,x,np.tan(th)))
        print(f"{prob[:26]:26s} n_start {ns:3d} n_axis {na}: x1 {x[1]:.4f} seg k {np.round(k[:3],4)} theta(0.025)-atan(0.025/R) {np.degrees(th25-np.arctan(0.025/2)):+.4f} deg ({time.time()-t:.0f}s)", flush=True)
