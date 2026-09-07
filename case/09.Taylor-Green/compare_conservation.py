#!/usr/bin/env python3
# KEEP explicit vs implicit(CFL掃引)を energy/entropy 保存の観点で比較。
# 保存的FVなので 質量 M=ΣρV・総エネルギー E=Σroe·V は機械精度保存のはず。
# KE が落ちた分は内部エネルギー Ei=E-K へ移り(数値散逸=加熱)、エントロピー S が増える、を確認する。
import glob, os, re
import numpy as np, h5py, yaml

def _cell_volume(f, path):
    """セル体積: res には既定 (output.level 1) で volume が無いので、mesh h5 (solverConfig の meshFileName) の CELLS/volume を読む。"""
    import os, yaml
    if "VALUE/volume" in f:
        return f["VALUE/volume"][:].astype(np.float64)
    rd = os.path.dirname(os.path.abspath(path))
    mesh = yaml.safe_load(open(os.path.join(rd, "solverConfig.yaml")))["mesh"]["meshFileName"]
    with h5py.File(os.path.join(rd, mesh), "r") as m:
        return m["CELLS/volume"][:].astype(np.float64)

gamma=1.4; cp=0.4; cv=cp/gamma
def step_of(p):
    m=re.search(r"res_(\d+)\.h5$",os.path.basename(p)); return int(m.group(1)) if m else -1
def integ(path):
    with h5py.File(path,"r") as f:
        V = _cell_volume(f, path); ro=f["VALUE/ro"][:].astype(np.float64)
        P=f["VALUE/P"][:].astype(np.float64); roe=f["VALUE/roe"][:].astype(np.float64)
        ux=f["VALUE/Ux"][:].astype(np.float64); uy=f["VALUE/Uy"][:].astype(np.float64); uz=f["VALUE/Uz"][:].astype(np.float64)
    M=np.sum(ro*V); E=np.sum(roe*V)
    K=np.sum(0.5*ro*(ux*ux+uy*uy+uz*uz)*V); Ei=E-K
    S=np.sum(ro*cv*np.log(P/np.power(ro,gamma))*V)
    return dict(M=M,E=E,K=K,Ei=Ei,S=S)
def first_last(d):
    fs=sorted([p for p in glob.glob(f"{d}/res_*.h5") if step_of(p)>=0 and "nan" not in p],key=step_of)
    return integ(fs[0]), integ(fs[-1])
RUNS=[("explicit  CFL.05","run_0011_cell_keep_expl_ref"),
("implicit  CFL.05","run_0012_cell_keep_impl_cfl005"),
("implicit  CFL.5 ","run_0014_cell_keep_impl_cfl05"),
("implicit  CFL1  ","run_0015_cell_keep_impl_cfl1"),
("implicit  CFL2  ","run_0016_cell_keep_impl_cfl2"),
("implicit  CFL4  ","run_0017_cell_keep_impl_cfl4"),
("implicit  CFL8  ","run_0018_cell_keep_impl_cfl8"),
("implicit  CFL16 ","run_0019_cell_keep_impl_cfl16")]
print(f"{'run':18s}{'dM/M0':>11}{'dE/E0':>11}{'dK/K0':>10}{'dEi/Ei0':>10}{'dS/|S0|':>11}")
for nm,d in RUNS:
    a,b=first_last(d)
    dM=(b['M']-a['M'])/a['M']; dE=(b['E']-a['E'])/a['E']
    dK=(b['K']-a['K'])/a['K']; dEi=(b['Ei']-a['Ei'])/a['Ei']; dS=(b['S']-a['S'])/abs(a['S'])
    print(f"{nm:18s}{dM:11.2e}{dE:11.2e}{dK:10.4f}{dEi:10.2e}{dS:11.2e}")
