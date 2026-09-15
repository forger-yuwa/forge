"""run_dir の最終 res_*.h5 から軸 (r=0) 上の全 VALUE/* を x 昇順の CSV に書き出す (h5 後処理)。
usage: python3 axis_csv_va.py <run_dir> [...]  → <run_dir>/axis_values.csv (列: x_m, x_over_rt, r_m, 全 VALUE 名 [ベクトル/3成分はそのまま])"""
import sys, json
from pathlib import Path
import numpy as np, h5py
for rd in map(Path, sys.argv[1:]):
    res=sorted(rd.glob("res_[0-9]*.h5"), key=lambda f:int("".join(c for c in f.stem if c.isdigit())))[-1]
    info=json.loads((rd/"prepare_info.json").read_text()); S=float(info["scale_m"])
    with h5py.File(rd/"nozzle.h5") as n: nc=n["MESH/COORD"][:].reshape(-1,3)
    ax=nc[:,1]<1e-9*max(nc[:,1].max(),1.0); idx=np.where(ax)[0][np.argsort(nc[ax,0])]
    cols={"x_m":nc[idx,0],"x_over_rt":nc[idx,0]/S,"r_m":nc[idx,1]}
    with h5py.File(res) as f:
        for k in sorted(f["VALUE"].keys()):
            v=f["VALUE/"+k][:]
            if v.ndim==1 and len(v)==len(nc): cols[k]=v[idx]
            elif v.ndim==2 and v.shape[0]==len(nc):
                for j in range(v.shape[1]): cols[f"{k}_{j}"]=v[idx,j]
    # 後処理 Tsat (forge と同じ Murphy-Koop 過冷却液 p_sat の反転)。dry run (凝縮 off, condTsat 無し) でも出す。
    # p_H2O = Y_H2O P (Y1 があれば TP split の H2O 質量分率 → モル分率換算; 無ければ設計組成 x_H2O=0.0462711)
    def psat_l(T):
        Tc=np.maximum(T,120.0)
        return np.exp(54.842763-6763.22/Tc-4.210*np.log(Tc)+0.000367*Tc+np.tanh(0.0415*(Tc-218.8))*(53.878-1331.22/Tc-9.44523*np.log(Tc)+0.014025*Tc))
    P=cols["P"]; T=cols["T"]
    # H2O の配列は種名で解決する (tools/forge_species.py; 5 種 full では Y0、split_h2o では Y1)。Y1 決め打ちはしない (plan cea-mole-fraction §4.4)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "solver_density_cuda" / "tools"))
    from forge_species import species_info
    yv = species_info(rd).get("vapor_array")
    if yv and yv in cols and "g_0" in cols: pv=cols["ro"]*(cols[yv]-cols["g_0"])*461.5*T   # carrier 式 (forge と同一)
    elif yv and yv in cols: pv=cols["ro"]*cols[yv]*461.5*T
    else: pv=0.0462711*P
    Ts=np.full_like(T,250.0)
    for _ in range(40):
        f=np.log(psat_l(Ts))-np.log(np.maximum(pv,1e-30)); dfdT=np.maximum((3.1485e6-2370.0*Ts).clip(2e6,3e6)/(461.5*Ts**2),1e-6)
        Ts=np.clip(Ts-np.clip(f/dfdT,-0.3*Ts,0.3*Ts),50,1000)
    Ts=np.where(pv>1e-6,Ts,0.0)
    cols["pH2O_post"]=pv; cols["Tsat_post"]=Ts; cols["subcool_post"]=Ts-T; cols["S_post"]=pv/psat_l(T)
    hdr=",".join(cols); arr=np.column_stack([cols[k] for k in cols])
    np.savetxt(rd/"axis_values.csv", arr, delimiter=",", header=f"{res.name}; axis nodes r=0 of {rd.name}; scale_m={S}\n"+hdr, comments="# " if False else "", fmt="%.9g")
    # 先頭行のコメント記号を整える (1 行目 # 情報, 2 行目 ヘッダ)
    txt=(rd/"axis_values.csv").read_text().splitlines(); txt[0]="# "+txt[0]; (rd/"axis_values.csv").write_text("\n".join(txt)+"\n")
    print(f"{rd.name}: {len(idx)} axis nodes, {len(cols)} columns -> axis_values.csv")
