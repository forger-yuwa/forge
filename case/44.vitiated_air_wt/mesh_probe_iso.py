#!/usr/bin/env python3
"""冷却壁 (Tw 300 K) 用ノズルメッシュの事前設計 (plan tooling-nozzle-isothermal-wall-chain §4.2 / §4.4)。

1. 断熱の生産 run (run_0107) の壁データから、各 x で「冷却壁にしたときの y1+ 倍率」
   f = sqrt(ρ_w,c τ_c / (ρ_w,a τ_a)) · μ_a/μ_c  (ρ_w,c = P_w/(R Tw), τ_c/τ_a ≈ 1.5 [case/48 B/A 実測], μ Sutherland)
   を見積り、目標 y1+_cold に要る第一セル高さ y1(x) を出す。
2. 候補 (wall_first_frac, wall_first_frac_throat, ni, throat_refine) でメッシュを生成し、AR/skew と y1(x) を表にする。
usage: design/.venv-opt/bin/python case/44.vitiated_air_wt/mesh_probe_iso.py
"""
import sys, glob, json, subprocess
from pathlib import Path
import numpy as np, h5py
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "design"))
from forge_design.probdef import load_problem
from forge_design.evaluate.runner_axismach import design_chain, _gam_or_gas
from forge_design.geometry.wall_axismach import PhysicalNozzleWall
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d

CASE = ROOT / "case/44.vitiated_air_wt"; REF = CASE / "run_0107_va_R2_LU6_Lc8_ns_ib_pass0"
prob = CASE / "problem_va_R2_LU6_Lc8_ns.yaml"
TW = 300.0; TAU_RATIO = 1.5

# ---- 1. 断熱 run の壁から冷却壁 y1+ 倍率
fn = sorted(glob.glob(str(REF / "res_[0-9]*.h5")), key=lambda s: int(s.split("_")[-1][:-3]))[-1]
f = h5py.File(fn, "r"); c = f["MESH/COORD"][:].reshape(-1, 3); wd = f["VALUE/wall_dist"][:]; vl = f["VALUE/vis_lam"][:]
ro = f["VALUE/ro"][:]; Ux = f["VALUE/Ux"][:]; Uy = f["VALUE/Uy"][:]; T = f["VALUE/T"][:]; P = f["VALUE/P"][:]
info = json.loads((REF / "prepare_info.json").read_text()); S = info["scale_m"]
wall = np.where(wd <= 0)[0]; xr = np.round(c[:, 0], 6)
def mu(T): return 1.716e-5 * (T / 273.0) ** 1.5 * 384.0 / (T + 111.0)
rows = []
for w in wall[np.argsort(c[wall, 0])][::10]:
    col = np.where(abs(xr - xr[w]) < 1e-9)[0]; col = col[wd[col] > 0]
    if len(col) == 0: continue
    j = col[np.argmin(wd[col])]; ut = np.hypot(Ux[j], Uy[j]); tau = vl[j] * ut / wd[j]
    yp_a = wd[j] * np.sqrt(ro[w] * tau) / vl[w]
    ro_c = P[w] / (287.0 * TW) if False else ro[w] * T[w] / TW
    fac = np.sqrt(ro_c * TAU_RATIO / ro[w]) * (vl[w] / mu(TW))
    rows.append((c[w, 0] / S, wd[j], T[w], yp_a, fac, yp_a * fac))
rows = np.array(rows)
print("x/rt   y1[um]  Tw_ad  y1+_ad  factor  y1+_cold(est)  y1 for y1+_cold=1 [um]")
for x in (-3, -1, -0.3, 0, 0.3, 1, 2, 4, 8, 12, 16, 20):
    k = np.argmin(abs(rows[:, 0] - x)); r = rows[k]
    print(f"{r[0]:6.2f} {r[1]*1e6:7.1f} {r[2]:6.0f} {r[3]:6.2f} {r[4]:6.2f} {r[5]:8.1f}   {r[1]*1e6/r[5]:7.2f}")
np.savetxt(CASE / "mesh_probe_iso_yplus.csv", rows, delimiter=",", header="x_rt,y1_m,Tw_ad,y1p_ad,factor,y1p_cold_est", comments="")

# ---- 2. 候補メッシュ
p = load_problem(prob); d = design_chain(p); scale = float(p.spec["r_throat"])
tbl = np.loadtxt(REF / "delta_r_initial.csv", delimiter=",", skiprows=1)
wallobj = PhysicalNozzleWall(d["wall"], d["wall_inv"], scale, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp,
                             offset="radial", delta_r_x=lambda x, _t=tbl: np.interp(x, _t[:, 0], _t[:, 1]))
import os as _os
_c = _os.environ.get("CANDS")
cands = eval(_c) if _c else [dict(ni=601, nj=97, wff=3.5e-5, wfft=None, tr=3.0)]
tmp = CASE / "_mesh_probe"; tmp.mkdir(exist_ok=True)
(tmp / "solverConfig.yaml").write_text('mesh: {meshFormat: "hdf5", discretization: "node", nodeWallDirichlet: 1, isAxisymmetric: 1, meshFileName: "m.h5", valueFileName: "m.h5"}\ngpu: 1\nsolver: "SLAU"\nphysProp: {isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1190.2, gamma: 1.32752}\ntime:\n  unsteady: 0\n  dualTime: 0\n  last: {control: 0, nStepOuter: 10}\n  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-9, dt_max: 1.0, detectNaN: 1}\n  outStepStart: 0\n  outStepInterval: 10\n  timeIntegration: 11\n  nStepInner: 5\nspace: {convMethod: 0, limiter: 0}\nturbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0}\ninitial: "uniform_p101325_u10"\n')
(tmp / "bcondConfig.yaml").write_text("inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 1140000.0, Tt: 1060.0, k: 1.0, omega: 18000.0}}\noutlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {Ps: 5235.0, Pt: 5235.0, Tt: 300.0}}\nwall:   {physID: 3, kind: wall, outputHDFflg: 1, ints: , floats: }\naxis:   {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n")
import os
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
for cd in cands:
    mp = Mesh2DParams(ni=cd["ni"], nj=cd["nj"], wall_first_frac=cd["wff"], throat_refine=cd["tr"], throat_width=cd.get("tw", 1.5), wall_first_frac_throat=cd["wfft"], wall_first_blend_x0=cd.get("bx0", 0.5), wall_first_blend_x1=cd.get("bx1", 6.0), wall_first_up_x0=cd.get("ux0"), wall_first_up_x1=cd.get("ux1"), local_center=cd.get("lc", 0.0), local_refine=cd.get("lr", 1.0), local_width=cd.get("lw", 0.75), scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wallobj, mp)
    write_msh41_2d(tmp / "m.msh", coords, quads, bedges)
    r = subprocess.run([str(ROOT / "solver_density_cuda/build/convertGmshToForge"), "m.msh", "m.h5"], cwd=tmp, env=ENV, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_mesh_quality.py"), str(tmp / "m.h5"), "--ar-max", str(cd.get("armax", 1000))], capture_output=True, text=True)
    ar = [l for l in q.stdout.splitlines() if "aspect" in l or "VERDICT" in l]
    # y1(x): 壁点列と隣接点
    X = coords[:, 0].reshape(cd["ni"], cd["nj"]); Rr = coords[:, 1].reshape(cd["ni"], cd["nj"])
    y1 = (Rr[:, -1] - Rr[:, -2]); xs = X[:, -1] / scale
    est = []
    for x in (-3, -1, 0, 1, 2, 4, 8, 12, 20):
        k = np.argmin(abs(xs - x)); kk = np.argmin(abs(rows[:, 0] - x))
        y1p_cold = rows[kk, 5] * y1[k] / rows[kk, 1]
        est.append(f"x{x:+d}: y1 {y1[k]*1e6:.1f}um y1+c {y1p_cold:.1f}")
    # AR>1000 の位置 (構造格子: quad の dx/dy 近似)
    dx = np.abs(np.diff(X, axis=0))[:, :-1]; dy = np.abs(np.diff(Rr, axis=1))[:-1, :]
    arq = np.maximum(dx, dy) / np.maximum(np.minimum(dx, dy), 1e-30)
    bad = arq > 1000.0
    xb = (X[:-1, :-1])[bad] / scale
    hist = np.histogram(xb, bins=[-20, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 12, 16, 30])[0] if bad.any() else []
    print(f"\n[ni {cd['ni']} nj {cd['nj']} wff {cd['wff']} wfft {cd['wfft']} tr {cd['tr']}] nodes {len(coords)}  " + " | ".join(ar))
    print("   AR>1000 x-hist bins[-20,-6,-4,-2,-1,0,1,2,4,6,8,12,16,30]:", list(hist), " max AR at x/rt =", float(X[:-1, :-1][np.unravel_index(np.argmax(arq), arq.shape)] / scale) if bad.any() else None)
    print("   " + " | ".join(est))
