#!/usr/bin/env python3
"""case/44 ノズル CPG × SU2 等温クロスチェック (plan tooling-nozzle-isothermal-wall-chain §4.4)。

forge (CPG, 素 SST) と SU2 (IDEAL_GAS γ 1.32752, SST V2003m) を同一メッシュ (nozzle.msh → nozzle.su2) で
断熱 / 等温 300 K の 2 対にする。SU2 は forge の収束場から warm start (restart_flow.csv を forge res から生成)。
  run_0112_va_cpg_fine_ad_plain      forge 断熱 (積分法初期壁 断熱, IC = run_0108 [同一メッシュ index] か run_0107 cross-mesh)
  run_0113_va_cpg_fine_iso300_plain  forge 300 K (同じ壁 = run_0112/delta_r_initial.csv, IC = run_0112)
  su2_va_cpg_ad / su2_va_cpg_iso300  SU2 (それぞれ forge 場から warm start)
usage: design/.venv-opt/bin/python run_iso_chain_cpg.py --stage forge_ad|forge_iso|su2_ad|su2_iso|all [--wait-for FILE:TOKEN]
"""
import argparse, csv, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "design"))
CASE = Path(__file__).resolve().parent
EULER = CASE / "run_0005_va_R2_LU6_Lc8"
P_AD = CASE / "problem_va_R2_LU6_Lc8_ns_fine_cpg.yaml"; P_ISO = CASE / "problem_va_R2_LU6_Lc8_ns_fine_cpg_iso300.yaml"
RA = CASE / "run_0112_va_cpg_fine_ad_plain"; RB = CASE / "run_0113_va_cpg_fine_iso300_plain"
SA = CASE / "su2_va_cpg_ad"; SB = CASE / "su2_va_cpg_iso300"
SU2_BIN = ROOT / ".external/su2/bin/SU2_CFD"
GAM, CP = 1.32752, 1190.2; R = CP * (GAM - 1) / GAM
PLAIN = {"dilatationCorrection: 2": "dilatationCorrection: 0", "katoLaunder: 1": "katoLaunder: 0",
         "turbulentPrandtl: 0.9}": "turbulentPrandtl: 0.9, sstOmegaProdFromPk: 0, sstSigmaBlend: 0, sstEnergyIncludesK: 0}"}


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def make_plain(run_dir):
    cfg = (run_dir / "solverConfig.yaml").read_text()
    for a, b in PLAIN.items():
        assert a in cfg, a
        cfg = cfg.replace(a, b)
    (run_dir / "solverConfig.yaml").write_text(cfg)


def run_forge(problem, run_dir, ic_from, **kw):
    from forge_design.evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    info = prepare_ns(problem, run_dir, ic_from=ic_from, euler_ref=EULER, **kw)
    make_plain(run_dir)
    info["sst"] = "plain (dilatation 0, katoLaunder 0, omegaProdFromPk 0, sigmaBlend 0, energyIncludesK 0)"
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    log(f"prepared {run_dir.name}: {info['mesh']} wall_thermal {info.get('wall_thermal')}")
    log(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1])
    t0 = time.time(); rc = run_staged_ns(run_dir, stages="full")
    log(f"{run_dir.name} rc={rc} {time.time()-t0:.0f}s")
    m = collect(problem, run_dir); (run_dir / "metrics.json").write_text(json.dumps(m, indent=1, default=str))
    log(f"{run_dir.name} verdict: {m.get('convergence_verdict')}")
    return rc


SU2_CFG = """% case/44 ノズル CPG クロスチェック — SU2 v8.5 RANS-SST (forge run_0112/0113 と同一メッシュ・物性)
SOLVER= RANS
KIND_TURB_MODEL= SST
SST_OPTIONS= V2003m
MATH_PROBLEM= DIRECT
SYSTEM_MEASUREMENTS= SI
AXISYMMETRIC= YES
RESTART_SOL= YES
RESTART_FILENAME= restart_flow
SOLUTION_FILENAME= restart_flow_in
READ_BINARY_RESTART= NO
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= {gam}
GAS_CONSTANT= {R:.4f}
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716e-5
MU_T_REF= 273.0
SUTHERLAND_CONSTANT= 111.0
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
PRANDTL_LAM= 0.72
PRANDTL_TURB= 0.9
REF_DIMENSIONALIZATION= DIMENSIONAL
INIT_OPTION= TD_CONDITIONS
FREESTREAM_PRESSURE= 1.0e6
FREESTREAM_TEMPERATURE= 1000.0
MACH_NUMBER= 0.15
FREESTREAM_TURBULENCEINTENSITY= 0.01
FREESTREAM_TURB2LAMVISCRATIO= 10.0
INLET_TYPE= TOTAL_CONDITIONS
MARKER_INLET= ( inlet, 1060.0, 1140000.0, 1.0, 0.0, 0.0 )
MARKER_OUTLET= ( outlet, 5235.0 )
{wall}
MARKER_SYM= ( axis )
MARKER_PLOTTING= ( wall )
MARKER_MONITORING= ( wall )
NUM_METHOD_GRAD= GREEN_GAUSS
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
VENKAT_LIMITER_COEFF= 0.05
TIME_DISCRE_FLOW= EULER_IMPLICIT
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_NUMBER= 2.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.5, 1.5, 1.0, 30.0, 0.001, 50 )
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1e-4
LINEAR_SOLVER_ITER= 10
ITER= {iters}
CONV_FIELD= RMS_DENSITY
CONV_RESIDUAL_MINVAL= -12
CONV_STARTITER= 500
MARKER_ANALYZE= ( outlet )
MARKER_ANALYZE_AVERAGE= MASSFLUX
MESH_FORMAT= SU2
MESH_FILENAME= nozzle.su2
TABULAR_FORMAT= CSV
OUTPUT_FILES= ( RESTART_ASCII, PARAVIEW, SURFACE_CSV )
VOLUME_OUTPUT= ( COORDINATES, SOLUTION, PRIMITIVE )
CONV_FILENAME= history
VOLUME_FILENAME= flow
SURFACE_FILENAME= surface_flow
OUTPUT_WRT_FREQ= 2000
SCREEN_OUTPUT= ( INNER_ITER, RMS_DENSITY, RMS_MOMENTUM-Y, RMS_ENERGY, RMS_TKE, RMS_DISSIPATION, CFL_NUMBER, AVG_MASSFLOW, TOTAL_HEATFLUX )
HISTORY_OUTPUT= ( ITER, RMS_RES, AERO_COEFF, HEAT, FLOW_COEFF )
"""


def write_su2_restart(forge_run, su2_dir, su2_mesh, header=None):
    """forge の最終 res_*.h5 (node) を SU2 restart_flow_in.csv (RESTART_ASCII) に変換する。
    節点は座標で対応付ける (KD-tree, 一致距離 < 1e-9 m を要求)。"""
    import h5py
    from scipy.spatial import cKDTree
    res = sorted(Path(forge_run).glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
    with h5py.File(res, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3); V = f["VALUE"]
        ro = V["ro"][:].astype(float); u = V["Ux"][:].astype(float); v = V["Uy"][:].astype(float)
        T = V["T"][:].astype(float); k = V["k"][:].astype(float); om = V["omega"][:].astype(float)
    # SU2 メッシュの節点順を読む
    pts = []
    with open(su2_mesh) as fh:
        npoin = None
        for line in fh:
            if line.startswith("NPOIN"):
                npoin = int(line.split("=")[1].split()[0]); break
        for _ in range(npoin):
            s = fh.readline().split(); pts.append((float(s[0]), float(s[1])))
    pts = np.array(pts)
    tree = cKDTree(c[:, :2]); dist, idx = tree.query(pts)
    # forge h5 の座標は float32 → 数 e-7 m のずれは正常。最小セル 0.5 µm より十分小さい 1e-6 を許容し、対応の一意性を確認
    if dist.max() > 1e-6 or len(set(idx.tolist())) != len(idx):
        raise RuntimeError(f"SU2/forge 節点が一致しない (max dist {dist.max():.3e}, unique {len(set(idx.tolist()))}/{len(idx)})")
    e = CP / GAM * T + 0.5 * (u ** 2 + v ** 2)          # e_int = cv T (CPG), E = e + ek
    hdr = header or ["PointID", "x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Turb_Kin_Energy", "Omega"]
    with open(Path(su2_dir) / "restart_flow_in.csv", "w", newline="") as fo:
        w = csv.writer(fo, quoting=csv.QUOTE_NONE); w.writerow(hdr)
        for i, j in enumerate(idx):
            w.writerow([i, pts[i, 0], pts[i, 1], ro[j], ro[j] * u[j], ro[j] * v[j], ro[j] * e[j], ro[j] * k[j], ro[j] * om[j]])
    log(f"restart written: {len(idx)} points from {res.name} (max coord mismatch {dist.max():.2e})")


def launch_su2(forge_run, su2_dir, wall, threads=5, iters=6000):
    su2_dir.mkdir(exist_ok=True)
    if not (su2_dir / "nozzle.su2").exists():
        subprocess.run(["gmsh", str(forge_run / "nozzle.msh"), "-o", str(su2_dir / "nozzle.su2"), "-format", "su2", "-save", "-v", "1"], check=True)
    write_su2_restart(forge_run, su2_dir, su2_dir / "nozzle.su2")
    (su2_dir / "nozzle.cfg").write_text(SU2_CFG.format(gam=GAM, R=R, wall=wall, iters=iters))
    subprocess.Popen(f"cd {su2_dir} && OMP_NUM_THREADS={threads} nohup {SU2_BIN} nozzle.cfg > su2_run.log 2>&1 &", shell=True)
    log(f"SU2 launched in {su2_dir.name} ({threads} threads, {iters} iters)")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", default="all"); ap.add_argument("--wait-for", default=None)
    ap.add_argument("--ic-ad", default=str(CASE / "run_0108_va_ns_fine_ad_samewall")); ap.add_argument("--su2-threads", type=int, default=5)
    a = ap.parse_args()
    if a.wait_for:
        f, tok = a.wait_for.split(":")
        while not (Path(f).exists() and tok in Path(f).read_text()):
            time.sleep(30)
        log("GPU free")
    if a.stage in ("forge_ad", "all"):
        ic = Path(a.ic_ad) if Path(a.ic_ad).exists() else CASE / "run_0107_va_R2_LU6_Lc8_ns_ib_pass0"
        run_forge(P_AD, RA, ic_from=ic, initializer={"model": "contur"})
    if a.stage in ("su2_ad", "all"):
        launch_su2(RA, SA, "MARKER_HEATFLUX= ( wall, 0.0 )", threads=a.su2_threads)
    if a.stage in ("forge_iso", "all"):
        run_forge(P_ISO, RB, ic_from=RA, delta_r_csv=str(RA / "delta_r_initial.csv"), offset="radial")
    if a.stage in ("su2_iso", "all"):
        launch_su2(RB, SB, "MARKER_ISOTHERMAL= ( wall, 300.0 )", threads=a.su2_threads)


if __name__ == "__main__":
    main()
