#!/usr/bin/env python3
"""case/57 遷移平板 (ERCOFTAC T3A / T3B) の forge run を生成して段階起動する。

usage:
  python3 gen_runs.py --run run_0001_t3a_sst --mesh plate_base --model sst
  python3 gen_runs.py --run run_0002_t3a_lm  --mesh plate_base --model lm --ic-from run_0001_t3a_sst [--reset-turb]
  --model sst|lm|lam   lm = SST + turbulence.transition: lm2009
  --case  t3a|t3b      入口乱れ (k, omega) の組
  --reset-turb         --ic-from の場の k/omega を入口値に戻す (出発場の感度用。T3A では戻しても戻さなくても同じ解: run_0005 / run_0011)
段階 (一様場から): soft (1 次, cfl 0.5, nStepInner 10) → mid (1 次, cfl 1) → ramp (2 次, cfl 1/2) → main (2 次, --cfl)。
--ic-from のときは ramp (2 次, cfl 1) → main。段間は interp_field.py (同一の 2D 平面メッシュなので index 同値)。
条件の根拠は plans/active/turbulence-transition-lm2009.md §6。
"""
import argparse, os, shutil, subprocess, sys
from pathlib import Path
import h5py, numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "solver_density_cuda" / "tools"
HERE = Path(__file__).resolve().parent
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""),
           FORGE_CUDA_BLOCKSIZE="256", FORGE_BIN=str(ROOT / "solver_density_cuda/.build-native/relwithdebinfo/forge"))
sys.path.insert(0, str(TOOLS))

GAM, CP = 1.4, 1004.5; R = CP * (GAM - 1) / GAM
P_INF, T_INF, U_INF, MU = 101325.0, 300.0, 69.44, 2.269e-4        # M 0.2, Re/m = 3.6e5 (T3A 5.4 m/s / 1.5e-5 と同じ)
RO_INF = P_INF / (R * T_INF)
PT, TT = 104190.0, 302.4                                          # SU2 cfg と同じ入口全量
# 入口 (k, omega) と粘性。U は M 0.2 に固定し、Re/m を粘性で実験に合わせる (k = 1.5 (Tu U)^2, omega = rho k / (R_mu mu))。
INLET = {"t3a": (7.87, 3401.0),      # Tu 3.3 %, mu_t/mu 12  (Langtry-Menter 2009 表 4: 5.4 m/s, Re/m 3.6e5)
         "t3b": (30.56, 2759.0)}     # Tu 6.5 %, mu_t/mu 100 (同: 9.4 m/s, rho 1.2, mu 1.8e-5 -> Re/m 6.27e5)
MU_CASE = {"t3a": 2.269e-4, "t3b": 1.3033e-4}

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 0, visc: {mu}, thermCond: {kth:.6g}, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, dt_min: 1e-10, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}}}
{turb}
output: {{level: {olevel}}}
initial: "uniform_p101325_u10"
"""
BC = """inlet:  {{physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {{Pt: {pt}, Tt: {tt}, k: {k}, omega: {om}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p}, Pt: {pt}, Tt: {tt}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
wall:   {{physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
sym:    {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def turb_line(model, k, om):
    if model == "lam":
        return 'turbulence: {model: "none"}'
    tr = ', transition: "lm2009"' if model == "lm" else ""
    # SU2 (SST V2003m) と揃える素の SST: dilatationCorrection 0, katoLaunder 0 (平板では A/B で無影響 [recommended-settings §2])
    return (f'turbulence: {{model: "sst", scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: 0, wallTreatmentSST: 0, '
            f'turbulentPrandtl: 0.9, kInit: {k}, omegaInit: {om}{tr}}}')


def cfg(model, k, om, nsteps, cfl, relax, conv, lim, ninner, outint, olevel=1):
    return CFG.format(mu=MU, kth=MU * CP / 0.72, cp=CP, gam=GAM, nsteps=nsteps, cfl=cfl, relax=relax, conv=conv, lim=lim,
                      ninner=ninner, outint=outint, turb=turb_line(model, k, om), olevel=olevel)


def patch_ic(h5, k, om):
    """一様流 + 壁近傍の速度ランプ (幅 2 mm)。"""
    with h5py.File(h5, "r+") as f:
        wd = f["/VALUE/wall_dist"][:].astype(float)
        u = U_INF * np.clip(wd / 2.0e-3, 0.0, 1.0)
        ro = np.full_like(wd, RO_INF)
        vals = {"ro": ro, "roUx": ro * u, "roUy": 0 * ro, "roUz": 0 * ro, "roe": P_INF / (GAM - 1) + 0.5 * ro * u * u, "roK": ro * k, "roOmega": ro * om}
        for n, v in vals.items():
            f["/VALUE/" + n][:] = v.astype(np.float32)


def reset_turb(h5, k, om):
    with h5py.File(h5, "r+") as f:
        ro = f["/VALUE/ro"][:].astype(float)
        f["/VALUE/roK"][:] = (ro * k).astype(np.float32); f["/VALUE/roOmega"][:] = (ro * om).astype(np.float32)
        for n in ("roGamma", "roReth"):
            if "/VALUE/" + n in f: del f["/VALUE/" + n]


def run_forge(rd):
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def stage(rd, man, text, nsteps, tag):
    (rd / "solverConfig.yaml").write_text(text)
    rc = run_forge(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if rc != 0 or not res or int(res[-1].stem.split("_")[1]) < nsteps:
        raise SystemExit(f"stage {tag} failed rc={rc} res={[r.name for r in res][-2:]}")
    subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]), str(rd / "mesh.h5")], env=ENV, check=True, capture_output=True, text=True)
    for f in rd.glob("res_*"): f.unlink()
    for f in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
        if (rd / f).exists(): shutil.move(str(rd / f), str(rd / f"{Path(f).stem}_{tag}{Path(f).suffix}"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--mesh", default="plate_base")
    ap.add_argument("--model", required=True, choices=["sst", "lm", "lam"]); ap.add_argument("--case", default="t3a", choices=list(INLET))
    ap.add_argument("--main-steps", type=int, default=20000); ap.add_argument("--cfl", type=float, default=5.0); ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--out-int", type=int, default=2000); ap.add_argument("--ic-from", default=None); ap.add_argument("--reset-turb", action="store_true")
    ap.add_argument("--k", type=float, default=None); ap.add_argument("--omega", type=float, default=None)
    ap.add_argument("--level", type=int, default=1); ap.add_argument("--dry", action="store_true"); ap.add_argument("--no-ramp", action="store_true")
    a = ap.parse_args()
    global MU
    k, om = INLET[a.case]; MU = MU_CASE[a.case]
    if a.k is not None: k = a.k
    if a.omega is not None: om = a.omega
    rd = HERE / a.run
    if rd.exists(): raise SystemExit(f"{rd} exists")
    rd.mkdir()
    shutil.copy(HERE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    main_cfg = cfg(a.model, k, om, a.main_steps, a.cfl, a.relax, 1, 2, 4, a.out_int, a.level)
    (rd / "solverConfig.yaml").write_text(main_cfg)      # interp_field が DST の隣の config を照合するので先に置く
    (rd / "bcondConfig.yaml").write_text(BC.format(pt=PT, tt=TT, p=P_INF, k=k, om=om))
    (rd / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    if a.ic_from:
        src = sorted((HERE / a.ic_from).glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
        subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(src), str(rd / "mesh.h5")], env=ENV, check=True, capture_output=True, text=True)
        (rd / "CONTINUED_FROM").write_text(str(src.relative_to(ROOT)) + ("  (k/omega reset to inlet values)" if a.reset_turb else "") + "\n")
        if a.reset_turb: reset_turb(rd / "mesh.h5", k, om)
    else:
        patch_ic(rd / "mesh.h5", k, om)
    if a.dry: return
    from stage_manifest import StageManifest
    man = StageManifest(rd)
    def st(text, n, tag):
        man.add(tag, text, (rd / "bcondConfig.yaml").read_text()); man.write()
        stage(rd, man, text, n, tag)
    if not a.ic_from:
        st(cfg(a.model, k, om, 2000, 0.5, a.relax, 0, 0, 10, 2000), 2000, "soft")
        st(cfg(a.model, k, om, 2000, 1.0, a.relax, 0, 0, 10, 2000), 2000, "mid")
    if not a.no_ramp:
        for i, cv in enumerate((1.0, 2.0) if not a.ic_from else (1.0,)):
            st(cfg(a.model, k, om, 2000, cv, a.relax, 1, 2, 4, 2000), 2000, f"ramp{i}_cfl{cv:g}")
    (rd / "solverConfig.yaml").write_text(main_cfg)
    man.add("main", main_cfg, (rd / "bcondConfig.yaml").read_text(), history="residual_history.csv"); man.write()
    rc = run_forge(rd)
    print("main rc", rc); print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-800:])


if __name__ == "__main__":
    main()
