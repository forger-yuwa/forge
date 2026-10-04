#!/usr/bin/env python3
"""V5 (#10d): 面流束ダンプの A/B を回す。plan convection-slau-wall-normal-chi §6 V5。"""
import shutil, subprocess, sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_runs as G

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh" / "fp_y1_12um.h5"
OUT = HERE / "_v5"
OUT.mkdir(exist_ok=True)

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: 1}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 0.2, cfl_pseudo: 0.2, implicitRelax: 1.0, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-10, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: 1
  timeIntegration: 11
  nStepInner: 1
space: {{convMethod: 0, limiter: 0, reconT: 0, lowMachThornber: 0, slauContactFloor: 0.0, slauWallNormalChi: {flag}}}
turbulence: {{model: "none"}}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""

def make(name, flag):
    rd = OUT / name
    if rd.exists(): shutil.rmtree(rd)
    rd.mkdir(parents=True)
    shutil.copy(MESH, rd / "mesh.h5")
    subprocess.run([sys.executable, str(HERE / "_v5_make_state.py"), str(rd / "mesh.h5")], check=True)
    (rd / "solverConfig.yaml").write_text(CFG.format(flag=flag))
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    (rd / "bcondConfig.yaml").write_text(G.BC.format(ro=G.RO_INF, u=G.U_INF, p=G.P_INF, t=G.T_INF,
                                                    k=G.K_INF, om=G.OM_INF, wallline=G.wall_line("adiabatic")))
    return rd

def run(rd, bs, dump):
    env = dict(G.ENV, FORGE_CUDA_BLOCKSIZE=str(bs), FORGE_DUMP_MASSFLUX=str(dump))
    r = subprocess.run([str(G.TOOLS / "run_case.sh"), str(rd)], env=env, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    ok = Path(dump).exists()
    print(f"{rd.name}: rc={r.returncode} dump={ok} ({Path(dump).stat().st_size if ok else 0} B)")
    if not ok:
        print((rd / "run_case_stdout.log").read_text()[-1200:]); raise SystemExit(1)

if __name__ == "__main__":
    for name, flag, bs in (("v5_flag0_a", 0, 128), ("v5_flag0_b", 0, 128),
                           ("v5_flag1", 1, 128), ("v5_flag1_bs256", 1, 256)):
        rd = make(name, flag)
        run(rd, bs, OUT / f"{name}.massflux.bin")
