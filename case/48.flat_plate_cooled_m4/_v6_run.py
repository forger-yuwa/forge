#!/usr/bin/env python3
"""V6 (#7): 周期・軸対称で面流束の A/B を回す。plan convection-slau-wall-normal-chi §6 V6。

V5 (`_v5_run.py`) と同じ 4 本 (flag0 x2 / flag1 / flag1 bs256) を 1 step ずつ。
**メッシュと bcond は外から渡す** (周期 = 低解像度チャネル、軸対称 = 最小直管)。
判定は `_v5_compare.py` を共用する (V5 と V6 は同じツール鎖)。

usage: _v6_run.py --mesh M.h5 --bcond B.yaml --out DIR [--axisymmetric] [--period dx,dy,dz]
"""
import argparse, shutil, subprocess, sys, os
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1{axi}, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True); ap.add_argument("--bcond", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--axisymmetric", action="store_true")
    ap.add_argument("--period", default="")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    axi = ", isAxisymmetric: 1" if a.axisymmetric else ""
    for name, flag, bs in (("flag0_a", 0, 128), ("flag0_b", 0, 128), ("flag1", 1, 128), ("flag1_bs256", 1, 256)):
        rd = out / name
        if rd.exists(): shutil.rmtree(rd)
        rd.mkdir(parents=True)
        shutil.copy(a.mesh, rd / "mesh.h5"); shutil.copy(a.bcond, rd / "bcondConfig.yaml")
        cmd = [sys.executable, str(HERE / "_v5_make_state.py"), str(rd / "mesh.h5")]
        if a.period: cmd += ["--period", a.period]
        if a.axisymmetric: cmd += ["--axisymmetric"]
        subprocess.run(cmd, check=True)
        (rd / "solverConfig.yaml").write_text(CFG.format(flag=flag, axi=axi))
        (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
        dump = out / f"{name}.massflux.bin"
        env = dict(ENV, FORGE_CUDA_BLOCKSIZE=str(bs), FORGE_DUMP_MASSFLUX=str(dump))
        r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=env, capture_output=True, text=True)
        (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
        ok = dump.exists()
        print(f"{name}: rc={r.returncode} dump={ok}")
        if not ok:
            print((rd / "run_case_stdout.log").read_text()[-1500:]); raise SystemExit(1)

if __name__ == "__main__":
    main()
