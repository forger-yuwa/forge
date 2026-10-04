#!/usr/bin/env python3
"""case/57 遷移平板メッシュ: .geo → .msh (forge) + .su2 (SU2、同一節点) → node 変換 h5 + 品質判定。

usage: python3 gen_mesh.py [--tag base] [--nx-plate 300] [--y1 20] [--r-y 1.1]
"""
import argparse, math, os, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
CONVERT = ROOT / "solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge"
TOOLS = ROOT / "solver_density_cuda/tools"
HERE = Path(__file__).resolve().parent; MESH = HERE / "mesh"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""), FORGE_CUDA_BLOCKSIZE="256")
CONV_CFG = """mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 0, visc: 2.269e-4, thermCond: 0.3165, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 10}
  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-9, dt_max: 1.0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 10
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 0, limiter: 0}
turbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0}
initial: "uniform_p101325_u10"
"""
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 1.1768, Ux: 69.44, Uy: 0.0, Uz: 0.0, Ps: 101325.0, k: 7.87, omega: 3401.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 101325.0, Pt: 104000.0, Tt: 302.4}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
wall:   {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
sym:    {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="base"); ap.add_argument("--L", type=float, default=1.5); ap.add_argument("--x-up", type=float, default=0.04)
    ap.add_argument("--H", type=float, default=0.3); ap.add_argument("--nx-plate", type=int, default=300); ap.add_argument("--r-plate", type=float, default=1.006)
    ap.add_argument("--nx-up", type=int, default=24); ap.add_argument("--r-up", type=float, default=1.10)
    ap.add_argument("--y1", type=float, default=20.0, help="第一セル高さ [µm]"); ap.add_argument("--r-y", type=float, default=1.1)
    a = ap.parse_args(); y1 = a.y1 * 1e-6
    ny = int(math.ceil(math.log(1.0 + a.H * (a.r_y - 1.0) / y1) / math.log(a.r_y))) + 1
    lo, hi = 1.001, 1.5
    for _ in range(100):
        r = 0.5 * (lo + hi); d1 = a.H * (r - 1.0) / (r ** (ny - 1) - 1.0)
        if d1 > y1: lo = r
        else: hi = r
    geo = (MESH / "plate.geo.tpl").read_text()
    for k, v in (("__X_UP__", a.x_up), ("__L__", a.L), ("__H__", a.H), ("__NX_UP__", a.nx_up + 1), ("__R_UP__", a.r_up),
                 ("__NX_PLATE__", a.nx_plate + 1), ("__R_PLATE__", a.r_plate), ("__NY__", ny), ("__R_Y__", f"{r:.8f}")):
        geo = geo.replace(k, str(v))
    tag = f"plate_{a.tag}"; (MESH / f"{tag}.geo").write_text(geo)
    dx1 = a.L * (a.r_plate - 1.0) / (a.r_plate ** a.nx_plate - 1.0); dxN = dx1 * a.r_plate ** (a.nx_plate - 1)
    print(f"[{tag}] ny={ny} r_y={r:.5f} y1={y1*1e6:.2f}um  dx {dx1*1e3:.2f}..{dxN*1e3:.2f} mm  AR(end)~{dxN/y1:.0f}")
    for fmt, ext in (("msh41", "msh"), ("su2", "su2")):
        subprocess.run(["gmsh", "-2", str(MESH / f"{tag}.geo"), "-o", str(MESH / f"{tag}.{ext}"), "-format", fmt, "-v", "1"], check=True)
    conv = MESH / "_conv"; conv.mkdir(exist_ok=True)
    (conv / "solverConfig.yaml").write_text(CONV_CFG); (conv / "bcondConfig.yaml").write_text(CONV_BC)
    rr = subprocess.run([str(CONVERT), str(MESH / f"{tag}.msh"), "m.h5"], cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / f"convert_{tag}.log").write_text(rr.stdout + rr.stderr)
    if not (conv / "m.h5").exists():
        print(rr.stdout[-1500:], rr.stderr[-1500:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{tag}.h5"))
    q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(MESH / f"{tag}.h5")], capture_output=True, text=True)
    (MESH / f"quality_{tag}.txt").write_text(q.stdout + q.stderr)
    print("   ", [l for l in q.stdout.splitlines() if "VERDICT" in l][:2])


if __name__ == "__main__":
    main()
