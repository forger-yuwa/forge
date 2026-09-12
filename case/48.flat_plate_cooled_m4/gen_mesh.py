#!/usr/bin/env python3
"""case/48 平板メッシュ生成: 第一セル高さ y1 [µm] ごとに .geo → .msh (msh4.1, forge) + .su2 (SU2) → node 変換 h5 + 品質判定。

usage: python3 gen_mesh.py [--y1 3 6 12 24] [--H 0.2] [--nx-plate 1000]
出力: mesh/fp_y1_<µm>um.{geo,msh,su2,h5} と mesh/quality_<µm>um.txt
node 変換は mesh/_conv/ の solverConfig.yaml (discretization: node, no-slip wall) を cwd にして行う (wall_dist は no-slip 壁から)。
"""
import argparse, math, os, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "solver_density_cuda" / "build"
TOOLS = ROOT / "solver_density_cuda" / "tools"
HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

CONV_CFG = """mesh: {meshFormat: "hdf5", discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}
gpu: 1
solver: "SLAU"
physProp: {isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {control: 0, nStepOuter: 10}
  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-9, dt_max: 1.0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 10
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 0, limiter: 0}
turbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0}
initial: "uniform_p101325_u10"
"""
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 0.062, Ux: 1413.0, Uy: 0.0, Uz: 0.0, Ps: 5037.4, k: 75.0, omega: 26000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 5037.4, Pt: 5037.4, Tt: 283.0}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
wall:   {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
sym:    {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
"""


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1  # 節点数 = セル数 + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--y1", type=float, nargs="+", default=[3.0, 6.0, 12.0, 24.0], help="第一セル高さ [µm]")
    ap.add_argument("--H", type=float, default=0.2)
    ap.add_argument("--nx-plate", type=int, default=1000)
    ap.add_argument("--nx-up", type=int, default=40)
    ap.add_argument("--r-plate", type=float, default=1.0025)
    ap.add_argument("--r-y", type=float, default=1.09)
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()
    tpl = (MESH / "flat_plate_cooled.geo.tpl").read_text()
    for y1um in a.y1:
        y1 = y1um * 1e-6
        # 等比列の第一セルを y1 に合わせる: 節点数 ny のとき Δ_1 = H (r-1)/(r^(ny-1) - 1)。ny を決めてから r を微調整
        ny = ny_for(y1, a.H, a.r_y)
        # r を二分法で y1 に合わせる
        lo, hi = 1.001, 1.5
        for _ in range(100):
            r = 0.5 * (lo + hi)
            d1 = a.H * (r - 1.0) / (r ** (ny - 1) - 1.0)
            if d1 > y1: lo = r
            else: hi = r
        tag = f"fp_y1_{y1um:g}um" + (f"_nx{a.nx_plate}" if a.nx_plate != 1000 else "")
        geo = tpl.replace("__H__", f"{a.H}").replace("__NX_UP__", str(a.nx_up)).replace("__NX_PLATE__", str(a.nx_plate + 1)) \
                 .replace("__R_PLATE__", f"{a.r_plate}").replace("__NY__", str(ny)).replace("__R_Y__", f"{r:.8f}")
        (MESH / f"{tag}.geo").write_text(geo)
        dx1 = 1.0 * (a.r_plate - 1.0) / (a.r_plate ** a.nx_plate - 1.0)
        dxN = dx1 * a.r_plate ** (a.nx_plate - 1)
        print(f"[{tag}] ny={ny} r_y={r:.5f} y1={a.H*(r-1)/(r**(ny-1)-1)*1e6:.3f}um  dx: {dx1*1e3:.3f}..{dxN*1e3:.3f} mm  AR_wall(end)~{dxN/y1:.0f}")
        subprocess.run(["gmsh", "-2", str(MESH / f"{tag}.geo"), "-o", str(MESH / f"{tag}.msh"), "-format", "msh41", "-v", "1"], check=True)
        subprocess.run(["gmsh", "-2", str(MESH / f"{tag}.geo"), "-o", str(MESH / f"{tag}.su2"), "-format", "su2", "-v", "1"], check=True)
        if a.no_convert:
            continue
        conv = MESH / "_conv"
        conv.mkdir(exist_ok=True)
        (conv / "solverConfig.yaml").write_text(CONV_CFG)
        (conv / "bcondConfig.yaml").write_text(CONV_BC)
        r = subprocess.run([str(BUILD / "convertGmshToForge"), str(MESH / f"{tag}.msh"), "m.h5"], cwd=conv, env=ENV,
                           capture_output=True, text=True)
        (conv / f"convert_{tag}.log").write_text(r.stdout + r.stderr)
        if not (conv / "m.h5").exists():
            print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit("convert failed")
        shutil.move(str(conv / "m.h5"), str(MESH / f"{tag}.h5"))
        for x in conv.glob("m.xmf"): x.unlink()
        q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(MESH / f"{tag}.h5")],
                           capture_output=True, text=True)
        (MESH / f"quality_{tag}.txt").write_text(q.stdout + q.stderr)
        print("   ", [l for l in q.stdout.splitlines() if "VERDICT" in l or "AR" in l][:3])


if __name__ == "__main__":
    main()
