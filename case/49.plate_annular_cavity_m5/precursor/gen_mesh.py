#!/usr/bin/env python3
"""case/49 前駆 2D 平板メッシュ生成 (.geo -> msh4.1 -> node 変換 h5 -> 品質判定)。

usage: python3 gen_mesh.py [--y1 6] [--H 0.05] [--x0 -0.05] [--x1 0.0] [--x2 0.40] [--nx-plate 800]
出力: mesh/fp_y1_<µm>um.{geo,msh,h5} と mesh/quality_<tag>.txt

node 変換は mesh/_conv/ の solverConfig.yaml (discretization: node, **no-slip wall**) を cwd にして行う
(SST の wall_dist は no-slip 壁から作られるため。procedures/recommended-settings.md §1.1)。
派生元: case/48 gen_mesh.py。
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
BUILD = ROOT / "solver_density_cuda" / "build"
TOOLS = ROOT / "solver_density_cuda" / "tools"
MESH = HERE / "mesh"
ENV = dict(os.environ,
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

sys.path.insert(0, str(CASE))
from setup import load as load_conditions  # noqa: E402

CONV_CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1, visc: {mu:.6g}, thermCond: {kcond:.6g}, thermCondMethod: 1, prandtlLam: {pr}, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: 10}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-9, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: 10
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 0, limiter: 0, pRef: {pref:.6g}}}
turbulence: {{model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0}}
initial: "uniform_p101325_u10"
"""
CONV_BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.6g}, Ux: {u:.6g}, Uy: 0.0, Uz: 0.0, Ps: {p:.6g}, k: {k:.6g}, omega: {om:.6g}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6g}, Pt: {p:.6g}, Tt: {t:.6g}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
wall:   {{physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
sym:    {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--y1", type=float, nargs="+", default=[6.0], help="第一セル高さ [µm]")
    ap.add_argument("--H", type=float, default=0.05, help="領域高さ [m]")
    ap.add_argument("--x0", type=float, default=-0.05, help="入口 x [m] (slip 助走の始端)")
    ap.add_argument("--x1", type=float, default=0.0, help="平板前縁 x [m]")
    ap.add_argument("--x2", type=float, default=0.40, help="出口 x [m]")
    ap.add_argument("--nx-plate", type=int, default=800)
    ap.add_argument("--nx-up", type=int, default=40)
    ap.add_argument("--r-plate", type=float, default=1.0035)
    ap.add_argument("--r-y", type=float, default=1.09)
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()

    d = load_conditions()
    MESH.mkdir(exist_ok=True)
    tpl = (MESH / "flat_plate.geo.tpl").read_text()
    for y1um in a.y1:
        y1 = y1um * 1e-6
        ny = ny_for(y1, a.H, a.r_y)
        lo, hi = 1.001, 1.5                      # 節点数 ny で第一セルを y1 に合わせる比を二分法で
        for _ in range(100):
            r = 0.5 * (lo + hi)
            d1 = a.H * (r - 1.0) / (r ** (ny - 1) - 1.0)
            if d1 > y1:
                lo = r
            else:
                hi = r
        tag = "fp_y1_%gum" % y1um
        geo = (tpl.replace("__X0__", repr(a.x0)).replace("__X1__", repr(a.x1))
                  .replace("__X2__", repr(a.x2)).replace("__H__", repr(a.H))
                  .replace("__NX_UP__", str(a.nx_up)).replace("__NX_PLATE__", str(a.nx_plate + 1))
                  .replace("__R_PLATE__", "%.6f" % a.r_plate)
                  .replace("__NY__", str(ny)).replace("__R_Y__", "%.8f" % r))
        (MESH / (tag + ".geo")).write_text(geo)
        dx1 = (a.x2 - a.x1) * (a.r_plate - 1.0) / (a.r_plate ** a.nx_plate - 1.0)
        dxN = dx1 * a.r_plate ** (a.nx_plate - 1)
        print("[%s] ny=%d r_y=%.5f y1=%.3f um  dx %.3f..%.3f mm  AR_wall(end)~%.0f"
              % (tag, ny, r, a.H * (r - 1) / (r ** (ny - 1) - 1) * 1e6, dx1 * 1e3, dxN * 1e3, dxN / y1))
        subprocess.run(["gmsh", "-2", str(MESH / (tag + ".geo")), "-o", str(MESH / (tag + ".msh")),
                        "-format", "msh41", "-v", "1"], check=True)
        if a.no_convert:
            continue
        conv = MESH / "_conv"
        conv.mkdir(exist_ok=True)
        (conv / "solverConfig.yaml").write_text(CONV_CFG.format(
            pr=d["prandtl_lam"], mu=d["mu_inf"], kcond=d["mu_inf"]*d["cp"]/d["prandtl_lam"], cp=d["cp"], gam=d["gamma"], pref=d["P_inf"]))
        (conv / "bcondConfig.yaml").write_text(CONV_BC.format(
            ro=d["ro_inf"], u=d["U_inf"], p=d["P_inf"], t=d["T_inf"], k=d["k_inf"], om=d["omega_inf"]))
        r2 = subprocess.run([str(BUILD / "convertGmshToForge"), str(MESH / (tag + ".msh")), "m.h5"],
                            cwd=conv, env=ENV, capture_output=True, text=True)
        (conv / ("convert_%s.log" % tag)).write_text(r2.stdout + r2.stderr)
        if not (conv / "m.h5").exists():
            print(r2.stdout[-2000:], r2.stderr[-2000:])
            sys.exit("convert failed")
        shutil.move(str(conv / "m.h5"), str(MESH / (tag + ".h5")))
        for x in conv.glob("m.xmf"):
            x.unlink()
        q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(MESH / (tag + ".h5"))],
                           capture_output=True, text=True)
        (MESH / ("quality_%s.txt" % tag)).write_text(q.stdout + q.stderr)
        print("   ", [l for l in q.stdout.splitlines() if "VERDICT" in l or "AR " in l][:3])


if __name__ == "__main__":
    main()
