#!/usr/bin/env python3
"""V1 (共役スラブ) の一様メッシュ生成: .geo -> msh4.1 -> forge h5 (node) -> 静止 IC をパッチ。

流体層 H [m] を ny 等分、幅 W を nx 等分した構造квад。壁近傍クラスタは**入れない**:
純伝導の緩和は最小セル幅で律速され、クラスタメッシュだと定常解に到達するのに 1e5 step 級かかる
(実測: 壁クラスタ d1=8.5e-5 のメッシュは 40000 step でも q が定常値の 6 倍)。

usage: python3 gen_mesh.py [--nx 20] [--ny 16] [--H 0.01] [--W 0.02]
出力: template/mesh.h5 (+ mesh/slab.geo, mesh/slab.msh)
"""
import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = ROOT / "solver_density_cuda" / ".build-native" / "relwithdebinfo"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

GEO = """// V1 共役スラブ: 一様構造 quad (physID 1 左 / 2 右 / 3 下=共役壁 / 4 上=等温壁)
W = {W}; H = {H}; nx = {nx}; ny = {ny};
Point(1) = {{0, 0, 0, 1}};  Point(2) = {{W, 0, 0, 1}};
Point(3) = {{W, H, 0, 1}};  Point(4) = {{0, H, 0, 1}};
Line(1) = {{1, 2}}; Line(2) = {{2, 3}}; Line(3) = {{3, 4}}; Line(4) = {{4, 1}};
Line Loop(1) = {{1, 2, 3, 4}}; Plane Surface(1) = {{1}};
Transfinite Line{{1, 3}} = nx + 1;
Transfinite Line{{2, 4}} = ny + 1;
Transfinite Surface{{1}}; Recombine Surface{{1}};
Physical Curve("side_left",  1) = {{4}};
Physical Curve("side_right", 2) = {{2}};
Physical Curve("wall_bot",   3) = {{1}};
Physical Curve("wall_top",   4) = {{3}};
Physical Surface("fluid",    5) = {{1}};
Mesh.MshFileVersion = 4.1;
"""

CONV_CFG = """mesh: {discretization: "node", meshFileName: "slab.h5", valueFileName: "slab.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 0, visc: 1.81e-5, thermCond: 0.0241, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 1}
  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-11, dt_max: 1.0}
  outStepStart: 0
  outStepInterval: 1
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 0, limiter: 0}
turbulence: {model: "none"}
initial: "uniform_p101325_u10"
"""

CONV_BC = """side_left:  {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }
side_right: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }
wall_bot:   {physID: 3, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 300.0}}
wall_top:   {physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 350.0}}
"""


def patch_static_ic(h5path: Path, p0: float, T0: float, cp: float, gamma: float):
    """静止・一様 (p0, T0) の IC を VALUE に書く。roe = ro*cv*T (ek=0)。"""
    R = cp * (gamma - 1.0) / gamma
    cv = cp / gamma
    ro = p0 / (R * T0)
    with h5py.File(h5path, "a") as f:
        n = f["CELLS/centCoords"].shape[0] // 3
        g = f.require_group("VALUE")
        for name, val in (("ro", ro), ("roUx", 0.0), ("roUy", 0.0), ("roUz", 0.0),
                          ("roe", ro * cv * T0)):
            if name in g:
                del g[name]
            g.create_dataset(name, data=np.full(n, val, dtype=np.float32))
    return ro


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=20)
    ap.add_argument("--ny", type=int, default=16)
    ap.add_argument("--H", type=float, default=0.01)
    ap.add_argument("--W", type=float, default=0.02)
    ap.add_argument("--p0", type=float, default=1013.25,
                    help="静止圧 [Pa]。低圧ほど熱拡散率 k/(rho cp) が大きく、擬似時間の緩和が速い "
                         "(解析解 q = k dT/dy は圧力に依らないので検証は成立する)")
    a = ap.parse_args()

    mdir = HERE / "mesh"; mdir.mkdir(exist_ok=True)
    (mdir / "slab.geo").write_text(GEO.format(W=a.W, H=a.H, nx=a.nx, ny=a.ny))
    subprocess.run(["gmsh", "-2", "slab.geo", "-o", "slab.msh", "-format", "msh41", "-v", "1"],
                   cwd=mdir, check=True, stdout=subprocess.DEVNULL)
    (mdir / "solverConfig.yaml").write_text(CONV_CFG)
    (mdir / "bcondConfig.yaml").write_text(CONV_BC)
    subprocess.run([str(BUILD / "convertGmshToForge"), "slab.msh", "slab.h5"], cwd=mdir, check=True,
                   env=ENV, stdout=subprocess.DEVNULL)
    out = mdir / "slab.h5"
    if not out.exists():
        sys.exit(f"converter did not produce {out}")
    ro = patch_static_ic(out, a.p0, 300.0, 1004.5, 1.4)
    tpl = HERE / "template"; tpl.mkdir(exist_ok=True)
    (tpl / "mesh.h5").write_bytes(out.read_bytes())
    print(f"[gen_mesh] {a.nx}x{a.ny} uniform quads, H={a.H} W={a.W}, p0={a.p0} Pa, static IC ro={ro:.6f} kg/m3 "
          f"-> {tpl/'mesh.h5'}")
    subprocess.run([sys.executable, str(ROOT / "solver_density_cuda" / "tools" / "check_mesh_quality.py"),
                    str(tpl / "mesh.h5")], env=ENV)


if __name__ == "__main__":
    main()
