#!/usr/bin/env python3
"""case/56 メッシュ生成 — 較正パネル相当の 2 次元平板 (平面 2D, node)。

Gate B で条件は固めた。ここは **forge 側で分母 $q_{FP}$ を作るための平板**。
構成: B1 助走 (slip) | B2 平板 (等温壁) | B3 出口バッファ (slip)。
平板を出口境界で切らないよう、下流に slip 延長を置く。

乱流 run はトリップが前縁から 13 cm なので、**x=0 をトリップ位置**に取り、
位置 I / II は x = 1.17-0.13 = 1.04 m / 1.88-0.13 = 1.75 m で評価する
(層流 run は実前縁なのでオフセット 0)。
"""
import argparse, json, math, os, shutil, subprocess, sys
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh"
BUILD = HERE.parents[1] / "solver_density_cuda" / "build"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")

CONV_CFG = """mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
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
turbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0, kInit: 100.0, omegaInit: 50000.0}
initial: "uniform_p101325_u10"
"""
# 壁は no-slip で変換する (wall_dist が no-slip 壁から作られるため)
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 0.0284, Ux: 2039.8, Uy: 0.0, Uz: 0.0, Ps: 1738.2, k: 100.0, omega: 50000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 1738.2, Pt: 1738.2, Tt: 206.9}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
plate:  {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
slip:   {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
"""


def geo_text(x_in, x_out, x_plate_end, H, ny, r_y, nx_up, nx_pl, nx_buf, bump_pl, r_up=1.06, r_buf=1.03):
    L = []; A = L.append
    A("// case/56 — TP-1187 較正パネル相当の 2D 平板 (平面 2D, node)。gen_mesh.py が生成。")
    A("Geometry.PointNumbers = 0;  lc = 0.05;")
    A(f"x_in = {x_in:.9f}; xle = 0.0; xpe = {x_plate_end:.9f}; x_out = {x_out:.9f}; H = {H:.9f};")
    for i, xs in enumerate(("x_in", "xle", "xpe", "x_out"), start=1):
        A(f"Point({i}) = {{{xs}, 0.0, 0.0, lc}};")
    for i, xs in enumerate(("x_in", "xle", "xpe", "x_out"), start=5):
        A(f"Point({i}) = {{{xs}, H, 0.0, lc}};")
    for i in range(1, 4):
        A(f"Line({i}) = {{{i}, {i+1}}};")            # 1..3 底辺
    for i in range(4, 7):
        A(f"Line({i}) = {{{i+1}, {i+2}}};")          # 4..6 上辺
    for i, (b, t) in enumerate(zip(range(1, 5), range(5, 9)), start=7):
        A(f"Line({i}) = {{{b}, {t}}};")              # 7..10 縦線
    A(f"Transfinite Line {{7, 8, 9, 10}} = {ny} Using Progression {r_y:.8f};")
    A(f"Transfinite Line {{1, 4}} = {nx_up} Using Progression {r_up};")
    A(f"Transfinite Line {{2, 5}} = {nx_pl} Using Bump {bump_pl};")
    A(f"Transfinite Line {{3, 6}} = {nx_buf} Using Progression {r_buf};")
    for k in range(1, 4):
        A(f"Curve Loop({k}) = {{{k}, {7+k}, -{3+k}, -{6+k}}};")
        A(f"Plane Surface({k}) = {{{k}}};  Transfinite Surface {{{k}}};  Recombine Surface {{{k}}};")
    A('Physical Curve("inlet",  1) = {7};')
    A('Physical Curve("outlet", 2) = {10};')
    A('Physical Curve("top",    3) = {4, 5, 6};')
    A('Physical Curve("plate",  4) = {2};')
    A('Physical Curve("slip",   5) = {1, 3};')
    A('Physical Surface("fluid", 8) = {1, 2, 3};')
    return "\n".join(L) + "\n"


def geo_text_ylist(x_in, x_out, x_plate_end, ys, nx_up, nx_pl, nx_buf, bump_pl):
    """壁法線の節点座標列 `ys` (0 から H まで昇順) をそのまま使う版。

    底辺 3 本を transfinite で割り (既定版と同じ x 分布)、`Extrude ... Layers` の
    累積高さで y 方向に押し出す。2026-09-26 の平板 A/B (plan #60、codex diagnose 2 回目) で、
    3D 格子の入口列と同じ y 配列を 2D 平板に与えるために足した。"""
    ys = np.asarray(ys, dtype=float)
    H = float(ys[-1])
    L = []; A = L.append
    A("// case/56 — 2D 平板 (y 節点列指定)。gen_mesh.py --y-file が生成。")
    A("Geometry.PointNumbers = 0;  lc = 0.05;")
    A(f"x_in = {x_in:.9f}; xle = 0.0; xpe = {x_plate_end:.9f}; x_out = {x_out:.9f};")
    for i, xs in enumerate(("x_in", "xle", "xpe", "x_out"), start=1):
        A(f"Point({i}) = {{{xs}, 0.0, 0.0, lc}};")
    for i in range(1, 4):
        A(f"Line({i}) = {{{i}, {i+1}}};")
    A(f"Transfinite Line {{1}} = {nx_up} Using Progression 1.06;")
    A(f"Transfinite Line {{2}} = {nx_pl} Using Bump {bump_pl};")
    A(f"Transfinite Line {{3}} = {nx_buf} Using Progression 1.03;")
    n = len(ys) - 1
    ones = ",".join(["1"] * n)
    hs = ",".join(f"{v / H:.12f}" for v in ys[1:])
    A(f"e[] = Extrude {{0, {H:.12f}, 0}} {{ Line{{1, 2, 3}}; Layers{{ {{{ones}}}, {{{hs}}} }}; Recombine; }};")
    # e[] は各線ごとに [上辺, 面, 側線(終点側), 側線(始点側, 向き負)] の 4 つ (gmsh 4.x で Printf 確認)。
    # 隣り合う線は側線を共有するので、入口 = 線 1 の始点側、出口 = 線 3 の終点側
    A('Physical Curve("inlet",  1) = {Abs(e[3])};')
    A('Physical Curve("outlet", 2) = {e[10]};')
    A('Physical Curve("top",    3) = {e[0], e[4], e[8]};')
    A('Physical Curve("plate",  4) = {2};')
    A('Physical Curve("slip",   5) = {1, 3};')
    A('Physical Surface("fluid", 8) = {e[1], e[5], e[9]};')
    return "\n".join(L) + "\n", H


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--y1", type=float, default=8.0, help="壁第一セル [µm]")
    ap.add_argument("--r-y", type=float, default=1.06)
    ap.add_argument("--H", type=float, default=0.35, help="前縁マッハ波 (8.2 deg) が出口から抜ける高さ")
    ap.add_argument("--x-in", type=float, default=-0.10)
    ap.add_argument("--x-plate-end", type=float, default=1.90)
    ap.add_argument("--x-out", type=float, default=2.10)
    ap.add_argument("--nx-up", type=int, default=41)
    ap.add_argument("--nx-plate", type=int, default=901)
    ap.add_argument("--nx-buf", type=int, default=41)
    ap.add_argument("--bump-plate", type=float, default=0.15)
    ap.add_argument("--r-up", type=float, default=1.06, help="助走区間の等比 (既定は従来値)")
    ap.add_argument("--r-buf", type=float, default=1.03, help="出口バッファの等比 (既定は従来値)")
    ap.add_argument("--tag", default="fp")
    ap.add_argument("--no-convert", action="store_true")
    ap.add_argument("--y-file", default=None,
                    help="壁法線の節点 y 座標列 (1 行 1 値、0 から上端まで昇順)。指定時は --y1/--r-y/--H を使わない")
    a = ap.parse_args()

    MESH.mkdir(exist_ok=True)
    if a.y_file:
        ys = np.loadtxt(a.y_file)
        if ys[0] != 0.0 or np.any(np.diff(ys) <= 0):
            sys.exit(f"{a.y_file}: 0 から始まる狭義単調増加の列でない")
        txt, H = geo_text_ylist(a.x_in, a.x_out, a.x_plate_end, ys,
                                a.nx_up, a.nx_plate, a.nx_buf, a.bump_plate)
        (MESH / f"{a.tag}.geo").write_text(txt)
        print(f"[{a.tag}] y 節点列 {a.y_file}: ny = {len(ys)}, y1 = {ys[1]*1e6:.3f} µm, H = {H*1e2:.1f} cm")
    else:
        make_progression(a)
    subprocess.run(["gmsh", "-2", str(MESH / f"{a.tag}.geo"), "-o", str(MESH / f"{a.tag}.msh"),
                    "-format", "msh41", "-v", "1"], check=True)
    convert(a)


def make_progression(a):
    y1 = a.y1 * 1e-6
    ny = ny_for(y1, a.H, a.r_y)
    lo, hi = 1.001, 1.5
    for _ in range(200):
        r = 0.5 * (lo + hi)
        d1 = a.H * (r - 1.0) / (r ** (ny - 1) - 1.0)
        if d1 > y1:
            lo = r
        else:
            hi = r
    txt = geo_text(a.x_in, a.x_out, a.x_plate_end, a.H, ny, r,
                   a.nx_up, a.nx_plate, a.nx_buf, a.bump_plate, a.r_up, a.r_buf)
    (MESH / f"{a.tag}.geo").write_text(txt)
    print(f"[{a.tag}] 平板 {a.x_plate_end*1e2:.0f} cm, H = {a.H*1e2:.0f} cm")
    print(f"        ny = {ny} (y1 = {a.H*(r-1)/(r**(ny-1)-1)*1e6:.3f} µm, r = {r:.5f}), "
          f"節点 ~{(a.nx_up+a.nx_plate+a.nx_buf-2)*ny/1000:.0f}k")


def convert(a):
    if a.no_convert:
        return
    conv = MESH / "_conv"; conv.mkdir(exist_ok=True)
    (conv / "solverConfig.yaml").write_text(CONV_CFG)
    (conv / "bcondConfig.yaml").write_text(CONV_BC)
    r2 = subprocess.run([str(BUILD / "convertGmshToForge"), str(MESH / f"{a.tag}.msh"), "m.h5"],
                        cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / f"convert_{a.tag}.log").write_text(r2.stdout + r2.stderr)
    if not (conv / "m.h5").exists():
        print(r2.stdout[-3000:], r2.stderr[-2000:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{a.tag}.h5"))
    q = subprocess.run([sys.executable, str(HERE.parents[1] / "solver_density_cuda" / "tools"
                                            / "check_mesh_quality.py"), str(MESH / f"{a.tag}.h5")],
                       capture_output=True, text=True)
    print([l for l in q.stdout.splitlines() if "VERDICT" in l] or q.stdout[-400:])


if __name__ == "__main__":
    main()
