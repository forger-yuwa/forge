#!/usr/bin/env python3
"""case/52 メッシュ生成 (case/49 条件の 2D 単一すきま) — 平面 2D 構造 (node 用)。T0 (すきま無し) と T1 (深キャビティ) を**同一ブロック構成**で作る。

幾何は geometry.json (NASA TN D-5908 図 1 / Table I–III) が正本。
  領域: x ∈ [x_in, x_out], y ∈ [0, H]。x<0 は slip 助走、x≥0 が平板 (前縁 x=0)。
  キャビティ: 前壁 x_r − w、後壁 x_r = 158.496 mm、深さ d = 20.32 mm。
  T0 は同じブロック分割のまま開口部を壁にする (分母診断を同一格子で取るため)。

usage: python3 gen_mesh.py --case T0            # すきま無し
       python3 gen_mesh.py --case T1 --wd 0.063 # w/d = 0.063 のキャビティ
"""
import argparse, json, math, os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = ROOT / "solver_density_cuda" / "build"
TOOLS = ROOT / "solver_density_cuda" / "tools"
MESH = HERE / "mesh"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

# 変換用の最小 config (wall_dist は no-slip 壁から作られるので、壁タグは最終形で変換する)
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
turbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0, kInit: 81.609, omegaInit: 50530.8}
initial: "uniform_p101325_u10"
"""
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 0.08803358, Ux: 1475.2116, Uy: 0.0, Uz: 0.0, Ps: 5474.7177, k: 81.609, omega: 50530.8}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 5474.7177, Pt: 5474.7177, Tt: 216.65}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
plate:  {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
runup:  {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
cavity: {physID: 6, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 500.0}}
"""


def geo_text(g, w, d, ny, r_y, nx_up, nx_pl, nx_gap, nx_dn, nx_buf, ny_cav,
             bump_pl, bump_dn, bump_cav, with_cavity, bump_gap=0.0):
    """B1 助走(slip) | B2 平板上流 | B3 開口上 | B4 平板下流 | B5 出口バッファ(slip)  (+ B6 キャビティ)

    平板は模型長 197 mm で終わり、その下流は slip 延長にする (出口境界を境界層で切らない)。
    """
    x_in, x_out, H = g["x_in"], g["x_out"], g["H"]
    xr = g["x_rear"]; xf = xr - w; xpe = g["x_plate_end"]
    L = []
    A = L.append
    A("// case/50 — NASA TN D-5908 深キャビティ (平面 2D, node)。gen_mesh.py が生成。")
    A("Geometry.PointNumbers = 0;  lc = 0.01;")
    A(f"x_in = {x_in:.9f}; xle = 0.0; xf = {xf:.9f}; xr = {xr:.9f}; xpe = {xpe:.9f}; x_out = {x_out:.9f};")
    A(f"H = {H:.9f}; d = {d:.9f};")
    names = ["x_in", "xle", "xf", "xr", "xpe", "x_out"]
    for i, xs in enumerate(names, start=1):          # 1..6 底辺
        A(f"Point({i}) = {{{xs}, 0.0, 0.0, lc}};")
    for i, xs in enumerate(names, start=7):          # 7..12 上辺
        A(f"Point({i}) = {{{xs}, H, 0.0, lc}};")
    for i in range(1, 6):                            # 1..5 底辺 (助走/平板上流/開口/平板下流/バッファ)
        A(f"Line({i}) = {{{i}, {i+1}}};")
    for i in range(6, 11):                           # 6..10 上辺
        A(f"Line({i}) = {{{i+1}, {i+2}}};")
    for i, (b, t) in enumerate(zip(range(1, 7), range(7, 13)), start=11):   # 11..16 縦線
        A(f"Line({i}) = {{{b}, {t}}};")
    A(f"Transfinite Line {{11, 12, 13, 14, 15, 16}} = {ny} Using Progression {r_y:.8f};")
    A(f"Transfinite Line {{1, 6}} = {nx_up} Using Progression 1.08;")
    A(f"Transfinite Line {{2, 7}} = {nx_pl} Using Bump {bump_pl};")
    A(f"Transfinite Line {{3, 8}} = {nx_gap}" + (f" Using Bump {bump_gap};" if bump_gap > 0 else ";"))
    A(f"Transfinite Line {{4, 9}} = {nx_dn} Using Bump {bump_dn};")
    A(f"Transfinite Line {{5, 10}} = {nx_buf} Using Progression 1.02;")
    for k in range(1, 6):
        lb, lr, lt, ll = k, 11 + k, 5 + k, 10 + k
        A(f"Curve Loop({k}) = {{{lb}, {lr}, -{lt}, -{ll}}};")
        A(f"Plane Surface({k}) = {{{k}}};  Transfinite Surface {{{k}}};  Recombine Surface {{{k}}};")
    if with_cavity:
        A(f"Point(13) = {{xf, -d, 0.0, lc}};")
        A(f"Point(14) = {{xr, -d, 0.0, lc}};")
        A("Line(17) = {13, 14};")   # 床
        A("Line(18) = {4, 14};")    # 後壁 (下向き)
        A("Line(19) = {3, 13};")    # 前壁 (下向き)
        A(f"Transfinite Line {{17}} = {nx_gap}" + (f" Using Bump {bump_gap};" if bump_gap > 0 else ";"))
        A(f"Transfinite Line {{18, 19}} = {ny_cav} Using Bump {bump_cav};")
        A("Curve Loop(6) = {19, 17, -18, -3};")
        A("Plane Surface(6) = {6};  Transfinite Surface {6};  Recombine Surface(6);")
    A('Physical Curve("inlet",  1) = {11};')
    A('Physical Curve("outlet", 2) = {16};')
    A('Physical Curve("top",    3) = {6, 7, 8, 9, 10};')
    if with_cavity:
        A('Physical Curve("plate",  4) = {2, 4};')
        A('Physical Curve("cavity", 6) = {17, 18, 19};')
        A('Physical Surface("fluid", 8) = {1, 2, 3, 4, 5, 6};')
    else:
        A('Physical Curve("plate",  4) = {2, 3, 4};')
        A('Physical Surface("fluid", 8) = {1, 2, 3, 4, 5};')
    A('Physical Curve("runup",  5) = {1, 5};')
    return "\n".join(L) + "\n"


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=["T0", "T1"], default="T0")
    ap.add_argument("--wd", type=float, default=0.05, help="w/d (T1 用。T0 でもブロック分割に使う)")
    ap.add_argument("--y1", type=float, default=4.0, help="壁第一セル [µm]")
    ap.add_argument("--r-y", type=float, default=1.06)
    ap.add_argument("--H", type=float, default=0.060)
    ap.add_argument("--x-in", type=float, default=-0.020)
    ap.add_argument("--x-out", type=float, default=0.280)
    ap.add_argument("--x-plate-end", type=float, default=0.250, help="模型長 (ここから先は slip バッファ)")
    ap.add_argument("--nx-buf", type=int, default=61)
    ap.add_argument("--nx-up", type=int, default=41)
    ap.add_argument("--nx-plate", type=int, default=601)
    ap.add_argument("--nx-gap", type=int, default=101)
    ap.add_argument("--nx-down", type=int, default=401)
    ap.add_argument("--ny-cav", type=int, default=301)
    ap.add_argument("--bump-plate", type=float, default=0.02)
    ap.add_argument("--bump-down", type=float, default=0.05)
    ap.add_argument("--bump-cav", type=float, default=0.02)
    ap.add_argument("--bump-gap", type=float, default=0.0,
                    help="開口幅方向を側壁に寄せる (すきま側壁の y1+ を下げる)。0 で一様")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()

    geom = json.loads((HERE / "geometry.json").read_text(encoding="utf-8"))
    d = geom["cavity"]["depth"] * 1e-3
    w = geom["cavity"]["widths"][f"{a.wd}"] * 1e-3
    g = dict(x_in=a.x_in, x_out=a.x_out, H=a.H, x_plate_end=a.x_plate_end,
             x_rear=geom["cavity"]["x_rear_wall_from_le"] * 1e-3)

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
    tag = a.tag or (f"{a.case.lower()}_wd{a.wd:g}_y{a.y1:g}um")
    MESH.mkdir(exist_ok=True)
    txt = geo_text(g, w, d, ny, r, a.nx_up, a.nx_plate, a.nx_gap, a.nx_down, a.nx_buf, a.ny_cav,
                   a.bump_plate, a.bump_down, a.bump_cav, with_cavity=(a.case == "T1"),
                   bump_gap=a.bump_gap)
    (MESH / f"{tag}.geo").write_text(txt)
    print(f"[{tag}] w = {w*1e3:.3f} mm, d = {d*1e3:.2f} mm, 開口 {a.nx_gap-1} セル ({w*1e3/(a.nx_gap-1)*1e3:.1f} µm/セル)")
    print(f"        ny = {ny} (y1 = {a.H*(r-1)/(r**(ny-1)-1)*1e6:.3f} µm, r = {r:.5f})")
    subprocess.run(["gmsh", "-2", str(MESH / f"{tag}.geo"), "-o", str(MESH / f"{tag}.msh"),
                    "-format", "msh41", "-v", "1"], check=True)
    if a.no_convert:
        return
    conv = MESH / "_conv"
    conv.mkdir(exist_ok=True)
    (conv / "solverConfig.yaml").write_text(CONV_CFG)
    (conv / "bcondConfig.yaml").write_text(CONV_BC)
    r2 = subprocess.run([str(BUILD / "convertGmshToForge"), str(MESH / f"{tag}.msh"), "m.h5"],
                        cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / f"convert_{tag}.log").write_text(r2.stdout + r2.stderr)
    if not (conv / "m.h5").exists():
        print(r2.stdout[-3000:], r2.stderr[-2000:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{tag}.h5"))
    for x in conv.glob("m.xmf"):
        x.unlink()
    q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(MESH / f"{tag}.h5")],
                       capture_output=True, text=True)
    (MESH / f"quality_{tag}.txt").write_text(q.stdout + q.stderr)
    print("   ", [l for l in q.stdout.splitlines() if "VERDICT" in l][:2])


if __name__ == "__main__":
    main()
