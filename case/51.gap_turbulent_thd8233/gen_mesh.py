#!/usr/bin/env python3
"""case/51 (TN D-8233) の 2D メッシュ生成 — 平板 + 単独横すきま、node 用の全四角構造格子。

`case/55/gen_mesh.py` のブロック構成を踏襲する (あちらは skew 0.000 / AR 495 で
`VERDICT: PASS` の実績)。違いは 2 点:

  1. **前縁を持たない**。流入 BL は δ* 12.8 cm・δ 22.9 cm とトンネル壁 BL 並に厚く、
     CFD で育てられない。`bl_derived.json` の再構成プロファイルを `inletProfile` で
     入口に与えるので、底面は入口から no-slip 平板にする (case/55 の slip 助走は無い)。
  2. **上面 slip をトンネル中心線に置く** (H = 787.4/2 mm)。δ = 229 mm に対して
     自由流コアが 165 mm 残る。

**`case/56/gen_mesh_gap4.py` の欠陥を持ち込まないこと**: あちらは VOUT ブロックの
対辺の点数が合っておらず (下辺 105 点 vs 上辺 25 点)、gmsh が transfinite を諦めて
非構造再結合に落ち、skew 0.992 で `VERDICT: FAIL` になっていた。ここでは全ブロックが
対辺同数の矩形 transfinite なので、生成後に必ず `check_mesh_quality.py` で確認する。

ブロック:  B1 平板 (調整) | B2 平板 (細分) | B3 開口上 | B4 平板下流 | B5 slip バッファ
           + B6 すきま (前壁 / 床 / 後壁)

usage: python3 gen_mesh.py                      # 既定 (W 2.29 mm, D 45.72 mm)
       python3 gen_mesh.py --nx-gap 51 --ny-cav 151 --tag t2p_coarse   # 格子感度の粗い方
"""
import argparse, json, math, os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = ROOT / "solver_density_cuda" / "build"
TOOLS = ROOT / "solver_density_cuda" / "tools"
MESH = HERE / "mesh"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

# 変換用の最小 config。wall_dist は no-slip 壁から作られるので壁タグは最終形で変換する
# (memory: sst-mesh-walldist-gotcha)。値は変換にしか使わない。
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
# Re' 1.47e6/m 系列 (derived.json の Gate A)。T_w = 300 K。
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 0.0037703, Ux: 1498.69, Uy: 0.0, Uz: 0.0, Ps: 56.96, k: 100.0, omega: 50000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 56.96, Pt: 56.96, Tt: 1112.02}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
plate:  {physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 300.0}}
buffer: {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
gap:    {physID: 6, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 300.0}}
"""


def geo_text(g, w, d, ny, r_y, n1, n2, n_gap, n_dn, n_buf, ny_cav,
             bump2, bump_dn, bump_cav, bump_gap, with_gap):
    x_in, x_sp, xr, xpe, x_out, H = (g["x_in"], g["x_split"], g["x_rear"],
                                     g["x_plate_end"], g["x_out"], g["H"])
    xf = xr - w
    L = []
    A = L.append
    A("// case/51 — NASA TN D-8233 平板 + 単独横すきま (平面 2D, node)。gen_mesh.py が生成。")
    A("Geometry.PointNumbers = 0;  lc = 0.01;")
    A(f"x_in = {x_in:.9f}; x_sp = {x_sp:.9f}; xf = {xf:.9f}; xr = {xr:.9f};")
    A(f"xpe = {xpe:.9f}; x_out = {x_out:.9f}; H = {H:.9f}; d = {d:.9f};")
    names = ["x_in", "x_sp", "xf", "xr", "xpe", "x_out"]
    for i, xs in enumerate(names, start=1):          # 1..6 底辺
        A(f"Point({i}) = {{{xs}, 0.0, 0.0, lc}};")
    for i, xs in enumerate(names, start=7):          # 7..12 上辺
        A(f"Point({i}) = {{{xs}, H, 0.0, lc}};")
    for i in range(1, 6):                            # 1..5 底辺
        A(f"Line({i}) = {{{i}, {i+1}}};")
    for i in range(6, 11):                           # 6..10 上辺
        A(f"Line({i}) = {{{i+1}, {i+2}}};")
    for i in range(11, 17):                          # 11..16 縦線
        A(f"Line({i}) = {{{i-10}, {i-4}}};")
    A(f"Transfinite Line {{11, 12, 13, 14, 15, 16}} = {ny} Using Progression {r_y:.8f};")
    A(f"Transfinite Line {{1, 6}} = {n1} Using Progression 1.03;")
    A(f"Transfinite Line {{2, 7}} = {n2} Using Bump {bump2};")
    A(f"Transfinite Line {{3, 8}} = {n_gap}" + (f" Using Bump {bump_gap};" if bump_gap > 0 else ";"))
    A(f"Transfinite Line {{4, 9}} = {n_dn} Using Bump {bump_dn};")
    A(f"Transfinite Line {{5, 10}} = {n_buf} Using Progression 1.03;")
    for k in range(1, 6):
        lb, lr, lt, ll = k, 11 + k, 5 + k, 10 + k
        A(f"Curve Loop({k}) = {{{lb}, {lr}, -{lt}, -{ll}}};")
        A(f"Plane Surface({k}) = {{{k}}};  Transfinite Surface {{{k}}};  Recombine Surface {{{k}}};")
    if with_gap:
        A("Point(13) = {xf, -d, 0.0, lc};")
        A("Point(14) = {xr, -d, 0.0, lc};")
        A("Line(17) = {13, 14};")   # 床
        A("Line(18) = {4, 14};")    # 後壁 (下向き)
        A("Line(19) = {3, 13};")    # 前壁 (下向き)
        A(f"Transfinite Line {{17}} = {n_gap}" + (f" Using Bump {bump_gap};" if bump_gap > 0 else ";"))
        A(f"Transfinite Line {{18, 19}} = {ny_cav} Using Bump {bump_cav};")
        A("Curve Loop(6) = {19, 17, -18, -3};")
        A("Plane Surface(6) = {6};  Transfinite Surface {6};  Recombine Surface(6);")
    A('Physical Curve("inlet",  1) = {11};')
    A('Physical Curve("outlet", 2) = {16};')
    A('Physical Curve("top",    3) = {6, 7, 8, 9, 10};')
    if with_gap:
        A('Physical Curve("plate",  4) = {1, 2, 4};')
        A('Physical Curve("gap",    6) = {17, 18, 19};')
        A('Physical Surface("fluid", 8) = {1, 2, 3, 4, 5, 6};')
    else:
        A('Physical Curve("plate",  4) = {1, 2, 3, 4};')
        A('Physical Surface("fluid", 8) = {1, 2, 3, 4, 5};')
    A('Physical Curve("buffer", 5) = {5};')
    return "\n".join(L) + "\n"


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=["G0", "T2p"], default="T2p",
                    help="G0 = すきま無し (分母診断を同一格子で取る) / T2p = すきま有り")
    ap.add_argument("--y1", type=float, default=30.0, help="壁第一セル [µm]。y1+ = 2376*y1[m] (u_tau 66.31, ro_w 6.61e-4, mu_w 1.85e-5)")
    ap.add_argument("--r-y", type=float, default=1.06)
    ap.add_argument("--nx1", type=int, default=81, help="B1 調整区間")
    ap.add_argument("--nx2", type=int, default=301, help="B2 細分区間 (開口直前)")
    ap.add_argument("--nx-gap", type=int, default=101, help="開口幅方向")
    ap.add_argument("--nx-down", type=int, default=401)
    ap.add_argument("--nx-buf", type=int, default=61)
    ap.add_argument("--ny-cav", type=int, default=301, help="すきま深さ方向")
    ap.add_argument("--bump2", type=float, default=0.02)
    ap.add_argument("--bump-down", type=float, default=0.05)
    ap.add_argument("--bump-cav", type=float, default=0.02)
    ap.add_argument("--bump-gap", type=float, default=0.0)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()

    geom = json.loads((HERE / "geometry.json").read_text(encoding="utf-8"))
    w = geom["gap"]["width"] * 1e-3
    d = geom["gap"]["depth"] * 1e-3
    H = geom["tunnel"]["half_height"] * 1e-3
    dom = geom["domain"]
    g = dict(x_in=dom["x_in"] * 1e-3, x_split=dom["x_split"] * 1e-3,
             x_rear=geom["gap"]["x_rear_wall"] * 1e-3,
             x_plate_end=dom["x_plate_end"] * 1e-3, x_out=dom["x_out"] * 1e-3, H=H)

    y1 = a.y1 * 1e-6
    ny = ny_for(y1, H, a.r_y)
    lo, hi = 1.001, 1.5
    for _ in range(200):                      # ny を固定して y1 を厳密に合わせる
        r = 0.5 * (lo + hi)
        d1 = H * (r - 1.0) / (r ** (ny - 1) - 1.0)
        if d1 > y1:
            lo = r
        else:
            hi = r
    y1_eff = H * (r - 1.0) / (r ** (ny - 1) - 1.0)
    tag = a.tag or f"{a.case.lower()}_y{a.y1:g}um_g{a.nx_gap}"
    MESH.mkdir(exist_ok=True)
    txt = geo_text(g, w, d, ny, r, a.nx1, a.nx2, a.nx_gap, a.nx_down, a.nx_buf, a.ny_cav,
                   a.bump2, a.bump_down, a.bump_cav, a.bump_gap, with_gap=(a.case == "T2p"))
    (MESH / f"{tag}.geo").write_text(txt)
    print(f"[{tag}] W = {w*1e3:.3f} mm, D = {d*1e3:.2f} mm (D/W {d/w:.2f}), "
          f"開口 {a.nx_gap-1} セル ({w/(a.nx_gap-1)*1e6:.1f} µm/セル)")
    print(f"        ny = {ny}, y1 = {y1_eff*1e6:.3f} µm (y1+ ≈ {2376*y1_eff:.3f}), r = {r:.5f}")
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
    print("\n".join(l for l in q.stdout.splitlines()
                    if any(k in l for k in ("aspect", "skew", "VERDICT", "cells:"))))


if __name__ == "__main__":
    main()
