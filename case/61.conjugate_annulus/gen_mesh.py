#!/usr/bin/env python3
r"""V-ax2 同心円環の流体メッシュ (子午面、node・全四角・一様)。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/accepted/boundary-cht-axisymmetric-fem2d.md) §6 V-ax2。

    y=r_b=10 mm ┌──────────────┐ wall_outer (physID 4) = 共役壁 (ints: {conjugate: 1})
                │  静止ガス     │
    y=r_a= 5 mm └──────────────┘ wall_inner (physID 3) = 等温 350 K (非連成)
               x=0            x=L_x=2 mm
               end_x0 (1, slip)  end_xL (2, slip)   軸 (y=0) は含まない

半径方向 `--nr` (既定 32、感度 16) × 軸方向 `--nax` (既定 4) の一様四角。
変換後に**静止・一様 IC** (p0 1013.25 Pa、T 325 K、u=0) を VALUE にパッチする (axcht.patch_static_ic)。

    python3 case/61.conjugate_annulus/gen_mesh.py                 # mesh/annulus_r32_x4.h5
    python3 case/61.conjugate_annulus/gen_mesh.py --nr 16         # mesh/annulus_r16_x4.h5
"""
from __future__ import annotations

import argparse

import axcht

RA, RB, LX = 5.0e-3, 10.0e-3, 2.0e-3

GEO = """// case/61 V-ax2 同心円環 (子午面)。gen_mesh.py が生成。編集しないこと。
Lx = {Lx}; ra = {ra}; rb = {rb};
Point(1) = {{0, ra, 0, 1}};  Point(2) = {{Lx, ra, 0, 1}};
Point(3) = {{Lx, rb, 0, 1}}; Point(4) = {{0, rb, 0, 1}};
Line(1) = {{1, 2}}; Line(2) = {{2, 3}}; Line(3) = {{3, 4}}; Line(4) = {{4, 1}};
Curve Loop(1) = {{1, 2, 3, 4}}; Plane Surface(1) = {{1}};
Transfinite Line{{1, 3}} = {nax} + 1;
Transfinite Line{{2, 4}} = {nr} + 1;
Transfinite Surface{{1}}; Recombine Surface{{1}};
Physical Curve("end_x0",     1) = {{4}};
Physical Curve("end_xL",     2) = {{2}};
Physical Curve("wall_inner", 3) = {{1}};
Physical Curve("wall_outer", 4) = {{3}};
Physical Surface("fluid",    5) = {{1}};
Mesh.MshFileVersion = 4.1;
"""


def conv_cfg():
    # 変換用。mesh ブロックは run の template/solverConfig.yaml と同じにする (nodeWallDirichlet は変換時にも読む)
    return (HERE_TEMPLATE / "solverConfig_dry.yaml").read_text()


def conv_bc():
    return (HERE_TEMPLATE / "bcondConfig_dry.yaml").read_text()


HERE_TEMPLATE = axcht.HERE / "template"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nr", type=int, default=32, help="半径方向のセル数 (登録: 32、感度 16)")
    ap.add_argument("--nax", type=int, default=4, help="軸方向のセル数 (登録: 4)")
    ap.add_argument("--name", default=None, help="既定 annulus_r<nr>_x<nax>")
    a = ap.parse_args()
    name = a.name or f"annulus_r{a.nr}_x{a.nax}"
    mdir = axcht.HERE / "mesh"
    h5 = axcht.gmsh_and_convert(mdir, name, GEO.format(Lx=LX, ra=RA, rb=RB, nax=a.nax, nr=a.nr),
                                conv_cfg(), conv_bc())
    ro = axcht.patch_static_ic(h5)
    print(f"[gen_mesh] {h5}: {a.nr} (r) x {a.nax} (x) 一様四角、dr {(RB-RA)/a.nr*1e3:.5f} mm, "
          f"dx {LX/a.nax*1e3:.4f} mm、静止 IC p0 {axcht.P0} Pa / T {axcht.T_IC} K / ro {ro:.6e}")
    q = axcht.mesh_quality(h5)
    print(q)
    (mdir / f"{name}.quality.txt").write_text(q)


if __name__ == "__main__":
    main()
