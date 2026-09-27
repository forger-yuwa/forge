#!/usr/bin/env python3
r"""V-ax2b 同軸円板の流体メッシュ (子午面、node・全四角)。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md) §6 V-ax2b。

    r=r_2=20 mm ┌──────────┐ side_r2 (physID 2, slip)
                │ 静止ガス  │
    wall_hot    │          │ wall_cj (physID 4) = 共役壁 x=H (固体は x∈[H, H+t])
    (physID 3)  │          │
    x=0, 350 K  │          │
    r=r_1= 5 mm └──────────┘ side_r1 (physID 1, slip)     軸は含まない
               x=0        x=H=5 mm

格子 (登録、plan レビュー M5): $x$ 方向 16 一様 × $r$ 方向 `--nr` (8 / 16 / 32 一様、または 32 で `--ratio 1.1`
= 等比 1.1、**間隔は $r$ が大きいほど広がる**)。

**軸に触れる負例メッシュ** (`--axis-probe`): $r\in[0, r_2]$、$r=0$ は `axis`。共役はせず `interfaceDiag` の恒等式
$q_{\rm eff}A^r_{\rm fluid}=Q_{f,\rm eff}$ が `axisRFloor` の床で成り立つかを 1 step で見る専用 (template/*_axisprobe.yaml)。

    python3 case/62.conjugate_disk/gen_mesh.py --nr 32                 # mesh/disk_r32_u.h5
    python3 case/62.conjugate_disk/gen_mesh.py --nr 32 --ratio 1.1     # mesh/disk_r32_g1p1.h5
    python3 case/62.conjugate_disk/gen_mesh.py --nr 16 --axis-probe    # mesh/disk_axisprobe_r16.h5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "61.conjugate_annulus"))
import axcht  # noqa: E402  (case/61 の共通部品)

HERE = Path(__file__).resolve().parent
H, R1, R2, NX = 5.0e-3, 5.0e-3, 20.0e-3, 16

GEO = """// case/62 V-ax2b 同軸円板 (子午面)。gen_mesh.py が生成。編集しないこと。
H = {H}; r1 = {r1}; r2 = {r2};
Point(1) = {{0, r1, 0, 1}};  Point(2) = {{H, r1, 0, 1}};
Point(3) = {{H, r2, 0, 1}};  Point(4) = {{0, r2, 0, 1}};
Line(1) = {{1, 2}};   // r = r1
Line(2) = {{2, 3}};   // x = H  (r が増える向き)
Line(3) = {{4, 3}};   // r = r2
Line(4) = {{1, 4}};   // x = 0  (r が増える向き)
Curve Loop(1) = {{1, 2, -3, -4}}; Plane Surface(1) = {{1}};
Transfinite Line{{1, 3}} = {nx} + 1;
Transfinite Line{{2, 4}} = {nr} + 1{prog};
Transfinite Surface{{1}}; Recombine Surface{{1}};
Physical Curve("{n1}",    1) = {{1}};
Physical Curve("side_r2", 2) = {{3}};
Physical Curve("wall_hot", 3) = {{4}};
Physical Curve("wall_cj",  4) = {{2}};
Physical Surface("fluid",  5) = {{1}};
Mesh.MshFileVersion = 4.1;
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nr", type=int, required=True, help="r 方向のセル数 (登録: 8 / 16 / 32)")
    ap.add_argument("--ratio", type=float, default=1.0, help="r 方向の等比 (登録: 1.0 一様 / 1.1 を N_r=32 で)")
    ap.add_argument("--axis-probe", action="store_true", help="r ∈ [0, r_2] の負例メッシュ (共役なし)")
    a = ap.parse_args()
    if a.axis_probe and a.ratio != 1.0:
        raise SystemExit("--axis-probe は一様格子だけ")
    r1 = 0.0 if a.axis_probe else R1
    prog = "" if a.ratio == 1.0 else f" Using Progression {a.ratio}"
    name = (f"disk_axisprobe_r{a.nr}" if a.axis_probe else
            f"disk_r{a.nr}_" + ("u" if a.ratio == 1.0 else f"g{a.ratio:g}".replace(".", "p")))
    tdir = HERE / "template"
    sfx = "_axisprobe" if a.axis_probe else "_dry"
    h5 = axcht.gmsh_and_convert(HERE / "mesh", name,
                                GEO.format(H=H, r1=r1, r2=R2, nx=NX, nr=a.nr, prog=prog,
                                           n1="axis" if a.axis_probe else "side_r1"),
                                (tdir / f"solverConfig{sfx}.yaml").read_text(),
                                (tdir / f"bcondConfig{sfx}.yaml").read_text())
    ro = axcht.patch_static_ic(h5)
    import h5py
    import numpy as np
    with h5py.File(h5, "r") as f:
        c = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
    rr = np.unique(np.round(c[:, 1], 9))
    dr = np.diff(rr)
    print(f"[gen_mesh] {h5}: x {NX} 一様 × r {a.nr} (比 {a.ratio:g})、r {rr[0]*1e3:.3f} .. {rr[-1]*1e3:.3f} mm、"
          f"dr {dr[0]*1e3:.5f} → {dr[-1]*1e3:.5f} mm (隣接比 {dr[1]/dr[0]:.4f})、静止 IC ro {ro:.6e}")
    if a.ratio != 1.0 and not (dr[-1] > dr[0]):
        raise SystemExit("非一様の向きが登録 (r が大きいほど広い) と逆")
    q = axcht.mesh_quality(h5)
    print(q)
    (HERE / "mesh" / f"{name}.quality.txt").write_text(q)


if __name__ == "__main__":
    main()
