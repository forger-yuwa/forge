#!/usr/bin/env python3
r"""V-ax2b 同軸円板の固体 ($x\in[H, H+t]$、$r\in[r_1, r_2]$) を作る。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md) §6 V-ax2b の格子
(事前固定、plan レビュー M5): 固体は **$x$ 方向 16 層 × $r$ 方向 $N_r$** (界面節点は流体の壁節点と一致)、
**四角形は左下→右上の対角**で 2 三角形 (`gen_solid_strip.py` の規則)。

`gen_solid_strip.py` (case/58) と同じ帯トポロジ — 行 j = 0 (背面 $x=H+t$) … 16 (界面 $x=H$)、
界面の並び i は **$y (=r)$ 降順** — で組むと、各セルの対角 (j,i)–(j+1,i+1) は $x$ も $r$ も減る向き =
**左下↔右上**になる (case/58 は固体が界面の −x 側、ここは +x 側なので見かけの向きが鏡映される)。

- 背面 $x=H+t$ = Robin ($h$=1e8、$T_c$=300 K) = `hole1`、界面 $x=H$ = `outer_edges`、$r$ 両端 = 未登録 = 断熱
- $k_s$ = 0.027 W/mK

    python3 case/62.conjugate_disk/gen_solid.py --wall case/62.conjugate_disk/run_0001_dry_r8u/res_wall_cj_4_1.h5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "61.conjugate_annulus"))
import axcht  # noqa: E402

HERE = Path(__file__).resolve().parent
H, T, KS, HR, TC = 5.0e-3, 5.0e-3, 0.027, 1.0e8, 300.0
NS = 16


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True, help="乾式 1 step の wall_cj 壁ダンプ")
    ap.add_argument("--ns", type=int, default=NS, help="x 方向の層数 (登録: 16。変えない)")
    ap.add_argument("--out", default=None, help="出力の基底名 (既定 mesh/solid_<流体メッシュ名>)")
    a = ap.parse_args()

    c = axcht.read_wall_coords(Path(a.wall))
    if not np.allclose(c[:, 0], H, rtol=0, atol=1e-8) or np.abs(c[:, 2]).max() > 0:
        raise SystemExit("壁ダンプが x = H 一定・z = 0 の直線でない (wall_cj を渡しているか?)")
    order = np.argsort(-c[:, 1])                    # r (= y) 降順
    iface = c[order, :2]
    if len(np.unique(iface[:, 1])) != len(iface):
        raise SystemExit("壁節点の r が重複している")
    nodes, tris, outer, back = axcht.solid_strip(iface, normal_axis=0, back_coord=H + T, ns=a.ns)
    if axcht.diag_direction(nodes, len(iface), a.ns) != "左下→右上":
        raise SystemExit("対角が登録 (左下→右上) と違う")
    if nodes[:, 1].min() <= 0.0:
        raise SystemExit("固体に r <= 0 の節点がある")
    nr = len(iface) - 1
    dr = -np.diff(iface[:, 1])
    run_mesh = Path(a.wall).resolve().parent / "RUN_INPUTS.txt"
    tag = f"r{nr}_" + ("u" if np.allclose(dr, dr[0], rtol=1e-5) else f"g{dr[0]/dr[1]:.3g}".replace(".", "p"))
    out = Path(a.out) if a.out else HERE / "mesh" / f"solid_disk_{tag}"
    meta = {"case": "case/62 V-ax2b 同軸円板の固体",
            "wall_dump": str(Path(a.wall)), "wall_dump_sha256": axcht.sha256(a.wall),
            "fluid_run_inputs": run_mesh.read_text() if run_mesh.exists() else None,
            "n_iface": int(len(iface)), "n_layers": int(a.ns), "n_radial": int(nr),
            "axial_x": {"from": H, "to": H + T, "layers": a.ns, "spacing": "uniform"},
            "radial": {"nodes_r_desc": iface[:, 1].tolist(), "dr_first_last": [float(dr[-1]), float(dr[0])],
                       "note": "流体の壁節点と 1 対 1 (乾式の壁ダンプ)。並びは r 降順"},
            "iface_to_wall_dump_index": order.tolist(),
            "iface_npz_index": [int(a.ns * len(iface) + i) for i in range(len(iface))],
            "k_solid": KS, "robin": {"edge": "x = H + t", "h": HR, "T_c": TC},
            "adiabatic": "r = r_1 と r = r_2 (未登録 = 自然境界)"}
    note = (f"case/62 V-ax2b 固体 x {H} .. {H+T} m ({a.ns} 層) × r {nr}。背面 Robin h={HR:g} / T_c={TC} K。"
            "r 両端は断熱。対角は左下→右上 (gen_solid_strip.py の規則、界面の並びは r 降順)")
    h5, meta, log = axcht.write_solid(out, nodes, tris, outer, back, KS, HR, TC, note, meta)
    print(log)
    print(f"[gen_solid] {h5}: {len(nodes)} 節点 / {len(tris)} 三角形 / 界面 {len(iface)} 点 / 対角 {meta['diagonal']}"
          f" / dr {dr.min()*1e3:.5f} .. {dr.max()*1e3:.5f} mm")


if __name__ == "__main__":
    main()
