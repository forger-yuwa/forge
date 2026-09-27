#!/usr/bin/env python3
r"""V-ax2 同心円環の固体殻 (r_b ≤ r ≤ r_c) を作る。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/accepted/boundary-cht-axisymmetric-fem2d.md) §6 V-ax2
「固体格子の登録」: `gen_solid_strip.py` の規則 — **$y$ を降順に並べ** (行 j=0 が外面 $r_c$、j=n_s が界面 $r_b$)、
**各セルで同じ対角 (左上→右下)**。界面節点は**乾式 1 step の壁ダンプ** (`res_wall_outer_4_1.h5`) から取り、
流体の壁節点と 1 対 1 (角節点の physID 所属を推測しない。case/58 と同じ)。界面の並び i は $x$ 昇順。

- 外面 $y=r_c$ = Robin ($h$=1e8、$T_c$=300 K) = `hole1`、界面 $y=r_b$ = `outer_edges`、$x$ 両端 = 未登録 = 断熱
- 半径方向 `--ns` 層 (一様、既定 16)。軸方向は界面節点で決まる (= 流体の軸方向分割)
- $k_s$ = 0.027 W/mK、固体の初期温度はソルバが共役壁の `Ts` (325 K) から作る (json の T_init は外部ループ用)

出力 (mesh/solid_s<ns>_x<nax>.*): `.npz` (帯)、`.json` (物性)、`.h5` (ソルバ入力、RCM)、
`.grid.json` (座標順・接続・分割数・界面 ↔ 壁ダンプ節点の対応・ハッシュ)。

    python3 case/61.conjugate_annulus/gen_solid.py --wall case/61.conjugate_annulus/run_0001_dry_r32/res_wall_outer_4_1.h5 --ns 16
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import axcht

RB, RC, KS, H, TC = 10.0e-3, 20.0e-3, 0.027, 1.0e8, 300.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True, help="乾式 1 step の wall_outer 壁ダンプ")
    ap.add_argument("--ns", type=int, default=16, help="半径方向の層数 (登録: 既定 16、感度 8/16/32)")
    ap.add_argument("--out", default=None, help="出力の基底名 (既定 mesh/solid_s<ns>_x<nax>)")
    a = ap.parse_args()

    c = axcht.read_wall_coords(Path(a.wall))
    if not np.allclose(c[:, 1], RB, rtol=0, atol=1e-8) or np.abs(c[:, 2]).max() > 0:
        raise SystemExit("壁ダンプが y = r_b 一定・z = 0 の直線でない (wall_outer を渡しているか?)")
    order = np.argsort(c[:, 0])                     # x 昇順
    iface = c[order, :2]
    if len(np.unique(iface[:, 0])) != len(iface):
        raise SystemExit("壁節点の x が重複している")
    nax = len(iface) - 1
    nodes, tris, outer, back = axcht.solid_strip(iface, normal_axis=1, back_coord=RC, ns=a.ns)
    # 登録: 行は y 降順 (j=0 が r_c)、対角は左上→右下
    ys = nodes.reshape(a.ns + 1, -1, 2)[:, 0, 1]
    if not np.all(np.diff(ys) < 0):
        raise SystemExit("行が y 降順になっていない")
    if axcht.diag_direction(nodes, len(iface), a.ns) != "左上→右下":
        raise SystemExit("対角が登録 (左上→右下) と違う")
    if nodes[:, 1].min() <= 0.0:
        raise SystemExit("固体に r <= 0 の節点がある (軸対称は r = y > 0)")

    out = Path(a.out) if a.out else axcht.HERE / "mesh" / f"solid_s{a.ns}_x{nax}"
    meta = {"case": "case/61 V-ax2 同心円環の固体殻",
            "wall_dump": str(Path(a.wall)), "wall_dump_sha256": axcht.sha256(a.wall),
            "n_iface": int(len(iface)), "n_layers": int(a.ns), "n_axial": int(nax),
            "radial": {"from": RB, "to": RC, "layers": a.ns, "spacing": "uniform"},
            "axial": {"nodes_x": iface[:, 0].tolist(), "note": "流体の壁節点と 1 対 1 (乾式の壁ダンプ)"},
            "iface_to_wall_dump_index": order.tolist(),
            "iface_npz_index": [int((a.ns) * len(iface) + i) for i in range(len(iface))],
            "k_solid": KS, "robin": {"edge": "y = r_c", "h": H, "T_c": TC},
            "adiabatic": "x = 0 と x = L_x (未登録 = 自然境界)"}
    note = (f"case/61 V-ax2 固体殻 r {RB} .. {RC} m、{a.ns} 層 × 軸 {nax}。外面 Robin h={H:g} / T_c={TC} K。"
            "x 両端は断熱。対角は左上→右下 (gen_solid_strip.py の規則)")
    h5, meta, log = axcht.write_solid(out, nodes, tris, outer, back, KS, H, TC, note, meta)
    print(log)
    print(f"[gen_solid] {h5}: {len(nodes)} 節点 / {len(tris)} 三角形 / 界面 {len(iface)} 点 / 対角 {meta['diagonal']}")
    print(f"  半径 dr {(RC-RB)/a.ns*1e3:.4f} mm × 軸 dx {np.diff(iface[:,0]).max()*1e3:.4f} mm"
          f" (縦横比 {np.diff(iface[:,0]).max()/((RC-RB)/a.ns):.3f})")


if __name__ == "__main__":
    main()
