#!/usr/bin/env python3
r"""case/63 Graetz の固体殻 (R ≤ r ≤ R + t、加熱区間 0 ≤ x ≤ L_heat だけ)。

plan [`boundary-cht-axisymmetric-graetz.md`](../../plans/accepted/boundary-cht-axisymmetric-graetz.md) §4.4。
case/61 `gen_solid.py` と同じ帯トポロジ (行 j=0 が外面 r = R+t、j = n_s が界面 r = R、界面の並び i は x 昇順、
対角は左上→右下)。界面節点は**乾式 1 step の壁ダンプ** `res_wall_heat_4_1.h5` の座標そのもの (流体の壁節点と 1 対 1)。

- $k_s$ = 100 W/mK、厚さ t = 0.5 mm、半径方向 8 層一様 (登録。格子感度でも固定。感度 `--ns 16` を N_r=32 で 1 本)、外面 Robin $h$ = 1e8・$T_c = T_{in} + \Delta T$
- 固体の x 両端は断熱 (未登録 = 自然境界)
- $T_c$ は固体 h5 に入るので ΔT ごとに作る: mesh/solid_<乾式 run 名>_s<層数>_dT<ΔT>.h5

    python3 case/63.graetz_cht/gen_solid.py --wall case/63.graetz_cht/run_0001_dry_r16/res_wall_heat_4_1.h5 --dT 10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import graetz_common as gc
from graetz_common import axcht

T_S, K_S, H_ROBIN, NS_REG = 0.5e-3, 100.0, 1.0e8, 8


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True, help="乾式 1 step の wall_heat 壁ダンプ")
    ap.add_argument("--dT", type=float, required=True, help="T_c − T_in [K]")
    ap.add_argument("--ns", type=int, default=NS_REG, help="半径方向の層数 (登録 8、感度 16)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    NS = a.ns
    c = axcht.read_wall_coords(Path(a.wall))
    if not np.allclose(c[:, 1], gc.R, rtol=0, atol=1e-9) or np.abs(c[:, 2]).max() > 0:
        raise SystemExit("壁ダンプが r = R 一定・z = 0 の直線でない (wall_heat を渡しているか?)")
    if abs(c[:, 0].min()) > 1e-9 or abs(c[:, 0].max() - gc.L_HEAT) > 1e-8:
        raise SystemExit("壁ダンプの x 範囲が加熱区間 [0, L_heat] でない")
    order = np.argsort(c[:, 0])
    iface = c[order, :2]
    if len(np.unique(iface[:, 0])) != len(iface):
        raise SystemExit("壁節点の x が重複している")
    nax = len(iface) - 1
    rc = float(iface[0, 1]) + T_S                       # 界面の r (float32 由来) + t
    nodes, tris, outer, back = axcht.solid_strip(iface, normal_axis=1, back_coord=rc, ns=NS)
    ys = nodes.reshape(NS + 1, -1, 2)[:, 0, 1]
    if not np.all(np.diff(ys) < 0):
        raise SystemExit("行が y 降順になっていない")
    if axcht.diag_direction(nodes, len(iface), NS) != "左上→右下":
        raise SystemExit("対角が登録 (左上→右下) と違う")
    Tc = gc.T_IN + a.dT
    tag = Path(a.wall).resolve().parent.name
    out = Path(a.out) if a.out else gc.HERE / "mesh" / f"solid_{tag}_s{NS}_dT{a.dT:g}"
    meta = {"case": "case/63 Graetz の固体殻", "wall_dump": str(Path(a.wall)), "wall_dump_sha256": axcht.sha256(a.wall),
            "n_iface": int(len(iface)), "n_layers": NS, "n_axial": int(nax),
            "radial": {"from": float(iface[0, 1]), "to": rc, "layers": NS, "spacing": "uniform"},
            "axial": {"nodes_x": iface[:, 0].tolist(), "note": "流体の壁節点と 1 対 1 (乾式の壁ダンプ)"},
            "iface_to_wall_dump_index": order.tolist(),
            "k_solid": K_S, "robin": {"edge": "r = R + t", "h": H_ROBIN, "T_c": Tc},
            "adiabatic": "x = 0 と x = L_heat (未登録 = 自然境界)"}
    note = (f"case/63 Graetz 固体殻 r {iface[0,1]:.9g} .. {rc:.9g} m、{NS} 層 × 軸 {nax}。k_s {K_S}、外面 Robin h={H_ROBIN:g} / "
            f"T_c={Tc} K。x 両端は断熱")
    h5, meta, log = axcht.write_solid(out, nodes, tris, outer, back, K_S, H_ROBIN, Tc, note, meta)
    print(log)
    print(f"[gen_solid] {h5}: {len(nodes)} 節点 / {len(tris)} 三角形 / 界面 {len(iface)} 点 / T_c {Tc} K")


if __name__ == "__main__":
    main()
