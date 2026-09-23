#!/usr/bin/env python3
r"""case/48 の冷却平板の固体を **帯メッシュ**にする (CHT V2 の case/48 段)。

plan boundary-conjugate-heat-transfer §6 V2。`local1d` は面内伝導を落とした極限なので、
**SU2 の 2D 固体ゾーンと同じ方程式**で比べるには面内伝導つきの固体が要る。ソルバ内 `shell2d` は
不採用にした (§4.4c) ので、**壁を厚さ $t$ で押し出した帯メッシュを `fem2d` に食わせる**。

固体 (case/48 の `conjugate:` と同じ):

    t = 1.0e-3 m,  k_s = 0.217 W/mK   (t/k_s = 4.608e-3 m2K/W),  背面 300 K

**背面の扱い**: `fem2d` は Robin 辺しか持てない (Dirichlet は未実装) ので、
$h_{\rm back}$ を大きく取って近似する。背面抵抗 $1/h$ が固体抵抗 $t/k_s$ に対して十分小さければ
Dirichlet と区別がつかない。既定 $h$=1e8 で $1/h/(t/k_s)$ = 2.2e-6。
**$h$ を 2 水準振って答えが動かないことを確かめること** (`--h-back` を 1e6 と 1e8 で比較)。

界面 (y=0) の節点は**流体の壁節点と 1 対 1**にする (壁ダンプから x を読む)。

使い方:
  python3 case/48.flat_plate_cooled_m4/gen_solid_strip.py \
      --wall <run>/res_wall_4_30000.h5 --nl 8 --out case/48.flat_plate_cooled_m4/mesh/solid_strip_nl8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True, help="流体の壁ダンプ (界面 x を読む)")
    ap.add_argument("--t", type=float, default=1.0e-3, help="固体の厚さ [m]")
    ap.add_argument("--k-solid", type=float, default=0.217, help="[W/mK]")
    ap.add_argument("--Tb", type=float, default=300.0, help="背面温度 [K]")
    ap.add_argument("--h-back", type=float, default=1.0e8,
                    help="背面 Robin の h [W/m2K]。Dirichlet の近似 (大きいほど厳密)")
    ap.add_argument("--nl", type=int, default=8, help="厚さ方向の層数")
    ap.add_argument("--out", required=True, help="出力の基底名 (.npz と .json を書く)")
    a = ap.parse_args()

    with h5py.File(a.wall, "r") as f:
        c = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
    x = np.sort(np.unique(np.round(c[:, 0], 12)))
    if len(x) != len(c):
        raise SystemExit(f"壁節点 {len(c)} に対し x の相異なる値が {len(x)} — 平面 2D でない")
    nx, nl = len(x), a.nl

    ys = np.linspace(-a.t, 0.0, nl + 1)
    nodes = np.array([[xx, yy] for yy in ys for xx in x], float)
    nid = lambda i, j: j * nx + i                      # noqa: E731

    tris = []
    for j in range(nl):
        for i in range(nx - 1):
            tris.append([nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)])
            tris.append([nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)])
    tris = np.array(tris, int)

    outer = np.array([[nid(i, nl), nid(i + 1, nl)] for i in range(nx - 1)], int)   # y=0 界面
    back = np.array([[nid(i, 0), nid(i + 1, 0)] for i in range(nx - 1)], int)      # y=-t 背面

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out.with_suffix(".npz"), nodes=nodes, tris=tris, outer_edges=outer, hole1=back)
    spec = {"_note": f"case/48 の固体帯 (厚さ {a.t} m, {nl} 層)。背面は h={a.h_back:g} の Robin で "
                     f"Dirichlet {a.Tb} K を近似 (1/h / (t/k_s) = {(1.0/a.h_back)/(a.t/a.k_solid):.3g})",
            "mesh_npz": str(out.with_suffix(".npz")),
            "k_solid": a.k_solid,
            "holes": [{"h": a.h_back, "T_c": a.Tb}],
            "T_init": 500.0}
    out.with_suffix(".json").write_text(json.dumps(spec, ensure_ascii=False, indent=2))

    print(f"[gen_solid_strip] {out.with_suffix('.npz')}")
    print(f"  節点 {len(nodes)} / 三角形 {len(tris)} / 界面辺 {len(outer)} / 背面辺 {len(back)}")
    print(f"  厚さ {a.t} m を {nl} 層,  x は壁節点と 1 対 1 ({nx} 点)")
    print(f"  t/k_s = {a.t/a.k_solid:.6e},  1/h = {1.0/a.h_back:.3e} "
          f"({(1.0/a.h_back)/(a.t/a.k_solid):.3g} 倍)")


if __name__ == "__main__":
    main()
