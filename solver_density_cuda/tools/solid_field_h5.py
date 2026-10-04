#!/usr/bin/env python3
r"""CHT の**固体内部の温度場**を HDF5 + XDMF に出す (ParaView 用)。

これまで連成が残していたのは界面の壁温 (`Tw_final.csv`) と履歴 (`cht_history.csv`) だけで、
**固体内部は `cht_loop.py` の中で消えていた**。移植後のソルバ内 `fem2d` は同じ内容を
`res_solid_<step>.h5` として自前で書く (plan boundary-conjugate-heat-transfer §4.6a) が、
本ツールは**既に収束している外部ループの run** からその場を再構成する。
移植の同値試験を目視で確かめる基準にもなる。

出力 (forge の `res_*.h5` と同じ XDMF 規約):

    MESH/COORD    節点座標 (x, y, 0)
    MESH/CONNE    三角形 (XDMF Mixed の要素コード 4)
    VALUE/T       節点温度 [K]
    VALUE/k_s     その温度での熱伝導率 [W/mK]
    VALUE/q_iface 界面節点がガス側から受け取った熱 [W/m] (非界面節点は 0)
    VALUE/q_hole  冷却孔 (Robin) が持ち去った熱 [W/m] (孔以外は 0)

`q_iface` は解から作った**反力** $(K_su-b_s)_{\rm iface}$ で、荷重を読み直したものではない。
したがって「固体が実際に受け取った熱」であり、G-cons の流体側 $\sum Q_f$ と突き合わせられる。

使い方:
  python3 solver_density_cuda/tools/solid_field_h5.py \
      --solid case/53.c3x_vane_cht/solid_published.json \
      --tw    case/53.c3x_vane_cht/run_0144_cht_published/Tw_final.csv \
      --out   case/53.c3x_vane_cht/run_0144_cht_published/res_solid.h5
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cht_loop import build_fem2d  # noqa: E402


def write_xmf(xmf: Path, h5name: str, n_nodes: int, conne_dim: int, n_elem: int,
              names: list[str], time_value: float = 0.0):
    with open(xmf, "w") as f:
        f.write("<?xml version='1.0' ?>\n<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n<Xdmf>\n  <Domain>\n")
        f.write("    <Grid GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n")
        f.write(f"    <Time TimeType='Single' Value='{time_value}'/> \n")
        f.write("      <Grid Name='solid'>\n")
        f.write(f"        <Topology Type='Mixed' NumberOfElements='{n_elem}'>\n")
        f.write(f"          <DataItem Format='HDF' DataType='Int' Dimensions='{conne_dim}'>\n")
        f.write(f"            {h5name}:MESH/CONNE\n          </DataItem>\n        </Topology>\n")
        f.write("        <Geometry Type='XYZ'>\n")
        f.write(f"          <DataItem Format='HDF' DataType='Float' Dimensions='{n_nodes*3}'>\n")
        f.write(f"            {h5name}:MESH/COORD\n          </DataItem>\n        </Geometry>\n")
        for nm in names:
            f.write(f"        <Attribute Name='{nm}' Center='Node' >\n")
            f.write(f"          <DataItem Format='HDF' DataType='Float' Dimensions='{n_nodes}'>\n")
            f.write(f"            {h5name}:VALUE/{nm}\n          </DataItem>\n        </Attribute>\n")
        f.write("      </Grid>\n    </Grid>\n  </Domain>\n</Xdmf>\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solid", required=True, help="固体 JSON (fem2d: mesh_npz / holes / k_table)")
    ap.add_argument("--tw", required=True, help="界面温度 CSV (x y z Ts。壁ダンプ順)")
    ap.add_argument("--out", required=True, help="出力 h5 (同名の .xmf も書く)")
    ap.add_argument("--time", type=float, default=0.0, help="XDMF に書く時刻")
    a = ap.parse_args()

    spec = json.loads(Path(a.solid).read_text())
    if "mesh_npz" not in spec:
        sys.exit("[solid_field_h5] この JSON は fem2d 用ではない (mesh_npz が無い)")

    # 区切りは `cht_loop` が CSV、ソルバ側 (`conjugate_Tw_*.csv`) が空白。両方受ける。
    _head = Path(a.tw).read_text().splitlines()[0]
    tw = np.loadtxt(a.tw, skiprows=1, delimiter="," if "," in _head else None)
    coords, Tw_wall = tw[:, :3], tw[:, 3]
    op, perm = build_fem2d(spec, coords)
    T_iface = Tw_wall[perm]                       # 壁ダンプ順 -> 界面節点順

    # $k_s(T)$ を内部温度と自己整合させる (plan §5.1 #33 と同型: 1 回組んで復元しただけでは
    # 組立てと復元で評価温度がずれ、内部残差が 0 にならない = 収支が 1 % 級で合わない)。
    u = None
    for it in range(20):
        op.assemble(T_iface)                      # 現在の op.u の k_s(T) で組む (内部消去の LU)
        u_new = op.recover_interior(T_iface)      # 全節点温度
        if u is not None and np.max(np.abs(u_new - u)) < 1e-10:
            u = u_new
            break
        u = u_new
    K, b = op.assemble_full(u)
    resid = np.asarray(K.dot(u)).ravel() - np.asarray(b).ravel()
    react = resid                                  # 反力 [W/m] (界面以外では内部残差 = 0 のはず)
    r_int = float(np.max(np.abs(resid[op._other]))) if len(op._other) else 0.0

    q_iface = np.zeros(op.N)
    q_iface[op.iface] = react[op.iface]

    # 孔 Robin が持ち去った熱を節点へ集中させる (辺の中点則)。
    q_hole = np.zeros(op.N)
    for (n0, n1, h, Tc) in op.robin_edges:
        L = float(np.linalg.norm(op.xy[n1] - op.xy[n0]))
        for g in (n0, n1):
            q_hole[g] += 0.5 * L * h * (u[g] - Tc)

    xy3 = np.zeros((op.N, 3))
    xy3[:, :2] = op.xy
    conne = np.concatenate([np.column_stack([np.full(len(op.tris), 4), op.tris]).ravel()])

    out = Path(a.out)
    with h5py.File(out, "w") as f:
        f.create_dataset("MESH/COORD", data=xy3.ravel().astype(np.float32))
        f.create_dataset("MESH/CONNE", data=conne.astype(np.int32))
        f.create_dataset("VALUE/T", data=u.astype(np.float32))
        f.create_dataset("VALUE/k_s", data=np.asarray(op.k_of(u), float).astype(np.float32))
        f.create_dataset("VALUE/q_iface", data=q_iface.astype(np.float32))
        f.create_dataset("VALUE/q_hole", data=q_hole.astype(np.float32))
    write_xmf(out.with_suffix(".xmf"), out.name, op.N, len(conne), len(op.tris),
              ["T", "k_s", "q_iface", "q_hole"], a.time)

    print(f"[solid_field_h5] {out}  nodes={op.N} tris={len(op.tris)} iface={len(op.iface)}")
    print(f"  T        {u.min():.2f} .. {u.max():.2f} K   (界面 {T_iface.min():.2f} .. {T_iface.max():.2f})")
    print(f"  k_s      {op.k_of(u).min():.3f} .. {op.k_of(u).max():.3f} W/mK")
    print(f"  界面から入った熱 sum(q_iface) = {q_iface.sum():.3f} W/m")
    print(f"  孔が持ち去った熱 sum(q_hole)  = {q_hole.sum():.3f} W/m")
    print(f"  収支の残り                    = {q_iface.sum() - q_hole.sum():.3e} W/m"
          f"  ({abs(q_iface.sum()-q_hole.sum())/max(1e-30,abs(q_iface.sum()))*100:.4f} %)")
    print(f"  内部残差 max|Ku-b|_interior   = {r_int:.3e} W/m  (Picard {it+1} 回)")


if __name__ == "__main__":
    main()
