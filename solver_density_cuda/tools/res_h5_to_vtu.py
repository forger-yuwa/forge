#!/usr/bin/env python3
"""forge の res_*.h5 / mesh.h5 を VTU (VTK UnstructuredGrid) に変換する。

forge の出力する `res_*.xmf` は XDMF "Mixed" topology を使うため、3D の
tet/prism 混在メッシュや prism(Wedge) を ParaView の XDMF リーダが開けない
ことがある。本ツールは `/MESH/CONNE` (XDMF コード + 節点) と `/MESH/COORD`、
`/VALUE/*` を読み、ParaView が確実に開ける VTU (XML, ascii) を書き出す。

注意:
- 多面体セル (forge の ieleType=0) を含むメッシュ (Fluent poly-hexcore 等) では
  forge 出力の CONNE 自体が壊れているため変換できない (tet/hex/prism/pyramid/
  tri/quad のみ対応)。
- 値が NaN/Inf (発散) の場合もそのまま書く。ParaView では色付けで NaN になるが
  形状は表示できる。

使い方:
    python3 res_h5_to_vtu.py res_0.h5            # -> res_0.vtu
    python3 res_h5_to_vtu.py res_0.h5 out.vtu
    python3 res_h5_to_vtu.py res_0.h5 --fields ro Ux Uy Uz P T
"""
import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("h5py が必要です: pip install h5py")

# XDMF topology code -> (節点数, VTK cell type)
# VTK: TRIANGLE=5 QUAD=9 TETRA=10 HEXAHEDRON=12 WEDGE=13 PYRAMID=14
XDMF2VTK = {
    4: (3, 5),    # triangle
    5: (4, 9),    # quad
    6: (4, 10),   # tetra
    7: (5, 14),   # pyramid
    8: (6, 13),   # wedge (prism) -- 節点順は gmsh/XDMF と VTK で一致
    9: (8, 12),   # hexahedron
}


def _count_cells(conne):
    """CONNE を走査してセル数を数える (nCells 不明時)。"""
    i = 0
    n = 0
    while i < len(conne):
        code = int(conne[i])
        if code not in XDMF2VTK:
            raise ValueError(f"未対応の XDMF コード {code} (多面体等)。VTU 化できない。")
        i += 1 + XDMF2VTK[code][0]
        n += 1
    return n


def parse_conne(conne, ncells):
    """XDMF Mixed CONNE を (connectivity, offsets, vtk_types) に変換。"""
    conn = []
    offsets = []
    vtk_types = []
    i = 0
    for _ in range(ncells):
        if i >= len(conne):
            raise ValueError("CONNE が NumberOfElements より短い (壊れている可能性)")
        code = int(conne[i]); i += 1
        if code not in XDMF2VTK:
            raise ValueError(f"未対応の XDMF コード {code} (多面体等)。"
                             f" このメッシュは VTU 化できない。")
        nn, vt = XDMF2VTK[code]
        conn.extend(int(x) for x in conne[i:i + nn]); i += nn
        offsets.append(len(conn))
        vtk_types.append(vt)
    return np.array(conn, dtype=np.int64), np.array(offsets, dtype=np.int64), \
        np.array(vtk_types, dtype=np.uint8)


def write_vtu(path, points, conn, offsets, types, cell_data, point_data=None):
    n_pts = len(points)
    n_cells = len(offsets)
    out = []
    out.append('<?xml version="1.0"?>')
    out.append('<VTKFile type="UnstructuredGrid" version="1.0" '
               'byte_order="LittleEndian">')
    out.append("  <UnstructuredGrid>")
    out.append(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">')

    # points
    out.append("      <Points>")
    out.append('        <DataArray type="Float64" NumberOfComponents="3" '
               'format="ascii">')
    out.append(_fmt(points.reshape(-1), 6))
    out.append("        </DataArray>")
    out.append("      </Points>")

    # cells
    out.append("      <Cells>")
    out.append('        <DataArray type="Int64" Name="connectivity" format="ascii">')
    out.append(_fmt(conn, 0))
    out.append("        </DataArray>")
    out.append('        <DataArray type="Int64" Name="offsets" format="ascii">')
    out.append(_fmt(offsets, 0))
    out.append("        </DataArray>")
    out.append('        <DataArray type="UInt8" Name="types" format="ascii">')
    out.append(_fmt(types, 0))
    out.append("        </DataArray>")
    out.append("      </Cells>")

    # cell data
    if cell_data:
        names = ",".join(cell_data.keys())
        out.append(f'      <CellData Scalars="{names}">')
        for name, arr in cell_data.items():
            out.append(f'        <DataArray type="Float32" Name="{name}" '
                       f'format="ascii">')
            out.append(_fmt(arr, 6))
            out.append("        </DataArray>")
        out.append("      </CellData>")

    # point data (node 離散化の節点値)
    if point_data:
        names = ",".join(point_data.keys())
        out.append(f'      <PointData Scalars="{names}">')
        for name, arr in point_data.items():
            out.append(f'        <DataArray type="Float32" Name="{name}" '
                       f'format="ascii">')
            out.append(_fmt(arr, 6))
            out.append("        </DataArray>")
        out.append("      </PointData>")

    out.append("    </Piece>")
    out.append("  </UnstructuredGrid>")
    out.append("</VTKFile>")
    Path(path).write_text("\n".join(out) + "\n")


def _fmt(arr, decimals):
    a = np.asarray(arr).reshape(-1)
    if decimals == 0:
        return "          " + " ".join(str(int(x)) for x in a)
    return "          " + " ".join(np.format_float_scientific(x, precision=decimals)
                                   if np.isfinite(x) else "nan" for x in a)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="forge res_*.h5 / mesh.h5")
    ap.add_argument("output", nargs="?", default=None, help="出力 .vtu (省略時は同名)")
    ap.add_argument("--fields", nargs="*", default=None,
                    help="出力する VALUE 名 (省略時は全て)")
    args = ap.parse_args()

    out = args.output or str(Path(args.input).with_suffix(".vtu"))
    with h5py.File(args.input, "r") as f:
        coord = np.array(f["MESH/COORD"], dtype=np.float64).reshape(-1, 3)
        conne = np.array(f["MESH/CONNE"])
        # nCells: 属性が無い res_*.h5 もあるので VALUE 配列長 → 属性 → CONNE の順で決める
        # node 離散化では VALUE 長 = **節点数** でセル数と違うので、VALUE 長を無条件に
        # nCells とみなしてはいけない (2026-09-06 修正: case/46 の node run が
        # 「CONNE が NumberOfElements より短い」で弾かれていた)。CONNE の走査を第一の根拠にする。
        ncells = None
        if "MESH" in f and "nCells" in f["MESH"].attrs:
            ncells = int(np.array(f["MESH"].attrs["nCells"]).ravel()[0])
        if ncells is None:
            ncells = _count_cells(conne)
        nnodes = coord.shape[0]
        conn, offsets, types = parse_conne(conne, ncells)

        # VALUE 長で cell データか point データかを振り分ける
        cell_data, point_data = {}, {}
        if "VALUE" in f:
            names = args.fields if args.fields else list(f["VALUE"].keys())
            for name in names:
                if name not in f["VALUE"]:
                    print(f"  警告: VALUE/{name} が無い")
                    continue
                arr = np.array(f["VALUE/" + name], dtype=np.float32)
                if arr.size == ncells:
                    cell_data[name] = arr
                elif arr.size == nnodes:      # node 離散化の節点値
                    point_data[name] = arr

    write_vtu(out, coord, conn, offsets, types, cell_data, point_data)
    print(f"-> {out}  (points={len(coord)}, cells={ncells}, "
          f"cell fields={len(cell_data)}, point fields={len(point_data)})")


if __name__ == "__main__":
    main()
