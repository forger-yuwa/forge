#!/usr/bin/env python3
"""Salome MED メッシュ -> forge 互換 gmsh 4.1 変換 (tet/prism/pyramid + 面グループ対応)。

使い方 (forge の mesh venv で):
    python3 med_to_msh41.py pintle_salome.med forge.msh

- meshio で MED を読み、面要素 (triangle/quad) の MED グループ名を physID にマップ、
  体要素 (tetra/wedge/pyramid) は fluid にまとめる。
- 節点順: meshio (VTK 準拠) の wedge/pyramid は gmsh の prism/pyramid と同一なので変換不要
  (solver_density_cuda/output/res_h5_to_vtu.py の対応表と同じ前提)。
- convertGmshToForge の必須要件 (gmsh4.1 / 1エンティティ=物理タグ1 / 体に3D物理グループ /
  外表面全カバー / 線形要素) を満たす形式で書く。physID は mesh_pintle.py と同一。
"""
import sys
from collections import defaultdict

import meshio
import numpy as np

# 既知グループ -> physID。MED に存在するものだけ書き出す (outlet はプルーム無し形状のみ、
# plume_* / base はプルームあり形状のみ)。physID は bcondConfig.yaml と一致させること。
PHYS = {"fluid": (3, 1), "wall": (2, 10), "inlet": (2, 11),
        "outlet": (2, 12), "symmetry": (2, 13), "wall_pintle": (2, 14),
        "plume_out": (2, 15), "plume_far": (2, 16), "base": (2, 17)}
SURF_ENT = {"wall": 1, "inlet": 2, "outlet": 3, "symmetry": 4, "wall_pintle": 5,
            "plume_out": 6, "plume_far": 7, "base": 8}
GMSH_TYPE = {"triangle": (2, 3), "quad": (3, 4),
             "tetra": (4, 4), "hexahedron": (5, 8), "wedge": (6, 6), "pyramid": (7, 5)}
SURF_TYPES = ("triangle", "quad")
VOL_TYPES = ("tetra", "wedge", "pyramid", "hexahedron")


def med_groups(mesh):
    """meshio MED: cell block ごとの int tag 配列 + tag->名前 から、block ごとの名前配列を返す。"""
    tags = mesh.cell_data.get("cell_tags")
    if tags is None:
        raise SystemExit("ERROR: MED に cell_tags が無い (グループ未定義?)")
    names = {}
    for tag, nm in mesh.cell_tags.items():
        names[tag] = nm[0] if nm else str(tag)
    return [np.array([names.get(int(t), "") for t in blk]) for blk in tags]


def main():
    src, dst = sys.argv[1], sys.argv[2]
    mesh = meshio.read(src)
    pts = mesh.points
    tagnames = med_groups(mesh)

    surf = defaultdict(lambda: defaultdict(list))   # name -> mtype -> [conn]
    vol = defaultdict(list)                         # mtype -> [conn]
    for blk, nmarr in zip(mesh.cells, tagnames):
        if blk.type in SURF_TYPES:
            for row, nm in zip(blk.data, nmarr):
                if nm in SURF_ENT:
                    surf[nm][blk.type].append(row)
        elif blk.type in VOL_TYPES:
            vol[blk.type].extend(blk.data)
    # MED に存在するグループのみ書く (最低限 inlet/wall/symmetry/wall_pintle は必須)
    need = {"inlet", "wall"} - set(surf)
    if need:
        raise SystemExit("ERROR: MED に必須面グループ %s が無い" % need)
    unknown = set(surf) - set(SURF_ENT)
    if unknown:
        raise SystemExit("ERROR: 未知の面グループ %s (PHYS/SURF_ENT に追加を)" % unknown)
    ents = {nm: et for nm, et in SURF_ENT.items() if nm in surf}

    nsurf = sum(len(v) for d in surf.values() for v in d.values())
    nvol = sum(len(v) for v in vol.values())
    nblocks = sum(len(d) for d in surf.values()) + len(vol)
    N = len(pts)
    print("nodes=%d surf=%d vol=%s" % (N, nsurf, {k: len(v) for k, v in vol.items()}))

    with open(dst, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        pnames = {nm: PHYS[nm] for nm in list(ents) + ["fluid"]}
        f.write("$PhysicalNames\n%d\n" % len(pnames))
        for nm, (d, pid) in pnames.items():
            f.write('%d %d "%s"\n' % (d, pid, nm))
        f.write("$EndPhysicalNames\n")
        f.write("$Entities\n0 0 %d 1\n" % len(ents))
        for nm, et in sorted(ents.items(), key=lambda kv: kv[1]):
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (et, PHYS[nm][1]))
        f.write("1 0 0 0 0 0 0 1 %d 0\n" % PHYS["fluid"][1])
        f.write("$EndEntities\n")
        f.write("$Nodes\n1 %d 1 %d\n" % (N, N))
        f.write("3 1 0 %d\n" % N)
        for i in range(1, N + 1):
            f.write("%d\n" % i)
        for x, y, z in pts:
            f.write("%.10g %.10g %.10g\n" % (x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d %d 1 %d\n" % (nblocks, nsurf + nvol, nsurf + nvol))
        tag = 1
        for nm, et in sorted(ents.items(), key=lambda kv: kv[1]):
            for mtype, rows in surf[nm].items():
                gt, nn = GMSH_TYPE[mtype]
                f.write("2 %d %d %d\n" % (et, gt, len(rows)))
                for r in rows:
                    f.write("%d %s\n" % (tag, " ".join(str(int(v) + 1) for v in r)))
                    tag += 1
        for mtype, rows in vol.items():
            gt, nn = GMSH_TYPE[mtype]
            f.write("3 1 %d %d\n" % (gt, len(rows)))
            for r in rows:
                f.write("%d %s\n" % (tag, " ".join(str(int(v) + 1) for v in r)))
                tag += 1
        f.write("$EndElements\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()
