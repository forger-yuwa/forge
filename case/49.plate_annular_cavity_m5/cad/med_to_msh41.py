#!/usr/bin/env python3
"""Salome MED -> forge 互換 gmsh 4.1 (tet/prism/pyramid + 面グループ) [case/49]

usage (mesh venv):  python3 med_to_msh41.py cavity_salome.med forge.msh

グループ -> physID は **../manifest.json の `phys_id`** を読む (bcondConfig と同じ単一ソース)。
派生元: case/37.pintle_nozzle/cad/med_to_msh41.py (physID のハードコードを manifest 化)。

convertGmshToForge の必須要件:
  gmsh 4.1 / 1 エンティティ = 物理タグ 1 / 体に 3D 物理グループ `fluid` / 外表面を全カバー / 線形要素
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import meshio
import numpy as np

HERE = Path(__file__).resolve().parent
MAN = json.loads((HERE.parent / "manifest.json").read_text())
PHYS = MAN["phys_id"]                      # 面グループ -> physID
FLUID_ID = MAN["fluid_id"]

GMSH_TYPE = {"triangle": (2, 3), "quad": (3, 4),
             "tetra": (4, 4), "hexahedron": (5, 8), "wedge": (6, 6), "pyramid": (7, 5)}
SURF_TYPES = ("triangle", "quad")
VOL_TYPES = ("tetra", "wedge", "pyramid", "hexahedron")


def med_groups(mesh):
    tags = mesh.cell_data.get("cell_tags")
    if tags is None:
        raise SystemExit("ERROR: MED に cell_tags が無い (グループ未定義?)")
    names = {t: (nm[0] if nm else str(t)) for t, nm in mesh.cell_tags.items()}
    return [np.array([names.get(int(t), "") for t in blk]) for blk in tags]


def main():
    src, dst = sys.argv[1], sys.argv[2]
    mesh = meshio.read(src)
    pts = mesh.points
    tagnames = med_groups(mesh)

    surf = defaultdict(lambda: defaultdict(list))
    vol = defaultdict(list)
    for blk, nmarr in zip(mesh.cells, tagnames):
        if blk.type in SURF_TYPES:
            for row, nm in zip(blk.data, nmarr):
                if nm in PHYS:
                    surf[nm][blk.type].append(row)
        elif blk.type in VOL_TYPES:
            vol[blk.type].extend(blk.data)
    missing = set(PHYS) - set(surf)
    if missing:
        raise SystemExit("ERROR: MED に面グループ %s が無い" % sorted(missing))
    ents = {nm: i + 1 for i, nm in enumerate(sorted(surf, key=lambda k: PHYS[k]))}

    nsurf = sum(len(v) for d in surf.values() for v in d.values())
    nvol = sum(len(v) for v in vol.values())
    nblocks = sum(len(d) for d in surf.values()) + len(vol)
    N = len(pts)
    print("nodes=%d surf=%d vol=%s" % (N, nsurf, {k: len(v) for k, v in vol.items()}))
    print("groups:", {k: (PHYS[k], sum(len(v) for v in surf[k].values())) for k in sorted(surf)})

    with open(dst, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n%d\n" % (len(ents) + 1))
        for nm in sorted(ents, key=lambda k: PHYS[k]):
            f.write('2 %d "%s"\n' % (PHYS[nm], nm))
        f.write('3 %d "fluid"\n' % FLUID_ID)
        f.write("$EndPhysicalNames\n")
        f.write("$Entities\n0 0 %d 1\n" % len(ents))
        for nm in sorted(ents, key=lambda k: ents[k]):
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (ents[nm], PHYS[nm]))
        f.write("1 0 0 0 0 0 0 1 %d 0\n" % FLUID_ID)
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
        for nm in sorted(ents, key=lambda k: ents[k]):
            for mtype, rows in surf[nm].items():
                gt, _ = GMSH_TYPE[mtype]
                f.write("2 %d %d %d\n" % (ents[nm], gt, len(rows)))
                for r in rows:
                    f.write("%d %s\n" % (tag, " ".join(str(int(v) + 1) for v in r)))
                    tag += 1
        for mtype, rows in vol.items():
            gt, _ = GMSH_TYPE[mtype]
            f.write("3 1 %d %d\n" % (gt, len(rows)))
            for r in rows:
                f.write("%d %s\n" % (tag, " ".join(str(int(v) + 1) for v in r)))
                tag += 1
        f.write("$EndElements\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()
