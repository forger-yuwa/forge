#!/usr/bin/env python3
"""直方体の構造化 hex メッシュ (gmsh 4.1) を書く。physID 1:xmin 2:xmax 3:ymin 4:ymax 5:zmin 6:zmax、体 7。
  python3 make_box_msh.py OUT.msh NX NY NZ LX LY LZ
(case/33 の make_wavy_node_msh.py と同じ書式。遠方境界 farfield の検証用、plan boundary-node-farfield-characteristic §6)
"""
import sys
import numpy as np


def write_hex_msh(out, nx, ny, nz, coord, bbox):
    """構造 hex 格子を書く。coord(i, j, k) -> (x, y, z)、bbox = (Lx, Ly, Lz) (Entities の包絡に使うだけ)"""
    Lx, Ly, Lz = bbox
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    hexes = [[nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
              nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
             for k in range(nz) for j in range(ny) for i in range(nx)]
    quads = {1: [], 2: [], 3: [], 4: [], 5: [], 6: []}
    for k in range(nz):
        for j in range(ny):
            quads[1].append([nid(0, j, k), nid(0, j, k + 1), nid(0, j + 1, k + 1), nid(0, j + 1, k)])
            quads[2].append([nid(nx, j, k), nid(nx, j + 1, k), nid(nx, j + 1, k + 1), nid(nx, j, k + 1)])
    for k in range(nz):
        for i in range(nx):
            quads[3].append([nid(i, 0, k), nid(i + 1, 0, k), nid(i + 1, 0, k + 1), nid(i, 0, k + 1)])
            quads[4].append([nid(i, ny, k), nid(i, ny, k + 1), nid(i + 1, ny, k + 1), nid(i + 1, ny, k)])
    for j in range(ny):
        for i in range(nx):
            quads[5].append([nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0)])
            quads[6].append([nid(i, j, nz), nid(i + 1, j, nz), nid(i + 1, j + 1, nz), nid(i, j + 1, nz)])
    nN = (nx + 1) * (ny + 1) * (nz + 1); nC = len(hexes); nQ = sum(len(v) for v in quads.values())
    names = {1: "xmin", 2: "xmax", 3: "ymin", 4: "ymax", 5: "zmin", 6: "zmax"}
    with open(out, "w") as w:
        w.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n7\n")
        for p in range(1, 7):
            w.write(f'2 {p} "{names[p]}"\n')
        w.write('3 7 "fluid"\n$EndPhysicalNames\n$Entities\n0 0 6 1\n')
        for p in range(1, 7):
            w.write(f"{p} 0 0 0 {Lx} {Ly} {Lz} 1 {p} 0\n")
        w.write(f"7 0 0 0 {Lx} {Ly} {Lz} 1 7 0\n$EndEntities\n")
        w.write(f"$Nodes\n1 {nN} 1 {nN}\n3 7 0 {nN}\n")
        for i in range(nN):
            w.write(f"{i + 1}\n")
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    x, y, z = coord(i, j, k)
                    w.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        w.write(f"$EndNodes\n$Elements\n7 {nQ + nC} 1 {nQ + nC}\n")
        eid = 1
        for p in range(1, 7):
            w.write(f"2 {p} 3 {len(quads[p])}\n")
            for q in quads[p]:
                w.write(f"{eid} " + " ".join(map(str, q)) + "\n"); eid += 1
        w.write(f"3 7 5 {nC}\n")
        for h in hexes:
            w.write(f"{eid} " + " ".join(map(str, h)) + "\n"); eid += 1
        w.write("$EndElements\n")
    print(f"wrote {out}: {nN} nodes, {nC} hex, {nQ} boundary quads")


def main():
    out = sys.argv[1]; nx, ny, nz = (int(v) for v in sys.argv[2:5]); Lx, Ly, Lz = (float(v) for v in sys.argv[5:8])
    xs, ys, zs = np.linspace(0, Lx, nx + 1), np.linspace(0, Ly, ny + 1), np.linspace(0, Lz, nz + 1)
    write_hex_msh(out, nx, ny, nz, lambda i, j, k: (xs[i], ys[j], zs[k]), (Lx, Ly, Lz))


if __name__ == "__main__":
    main()
