#!/usr/bin/env python3
r"""case/48 の固体帯を **SU2 の固体ゾーン** (.su2) にする (CHT V2 の case/48 段)。

plan boundary-conjugate-heat-transfer §6 V2。流体ゾーンは既存の `su2_B_tw300/mesh.su2` を使い、
その `wall` マーカの $x$ を読んで**界面節点を 1 対 1 に合わせた**固体帯を作る。

固体 (forge の `conjugate:` と同じ):  t=1.0e-3 m, k_s=0.217 W/mK, 背面 300 K (Dirichlet)

SU2 側は背面を `MARKER_ISOTHERMAL` で厳密な Dirichlet にできる (forge の `fem2d` は Robin しか
持てないので h=1e8 で近似している。$1/h$ は $t/k_s$ の 2.2e-6 倍なので差は出ない)。

使い方:
  python3 case/48.flat_plate_cooled_m4/gen_su2_solid.py \
      --fluid-mesh <...>/su2_B_tw300/mesh.su2 --nl 16 --out <...>/su2_cht/solid.su2
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def wall_x(mesh: Path, tag: str) -> np.ndarray:
    coords = None
    ids: set[int] = set()
    with open(mesh) as f:
        it = iter(f)
        for line in it:
            t = line.strip()
            if t.startswith("NPOIN="):
                n = int(t.split("=")[1].split()[0])
                coords = np.empty((n, 2))
                for i in range(n):
                    w = next(it).split()
                    coords[i] = (float(w[0]), float(w[1]))
            elif t.startswith("MARKER_TAG=") and t.split("=")[1].strip() == tag:
                ne = int(next(it).split("=")[1])
                for _ in range(ne):
                    w = next(it).split()
                    ids.update(int(v) for v in w[1:3])
    if coords is None or not ids:
        raise SystemExit(f"{mesh}: マーカ {tag} が見つからない")
    return np.sort(coords[sorted(ids)][:, 0])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fluid-mesh", required=True)
    ap.add_argument("--wall-tag", default="wall")
    ap.add_argument("--t", type=float, default=1.0e-3)
    ap.add_argument("--nl", type=int, default=16)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    x = wall_x(Path(a.fluid_mesh), a.wall_tag)
    nx, nl = len(x), a.nl
    ys = np.linspace(-a.t, 0.0, nl + 1)
    nid = lambda i, j: j * nx + i                      # noqa: E731

    lines = ["NDIME= 2", f"NELEM= {(nx-1)*nl}"]
    e = 0
    for j in range(nl):
        for i in range(nx - 1):
            lines.append(f"9 {nid(i,j)} {nid(i+1,j)} {nid(i+1,j+1)} {nid(i,j+1)} {e}")
            e += 1
    lines.append(f"NPOIN= {nx*(nl+1)}")
    for j in range(nl + 1):
        for i in range(nx):
            lines.append(f"{x[i]:.16e} {ys[j]:.16e} {nid(i,j)}")

    marks = [
        ("solid_back", [(nid(i, 0), nid(i + 1, 0)) for i in range(nx - 1)]),
        ("solid_interface", [(nid(i, nl), nid(i + 1, nl)) for i in range(nx - 1)]),
        ("solid_left", [(nid(0, j), nid(0, j + 1)) for j in range(nl)]),
        ("solid_right", [(nid(nx - 1, j), nid(nx - 1, j + 1)) for j in range(nl)]),
    ]
    lines.append(f"NMARK= {len(marks)}")
    for name, edges in marks:
        lines.append(f"MARKER_TAG= {name}")
        lines.append(f"MARKER_ELEMS= {len(edges)}")
        for p, q in edges:
            lines.append(f"3 {p} {q}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(f"[gen_su2_solid] {out}")
    print(f"  節点 {nx*(nl+1)} / quad {(nx-1)*nl}   界面 {nx} 点 (流体の {a.wall_tag} と 1 対 1)")
    print(f"  厚さ {a.t} m を {nl} 層,  x {x.min():.4f} .. {x.max():.4f}")


if __name__ == "__main__":
    main()
