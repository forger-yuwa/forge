#!/usr/bin/env python3
r"""C3X / Mark II 翼断面の**固体メッシュ**を作る (CHT の固体側 = `solid_fem2d.py` の入力)。

- 外形は `ref/vane_*_profile.csv` (表の点 + 前縁/後縁の円弧を復元した閉輪郭) を**弧長等間隔で再標本化**。
  各区間を Line + `Transfinite Curve = 2` にするので、**節点は指定した点にだけ**置かれる
  (流体メッシュと界面節点を 1 対 1 に合わせるための前提)。
- 冷却孔は `ref/cooling_holes_*.csv` の中心と直径から円を切り抜き、**孔ごとに物理曲線**にする
  (`fem2d` の Robin 辺を孔単位で持つため)。
- 出力: `mesh/solid_<vane>.npz` (nodes, tris, outer[], holes[[...]]) と品質のログ。

usage: python3 case/53.c3x_vane_cht/tools/gen_solid_mesh.py [--vane c3x|markii] [--n-outer 240] [--lc 0.12]
"""
import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REF = {"c3x": ROOT / "case/53.c3x_vane_cht/ref",
       "markii": ROOT / "case/54.markii_vane_cht/ref"}
OUT = {"c3x": ROOT / "case/53.c3x_vane_cht/mesh",
       "markii": ROOT / "case/54.markii_vane_cht/mesh"}


def load_profile(vane):
    f = REF[vane] / (f"vane_{vane}_profile.csv" if vane == "c3x" else "vane_markii_profile.csv")
    pts = []
    for line in open(f):
        if line.startswith("#") or line.startswith("i,"):
            continue
        _, x, y = line.split(",")
        pts.append((float(x), float(y)))
    return np.array(pts)


def load_holes(vane):
    f = REF[vane] / (f"cooling_holes_{vane}.csv")
    rows = []
    for line in open(f):
        if line.startswith("#") or line.startswith("hole"):
            continue
        p = line.split(",")
        rows.append((float(p[5]), float(p[6]), float(p[3])))    # x_cm, y_cm, D_cm
    return rows


def resample_closed(P, n, curv_weight=6.0, kappa_ref=2.0):
    """閉曲線を**曲率重み付き**で n 点に再標本化する (流体・固体で共通。界面節点を一致させる)。

    弧長等間隔だと**後縁 (R=1.7 mm) に 4 点しか載らず**、2 µm の境界層セルと組み合わさって
    step 6 で圧力が床に張り付いて発散した (実測 2026-09-20)。曲率 κ [1/cm] に応じて
    密度 ∝ 1 + curv_weight * min(κ/kappa_ref, 1) で詰める。
    """
    Q = np.vstack([P, P[:1]])
    seg = np.hypot(np.diff(Q[:, 0]), np.diff(Q[:, 1]))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    # 曲率: 3 点円の逆半径
    m = len(P)
    kap = np.zeros(m)
    for i in range(m):
        a, b, c = P[(i - 1) % m], P[i], P[(i + 1) % m]
        ab = np.linalg.norm(b - a); bc = np.linalg.norm(c - b); ca = np.linalg.norm(a - c)
        area2 = abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
        kap[i] = 0.0 if ab * bc * ca < 1e-12 else 2.0 * area2 / (ab * bc * ca)
    dens = 1.0 + curv_weight * np.minimum(kap / kappa_ref, 1.0)
    dens_s = np.concatenate([dens, dens[:1]])
    w = 0.5 * (dens_s[:-1] + dens_s[1:]) * seg          # 重み付き弧長
    W = np.concatenate([[0.0], np.cumsum(w)])
    t = np.linspace(0.0, W[-1], n, endpoint=False)
    su = np.interp(t, W, s)
    x = np.interp(su, s, Q[:, 0]); y = np.interp(su, s, Q[:, 1])
    return np.column_stack([x, y]), total


def write_geo(path, outer, holes, lc):
    L = []
    L.append("// C3X/Mark II 固体断面 (単位 cm を m に直して出力)")
    L.append(f"lc = {lc};")
    n = len(outer)
    for i, (x, y) in enumerate(outer, 1):
        L.append(f"Point({i}) = {{{x/100:.8f}, {y/100:.8f}, 0, lc}};")
    for i in range(1, n + 1):
        j = i % n + 1
        L.append(f"Line({i}) = {{{i}, {j}}};")
    L.append("Curve Loop(1) = {" + ",".join(str(i) for i in range(1, n + 1)) + "};")
    pid = n + 1
    loops = []
    hole_curves = []
    for k, (cx, cy, D) in enumerate(holes, 1):
        r = 0.5 * D
        c = pid; pid += 1
        L.append(f"Point({c}) = {{{cx/100:.8f}, {cy/100:.8f}, 0, lc}};")
        ptr = []
        for a in range(4):
            th = a * math.pi / 2
            L.append(f"Point({pid}) = {{{(cx + r*math.cos(th))/100:.8f}, "
                     f"{(cy + r*math.sin(th))/100:.8f}, 0, lc}};")
            ptr.append(pid); pid += 1
        cid0 = pid
        for a in range(4):
            L.append(f"Circle({pid}) = {{{ptr[a]}, {c}, {ptr[(a+1) % 4]}}};")
            pid += 1
        L.append(f"Curve Loop({100+k}) = {{{cid0},{cid0+1},{cid0+2},{cid0+3}}};")
        loops.append(100 + k)
        hole_curves.append([cid0, cid0 + 1, cid0 + 2, cid0 + 3])
        L.append(f"Transfinite Curve {{{cid0},{cid0+1},{cid0+2},{cid0+3}}} = 9;")
    L.append("Plane Surface(1) = {1," + ",".join(str(x) for x in loops) + "};")
    L.append("Transfinite Curve {" + ",".join(str(i) for i in range(1, n + 1)) + "} = 2;")
    L.append('Physical Curve("outer", 1) = {' + ",".join(str(i) for i in range(1, n + 1)) + "};")
    for k, cs in enumerate(hole_curves, 1):
        L.append(f'Physical Curve("hole{k}", {10+k}) = {{{",".join(str(c) for c in cs)}}};')
    L.append('Physical Surface("solid", 2) = {1};')
    L.append("Mesh.MshFileVersion = 4.1;")
    path.write_text("\n".join(L) + "\n")


def parse_msh41(path):
    """msh4.1 から nodes / 三角形 / physical curve ごとの辺 を取り出す最小パーサ。"""
    txt = path.read_text().splitlines()
    i = 0
    nodes = {}
    tris = []
    lines_by_phys = {}
    ent_curve_phys = {}
    while i < len(txt):
        ln = txt[i].strip()
        if ln == "$Entities":
            i += 1
            npnt, ncur, nsur, nvol = map(int, txt[i].split())
            i += 1
            for _ in range(npnt):
                i += 1
            for _ in range(ncur):
                f = txt[i].split()
                tag = int(f[0]); nphys = int(f[7])
                phys = [int(f[8 + k]) for k in range(nphys)]
                ent_curve_phys[tag] = phys
                i += 1
            while i < len(txt) and txt[i].strip() != "$EndEntities":
                i += 1
        elif ln == "$Nodes":
            i += 1
            nblk, nn, _, _ = map(int, txt[i].split()); i += 1
            for _ in range(nblk):
                dim, tag, param, num = map(int, txt[i].split()); i += 1
                ids = [int(txt[i + k]) for k in range(num)]
                i += num
                for k in range(num):
                    x, y, z = map(float, txt[i + k].split())
                    nodes[ids[k]] = (x, y)
                i += num
        elif ln == "$Elements":
            i += 1
            nblk, ne, _, _ = map(int, txt[i].split()); i += 1
            for _ in range(nblk):
                dim, tag, etype, num = map(int, txt[i].split()); i += 1
                for _ in range(num):
                    f = list(map(int, txt[i].split())); i += 1
                    if etype == 2:                       # 3-node triangle
                        tris.append(f[1:4])
                    elif etype == 1:                     # 2-node line
                        for p in ent_curve_phys.get(tag, []):
                            lines_by_phys.setdefault(p, []).append(f[1:3])
        i += 1
    ids = sorted(nodes)
    remap = {g: k for k, g in enumerate(ids)}
    xy = np.array([nodes[g] for g in ids])
    T = np.array([[remap[a] for a in t] for t in tris], int)
    P = {p: np.array([[remap[a], remap[b]] for a, b in e], int) for p, e in lines_by_phys.items()}
    return xy, T, P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="c3x", choices=["c3x", "markii"])
    ap.add_argument("--n-outer", type=int, default=240, help="外周の節点数 (流体メッシュと共有する)")
    ap.add_argument("--lc", type=float, default=0.0012, help="代表要素寸法 [m]")
    a = ap.parse_args()

    P = load_profile(a.vane)
    holes = load_holes(a.vane)
    outer, arc = resample_closed(P, a.n_outer)
    print(f"[{a.vane}] profile arc length {arc:.3f} cm -> {a.n_outer} nodes "
          f"(spacing {arc/a.n_outer*10:.2f} mm)")

    OUT[a.vane].mkdir(parents=True, exist_ok=True)
    geo = OUT[a.vane] / f"solid_{a.vane}.geo"
    msh = OUT[a.vane] / f"solid_{a.vane}.msh"
    write_geo(geo, outer, holes, a.lc)
    subprocess.run(["gmsh", "-2", str(geo), "-o", str(msh), "-format", "msh41", "-v", "1"], check=True)
    xy, tris, phys = parse_msh41(msh)
    outer_edges = phys.get(1, np.zeros((0, 2), int))
    outer_nodes = sorted(set(outer_edges.ravel().tolist()))
    hole_edges = {k: phys.get(10 + k, np.zeros((0, 2), int)) for k in range(1, len(holes) + 1)}
    print(f"[{a.vane}] mesh: {len(xy)} nodes, {len(tris)} triangles, "
          f"outer {len(outer_nodes)} nodes / {len(outer_edges)} edges, "
          f"holes {[len(v) for v in hole_edges.values()]} edges")
    if len(outer_nodes) != a.n_outer:
        sys.exit(f"[{a.vane}] 外周節点が {len(outer_nodes)} で指定 {a.n_outer} と違う "
                 "(Transfinite が効いていない)")
    npz = OUT[a.vane] / f"solid_{a.vane}.npz"
    np.savez(npz, nodes=xy, tris=tris, outer_edges=outer_edges,
             **{f"hole{k}": v for k, v in hole_edges.items()})
    print(f"  -> {npz.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
