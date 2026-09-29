#!/usr/bin/env python3
r"""case/65 共役平板の流体メッシュ (平面 2D、node・全四角・構造)。

    y=H  ┌──────────────────────── top (3、slip) ────────────────────────┐
         │ inlet (1)                                                     │ outlet (2)
    y=0  └─ slip_up (4) ─┬──────── plate (5、共役壁) ────────┬─ slip_down (6) ─┘
        x=−L/2          x=0 (前縁)                         x=L (後縁)       x=3L/2

`--ny` (16/32/64): y 方向は N=32 の点列 (壁で 5 µm から等比) を基準に入れ子で細分。x 方向も同様 (前縁・後縁で細かく)。
点列は明示 (点を y に押し出し → 線を x に押し出し) なので、3 水準が完全に入れ子になる。

    python3 gen_mesh.py --ny 32
"""
from __future__ import annotations
import argparse
import numpy as np
import plate_common as pc
from plate_common import axcht


def ratio_for(L, n, h0):
    lo, hi = 1.0 + 1e-12, 3.0
    for _ in range(200):
        q = 0.5 * (lo + hi); s = h0 * (q ** n - 1) / (q - 1)
        lo, hi = (q, hi) if s < L else (lo, q)
    return 0.5 * (lo + hi)


def geo_n(L, n, h0):
    q = ratio_for(L, n, h0); w = h0 * q ** np.arange(n); return w * (L / w.sum())


def base():
    """N=32: y は壁 5 µm から 32 セル。x は上流 20 (前縁へ細かく)、板 80 (両端で細かく)、下流 20 (後縁から)。"""
    y = geo_n(pc.H_TOP, 32, 5e-6)
    h0 = 0.004 * pc.L
    up = geo_n(pc.L_UP, 20, h0)[::-1]
    half = geo_n(pc.L / 2, 40, h0); plate = np.concatenate([half, half[::-1]])
    down = geo_n(pc.L_DOWN, 20, h0)
    return y, up, plate, down


def refine(w, f):
    if f == 1: return w
    if f == 2: return np.repeat(w / 2, 2)
    if f == 0.5: return w.reshape(-1, 2).sum(axis=1)
    raise SystemExit("--ny は 16 / 32 / 64")


def frac(w):
    fr = np.cumsum(w) / w.sum(); fr[-1] = 1.0; return ",".join(f"{v:.15g}" for v in fr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ny", type=int, default=32)
    a = ap.parse_args()
    f = a.ny / 32.0
    y, up, plate, down = (refine(w, f) for w in base())
    geo = f"""// case/65 共役平板。gen_mesh.py が生成。
Point(1) = {{{-pc.L_UP}, 0, 0, 1}};
v[] = Extrude {{0, {pc.H_TOP}, 0}} {{ Point{{1}}; Layers{{ {{{",".join(["1"]*len(y))}}}, {{{frac(y)}}} }}; }};
"""
    src = "v[1]"
    for k, (w, L) in enumerate(((up, pc.L_UP), (plate, pc.L), (down, pc.L_DOWN)), start=1):
        geo += f"o{k}[] = Extrude {{{L}, 0, 0}} {{ Line{{{src}}}; Layers{{ {{{','.join(['1']*len(w))}}}, {{{frac(w)}}} }}; Recombine; }};\n"
        src = f"o{k}[0]"
    geo += """// 源の線は y=0 → y=H。[2] = 終点側 (上 y=H)、[3] = 始点側 (下 y=0) (case/63 と同じ規則。乾式の壁ダンプで確認する)
Physical Curve("inlet", 1) = {Abs(v[1])};
Physical Curve("outlet", 2) = {Abs(o3[0])};
Physical Curve("top", 3) = {Abs(o1[2]), Abs(o2[2]), Abs(o3[2])};
Physical Curve("slip_up", 4) = {Abs(o1[3])};
Physical Curve("plate", 5) = {Abs(o2[3])};
Physical Curve("slip_down", 6) = {Abs(o3[3])};
Physical Surface("fluid", 7) = {o1[1], o2[1], o3[1]};
Mesh.MshFileVersion = 4.1;
"""
    tdir = pc.HERE / "template"
    name = f"plate_n{a.ny}"
    h5 = axcht.gmsh_and_convert(pc.HERE / "mesh", name, geo, (tdir / "solverConfig_dry.yaml").read_text(), (tdir / "bcondConfig_dry.yaml").read_text())
    allx = np.concatenate([up, plate, down])
    print(f"[gen_mesh] {h5}: y {len(y)} セル (壁 {y[0]*1e6:.2f} µm)、x {len(up)}+{len(plate)}+{len(down)} (最小 {allx.min()*1e6:.1f} µm)、"
          f"AR 最大 {max(allx.max()/y[0], allx.min()/y[0], y.max()/allx.min()):.0f}")
    q = axcht.mesh_quality(h5); print(q.splitlines()[0]); (pc.HERE / "mesh" / f"{name}.quality.txt").write_text(q)


if __name__ == "__main__":
    main()
