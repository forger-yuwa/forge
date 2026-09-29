#!/usr/bin/env python3
r"""case/64 の流体メッシュ (子午面、node・全四角・構造)。壁は全長が 1 つの共役壁 (physID 3)。

    r=R  ┌──────────────── wall (3、全長が共役) ─────────────────┐
         │ inlet (1)                                              │ outlet (2)
    r=0  └──────────────────── axis (4) ────────────────────────────┘
        x=−80R              x=0 (加熱開始)   x=10R (加熱終了)   x=20R

半径 `--nr` 一様 (16/32/64)。軸方向は N_r=32 の点列を基準に入れ子で細分 (64 は 2 等分、16 は 2 つずつ併合)。
点列: 上流は x=0 に向かって細かく (0.02 mm)、加熱区間は両端 (Robin の端) で細かく中央で粗く、下流は x=10R から。

    python3 case/64.conjugate_pipe_wall/gen_mesh.py --nr 32
"""
from __future__ import annotations
import argparse
import numpy as np
import pipe_common as pc
from pipe_common import axcht

GEO_HEAD = """// case/64 厚肉管 (子午面)。gen_mesh.py が生成。
Point(1) = {{{xa}, 0, 0, 1}}; Point(2) = {{{xa}, {R}, 0, 1}};
Line(1) = {{1, 2}};
Transfinite Line{{1}} = {nr} + 1;
"""
EXTRUDE = "o{k}[] = Extrude {{{L}, 0, 0}} {{ Line{{{src}}}; Layers{{ {{{ones}}}, {{{fr}}} }}; Recombine; }};\n"
GEO_TAIL = """// Extrude の戻り値は向き付き: Abs()。[2] = 点 2 側 (壁 r=R)、[3] = 点 1 側 (軸)。case/63 で確認済み
Physical Curve("inlet",  1) = {{1}};
Physical Curve("outlet", 2) = {{Abs(o3[0])}};
Physical Curve("wall",   3) = {{Abs(o1[2]), Abs(o2[2]), Abs(o3[2])}};
Physical Curve("axis",   4) = {{Abs(o1[3]), Abs(o2[3]), Abs(o3[3])}};
Physical Surface("fluid", 5) = {{o1[1], o2[1], o3[1]}};
Mesh.MshFileVersion = 4.1;
"""


def ratio_for(L, n, h0):
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi); s = h0 * (q ** n - 1) / (q - 1)
        lo, hi = (q, hi) if s < L else (lo, q)
    return 0.5 * (lo + hi)


def capped(h0, q, hmax, L):
    """幅 h0 から等比 q で hmax まで伸ばし、残りを幅 ≈hmax の一様で埋める (総数は偶数)。
    等比のまま hmax に届く前に区間を埋めてしまう場合は、等比数列だけで区間ちょうどに合わせる
    (初版は最後の要素が負になり、下流 10R で変換が壊れた)。"""
    w = []
    while w == [] or (w[-1] * q <= hmax and sum(w) + w[-1] * q <= L):
        w.append(h0 * q ** len(w))
    rest = L - sum(w)
    if rest < hmax:                                        # 等比だけで埋める
        n = len(w) + 1
        if n % 2: n += 1
        qq = ratio_for(L, n, h0)
        return h0 * qq ** np.arange(n)
    n2 = max(1, round(rest / hmax))
    if (len(w) + n2) % 2: n2 += 1
    return np.array(w + [rest / n2] * n2)


def base_widths():
    """N_r=32 の区間幅 (各区間の数は偶数)。予熱が上流へ数十 mm 届く (A2) ので上下流は 1 mm で頭打ち、加熱区間は両端で細かく中央 0.25 mm。"""
    h0 = 0.02e-3
    up = capped(h0, 1.08, 1.0e-3, pc.L_UP)[::-1]
    half = capped(h0, 1.08, 0.25e-3, pc.L_HEAT / 2)
    heat = np.concatenate([half, half[::-1]])
    down = capped(h0, 1.08, 1.0e-3, pc.L_DOWN)
    return up, heat, down


def refine(w, f):
    if f == 1: return w
    if f == 2: return np.repeat(w / 2, 2)
    if f == 0.5:
        assert len(w) % 2 == 0; return w.reshape(-1, 2).sum(axis=1)
    raise SystemExit("--nr は 16 / 32 / 64")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nr", type=int, default=32)
    a = ap.parse_args()
    f = a.nr / 32.0
    segs = [refine(w, f) for w in base_widths()]
    name = f"pipe_r{a.nr}"
    geo = GEO_HEAD.format(xa=-pc.L_UP, R=pc.R, nr=a.nr); src = "1"
    for k, (w, L) in enumerate(zip(segs, (pc.L_UP, pc.L_HEAT, pc.L_DOWN)), start=1):
        fr = np.cumsum(w) / w.sum(); fr[-1] = 1.0
        geo += EXTRUDE.format(k=k, L=L, src=src, ones=",".join(["1"] * len(w)), fr=",".join(f"{v:.15g}" for v in fr))
        src = f"o{k}[0]"
    geo += GEO_TAIL.format()
    tdir = pc.HERE / "template"
    h5 = axcht.gmsh_and_convert(pc.HERE / "mesh", name, geo, (tdir / "solverConfig_dry.yaml").read_text(),
                                (tdir / "bcondConfig_dry.yaml").read_text())
    dr = pc.R / a.nr
    allw = np.concatenate(segs)
    print(f"[gen_mesh] {h5}: nr {a.nr} × 軸 {len(segs[0])}+{len(segs[1])}+{len(segs[2])}、最小幅 {allw.min()*1e6:.2f} um、"
          f"最大幅 {allw.max()*1e3:.3f} mm (AR {allw.max()/dr:.1f})、隣接比の最大 {max(np.max(np.maximum(w[1:]/w[:-1], w[:-1]/w[1:])) for w in segs):.4f}")
    q = axcht.mesh_quality(h5); print(q.splitlines()[0]); (pc.HERE / "mesh" / f"{name}.quality.txt").write_text(q)


if __name__ == "__main__":
    main()
