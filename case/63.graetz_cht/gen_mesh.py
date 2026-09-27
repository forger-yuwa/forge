#!/usr/bin/env python3
r"""case/63 Graetz の流体メッシュ (子午面、node・全四角・構造)。

plan [`boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md) §4.1・§5 #4。

    r=R ┌─ wall_up (3) ─┬──── wall_heat (4, 共役 / 等温) ────┬─ wall_down (5) ─┐
        │ inlet (1)    │                                     │                 │ outlet (2)
    r=0 └──────────────┴──────────── axis (6) ───────────────┴─────────────────┘
       x=−L_up        x=0                                 x=L_heat          x=L_heat+L_down

- 半径方向 `--nr` 一様 (登録 32、感度 16/64)。軸方向は加熱開始点 x=0 に向けて両側から等比で細かくする。
  **軸方向も N_r に比例して入れ子で細分** (格子感度で 2 方向を同時に refine する): N_r = 32 の点列
  (上流 48 セル・x=0 の幅 0.02 mm、加熱は 0.02 mm から等比 1.05 で 1 mm まで伸ばして一様、下流 1 mm × 10) を基準に、
  N_r = 64 は各区間を 2 等分、N_r = 16 は隣り合う 2 区間を併合する。
- 変換後に IC (Poiseuille + 線形圧力 + T_in) をパッチする (`make_run.py` 側。ここでは変換だけ)。

    python3 case/63.graetz_cht/gen_mesh.py --nr 32      # mesh/graetz_r32.h5
"""
from __future__ import annotations

import argparse

import numpy as np

import graetz_common as gc
from graetz_common import axcht

GEO_HEAD = """// case/63 Graetz 管 (子午面)。gen_mesh.py が生成。編集しないこと。
// 入口の半径線を x 方向に 3 回押し出す (上流・加熱・下流)。層の位置は gen_mesh.py が明示する (入れ子の細分)。
Point(1) = {{{xa}, 0, 0, 1}}; Point(2) = {{{xa}, {R}, 0, 1}};
Line(1) = {{1, 2}};
Transfinite Line{{1}} = {nr} + 1;
"""
EXTRUDE = """o{k}[] = Extrude {{{L}, 0, 0}} {{ Line{{{src}}}; Layers{{ {{{ones}}}, {{{fr}}} }}; Recombine; }};
"""
GEO_TAIL = """// Extrude の戻り値は向き付き (負号あり) なので Abs() を取る (負のまま渡すと変換器が physID を負で読む)
// o{{k}}[0] = 押し出し先の半径線、[1] = 面、[2] = 点 2 側 (壁 r=R)、[3] = 点 1 側 (軸 r=0)
// (2026-09-27 乾式の壁ダンプ座標で確認: 逆に書いた初版は壁が y=0 に出た。確認は乾式 1 step の壁ダンプで行う: check_dry.py)
Physical Curve("inlet",     1) = {{1}};
Physical Curve("outlet",    2) = {{Abs(o3[0])}};
Physical Curve("wall_up",   3) = {{Abs(o1[2])}};
Physical Curve("wall_heat", 4) = {{Abs(o2[2])}};
Physical Curve("wall_down", 5) = {{Abs(o3[2])}};
Physical Curve("axis",      6) = {{Abs(o1[3]), Abs(o2[3]), Abs(o3[3])}};
Physical Surface("fluid",   7) = {{o1[1], o2[1], o3[1]}};
Mesh.MshFileVersion = 4.1;
"""


def ratio_for(L: float, n: int, h0: float) -> float:
    """h0 (1 + q + … + q^{n-1}) = L を満たす等比 q (>1)。"""
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        s = h0 * (q ** n - 1.0) / (q - 1.0)
        lo, hi = (q, hi) if s < L else (lo, q)
    return 0.5 * (lo + hi)


def geometric(h0: float, q: float, hmax: float, L: float):
    """幅 h0 から等比 q で hmax まで伸ばし、残りを幅 ≈hmax の一様で埋めた区間幅の列 (合計 L)。"""
    w = []
    while (not w or w[-1] * q <= hmax) and sum(w) < L:
        w.append(h0 * q ** len(w))
    rest = L - sum(w)
    n2 = max(1, round(rest / hmax))
    if (len(w) + n2) % 2:                               # 総数を偶数に (N_r = 16 で 2 つずつ併合するため)
        n2 += 1
    return np.array(w + [rest / n2] * n2)


def base_widths():
    """f = 1 (N_r = 32) の軸方向の区間幅。各区間のセル数は偶数 (N_r = 16 で 2 つずつ併合するため)。"""
    h0 = 0.02e-3
    heat = geometric(h0, 1.05, 1.0e-3, gc.L_HEAT)
    nu = 48
    qu = ratio_for(gc.L_UP, nu, h0)
    up = (h0 * qu ** np.arange(nu))[::-1]               # x=0 に向かって細かく
    down = np.full(10, gc.L_DOWN / 10)
    return up, heat, down


def refine(w: np.ndarray, f: float) -> np.ndarray:
    """f = 2 は各区間を 2 等分、f = 0.5 は隣り合う 2 区間を併合、f = 1 はそのまま (入れ子の細分)。"""
    if f == 1:
        return w
    if f == 2:
        return np.repeat(w / 2, 2)
    if f == 0.5:
        assert len(w) % 2 == 0
        return w.reshape(-1, 2).sum(axis=1)
    raise SystemExit("--nr は 16 / 32 / 64 のみ (入れ子の細分)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nr", type=int, default=32)
    ap.add_argument("--name", default=None)
    a = ap.parse_args()
    f = a.nr / 32.0
    segs = [refine(w, f) for w in base_widths()]
    name = a.name or f"graetz_r{a.nr}"
    geo = GEO_HEAD.format(xa=-gc.L_UP, R=gc.R, nr=a.nr)
    src = "1"
    for k, (w, L) in enumerate(zip(segs, (gc.L_UP, gc.L_HEAT, gc.L_DOWN)), start=1):
        fr = np.cumsum(w) / w.sum()
        fr[-1] = 1.0
        geo += EXTRUDE.format(k=k, L=L, src=src, ones=",".join(["1"] * len(w)),
                              fr=",".join(f"{v:.15g}" for v in fr))
        src = f"o{k}[0]"
    geo += GEO_TAIL.format()
    tdir = gc.HERE / "template"
    h5 = axcht.gmsh_and_convert(gc.HERE / "mesh", name, geo,
                                (tdir / "solverConfig_dry.yaml").read_text(), (tdir / "bcondConfig_dry.yaml").read_text())
    dr = gc.R / a.nr
    up, heat, down = segs
    print(f"[gen_mesh] {h5}: nr {a.nr} (dr {dr*1e6:.3f} um) × 軸 {len(up)}+{len(heat)}+{len(down)} セル、"
          f"x=0 の幅 上流側 {up[-1]*1e6:.2f} / 加熱側 {heat[0]*1e6:.2f} um、隣接比の最大 "
          f"{max(np.max(np.maximum(w[1:]/w[:-1], w[:-1]/w[1:])) for w in segs):.4f}、"
          f"最大幅 {max(w.max() for w in segs)*1e3:.4f} mm (AR {max(w.max() for w in segs)/dr:.1f})")
    np.savez(gc.HERE / "mesh" / f"{name}.xwidths.npz", up=up, heat=heat, down=down)
    q = axcht.mesh_quality(h5)
    print(q.splitlines()[0])
    (gc.HERE / "mesh" / f"{name}.quality.txt").write_text(q)


if __name__ == "__main__":
    main()
