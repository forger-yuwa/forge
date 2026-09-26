#!/usr/bin/env python3
r"""V6′ の固体帯 (スロット前壁の外側 $t$) を作る。

plan [`boundary-conjugate-heat-transfer.md`](../../plans/accepted/boundary-conjugate-heat-transfer.md) §6 V6′。
`case/48.flat_plate_cooled_m4/gen_solid_strip.py` と**同じ帯トポロジ**:

- **界面** = `outer_edges` (流体の `slot_front` と 1 対 1。**壁ダンプから $y$ を読む** —
  角節点 (0,0)/(0,−D) がどの physID に属するかを変換器の規則から推測しない)
- **背面** = `hole1` の **Robin** ($h$ 大で Dirichlet を近似)。`fem2d` は Dirichlet を持たない
- **両端 (上端 $y$=0 と下端 $y$=−D)** = **未登録 = 自然境界 = 断熱**

固体は $x\in[-t,0]$ に置く (界面が $x$=0、背面が $x$=−t)。流体のスロットは $x\in[0,W]$ なので重ならない。

## 設計値 (plan §6 V6′)

    t = 1.0 mm,  k_s = 0.05 W/mK  (R_s = 0.0200 m²K/W),  背面 T_c = 300 K (冷却剤),  h = 1e8
    向かい側 slot_back = 固定 500 K
    → 漸近解 q* = (500-300)/(1/h + t/k_s + W/k_f) = 4709.0 W/m²
      T_w1* = 300 + q*(1/h + t/k_s) = 394.18 K   (固体降下 94.18 K、0.5 % ゲート = 0.471 K)

**向きに注意** (2026-09-26 実測): 当初は背面 800 K でガスより熱い固体にしたが、ソルバの安全停止が
step 0 で拒否した。判定は固体温度を **[min(T_c)-20, 流体の最大全温+20]** に制限しており
(`conjugateWall.cpp`:793-795)、**「冷却剤が冷側・ガスが熱側」を前提にしている**。
背面を加熱してガスより熱い固体にする構成は、いまのソルバでは通らない。

**$k_s$ は 0.217 ではなく 0.05** — 0.217 だと $R_s/R_f$=0.205 で固体降下が 68 K しかなく、結合の感度が薄い。

使い方:
  python3 case/58.conjugate_slot/gen_solid_strip.py \
      --wall case/58.conjugate_slot/run_0001_dry/res_slot_front_5_1.h5 --nl 16
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True, help="流体の slot_front 壁ダンプ (界面 y を読む)")
    ap.add_argument("--t", type=float, default=1.0e-3, help="固体の厚さ [m]")
    ap.add_argument("--k-solid", type=float, default=0.05, help="[W/mK] (0.217 ではない — §6 V6′)")
    ap.add_argument("--Tc", type=float, default=300.0,
                    help="背面 Robin の温度 [K] = 冷却剤。**ガスより冷たくすること** — ソルバの安全停止は\n"
                         "固体温度を [min(T_c)-20, 流体の最大全温+20] に制限する (conjugateWall.cpp:793-795)")
    ap.add_argument("--h-back", type=float, default=1.0e8,
                    help="背面 Robin の h [W/m2K]。Dirichlet の近似 (大きいほど厳密)")
    ap.add_argument("--nl", type=int, default=16, help="厚さ方向の層数")
    ap.add_argument("--T-init", type=float, default=350.0, help="固体の初期温度 [K]")
    ap.add_argument("--Tw2", type=float, default=500.0,
                    help="向かい側 (slot_back) の固定壁温 [K]。**表示する漸近解にだけ使う** (固体 h5 には入らない)")
    ap.add_argument("--kf", type=float, default=0.0445, help="流体の熱伝導率 [W/mK] (同上)")
    ap.add_argument("--W", type=float, default=1.0e-3, help="スロット幅 [m] (同上)")
    ap.add_argument("--out", default=None, help="出力の基底名 (既定 mesh/solid_front_nl<nl>)")
    a = ap.parse_args()

    with h5py.File(a.wall, "r") as f:
        c = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
    if not np.allclose(c[:, 0], c[0, 0]) or not np.allclose(c[:, 2], 0.0):
        raise SystemExit("壁ダンプが x 一定・z=0 の平面でない (slot_front を渡しているか?)")
    y = np.sort(np.unique(np.round(c[:, 1], 12)))[::-1]      # 開口 (y=0) → 底 (y=-D)
    if len(y) != len(c):
        raise SystemExit(f"壁節点 {len(c)} に対し y の相異なる値が {len(y)} — 平面でない")
    ny, nl = len(y), a.nl

    # 固体は x ∈ [-t, 0]。界面は x=0 (列 j=nl)、背面は x=-t (列 j=0)。
    xs = np.linspace(-a.t, 0.0, nl + 1)
    nodes = np.array([[xx, yy] for xx in xs for yy in y], float)
    nid = lambda j, i: j * ny + i                            # noqa: E731  (j=厚さ, i=深さ)

    tris = []
    for j in range(nl):
        for i in range(ny - 1):
            tris.append([nid(j, i), nid(j + 1, i), nid(j + 1, i + 1)])
            tris.append([nid(j, i), nid(j + 1, i + 1), nid(j, i + 1)])
    tris = np.array(tris, int)

    outer = np.array([[nid(nl, i), nid(nl, i + 1)] for i in range(ny - 1)], int)   # x=0  界面
    back = np.array([[nid(0, i), nid(0, i + 1)] for i in range(ny - 1)], int)      # x=-t 背面

    out = Path(a.out) if a.out else Path(a.wall).resolve().parents[1] / "mesh" / f"solid_front_nl{nl}"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out.with_suffix(".npz"), nodes=nodes, tris=tris, outer_edges=outer, hole1=back)

    Rh, Rs = 1.0 / a.h_back, a.t / a.k_solid
    spec = {"_note": (f"case/58 V6′ スロット前壁の固体帯 (厚さ {a.t} m, {nl} 層)。"
                      f"背面は h={a.h_back:g} の Robin で Dirichlet {a.Tc} K を近似 "
                      f"(1/h / (t/k_s) = {Rh/Rs:.3g})。両端は未登録 = 断熱。"),
            "mesh_npz": str(out.with_suffix(".npz")),
            "k_solid": a.k_solid,
            "holes": [{"h": a.h_back, "T_c": a.Tc}],
            "T_init": a.T_init}
    out.with_suffix(".json").write_text(json.dumps(spec, ensure_ascii=False, indent=2))

    kf, W, Tw2 = a.kf, a.W, a.Tw2
    q = (Tw2 - a.Tc) / (Rh + Rs + W / kf)   # 固体に入る向きを正 (ガス側 T_w2 が熱側)
    print(f"[gen_solid_strip] {out.with_suffix('.npz')}")
    print(f"  節点 {len(nodes)} / 三角形 {len(tris)} / 界面辺 {len(outer)} / 背面辺 {len(back)}")
    print(f"  厚さ {a.t*1e3:.3f} mm を {nl} 層、深さ方向は壁節点と 1 対 1 ({ny} 点, y 0 .. {y[-1]:.4f} m)")
    print(f"  t/k_s = {Rs:.6e},  1/h = {Rh:.3e} ({Rh/Rs:.3g} 倍)")
    print(f"  漸近解 (k_f={kf}, W={W*1e3:.1f} mm, T_w2={Tw2:.0f} K):")
    print(f"    q* = {q:.1f} W/m2 (固体に入る向きを正),  T_w1* = {a.Tc + q*(Rh+Rs):.2f} K,"
          f"  固体の温度上昇 {q*(Rh+Rs):.2f} K")
    print(f"    0.5 % ゲート = {0.005*q*(Rh+Rs):.3f} K,  q の 0.5 % = {0.005*q:.2f} W/m2")


if __name__ == "__main__":
    main()
