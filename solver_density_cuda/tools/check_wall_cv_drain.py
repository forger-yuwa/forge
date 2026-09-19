#!/usr/bin/env python3
"""node の**速度ピン壁ノード**で、対流再構成が制御体積からどれだけ質量を抜いているかを測る。

背景 (2026-09-19, case/46 R4e): node (median-dual) の壁ノードは

- 状態が `u = 0` に固定され (`nodeWallDirichlet`、`cuda_forge/nodeWallDirichlet_d.cu`)、
  運動量の残差も 0 に射影される。**連続の式の残差 `res_ro` は固定されない**。
- 壁の境界半割面の質量流束は bvar 状態のみで作られ、no-slip なら厳密に 0
  (`convectiveFlux_boundary_d.inc.cuh`)。**つまり壁 CV に対流で質量が入る道は内部面しかない**。
- その内部面で 2 次 MUSCL は `phi_face = phi_node + psi * (grad . d)` を使う
  (`convectiveFlux_common_d.cuh`)。`phi_node` (速度) は 0 でも **grad が大きければ面に速度が乗る**。

1 次なら `phi_face = phi_node = 0` なので壁 CV は対流で質量を失えない。2 次はその性質を壊す。
流れに正対する薄い壁 (機体ベース等) では全ての内部面が外向きになり、CV が単調に空になって発散する。

本ツールはその「抜ける速さ」を **1/s** (= 正味外向き質量流束 / CV 質量) で測る。

**限界 (必ず読むこと)**: SLAU/Roe の風上化を含まない見積りである。再構成した左状態から
`rho_L * u_L . S` を素直に積むだけなので、実際の流束とは一致しない。**これはふるいであって予測器ではない。**

**偽陽性の実例** (case/46 run_0208): 本ツールが最大 (9.4e5 /s) と出したカウル後縁のノード 9179 は、
実際には密度が 0.0253 -> 0.0285 と**上がって**おり漏れていない。一方 3 位のベース壁ノード 62104 (1.5e5 /s) は
27 step で 0.00358 -> 0.000027 と **132 倍**薄くなり破綻した。**判定は必ず `ro` の時系列で確認すること**
(`res_N.h5` を複数出して当該ノードの `ro` を追う)。

また収束解では正味流束は定義上ほぼ 0 になるので、**収束解同士でないと公平な比較にならない**
(1 次収束場に 2 次を当てた直後の値と、2 次収束解の値を並べない)。

使い方:
    python3 check_wall_cv_drain.py MESH.h5 RES.h5 [--top N] [--threshold 1e3]

RES は `output.level: 2` (勾配とリミッタが要る) で出したもの。
"""
from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_lsq_gradient import parse_plane_cells  # noqa: E402

NEED = ("ro", "Ux", "Uy", "Uz", "dUxdx", "dUxdy", "dUxdz", "dUydx", "dUydy", "dUydz",
        "dUzdx", "dUzdy", "dUzdz", "drodx", "drody", "drodz",
        "limiter_Ux", "limiter_Uy", "limiter_Uz", "limiter_ro")


def drain_rates(mesh_h5, res_h5):
    """(節点 index, 抜ける速さ [1/s], 正味外向き質量流束 [kg/s], CV 質量 [kg]) を返す。"""
    with h5py.File(mesh_h5, "r") as f:
        cc = np.asarray(f["CELLS/centCoords"][()]).reshape(-1, 3).astype(np.float64)
        st = np.asarray(f["PLANES/STRUCT"][()])
        nfaces = np.asarray(f["PLANES/surfArea"][()]).shape[0]
        sv = np.asarray(f["PLANES/surfVect"][()]).reshape(-1, 3).astype(np.float64)
        vol = np.asarray(f["CELLS/volume"][()]).ravel().astype(np.float64)
    ic0, ic1 = parse_plane_cells(st, nfaces)
    with h5py.File(res_h5, "r") as f:
        V = f["VALUE"]
        missing = [k for k in NEED if k not in V]
        if missing:
            raise SystemExit(f"{res_h5}: 変数が足りない {missing}  → `output: {{level: 2}}` で出し直すこと")
        q = {k: np.asarray(V[k][()]).ravel().astype(np.float64) for k in NEED}
    # 速度 3 成分が厳密 0 のノード = Dirichlet でピンされた壁ノード
    pin = np.where((np.abs(q["Ux"]) < 1e-9) & (np.abs(q["Uy"]) < 1e-9) & (np.abs(q["Uz"]) < 1e-9))[0]
    sel = ic1 >= 0
    i0, i1, S = ic0[sel], ic1[sel], sv[sel]
    d0 = 0.5 * (cc[i1] - cc[i0])                     # node: 再構成の目標点はエッジ中点
    mdot = np.zeros(len(q["ro"]))

    def rec(i, dv, base, gk, lk):
        return q[base][i] + q[lk][i] * (q[gk[0]][i] * dv[:, 0] + q[gk[1]][i] * dv[:, 1] + q[gk[2]][i] * dv[:, 2])

    for i, dv, sgn in ((i0, d0, 1.0), (i1, -d0, -1.0)):
        u = rec(i, dv, "Ux", ("dUxdx", "dUxdy", "dUxdz"), "limiter_Ux")
        v = rec(i, dv, "Uy", ("dUydx", "dUydy", "dUydz"), "limiter_Uy")
        w = rec(i, dv, "Uz", ("dUzdx", "dUzdy", "dUzdz"), "limiter_Uz")
        r = rec(i, dv, "ro", ("drodx", "drody", "drodz"), "limiter_ro")
        np.add.at(mdot, i, sgn * r * (u * S[:, 0] + v * S[:, 1] + w * S[:, 2]))
    mass = q["ro"][pin] * vol[pin]
    return pin, mdot[pin] / np.maximum(mass, 1e-300), mdot[pin], mass


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="node 速度ピン壁ノードの対流による質量漏れを測る")
    ap.add_argument("mesh"); ap.add_argument("res")
    ap.add_argument("--top", type=int, default=8, help="漏れの大きい順に表示する節点数")
    ap.add_argument("--threshold", type=float, default=1.0e3,
                    help="この速さ [1/s] を超える節点を数える (既定 1e3 = 1 ms で空になる)")
    a = ap.parse_args(argv)
    pin, rate, mdot, mass = drain_rates(a.mesh, a.res)
    if len(pin) == 0:
        print("速度ピン壁ノードが無い (nodeWallDirichlet が 0 か、node でない)")
        return 0
    with h5py.File(a.res, "r") as f:
        xyz = np.asarray(f["MESH/COORD"][()]).reshape(-1, 3)
    n_bad = int((rate > a.threshold).sum())
    print(f"{a.res}: 速度ピン壁ノード {len(pin)}")
    print(f"  抜ける速さ/質量 [1/s]: 中央値 {np.median(rate):10.3e}  p90 {np.percentile(rate, 90):10.3e}  最大 {rate.max():10.3e}")
    print(f"  > {a.threshold:g} /s の節点: {n_bad} / {len(pin)}")
    order = pin[np.argsort(-rate)][:a.top]
    ranks = np.argsort(-rate)[:a.top]
    print(f"  {'節点':>8s} {'速さ [1/s]':>12s} {'流束 [kg/s]':>13s} {'CV 質量 [kg]':>13s}   座標 (x, y, z) [m]")
    for n, k in zip(order, ranks):
        print(f"  {n:8d} {rate[k]:12.3e} {mdot[k]:13.3e} {mass[k]:13.3e}   "
              f"({xyz[n,0]:.5f}, {xyz[n,1]:.5f}, {xyz[n,2]:.5f})")
    print()
    print("VERDICT: " + ("OK (1 ms 以内に空になる壁ノードは無い)" if n_bad == 0
                         else f"DRAINING ({n_bad} 個の壁 CV が {1e3/a.threshold:g} ms 以内に空になる勢いで漏れている)"))
    print("注意: これはふるいであって予測器ではない (偽陽性あり)。上位ノードの `ro` 時系列で必ず確認すること。")
    print("      収束解同士でないと公平に比較できない。詳細は docstring。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
