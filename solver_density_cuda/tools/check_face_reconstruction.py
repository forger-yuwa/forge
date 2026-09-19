#!/usr/bin/env python3
"""再構成した面状態が物理的か (P>0, rho>0) を、**実行時と同じ座標**で検査する。

SU2 は MUSCL 再構成の後に `P<0 || rho<0 || 音速^2<0` を検査し、該当エッジを 20 反復 1 次に落とす
(`CEulerSolver.cpp` の `bad_recon` / `CFlowVariable.hpp` の `UpdateNonPhysicalEdgeCounter`)。
**forge にはこの検査が無い**ので、負の面状態がそのまま流束へ入る。本ツールはその有無を後処理で数える。

再構成は実装と同じ式:

    phi_face = phi_node + psi * (grad . d)          (`convectiveFlux_common_d.cuh` interp_MUSCL_2nd)
    d        = ±0.5 * (x_neighbour - x_node)        (node は常にエッジ中点。`convectiveFlux_d.cu` の g_reconEdgeMid)

**座標に注意 (2026-09-19 の事故)**: `CELLS/centCoords` は**実行時にノード座標へ置換される**
(`main.cpp` の `nodeValueAtNode` -> `mesh.cpp`。`forge_run.log` に `centCoords <- node coords ... max centroid shift`
が出る)。メッシュ HDF5 の `centCoords` をそのまま使うと**別の点で再構成したことになり結論が変わる**
(case/46 run_0208 で `uL/(0.5 u_nb)` が 0.70 と 1.00、負圧の有無まで反転した)。
本ツールは `res_*.h5` の `MESH/COORD` (= ノード座標) を使う。

**この検査が見ないもの**: SLAU の質量流束には圧力差の散逸項
`-chi/c_diss * (P_R - P_L)` が入る (`convectiveFlux_slau_d.inc.cuh`)。
**両側の速度が 0 でも圧力差だけで質量が動く**ので、「1 次なら壁 CV は質量を失えない」は成り立たない。
面流束そのものを論じたいなら SLAU の式を丸ごと評価すること。

使い方:
    python3 check_face_reconstruction.py RUN_DIR [--mesh sern.h5] [--steps 1,5,10] [--top 5]

`res_*.h5` は `output.level: 2` (勾配とリミッタが要る) で出したもの。
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_lsq_gradient import parse_plane_cells  # noqa: E402

NEED = ("P", "dPdx", "dPdy", "dPdz", "limiter_P", "ro", "drodx", "drody", "drodz", "limiter_ro")


def scan_step(res_path, i0, i1):
    """1 つの res から (非物理面数, 最小再構成 P, 最小の面 (ic0, ic1), その座標) を返す。"""
    with h5py.File(res_path, "r") as f:
        xyz = np.asarray(f["MESH/COORD"][()]).reshape(-1, 3).astype(np.float64)
        V = f["VALUE"]
        missing = [k for k in NEED if k not in V]
        if missing:
            raise SystemExit(f"{res_path}: 変数不足 {missing}  → `output: {{level: 2}}` で出し直す")
        q = {k: np.asarray(V[k][()]).ravel().astype(np.float64) for k in NEED}
    d0 = 0.5 * (xyz[i1] - xyz[i0])

    def rec(i, dv, base, gk, lk):
        return q[base][i] + q[lk][i] * (q[gk[0]][i] * dv[:, 0] + q[gk[1]][i] * dv[:, 1] + q[gk[2]][i] * dv[:, 2])

    PL = rec(i0, d0, "P", ("dPdx", "dPdy", "dPdz"), "limiter_P")
    PR = rec(i1, -d0, "P", ("dPdx", "dPdy", "dPdz"), "limiter_P")
    RL = rec(i0, d0, "ro", ("drodx", "drody", "drodz"), "limiter_ro")
    RR = rec(i1, -d0, "ro", ("drodx", "drody", "drodz"), "limiter_ro")
    bad = (PL <= 0) | (PR <= 0) | (RL <= 0) | (RR <= 0)
    pmin = np.minimum(PL, PR)
    k = int(np.argmin(pmin))
    return int(bad.sum()), float(pmin[k]), (int(i0[k]), int(i1[k])), xyz[i0[k]]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="再構成した面状態の物理性を検査する (SU2 の bad_recon 相当)")
    ap.add_argument("run_dir")
    ap.add_argument("--mesh", default=None, help="メッシュ h5 (既定: run_dir 内の res_ でない唯一の h5)")
    ap.add_argument("--steps", default=None, help="カンマ区切り (既定: 全 res_*.h5)")
    a = ap.parse_args(argv)
    rd = a.run_dir
    mesh = a.mesh or next(p for p in sorted(glob.glob(f"{rd}/*.h5")) if "/res_" not in p)
    with h5py.File(mesh, "r") as f:
        st = np.asarray(f["PLANES/STRUCT"][()])
        nfaces = np.asarray(f["PLANES/surfArea"][()]).shape[0]
    ic0, ic1 = parse_plane_cells(st, nfaces)
    sel = ic1 >= 0
    i0, i1 = ic0[sel], ic1[sel]
    if a.steps:
        steps = [int(x) for x in a.steps.split(",")]
    else:
        steps = sorted(int(os.path.basename(p)[4:-3]) for p in glob.glob(f"{rd}/res_[0-9]*.h5"))
    print(f"{rd}: 内部面 {int(sel.sum())} / 全面 {nfaces}  (座標はノード座標 = 実行時と同じ)")
    print(f"  {'step':>6s} {'非物理な面':>10s} {'再構成 P の最小 [Pa]':>20s}   最小の面と座標")
    worst_any = 0
    for s in steps:
        p = f"{rd}/res_{s}.h5"
        if not os.path.exists(p):
            continue
        nbad, pmin, (ja, jb), xy = scan_step(p, i0, i1)
        worst_any = max(worst_any, nbad)
        print(f"  {s:6d} {nbad:10d} {pmin:20.2f}   {ja}->{jb} at ({xy[0]:.5f}, {xy[1]:.5f}, {xy[2]:.5f})")
    print()
    print("VERDICT: " + ("OK (非物理な再構成なし)" if worst_any == 0
                         else f"NON-PHYSICAL ({worst_any} 面まで。SU2 ならこのエッジを 1 次へ落とす)"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
