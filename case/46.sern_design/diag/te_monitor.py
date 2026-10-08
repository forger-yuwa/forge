#!/usr/bin/env python3
"""後縁まわりの温度の監視 (plan convection-zero-thickness-edge-reconstruction §6.0・§6.2 #2)。

全場のスナップショット res_<n>.h5 1 枚から、事前に固定した物理領域の最低温度・低温域・床近傍の節点数を出し、
CSV に 1 行追記する (読むだけ。スナップショットの削除は呼び出し側)。

領域 (流れを参照せず固定、§6.2 の定義。H = prepare_info の H_m):
  R_TE: x ∈ [L_cowl − 0.05H, L_cowl + 0.3H]、|y − y_te| ≤ 0.2H、全 z
  R_SE: 露出したカウル側端の線分 (P0–P1) から距離 0.2H 以内
  RET : 変形の復帰区間 x ∈ [L_cowl, L_cowl + L_b]、全 y・全 z (L_b は B の te_wake_blend_H。A′ でも同じ領域を見る)
  ALL : 全域
各領域で: T の最小とその節点・座標、T < 150 K の節点数、T < 200 K の節点数と体積 (CELLS/volume の和)。
全域で: T ≤ 床 + 1 K (既定 51 K) の節点数 (§6.0: 床補正のカウンタが無いときの代用)。

使い方: te_monitor.py <run_dir> <res_n.h5> [--csv TE_MONITOR.csv] [--lb-h 1.0]
"""
import argparse
import csv
import json
import os
import re
import sys

import h5py
import numpy as np

# 露出したカウル側端 (壁終端の幾何点検 L1、g3 = run_1055 の格子): x 0.08 → 0.12、z = W/2
SE_P0 = np.array([0.08, -0.00699909, 0.1])
SE_P1 = np.array([0.12, -0.0104986, 0.1])
Y_TE = -0.0104986          # 後縁の節点の y (同上、L2)


def seg_dist(p, a, b):
    ab = b - a
    t = np.clip(((p - a) @ ab) / (ab @ ab), 0.0, 1.0)
    q = a + t[:, None] * ab
    return np.linalg.norm(p - q, axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("res")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--lb-h", type=float, default=1.0)
    ap.add_argument("--t-floor", type=float, default=50.0)
    a = ap.parse_args()

    info = json.load(open(os.path.join(a.run_dir, "prepare_info.json")))
    H = float(info["H_m"])
    L_cowl = float(info["design"].get("L_cowl", 1.2)) * H if "L_cowl" in info.get("design", {}) else 1.2 * H
    mesh = os.path.join(a.run_dir, "sern.h5")
    with h5py.File(mesh, "r") as f:
        xyz = f["MESH/COORD"][...].astype(np.float64).reshape(-1, 3)
        vol = f["CELLS/volume"][...].astype(np.float64)
    with h5py.File(a.res, "r") as f:
        T = f["VALUE/T"][...].astype(np.float64)
    if T.shape[0] != xyz.shape[0]:
        sys.exit(f"節点数が合わない: T {T.shape[0]} / 格子 {xyz.shape[0]}")
    m = re.search(r"res_(\d+)\.h5$", a.res)
    step = int(m.group(1)) if m else -1

    x, y = xyz[:, 0], xyz[:, 1]
    regions = {
        "R_TE": (x >= L_cowl - 0.05 * H) & (x <= L_cowl + 0.3 * H) & (np.abs(y - Y_TE) <= 0.2 * H),
        "R_SE": seg_dist(xyz, SE_P0, SE_P1) <= 0.2 * H,
        "RET": (x >= L_cowl) & (x <= L_cowl + a.lb_h * H),
        "ALL": np.ones_like(x, dtype=bool),
    }
    row = {"step": step, "nonfinite_T": int(np.count_nonzero(~np.isfinite(T))),
           "n_floor_plus1": int(np.count_nonzero(T <= a.t_floor + 1.0))}
    for k, msk in regions.items():
        idx = np.flatnonzero(msk & np.isfinite(T))
        if idx.size == 0:
            row.update({f"{k}_n": 0})
            continue
        j = idx[np.argmin(T[idx])]
        row.update({f"{k}_n": int(idx.size), f"{k}_Tmin": float(T[j]), f"{k}_Tmin_node": int(j),
                    f"{k}_Tmin_x": float(xyz[j, 0]), f"{k}_Tmin_y": float(xyz[j, 1]), f"{k}_Tmin_z": float(xyz[j, 2]),
                    f"{k}_n_lt150": int(np.count_nonzero(T[idx] < 150.0)),
                    f"{k}_n_lt200": int(np.count_nonzero(T[idx] < 200.0)),
                    f"{k}_vol_lt200": float(vol[idx][T[idx] < 200.0].sum())})
    out = a.csv or os.path.join(a.run_dir, "TE_MONITOR.csv")
    new = not os.path.exists(out)
    with open(out, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)
    print(" ".join(f"{k}={v}" for k, v in row.items() if k in ("step", "n_floor_plus1") or k.endswith("Tmin") or k.endswith("lt150")))


if __name__ == "__main__":
    main()
