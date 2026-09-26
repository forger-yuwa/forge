#!/usr/bin/env python3
"""衝突ガスの**全温**を経路に沿って監査する (D1 不合格の性格づけ用)。

`diagnostician` (2026-09-26) の第 1 仮説: D1 の不合格は h でも格子でもなく
**前向き壁に届くガスが冷たい** (駆動温度差 T_t − T_w が実測に必要な値の 1/5〜1/30)。
それを確かめるには q でなく T_t を見る必要がある。

AGENTS.md の「出力と後処理の原則」どおり **T_t は `VALUE/h0` から** 作る
(`solver_density_cuda/tools/total_quantities.py` の `total_state`)。
スクリプト側で T + u²/2c_p を組まない。

    python3 tools/tt_audit.py --run run_0052_ab3d_open_long [--res res_60000.h5]
"""
import argparse, sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
ROOT = CASE.parents[1]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
from total_quantities import total_state                      # noqa: E402

WALL_TC = [(92, 0.25), (91, 0.51), (90, 0.76), (89, 1.52), (88, 2.54), (87, 3.81)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--res", default=None)
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--zband", type=float, default=6.0e-3)
    ap.add_argument("--off", type=float, default=0.3e-3, help="壁から手前へのオフセット [m]")
    ap.add_argument("--tw", type=float, default=300.0)
    a = ap.parse_args()
    rd = CASE / a.run

    st = total_state(str(rd), str(rd / a.res) if a.res else None)
    T0 = np.asarray(st["T0"], dtype=float)
    with h5py.File(rd / "mesh.h5") as m:
        c = np.asarray(m["MESH/COORD"], dtype=float).reshape(-1, 3)
    print(f"[{a.run}] T0 の作り方: {st['method']}  (k 込み: {st['includes_k']})  節点 {len(T0)}")
    xd = 0.5 * a.w

    def near(mask, label):
        if mask.sum() == 0:
            print(f"  {label}: 該当節点なし"); return
        print(f"  {label}: 節点 {mask.sum()}  T0 中央 {np.median(T0[mask]):7.1f} K  "
              f"最大 {T0[mask].max():7.1f} K")

    # (1) 入口 BL の高さ別 T0 (x = 入口から少し下流)
    x0 = c[:, 0].min() + 5e-3
    col = np.abs(c[:, 0] - x0) < 2e-3
    print("\n  == 入口 BL (x = 入口 +5 mm) の高さ別 T0 ==")
    for y in (0.1e-3, 0.3e-3, 1.0e-3, 1.8e-3, 5.5e-3, 2.0e-2):
        m = col & (np.abs(c[:, 1] - y) < max(0.15 * y, 3e-5))
        if m.sum():
            print(f"    y = {y*1e3:6.2f} mm : T0 = {np.median(T0[m]):7.1f} K  ({m.sum()} 節点)")

    # (2) 前向き壁の手前 --off のガス (熱電対の深さ)
    print(f"\n  == 前向き壁 (x = W/2) の {a.off*1e3:.1f} mm 手前、z<= {a.zband*1e3:.0f} mm ==")
    print("    TC  深さ[cm]    T0 [K]   T0-Tw [K]")
    for tc, d in WALL_TC:
        m = ((np.abs(c[:, 0] - (xd - a.off)) < 0.5 * a.off)
             & (np.abs(c[:, 1] + d * 1e-2) < 4e-4) & (c[:, 2] <= a.zband))
        if m.sum() == 0:
            print(f"    {tc:3d} {d:7.2f}   (節点なし)"); continue
        t0 = float(np.median(T0[m]))
        print(f"    {tc:3d} {d:7.2f}  {t0:8.1f}  {t0 - a.tw:9.1f}   ({m.sum()} 節点)")

    # (3) 交差部 (|x| < W/2+r) の深さ別 T0
    print(f"\n  == 交差部 (|x| < W/2+r, z <= {a.zband*1e3:.0f} mm) の深さ別 T0 ==")
    cross = (np.abs(c[:, 0]) < 0.5 * a.w + a.r) & (c[:, 2] <= a.zband)
    for d in (0.25, 0.55, 0.76, 1.5, 2.5, 4.0):
        m = cross & (np.abs(c[:, 1] + d * 1e-2) < 4e-4)
        if m.sum():
            print(f"    深さ {d:5.2f} cm : T0 中央 {np.median(T0[m]):7.1f} K  "
                  f"最大 {T0[m].max():7.1f} K  ({m.sum()} 節点)")

    # (4) 上流縦すきま内 (トレンチ) の深さ別 T0
    print(f"\n  == 上流縦すきま (x < -W/2-r, y < 0, z <= W/2) の深さ別 T0 ==")
    tr = (c[:, 0] < -0.5 * a.w - a.r) & (c[:, 1] < 0) & (c[:, 2] <= 0.5 * a.w)
    for d in (0.1, 0.3, 1.0, 2.5, 5.0):
        m = tr & (np.abs(c[:, 1] + d * 1e-2) < 5e-4)
        if m.sum():
            print(f"    深さ {d:5.2f} cm : T0 中央 {np.median(T0[m]):7.1f} K  "
                  f"最大 {T0[m].max():7.1f} K  ({m.sum()} 節点)")


if __name__ == "__main__":
    main()
