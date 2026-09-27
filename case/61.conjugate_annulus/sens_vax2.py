#!/usr/bin/env python3
r"""V-ax2 の格子感度 — **共通位置**の界面温度の差 (plan §6 V-ax2「感度」)。

登録: 細分化の系列 (粗 → 細) を渡し、**最後の細分化** (最後の 2 本) で、共通位置 (両 run に同じ $x$ の壁節点がある位置)
の界面温度 `Ts` の差の最大が **0.023581 K (= 0.1 % of 降下 23.581 K) 以下**。途中の組の差は参考として出す。

  固体: 半径方向 8 / 16 / 32 層 × 軸方向を同率 (**系列の組み方は保留中** — README「保留」)
  流体: 半径方向 16 / 32 (固体は 16 層 × 軸 4)

使い方 (粗 → 細の順に並べる):
  python3 case/61.conjugate_annulus/sens_vax2.py <run_coarse> [<run_mid>] <run_fine> [--step N]
各 run の最終 step (または --step) の壁ダンプ `res_wall_outer_4_<step>.h5` を読む。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import axcht
from eval_vax2 import NAME_CJ, PID_CJ, analytic


def wall_T(run: Path, step):
    st = step if step is not None else axcht.last_step(run)
    d = axcht.wall_dump(run, NAME_CJ, PID_CJ, st)
    return st, d["xyz"][:, 0], d["Ts"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="粗 → 細の順")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--tol-K", type=float, default=0.023581, help="最後の細分化の許容 [K] (登録値)")
    a = ap.parse_args()
    if len(a.runs) < 2:
        raise SystemExit("run を 2 本以上 (粗 → 細) 渡すこと")
    an = analytic()
    data = [(Path(r),) + wall_T(Path(r), a.step) for r in a.runs]
    print(f"解析 T_w* = {an['Tw']:.5f} K、許容 (最後の細分化) {a.tol_K} K")
    last = None
    for (r0, s0, x0, T0), (r1, s1, x1, T1) in zip(data[:-1], data[1:]):
        k0 = np.round(x0, 9)
        k1 = np.round(x1, 9)
        common = np.intersect1d(k0, k1)
        if len(common) == 0:
            raise SystemExit(f"{r0} と {r1} に共通位置が無い")
        i0 = np.array([np.where(k0 == c)[0][0] for c in common])
        i1 = np.array([np.where(k1 == c)[0][0] for c in common])
        d = np.abs(T1[i1] - T0[i0])
        last = d.max()
        print(f"{r0.name} (step {s0}) → {r1.name} (step {s1}): 共通 {len(common)} 位置 "
              f"(x = {', '.join(f'{c*1e3:.3f}' for c in common)} mm)  max|ΔT_w| = {d.max():.6f} K"
              f"  [解析との差 粗 {np.abs(T0[i0]-an['Tw']).max():.5f} / 細 {np.abs(T1[i1]-an['Tw']).max():.5f} K]")
    ok = last <= a.tol_K
    print(f"\n最後の細分化 max|ΔT_w| = {last:.6f} K (許容 {a.tol_K})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
