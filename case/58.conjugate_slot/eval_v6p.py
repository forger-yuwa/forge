#!/usr/bin/env python3
r"""V6′ の帯を決めて合否を当てる (plan boundary-conjugate-heat-transfer §6 V6′)。

**帯を決めてから 1 回だけゲートを当てる**のが登録された手順なので、このスクリプトは
① 帯の決定 → ② 合否 の順に出す。**閾値と許容は引数の既定値に固定してあり、結果を見て動かさない**。

## 帯の決め方 (2026-09-26 に 1/4 へ締めた登録値)

深さ行ごとに

    ① Pe_j = max_x|u| W / alpha_j          <= 2.5e-3
    ② T(x) の線形フィット残差 / (T_w1-T_w2) <= 2.5e-4
    ③ 固体帯内の |dT/dy| / |dT/dx|          <= 2.5e-3

を満たす**最長連続区間**、**底側 2 行を除外**。**帯が 5W 未満なら FAIL**。

旧登録 (1e-2 / 1e-3 / 1e-2) は**ゲート 0.5 % に対して 12 倍緩かった**: 2 次のずれ c·x(W-x) の
最小二乗線形残差は cW²/6、両壁の勾配差は 2cW なので **勾配差/ΔT = 12 × 残差比**。

## 合否 (すべて帯内 max。帯平均も併記)

    (a) |T_w1 - T_w1*| / (固体の温度上昇)      <= 0.5 %
    (b) |q_w1 - q*| / q*                      <= 0.5 %
    (c) |q_w2 - q*| / q*                      <= 0.5 %
    (d) G-cons |q_w1 - q_w2| / q*             <= 0.5 %   ← 流体柱の 1 次元性 (固体に依らない検査ではない)
    (e) |q_iface/L - q_eff| / q*              <= 0.1 %   ← **連成の保存性**
    (f) |固体 T(x=0) - 流体 Ts|                <= 1e-3 K ← **連成の温度連続**

使い方:
  python3 case/58.conjugate_slot/eval_v6p.py <run_dir> [--step 100000]
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import h5py
import numpy as np

# --- 問題の定義 (固体 h5 と config と一致させること) ---
TC, TW2 = 300.0, 500.0          # 冷却剤 (背面 Robin) / 後壁 (固定等温・熱側)
KF, KS, T_S, W, H_BACK = 0.0445, 0.05, 1.0e-3, 1.0e-3, 1.0e8
CP = 1004.5


def analytic():
    Rh, Rs, Rf = 1.0 / H_BACK, T_S / KS, W / KF
    q = (TW2 - TC) / (Rh + Rs + Rf)
    return q, TC + q * (Rh + Rs), q * (Rh + Rs)      # q*, T_w1*, 固体の温度上昇


def last_step(run):
    st = [int(re.search(r"res_(\d+)\.h5$", p).group(1))
          for p in glob.glob(os.path.join(run, "res_[0-9]*.h5"))]
    return max(st)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--pe", type=float, default=2.5e-3, help="帯の Pe 上限 (登録値。動かさない)")
    ap.add_argument("--lin", type=float, default=2.5e-4, help="帯の線形残差上限 (同)")
    ap.add_argument("--rat", type=float, default=2.5e-3, help="帯の勾配比上限 (同)")
    ap.add_argument("--tol", type=float, default=0.5, help="(a)-(d) の許容 [%%]")
    ap.add_argument("--tol-cons", type=float, default=0.1, help="(e) の許容 [%%]")
    ap.add_argument("--tol-tc", type=float, default=1.0e-3, help="(f) の許容 [K]")
    a = ap.parse_args()
    run = a.run.rstrip("/")
    step = a.step if a.step else last_step(run)
    qs, Tw1s, drop = analytic()

    # ---- 流体場: スロット内の深さ行ごとに Pe と線形残差
    with h5py.File(f"{run}/res_{step}.h5", "r") as f:
        c = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
        V = f["VALUE"]
        T = np.asarray(V["T"][:], float); ro = np.asarray(V["ro"][:], float)
        U = np.hypot(np.asarray(V["Ux"][:], float), np.asarray(V["Uy"][:], float))
    m = (c[:, 0] >= -1e-12) & (c[:, 0] <= W + 1e-12) & (c[:, 1] < 1e-12)
    x, y, Tf, rof, Uf = c[m, 0], c[m, 1], T[m], ro[m], U[m]
    rows = np.unique(np.round(y, 9))[::-1]

    # ---- 固体: 界面列の縦勾配 / 厚さ方向勾配 と、界面温度・界面熱流束
    with h5py.File(f"{run}/res_solid_5_{step}.h5", "r") as g:
        sc = np.asarray(g["MESH/COORD"][:], float).reshape(-1, 3)
        Ts_sol = np.asarray(g["VALUE/T"][:], float)
        q_if = np.asarray(g["VALUE/q_iface"][:], float)
    ys = np.unique(np.round(sc[:, 1], 9))[::-1]
    gx = np.empty(len(ys)); tif = np.empty(len(ys)); qif = np.empty(len(ys))
    for i, yr in enumerate(ys):
        s = np.abs(sc[:, 1] - yr) < 1e-9
        xs, ts, qs_ = sc[s, 0], Ts_sol[s], q_if[s]
        o = np.argsort(xs)
        gx[i] = abs(np.polyfit(xs[o], ts[o], 1)[0])
        tif[i] = ts[o][-1]                     # x=0 = 界面
        qif[i] = qs_[o][-1]
    gy = np.abs(np.gradient(tif, ys))
    rat_i = np.interp(rows[::-1], ys[::-1], (gy / np.maximum(gx, 1e-30))[::-1])[::-1]
    tif_i = np.interp(rows[::-1], ys[::-1], tif[::-1])[::-1]
    qif_i = np.interp(rows[::-1], ys[::-1], qif[::-1])[::-1]

    Pe = np.empty(len(rows)); lin = np.empty(len(rows))
    for i, yr in enumerate(rows):
        s = np.abs(y - yr) < 5e-10
        xs, Ts_, ros, Us = x[s], Tf[s], rof[s], Uf[s]
        o = np.argsort(xs); xs, Ts_, ros, Us = xs[o], Ts_[o], ros[o], Us[o]
        Pe[i] = Us.max() * W / (KF / (ros.mean() * CP))
        A = np.polyfit(xs, Ts_, 1)
        lin[i] = np.max(np.abs(Ts_ - np.polyval(A, xs))) / max(abs(Ts_[0] - Ts_[-1]), 1e-9)

    ok = (Pe <= a.pe) & (lin <= a.lin) & (rat_i <= a.rat)
    idx = np.where(ok)[0]
    if len(idx) == 0:
        print("VERDICT: FAIL (帯の条件を満たす行が無い)"); return 1
    segs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    seg = max(segs, key=len)[:-2]                      # 底側 2 行を除外
    span = (rows[seg[0]] - rows[seg[-1]]) / W
    band = rows[seg]
    print(f"=== V6′ 評価: {run}  step {step} ===")
    print(f"解析解: q* = {qs:.1f} W/m2,  T_w1* = {Tw1s:.3f} K,  固体の温度上昇 {drop:.2f} K "
          f"(0.5 % = {0.005*drop:.3f} K)")
    print(f"帯 (Pe<={a.pe:.1e} / 線形残差<={a.lin:.1e} / 勾配比<={a.rat:.1e}、底 2 行除外): "
          f"{len(band)} 行  y {band.max()*1e3:.3f} .. {band.min()*1e3:.3f} mm = {span:.1f} W"
          f"  -> {'OK' if span >= 5 else '**FAIL (5W 未満)**'}")

    # ---- 壁ダンプ
    def wall(stem):
        with h5py.File(f"{run}/{stem}_{step}.h5", "r") as f:
            cc = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
            return cc[:, 1], np.asarray(f["VALUE/Ts"][:], float), \
                   np.abs(np.asarray(f["VALUE/iface_q_eff"][:], float))
    y1, Ts1, q1 = wall("res_slot_front_5")
    y2, _, q2 = wall("res_slot_back_6")
    i1 = np.array([np.argmin(np.abs(y1 - b)) for b in band])
    i2 = np.array([np.argmin(np.abs(y2 - b)) for b in band])
    ib = np.array([np.argmin(np.abs(rows - b)) for b in band])
    T1, Q1, Q2 = Ts1[i1], q1[i1], q2[i2]
    # 固体側の界面熱流束 [W/m2] = q_iface / 集中辺長。q_iface は節点あたり [W/m] なので
    # 隣接辺長の半和で割る。ここでは帯内で滑らかなので中心差分の辺長を使う。
    L = np.abs(np.gradient(rows))[ib]
    Qsol = np.abs(qif_i[ib]) / np.maximum(L, 1e-30)

    crit = [("(a) T_w1 誤差 [% of 固体上昇]", np.abs(T1 - Tw1s) / drop * 100, a.tol),
            ("(b) q_w1 誤差 [% of q*]", np.abs(Q1 - qs) / qs * 100, a.tol),
            ("(c) q_w2 誤差 [% of q*]", np.abs(Q2 - qs) / qs * 100, a.tol),
            ("(d) G-cons |q_w1-q_w2| [% of q*]", np.abs(Q1 - Q2) / qs * 100, a.tol),
            ("(e) 連成の保存 |q_sol-q_eff| [% of q*]", np.abs(Qsol - Q1) / qs * 100, a.tol_cons),
            ("(f) 連成の T 連続 [K]", np.abs(tif_i[ib] - T1), a.tol_tc)]
    print(f"\n{'量':<40}{'帯内 max':>11}{'帯平均':>11}{'許容':>9}  判定")
    bad = []
    for nm, v, tol in crit:
        p = v.max() <= tol
        if not p:
            bad.append(nm)
        print(f"{nm:<40}{v.max():11.4f}{v.mean():11.4f}{tol:9.3f}  {'PASS' if p else '**FAIL**'}")
    print(f"\n参考: T_w1 帯内 {T1.min():.4f}..{T1.max():.4f} K   "
          f"q_w1 {Q1.min():.1f}..{Q1.max():.1f}   q_w2 {Q2.min():.1f}..{Q2.max():.1f} W/m2")
    v = "PASS" if (not bad and span >= 5) else "FAIL"
    print(f"\nVERDICT: {v}" + (f"  (外れた量: {', '.join(bad)})" if bad else ""))
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
