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
    (e) |Q_sol - Q_f| / A_i / q*              <= 0.1 %   ← **連成の保存性**
        Q_sol = (K_s u - b_s)_iface は**固体の物理作用素**が界面節点で受け持つ熱 [W/m]、
        Q_f = 壁ダンプの `iface_Qf_eff` (流体が渡す積分済み荷重)、A_i = 固体側の集中辺長。
        **符号を保って**比べる。旧版は固体ダンプの `q_iface` を使っていたが、これは
        `conjugateWall.cpp` で $Q_f$ をコピーしただけの量なので、渡した荷重を渡した荷重と
        比べていた (2026-09-26 codex diagnose、plan §5.1 #100/#101)。
    (f) |固体 T(x=0) - 流体 Ts|                <= 1e-3 K ← **連成の温度連続**

使い方:
  python3 case/58.conjugate_slot/eval_v6p.py <run_dir> [--step 100000]

副産物 (plan §5.1 #101 ②③):
  <run>/v6p_band.json        決めた帯 (y 上端・下端)。G-if の帯判定
                             (`check_cht_interface.py --band-y`) にこの値を**そのまま**渡す
  <run>/v6p_band_series.csv  節点ログ (`conjugate.node_log: 1`) があれば、帯内の毎更新の
                             T_w-300 と Q_f/A を帯平均 (`Tw_m300_mean`,`q_mean`) と節点ごと
                             (`Tw_m300_<i>`,`q_<i>`) で書く → `check_quasisteady.py --series-csv`
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "solver_density_cuda", "tools"))
from solid_fem2d import Fem2DOperator  # noqa: E402

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
    ys = np.unique(np.round(sc[:, 1], 9))[::-1]
    gx = np.empty(len(ys)); tif = np.empty(len(ys))
    for i, yr in enumerate(ys):
        s = np.abs(sc[:, 1] - yr) < 1e-9
        xs, ts = sc[s, 0], Ts_sol[s]
        o = np.argsort(xs)
        gx[i] = abs(np.polyfit(xs[o], ts[o], 1)[0])
        tif[i] = ts[o][-1]                     # x=0 = 界面
    gy = np.abs(np.gradient(tif, ys))
    rat_i = np.interp(rows[::-1], ys[::-1], (gy / np.maximum(gx, 1e-30))[::-1])[::-1]
    tif_i = np.interp(rows[::-1], ys[::-1], tif[::-1])[::-1]

    # ---- (e) 用: 固体の物理作用素が界面で受け持つ熱 Q_sol = (K_s u - b_s)_iface [W/m]
    # 固体は run が実際に読んだ `solid.h5` から組む (固体ダンプの節点順は solid.h5 と同一。下で検査)。
    with h5py.File(f"{run}/solid.h5", "r") as g:
        sxy = np.asarray(g["MESH/COORD"][:], float)
        ifn = np.asarray(g["IFACE/NODES"][:], int)
        robin = [(int(e[0]), int(e[1]), float(h), float(t)) for e, h, t in
                 zip(g["ROBIN/EDGES"][:], g["ROBIN/H"][:], g["ROBIN/TC"][:])]
        kT, kV = np.asarray(g["SOLID/K_T"][:], float), np.asarray(g["SOLID/K_V"][:], float)
        op = Fem2DOperator(sxy, np.asarray(g["MESH/TRIS"][:], int), ifn,
                           np.asarray(g["IFACE/EDGES"][:], int), robin,
                           float(kV[0]) if len(kV) == 1 else (kT, kV))
    if np.abs(sxy - sc[:, :2]).max() > 1e-12:
        raise SystemExit("固体ダンプの節点順が solid.h5 と一致しない (Q_sol を組めない)")
    K, bs = op.assemble_full(Ts_sol)
    Qsol_if = (K @ Ts_sol - bs)[ifn]          # [W/m]、界面節点順
    A_if = op.area                            # 集中辺長 [m]
    y_if = sxy[ifn, 1]

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
            Qf = np.asarray(f["VALUE/iface_Qf_eff"][:], float) if "iface_Qf_eff" in f["VALUE"] else None
            return cc[:, 1], np.asarray(f["VALUE/Ts"][:], float), \
                   np.abs(np.asarray(f["VALUE/iface_q_eff"][:], float)), Qf
    y1, Ts1, q1, Qf1 = wall("res_slot_front_5")
    y2, _, q2, _ = wall("res_slot_back_6")
    i1 = np.array([np.argmin(np.abs(y1 - b)) for b in band])
    i2 = np.array([np.argmin(np.abs(y2 - b)) for b in band])
    ib = np.array([np.argmin(np.abs(rows - b)) for b in band])
    T1, Q1, Q2 = Ts1[i1], q1[i1], q2[i2]
    # (e): 界面節点 ↔ 壁節点は座標で 1 対 1 (ソルバも 1e-7 m 一致を要求している)
    if Qf1 is None:
        raise SystemExit("壁ダンプに iface_Qf_eff が無い (flux: q_eff の run でない)")
    jw = np.array([np.argmin(np.abs(y1 - yy)) for yy in y_if])
    if np.abs(y1[jw] - y_if).max() > 1e-9:
        raise SystemExit("界面節点と壁節点が座標で対応しない")
    jb = np.array([np.argmin(np.abs(y_if - b)) for b in band])
    Econs = np.abs(Qsol_if[jb] - Qf1[jw][jb]) / A_if[jb]          # [W/m2]、符号を保った差

    crit = [("(a) T_w1 誤差 [% of 固体上昇]", np.abs(T1 - Tw1s) / drop * 100, a.tol),
            ("(b) q_w1 誤差 [% of q*]", np.abs(Q1 - qs) / qs * 100, a.tol),
            ("(c) q_w2 誤差 [% of q*]", np.abs(Q2 - qs) / qs * 100, a.tol),
            ("(d) G-cons |q_w1-q_w2| [% of q*]", np.abs(Q1 - Q2) / qs * 100, a.tol),
            ("(e) 連成の保存 |Q_sol-Q_f|/A [% of q*]", Econs / qs * 100, a.tol_cons),
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

    # ---- 副産物: 帯 (G-if の帯判定に渡す) と帯内の毎更新系列 (準定常の判定に渡す)
    ytop, ybot = float(band.max()), float(band.min())
    with open(f"{run}/v6p_band.json", "w") as f:
        json.dump({"step": step, "ytop": ytop, "ybot": ybot, "rows": int(len(band)),
                   "span_W": float(span)}, f, indent=1)
    print(f"\n帯を {run}/v6p_band.json に書いた: --band-y {ytop:.9g} {ybot:.9g}")
    fl, fn = f"{run}/conjugate_iface_log_5.csv", f"{run}/conjugate_iface_nodes_5.csv"
    if os.path.exists(fl) and os.path.exists(fn):
        nd = np.loadtxt(fn, delimiter=",", skiprows=1, ndmin=2)
        lg = np.loadtxt(fl, delimiter=",", skiprows=1, ndmin=2)
        yN, AN = nd[:, 3], nd[:, 4]
        inb = np.where((yN >= ybot - 1e-9) & (yN <= ytop + 1e-9))[0]
        ups = np.unique(lg[:, 0].astype(int))
        steps = np.empty(len(ups), int)
        Tw = np.empty((len(ups), len(yN))); Q = np.empty_like(Tw)
        for k, u in enumerate(ups):
            blk = lg[lg[:, 0].astype(int) == u]
            idx = blk[:, 2].astype(int)
            Tw[k, idx] = blk[:, 6]; Q[k, idx] = blk[:, 4] / AN[idx]
            steps[k] = int(blk[0, 1])
        cols = ["step", "Tw_m300_mean", "q_mean"] + [f"Tw_m300_{i}" for i in inb] + [f"q_{i}" for i in inb]
        data = np.column_stack([steps, Tw[:, inb].mean(1) - TC, Q[:, inb].mean(1),
                                Tw[:, inb] - TC, Q[:, inb]])
        np.savetxt(f"{run}/v6p_band_series.csv", data, delimiter=",", header=",".join(cols),
                   comments="", fmt=["%d"] + ["%.10e"] * (data.shape[1] - 1))
        print(f"帯内の毎更新系列を {run}/v6p_band_series.csv に書いた ({len(ups)} 更新 × {len(inb)} 節点)")
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
