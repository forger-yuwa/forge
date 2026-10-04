#!/usr/bin/env python3
"""#14 (codex result-3 M3): SERN 2D 衝撃足の**各 run の準定常**を判定するための系列を作る。
plan convection-slau-wall-normal-chi §5.1 #14。差ノルム (v3sern_series.py) を評価する前に、
flag0 / flag1 **それぞれの場**が衝撃足付近で定常かを見る (両側が一緒に動くと差は小さくても各場は未定常)。

run ごとに書く列 (判定には run の**量そのもの**を当てる。偏差の系列に相対閾値を当てない — #10b の operand 取り違え):
  xfoot        = その run の max|dp_w/dx| の位置 [m] (絶対位置)
  pwin_mean    = 固定窓 [xf_ref − 2t, xf_ref + 5t] の p_w 平均 [Pa] (xf_ref = flag0 の末尾 dump 群の x_foot 中央値。両 run 共通)
  p_m1t / p_p1t / p_p3t = 窓内 3 点 (xf_ref − t, + t, + 3t) の p_w [Pa] (ランプ節点の線形補間)
  selfL2_pct   = 窓内 p_w の「各 dump vs 自 run 末尾平均」相対 L2 [%] (参考列。偏差なので --series-csv の判定には当てない)

usage: v3sern_foot_perrun.py RUN0 RUN1 --out-prefix PFX [--ramp-id 4] [--t-cowl 0.002] [--tail 0.4]
"""
import argparse, sys
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v3sern_series import dumps, ramp_geom

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run0"); ap.add_argument("run1")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ramp-id", type=int, default=4)
    ap.add_argument("--t-cowl", type=float, default=0.002)
    ap.add_argument("--tail", type=float, default=0.4)
    ap.add_argument("--mesh", default=None)
    a = ap.parse_args()
    mesh = a.mesh or str(Path(a.run0) / "sern.h5")
    nodes, x, _, _, _ = ramp_geom(mesh, a.ramp_id)
    t = a.t_cowl
    runs = {"flag0": a.run0, "flag1": a.run1}
    prof, xf = {}, {}
    for tag, rd in runs.items():
        d = dumps(rd); steps = sorted(s for s in d if s > 0)
        P = []
        for s in steps:
            with h5py.File(d[s]) as g:
                P.append(np.asarray(g["/VALUE/P"], np.float64)[nodes])
        P = np.array(P)
        prof[tag] = (np.array(steps), P)
        xf[tag] = np.array([x[int(np.argmax(np.abs(np.gradient(p, x))))] for p in P])
    s0, _ = prof["flag0"]; ntail = max(3, int(round(len(s0) * a.tail)))
    xf_ref = float(np.median(xf["flag0"][-ntail:]))
    win = (x >= xf_ref - 2 * t) & (x <= xf_ref + 5 * t)
    j = int(np.argmin(np.abs(x - xf_ref)))
    dx_local = float(np.median(np.diff(x[max(j - 3, 0): j + 4])))
    print(f"# ramp ノード {len(nodes)}  xf_ref {xf_ref:.6f} m  窓内ノード {int(win.sum())}  局所格子間隔 {dx_local:.3e} m")
    for tag in runs:
        steps, P = prof[tag]
        pm = P[-ntail:].mean(axis=0)
        rows = []
        for k, s in enumerate(steps):
            p = P[k]
            l2 = 100.0 * np.sqrt(np.mean((p[win] - pm[win]) ** 2)) / np.sqrt(np.mean(pm[win] ** 2))
            pts = [np.interp(xf_ref + c * t, x, p) for c in (-1, 1, 3)]
            rows.append((s, xf[tag][k], p[win].mean(), *pts, l2))
        out = f"{a.out_prefix}_{tag}.csv"
        with open(out, "w") as f:
            f.write("step,xfoot,pwin_mean,p_m1t,p_p1t,p_p3t,selfL2_pct\n")
            for r in rows:
                f.write(f"{r[0]},{r[1]:.9f},{r[2]:.6f},{r[3]:.6f},{r[4]:.6f},{r[5]:.6f},{r[6]:.6f}\n")
        tx = xf[tag][-ntail:]
        print(f"# {tag}: {len(steps)} dump, 末尾 {ntail} の x_foot span {tx.max() - tx.min():.3e} m "
              f"(= 格子間隔の {(tx.max() - tx.min()) / dx_local:.2f} 倍)、書き出し {out}")
