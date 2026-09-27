#!/usr/bin/env python3
r"""非連成の静止保持の判定 (plan `axisymmetric-graded-grid-static-gas.md` §5.1 #3、事前登録 2026-09-27)。

run 直下の `res_<step>.h5` (場) と `res_wall_hot_3_<step>.h5` / `res_wall_cj_4_<step>.h5` (両壁) から、
各スナップショットの

    Umax      : max|U| [m/s]
    Uymax     : max|U_y| [m/s]
    dPrel     : (P_max − P_min) / <P>_V   (<P>_V は双対体積 × 双対重心 r の重み平均)
    checker   : 冷却壁 (x=H) の iface_q_eff の市松振幅 [W/m²] (内部壁節点の値と左右隣接からの r 線形補間値の差の max)
    qerr_hot / qerr_cj : 全壁節点の |q_eff| と伝導基準 120.5 W/m² の差の max [% of 基準]

を系列 CSV (`static_hold_series.csv`) に書き、最終比較窓 (step 10000–20000) で

    max|U| ≤ 1e-3 m/s、全壁熱流束誤差 ≤ 0.5 %、相対圧力差 ≤ 1e-4

を判定する。残差 (`check_convergence.py`) と準定常 (`check_quasisteady.py --series-csv --tail 0.5
--drift 0.001 --osc 0.001`) は別ツールで回す (このスクリプトは準定常のコマンドを表示する)。

    python3 case/62.conjugate_disk/eval_static_hold.py <run> [--win0 10000 --win1 20000]
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import h5py
import numpy as np

Q_REF = 0.0241 * (350.0 - 325.0) / 0.005      # 120.5 W/m²


def checker(r, q):
    o = np.argsort(r); r, q = r[o], q[o]
    return max(abs(q[i] - (q[i - 1] + (q[i + 1] - q[i - 1]) * (r[i] - r[i - 1]) / (r[i + 1] - r[i - 1])))
               for i in range(1, len(r) - 1))


def wall(run, stem, st):
    f = f"{run}/{stem}_{st}.h5"
    if not os.path.exists(f):
        return None
    with h5py.File(f, "r") as w:
        return (np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)[:, 1],
                np.asarray(w["VALUE/iface_q_eff"][:], float))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--win0", type=int, default=10000)
    ap.add_argument("--win1", type=int, default=20000)
    a = ap.parse_args()
    run = a.run.rstrip("/")
    with h5py.File(f"{run}/mesh.h5", "r") as m:
        V = np.asarray(m["CELLS/volume"][:], float) * np.asarray(m["CELLS/centCoords"][:], float).reshape(-1, 3)[:, 1]
    rows = []
    for f in glob.glob(f"{run}/res_[0-9]*.h5"):
        st = int(re.search(r"res_(\d+)\.h5$", f).group(1))
        with h5py.File(f, "r") as h:
            Vv = h["VALUE"]
            Ux, Uy, P = (np.asarray(Vv[k][:], float) for k in ("Ux", "Uy", "P"))
        if not (np.isfinite(Ux).all() and np.isfinite(Uy).all() and np.isfinite(P).all()):
            raise SystemExit(f"{f}: 非有限値")
        hot, cj = wall(run, "res_wall_hot_3", st), wall(run, "res_wall_cj_4", st)
        if hot is None or cj is None:
            continue
        rows.append((st, float(np.hypot(Ux, Uy).max()), float(np.abs(Uy).max()),
                     float((P.max() - P.min()) / ((P * V).sum() / V.sum())),
                     float(checker(*cj)),
                     float(np.abs(np.abs(hot[1]) - Q_REF).max() / Q_REF * 100),
                     float(np.abs(np.abs(cj[1]) - Q_REF).max() / Q_REF * 100)))
    rows.sort()
    X = np.array(rows)
    cols = ["step", "Umax", "Uymax", "dPrel", "checker", "qerr_hot", "qerr_cj"]
    out = f"{run}/static_hold_series.csv"
    np.savetxt(out, X, delimiter=",", header=",".join(cols), comments="", fmt=["%d"] + ["%.10e"] * 6)
    w = X[(X[:, 0] >= a.win0) & (X[:, 0] <= a.win1)]
    if len(w) == 0:
        raise SystemExit(f"比較窓 {a.win0}–{a.win1} にスナップショットが無い")
    crit = [("max|U| [m/s]", w[:, 1].max(), 1e-3),
            ("全壁熱流束誤差 hot [% of 120.5]", w[:, 5].max(), 0.5),
            ("全壁熱流束誤差 cj [% of 120.5]", w[:, 6].max(), 0.5),
            ("相対圧力差 (Pmax−Pmin)/<P>", w[:, 3].max(), 1e-4)]
    print(f"=== 静止保持: {run}  ({len(X)} スナップショット、窓 {a.win0}–{a.win1} に {len(w)}) ===")
    bad = []
    for nm, v, tol in crit:
        ok = v <= tol
        bad += [] if ok else [nm]
        print(f"  {'PASS' if ok else 'FAIL'}  {nm:<34} max {v:.4e}  (許容 {tol:g})")
    print(f"  参考: max|U_y| {w[:, 2].max():.4e} m/s、市松 {w[:, 4].max():.4e} W/m²")
    print(f"\n準定常: python3 solver_density_cuda/tools/check_quasisteady.py --series-csv {out} "
          f"--series-cols Umax,Uymax,dPrel,checker,qerr_hot,qerr_cj --tail 0.5 --drift 0.001 --osc 0.001")
    print(f"\nVERDICT (閾値のみ): {'PASS' if not bad else 'FAIL'}" + (f"  ({', '.join(bad)})" if bad else ""))
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
