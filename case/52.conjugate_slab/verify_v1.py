#!/usr/bin/env python3
r"""V1 の合否判定: 1 次元純伝導の共役解を解析解と照合する。

問題: 静止流体層 (厚さ $H$, 熱伝導率 $k_f$ 一定) の上面を $T_\infty$ に固定し、下面を
固体 (全抵抗 $R_s = t/k_s + R_{back}$、背面 $T_b$) と連成させる。直列抵抗なので

$$R_f = \frac{H}{k_f},\qquad q = \frac{T_\infty - T_b}{R_f + R_s},\qquad T_w = T_b + q R_s
   \;=\; \frac{T_\infty/R_f + T_b/R_s}{1/R_f + 1/R_s}$$

合格ライン (plan §6 V1): $T_w$ が解析解と **0.5 % 以内** (温度上昇 $T_w-T_b$ に対する相対)、
**両側の $q$ 不一致が 0.1 % 以内**。

usage: python3 verify_v1.py run_0001_v1_slab [--H 0.01] [--kf 0.0241] [--Tinf 350] [--flux q_compact]
"""
import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda" / "tools"))
from solid_shell import ShellOperator, SolidModel          # noqa: E402
from cht_loop import read_wall_dump, latest_wall_dump      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--H", type=float, default=0.01, help="流体層の厚さ [m]")
    ap.add_argument("--kf", type=float, default=0.0241, help="流体の熱伝導率 [W/mK]")
    ap.add_argument("--Tinf", type=float, default=350.0, help="上壁温度 [K]")
    ap.add_argument("--phys-id", type=int, default=3)
    ap.add_argument("--phys-name", default="wall_bot")
    ap.add_argument("--flux", default="q_compact")
    ap.add_argument("--tol-Tw", type=float, default=0.5, help="[%] 温度上昇に対する相対許容")
    ap.add_argument("--tol-q", type=float, default=0.1, help="[%] 両側 q の不一致許容")
    a = ap.parse_args()

    run = Path(a.run_dir)
    model = SolidModel.from_json(str(run / "solid.json"))
    its = sorted(glob.glob(str(run / "it_*")))
    if not its:
        sys.exit(f"no iterations in {run}")
    last = Path(its[-1])

    coords, faces, vals = read_wall_dump(latest_wall_dump(last, a.phys_name, a.phys_id))
    op = ShellOperator(coords, faces, model)
    Tw = np.array(vals["Ts"], float)                 # 場に課された壁温 (bvar の再出力)
    q_f = np.array(vals["iface_" + a.flux], float)   # 流体側 [W/m2] 固体向き正

    Rs = float(model.R_tot(Tw).mean())
    Rf = a.H / a.kf
    q_ex = (a.Tinf - model.T_b) / (Rf + Rs)
    Tw_ex = model.T_b + q_ex * Rs

    q_s = (Tw - model.T_b) / Rs                      # 固体側の熱流束 [W/m2]
    dT = float(np.mean(Tw)) - model.T_b
    err_Tw = 100.0 * abs(float(np.mean(Tw)) - Tw_ex) / (Tw_ex - model.T_b)
    err_q = 100.0 * abs(float(np.mean(q_f)) - float(np.mean(q_s))) / abs(float(np.mean(q_s)))
    err_qex = 100.0 * abs(float(np.mean(q_f)) - q_ex) / q_ex

    print(f"case      : H={a.H} m, k_f={a.kf} W/mK -> R_f={Rf:.6f} ; R_s={Rs:.6f} m2K/W "
          f"(t/k_s + R_back), T_inf={a.Tinf} K, T_b={model.T_b} K")
    print(f"analytic  : q = {q_ex:.4f} W/m2 , T_w = {Tw_ex:.4f} K")
    print(f"forge+shell: q_fluid = {np.mean(q_f):.4f} W/m2 (spread {q_f.max()-q_f.min():.2e}) , "
          f"q_solid = {np.mean(q_s):.4f} W/m2 , T_w = {np.mean(Tw):.4f} K "
          f"(spread {Tw.max()-Tw.min():.2e} K)")
    print(f"errors    : T_w {err_Tw:.4f} % of the {dT:.3f} K rise | "
          f"two-sided q mismatch {err_q:.4f} % | q vs analytic {err_qex:.4f} %")
    print(f"source    : {last}/  (flux definition: iface_{a.flux})")

    ok = (err_Tw <= a.tol_Tw) and (err_q <= a.tol_q)
    print(f"\nVERDICT: {'PASS' if ok else 'FAIL'} "
          f"(T_w <= {a.tol_Tw} % : {'ok' if err_Tw <= a.tol_Tw else 'NG'} ; "
          f"two-sided q <= {a.tol_q} % : {'ok' if err_q <= a.tol_q else 'NG'})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
