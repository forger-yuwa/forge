#!/usr/bin/env python3
"""S2 case/40 (軸対称ベルノズル) の η_CF・ṁ・壁温 (plan gradient-scalar-lsq-unification §6 S2 表)。

  python3 s2_eval_case40.py --gg RUN [RUN_ext] --lsq RUN [RUN_ext] --out-prefix PREFIX

- η_CF は `design/forge_design/metrics/extract.py` の `thrust_metrics` (README:172 と同じ抽出。Pt 4e6、背圧 2e4、r_t 0.01、出口 physID 2)。
- ṁ = 出口面の ∫ρu 2πr dr。壁温 = 壁 (physID 3) の節点の VALUE/T。
- 複数 run を渡すと step を通算して 1 系列にする (初段 + 延長)。系列 CSV を `PREFIX_{gg,lsq}.csv` に書く
  (`check_quasisteady.py --series-csv` 用)。最終スナップショットで双子の差 (η_CF、ṁ 相対、壁温 L∞) を出す。
"""
import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "design"))
from forge_design.metrics.extract import thrust_metrics  # noqa: E402

GAMMA, PT, PA, RT = 1.4, 4.0e6, 2.0e4, 0.01


def snaps(run):
    out = []
    for p in glob.glob(os.path.join(run, "res_*.h5")):
        m = re.search(r"/res_(\d+)\.h5$", p)
        if m:
            out.append((int(m.group(1)), p))
    return sorted(out)


def mdot_and_wall(mesh, res):
    with h5py.File(mesh, "r") as nz:
        ip = nz["/BCONDS/2/iPlanes"][:]; ic = nz["/BCONDS/2/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        dA = 2.0 * np.pi * pc[ip, 1] * nz["/PLANES/surfArea"][:][ip]
        iw = np.unique(nz["/BCONDS/3/iCells"][:])
    with h5py.File(res, "r") as f:
        md = float(np.sum(f["/VALUE/roUx"][:][ic] * dA))
        Tw = f["/VALUE/T"][:][iw].astype(np.float64)
    return md, Tw


def series(runs):
    rows, off, last_w = [], 0, None
    for k, run in enumerate(runs):
        mesh = os.path.join(run, "nozzle.h5")
        ss = snaps(run)
        for st, p in ss:
            if k > 0 and st == 0:
                continue   # 延長の res_0 は前段の最終と同じ場
            tm = thrust_metrics(mesh, p, GAMMA, PT, PA, RT)
            md, Tw = mdot_and_wall(mesh, p)
            rows.append((off + st, tm["eta_cf"], md, float(Tw.max()), float(Tw.mean())))
            last_w = Tw
        off += ss[-1][0]
    return rows, last_w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gg", nargs="+", required=True)
    ap.add_argument("--lsq", nargs="+", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    res = {}
    for tag, runs in (("gg", a.gg), ("lsq", a.lsq)):
        rows, Tw = series(runs)
        res[tag] = (rows, Tw)
        with open(f"{a.out_prefix}_{tag}.csv", "w") as f:
            f.write("step,eta_cf,mdot,Tw_max,Tw_mean\n")
            for r in rows:
                f.write(",".join(f"{x:.9g}" for x in r) + "\n")
    print("| 双子 | step | η_CF | ṁ [kg/s] | 壁温 max [K] | 壁温 平均 [K] |")
    print("| --- | --- | --- | --- | --- | --- |")
    for tag in ("gg", "lsq"):
        r = res[tag][0][-1]
        print(f"| {tag} | {r[0]} | {r[1]:.5f} | {r[2]:.5f} | {r[3]:.1f} | {r[4]:.1f} |")
    g, l = res["gg"][0][-1], res["lsq"][0][-1]
    dT = np.abs(res["lsq"][1] - res["gg"][1]).max()
    print(f"\nlsq − gg: Δη_CF = {l[1] - g[1]:+.5f} ({(l[1] - g[1]) / g[1] * 100:+.3f} %)、"
          f"Δṁ/ṁ = {(l[2] - g[2]) / g[2] * 100:+.3f} %、壁温 L∞ = {dT:.2f} K")
    print("上限 (§6): η_CF ≤ 0.1 %、壁温 L∞ ≤ 15 K、ṁ ≤ 0.2 %。ゲート: η_CF ∈ 0.978 ± 0.003、ṁ ∈ 1.29–1.30 kg/s")


if __name__ == "__main__":
    main()
