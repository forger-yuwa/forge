#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B1-a (case/40 軸対称ノズル): chi 明示 0 と省略 (auto 1) の比較。

量と許容 (plan §6 B1-a 表、測る前に固定):
  ① η_CF と ṁ (`s2_eval_case40.py` = `thrust_metrics`、README:172 と同じ抽出): 末尾 40 % の差区間 (usage-rule §4.2 の式) が ±0.1 % (相対) 以内
  ② 輪郭壁 p/p0 (p0 = 入口全圧 4 MPa、壁 physID 3 の節点): 最終スナップショットの L∞ ≤ 0.5 %
  ③ 等温壁 q_w: **本ケースの壁は断熱 (`kind: wall`) なので対象外** (記録)
  準定常: η_CF・ṁ `--drift 0.0002 --osc 0.0005`、壁 3 点 p/p0 (壁 x 範囲の 25/50/75 %) `--drift 0.001 --osc 0.0025`
  追加の成立条件: 壁∩軸ノード数を記録

  python3 chidef_b1a.py RUN_FLAG0 RUN_OMIT
"""
import csv
import glob
import os
import re
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "case", "09.Taylor-Green", "_g0_lsq_seam"))
import s2_eval_case40 as E  # noqa: E402

P0 = 4.0e6
TAIL = 0.4


def wall_nodes(mesh):
    with h5py.File(mesh, "r") as f:
        iw = np.unique(f["/BCONDS/3/iCells"][:])
        ia = np.unique(f["/BCONDS/4/iCells"][:]) if "/BCONDS/4" in f else np.array([], int)
        c = f["/MESH/COORD"][:].reshape(-1, 3)
    return iw, np.intersect1d(iw, ia), c


def snaps(run):
    return sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(os.path.join(run, "res_*.h5"))
                  if re.search(r"/res_\d+\.h5$", p))


def wall_series(run, iw, c):
    order = iw[np.argsort(c[iw, 0])]
    xs = c[order, 0]
    probes = [order[np.argmin(np.abs(xs - (xs.min() + q * (xs.max() - xs.min()))))] for q in (0.25, 0.5, 0.75)]
    rows, last = [], None
    for st, p in snaps(run):
        with h5py.File(p, "r") as f:
            P = f["/VALUE/P"][:].astype(np.float64)
        rows.append([st] + [P[n] / P0 for n in probes])
        last = P[order] / P0
    return rows, last, xs


def interval(v0, v1):
    k0 = max(3, int(np.ceil(TAIL * len(v0)))); k1 = max(3, int(np.ceil(TAIL * len(v1))))
    t0, t1 = np.array(v0[-k0:]), np.array(v1[-k1:])
    return t1.min() - t0.max(), t1.max() - t0.min(), t0.mean()


def qs(csvp, cols, drift, osc):
    r = subprocess.run([sys.executable, os.path.join(REPO, "solver_density_cuda", "tools", "check_quasisteady.py"), "--series-csv", csvp,
                        "--series-cols", cols, "--drift", str(drift), "--osc", str(osc)], capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if l.startswith("===") or l.startswith("  ")]


def main():
    r0, r1 = sys.argv[1], sys.argv[2]
    mesh = os.path.join(r0, "nozzle.h5")
    iw, iwa, c = wall_nodes(mesh)
    print(f"# B1-a 比較: flag0 = {r0} / 省略 = {r1}")
    print(f"追加の成立条件 (記録): 壁ノード {len(iw)}、壁∩軸ノード {len(iwa)}")
    ser = {}
    for tag, r in (("flag0", r0), ("omit", r1)):
        rows, _ = E.series([r])
        p = os.path.join(r, "chidef_b1a_eta.csv")
        with open(p, "w") as f:
            f.write("step,eta_cf,mdot\n")
            for s in rows:
                f.write(f"{s[0]},{s[1]:.9g},{s[2]:.9g}\n")
        ser[tag] = rows
        for l in qs(p, "eta_cf,mdot", 0.0002, 0.0005):
            print(f"[{tag}] {l}")
        wr, last, xs = wall_series(r, iw, c)
        pw = os.path.join(r, "chidef_b1a_wall.csv")
        with open(pw, "w") as f:
            f.write("step,pp0_25,pp0_50,pp0_75\n")
            for s in wr:
                f.write(",".join(f"{v:.9g}" for v in s) + "\n")
        for l in qs(pw, "pp0_25,pp0_50,pp0_75", 0.001, 0.0025):
            print(f"[{tag}] {l}")
        ser[tag + "_wall"] = last
    ok = True
    for j, nm in ((1, "η_CF"), (2, "ṁ")):
        lo, hi, m0 = interval([s[j] for s in ser["flag0"]], [s[j] for s in ser["omit"]])
        rl, rh = lo / m0, hi / m0
        good = -0.001 <= rl and rh <= 0.001
        ok &= good
        print(f"① {nm}: 差区間 [{rl * 100:+.4f}, {rh * 100:+.4f}] % (許容 ±0.1 %) → {'帯内' if good else '帯外/判定不能'}")
    linf = np.max(np.abs(ser["omit_wall"] - ser["flag0_wall"]) / np.abs(ser["flag0_wall"]))
    good = linf <= 0.005
    ok &= good
    print(f"② 輪郭壁 p/p0 L∞ = {linf * 100:.4f} % (許容 0.5 %) → {'PASS' if good else 'FAIL'}")
    print("③ 等温壁 q_w: 本ケースは断熱壁のため対象外")
    print(f"VERDICT B1-a 量: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
