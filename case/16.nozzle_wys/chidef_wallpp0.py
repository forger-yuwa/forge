#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B1 (iii) と B1-c (case/16): chi 明示 0 と省略 (auto 1) の壁圧比較。

前 plan V3 の式をそのまま使う (測る前に固定済み): p0 = 59070 Pa、輪郭壁 = `extract_wall_pp0.py` と同じ抽出
(wall_dist ≤ 0 ∧ y > 0 の各 x 列で y 最大のノード、座標は `MESH/COORD`)、x ∈ [10, 94] mm の L∞ (相対) ≤ 0.5 %、
3 点 (16.4 / 45.6 / 85 mm) の系列を `check_quasisteady --series-csv --drift 0.001 --osc 0.0025` で判定。
`--cond` のとき B1-c の onset (中心線 g が初めて 1e-3 を超える x、線形補間) と g_exit (報告のみ) も出し、
onset の差 ≤ 中心線ノード間隔 1 つで判定、onset 系列も同じ閾値で準定常判定。

  python3 chidef_wallpp0.py RUN_FLAG0 RUN_OMIT [--cond]
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
P0 = 59070.0
XP = (16.4, 45.6, 85.0)


def snaps(run):
    return sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(os.path.join(run, "res_*.h5"))
                  if re.search(r"/res_\d+\.h5$", p))


def lines(p):
    with h5py.File(p, "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3); P = np.array(f["VALUE/P"]); wd = np.array(f["VALUE/wall_dist"])
        g = np.array(f["VALUE/g_0"]) if "VALUE/g_0" in f else None
    xs = np.round(c[:, 0], 7)
    wx, wp, cx, cg = [], [], [], []
    for xv in np.unique(xs):
        k = np.where((xs == xv) & (wd <= 0) & (c[:, 1] > 0))[0]
        if len(k):
            j = k[np.argmax(c[k, 1])]; wx.append(c[j, 0] * 1e3); wp.append(P[j] / P0)
        if g is not None:
            k = np.where(xs == xv)[0]; j = k[np.argmin(np.abs(c[k, 1]))]; cx.append(c[j, 0] * 1e3); cg.append(g[j])
    return np.array(wx), np.array(wp), np.array(cx), np.array(cg)


def onset(cx, cg, thr=1e-3):
    o = np.argsort(cx); cx, cg = cx[o], cg[o]
    k = np.where(cg > thr)[0]
    if not len(k) or k[0] == 0:
        return float("nan"), float("nan")
    i = k[0]
    x = cx[i - 1] + (thr - cg[i - 1]) * (cx[i] - cx[i - 1]) / (cg[i] - cg[i - 1])
    return x, cx[i] - cx[i - 1]


def qs(csvp, cols):
    r = subprocess.run([sys.executable, os.path.join(REPO, "solver_density_cuda", "tools", "check_quasisteady.py"), "--series-csv", csvp,
                        "--series-cols", cols, "--drift", "0.001", "--osc", "0.0025"], capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if l.startswith("===") or l.startswith("  ")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("flag0"); ap.add_argument("omit"); ap.add_argument("--cond", action="store_true")
    a = ap.parse_args()
    print(f"# case/16 壁圧比較: flag0 = {a.flag0} / 省略 = {a.omit}{' (凝縮: onset も)' if a.cond else ''}")
    fin = {}
    for tag, r in (("flag0", a.flag0), ("omit", a.omit)):
        rows = []
        for st, p in snaps(r):
            wx, wp, cx, cg = lines(p)
            row = [st] + [float(np.interp(x, wx, wp)) for x in XP]
            if a.cond:
                row.append(onset(cx, cg)[0])
            rows.append(row)
            fin[tag] = (wx, wp, cx, cg)
        pth = os.path.join(r, "chidef_wallpp0.csv")
        with open(pth, "w") as f:
            f.write("step,pp0_16p4,pp0_45p6,pp0_85" + (",onset_mm" if a.cond else "") + "\n")
            for row in rows:
                f.write(",".join(f"{v:.9g}" for v in row) + "\n")
        for l in qs(pth, "pp0_16p4,pp0_45p6,pp0_85" + (",onset_mm" if a.cond else "")):
            print(f"[{tag}] {l}")
    (x0, p0, cx0, cg0), (x1, p1, cx1, cg1) = fin["flag0"], fin["omit"]
    m = (x0 >= 10) & (x0 <= 94)
    linf = float(np.max(np.abs(np.interp(x0[m], x1, p1) - p0[m]) / np.abs(p0[m])))
    ok = linf <= 0.005
    print(f"輪郭壁 p/p0 L∞ (x 10–94 mm、最終) = {linf * 100:.4f} % (許容 0.5 %) → {'PASS' if linf <= 0.005 else 'FAIL'}")
    print("3 点 (最終) flag0 / 省略: " + ", ".join(f"{x} mm {np.interp(x, x0, p0):.5f}/{np.interp(x, x1, p1):.5f}" for x in XP))
    if a.cond:
        o0, dx = onset(cx0, cg0); o1, _ = onset(cx1, cg1)
        good = np.isfinite(o0) and np.isfinite(o1) and abs(o1 - o0) <= dx
        ok &= good
        print(f"onset = {o0:.3f} / {o1:.3f} mm、差 {abs(o1 - o0):.4f} (許容 = 中心線ノード間隔 {dx:.4f} mm) → {'PASS' if good else 'FAIL'}")
        o = np.argsort(cx0)
        print(f"g_exit (報告のみ) = {cg0[o][-1]:.5f} / {cg1[np.argsort(cx1)][-1]:.5f}")
    print(f"VERDICT 量: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
