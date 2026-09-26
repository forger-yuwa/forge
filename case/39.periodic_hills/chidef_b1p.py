#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B1-p (case/39 周期丘): chi 明示 0 と省略 (auto 1) の比較量を出す。
**結果を見る前に commit する** (plan §6「新規に書く抽出器は case/39 の Cf / x_r / seam 比だけ」)。

量 (plan §6 B1-p、許容も同表):
  ① 下壁 C_f(x) (`res_ylo_3_*.h5` の twall を壁接線へ射影、z は一意 DOF 平均)。½ρ_b U_b² は**両 run 共通に flag 0 の値**。
     比較は最終スナップショットの相対 L2 = ‖Cf_1 − Cf_0‖₂ / ‖Cf_0‖₂ (x 刻みは共通の壁ノード列)。許容 ≤ 1 %
  ② 再付着点 x_r (下壁 C_f の負→正ゼロ交差、x/h ∈ [1, 8])。許容 ≤ 下壁の局所ノード間隔 1 つ
  ③ 丘頂 (x = 0) の z 継ぎ目列 (z = 0) と隣接内部列 (z = Δz) の質量平均バルク速度の比 r = U_b,seam / U_b,adj。
     比較は |r_1 − r_0| / r_0 ≤ 0.5 % (case/16 V3 の値を流用 — 出典は流用と明記)
各 run の系列 (Cf 3 位置・x_r・r) は `chidef_b1p_series.csv` に書き、`check_quasisteady --series-csv --drift 0.002 --osc 0.005` で判定する。

  python3 chidef_b1p.py series RUN            # RUN/chidef_b1p_series.csv を作る
  python3 chidef_b1p.py compare RUN0 RUN1     # 最終スナップショットで ①②③
"""
import argparse
import csv
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r1_extract as R1  # noqa: E402

TOL = 1.0e-6 * R1.H


def seam_bulk_ratio(f):
    """丘頂 x = 0 の列で、z = 0 (継ぎ目) と z = Δz (最初の内部) の質量平均バルク速度の比。"""
    c = f["MESH/COORD"][:].astype(np.float64).reshape(-1, 3)
    ro = f["VALUE/ro"][:].astype(np.float64); ux = f["VALUE/Ux"][:].astype(np.float64)
    sel = np.abs(c[:, 0]) < TOL
    zs = np.unique(np.round(c[sel, 2], 9))
    z0, z1 = zs[0], zs[1]
    out = []
    for zv in (z0, z1):
        m = sel & (np.abs(c[:, 2] - zv) < TOL)
        o = np.argsort(c[m, 1]); y = c[m, 1][o]
        out.append(np.trapezoid(ro[m][o] * ux[m][o], y) / np.trapezoid(ro[m][o], y))
    return out[0] / out[1]


def snap(run, st, bulk_from=None):
    pv = os.path.join(run, f"res_{st}.h5"); pw = os.path.join(run, f"res_ylo_3_{st}.h5")
    with h5py.File(pv, "r") as f, h5py.File(pw, "r") as fw:
        rho_b, u_b = bulk_from if bulk_from else R1.bulk(f, TOL)
        xh, cf = R1.cf_line(fw, rho_b, u_b, TOL)
        r = seam_bulk_ratio(f)
    return xh, cf, R1.reattach(xh, cf), r, (rho_b, u_b)


def cmd_series(run):
    rows = []
    for st in R1.steps_of(run):
        if not os.path.exists(os.path.join(run, f"res_ylo_3_{st}.h5")):
            continue
        xh, cf, xr, r, _ = snap(run, st)
        rows.append([st] + [float(np.interp(x0, xh, cf)) for x0 in R1.XS_CF] + [xr, r])
    p = os.path.join(run, "chidef_b1p_series.csv")
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["step", "Cf_x05", "Cf_x2", "Cf_x6", "xr_h", "seam_ub_ratio"]); w.writerows(rows)
    print("wrote", p, len(rows), "rows")


def cmd_compare(run0, run1):
    st = min(max(R1.steps_of(run0)), max(R1.steps_of(run1)))
    x0, cf0, xr0, r0, b0 = snap(run0, st)
    x1, cf1, xr1, r1, _ = snap(run1, st, bulk_from=b0)   # ½ρU² は flag 0 の値を共通に使う
    if len(x0) != len(x1) or np.max(np.abs(x0 - x1)) > 1e-9:
        sys.exit("下壁のノード列が一致しない")
    l2 = np.linalg.norm(cf1 - cf0) / np.linalg.norm(cf0)
    i = int(np.argmin(np.abs(x0 - xr0))) if np.isfinite(xr0) else 0
    dx_local = float(np.max(np.diff(x0)[max(i - 1, 0):i + 1])) if np.isfinite(xr0) else float("nan")
    dxr = abs(xr1 - xr0)
    dr = abs(r1 - r0) / abs(r0)
    ok = (l2 <= 0.01) and (np.isfinite(dxr) and dxr <= dx_local) and (dr <= 0.005)
    print(f"# B1-p 比較 @ step {st}: flag0 = {run0} / 省略 = {run1}")
    print(f"① 下壁 Cf(x) 相対 L2 = {l2 * 100:.4f} % (許容 1 %)")
    print(f"② x_r/h = {xr0:.4f} / {xr1:.4f}、差 {dxr:.4f} (許容 = 局所ノード間隔 {dx_local:.4f} h)")
    print(f"③ 継ぎ目/隣接 バルク速度比 = {r0:.6f} / {r1:.6f}、相対差 {dr * 100:.4f} % (許容 0.5 %、case/16 V3 からの流用)")
    print(f"VERDICT B1-p 量: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("series", "compare")); ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    cmd_series(a.runs[0]) if a.cmd == "series" else cmd_compare(a.runs[0], a.runs[1])
