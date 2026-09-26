#!/usr/bin/env python3
"""R1 (plan boundary-node-periodic-gradient-fix §6 R1) の抽出: snapshot 時系列の CSV を作る。

    python3 r1_extract.py RUN_DIR [--out r1_series.csv] [--min-step 0]

出力 `RUN_DIR/r1_series.csv` の列は §6 R1 で固定したとおり:

    step,Cf_x05,Cf_x2,Cf_x6,xr_h,r_gradu,r_gradk,r_gradw,dF1_inf

副出力 `RUN_DIR/r1_bulk.csv` (step,rho_b,U_b,Re_b) は正規化に使ったバルク量の記録。
判定は `check_quasisteady.py --series-csv RUN_DIR/r1_series.csv --series-cols ... --drift 0.002 --osc 0.005 --tail 0.4`。

抽出規則 (§6 R1 の文言をそのまま実装。**解釈を加えた箇所は [解釈] と書く**):

- $C_f=\tau_{w,t}/(\tfrac12\rho_bU_b^2)$。$\tau_{w,t}$ は下壁 (`res_ylo_3_<step>.h5`) の `twall` を
  +x 向きの壁接線 $\mathbf t=(1,y_w'(x),0)/\sqrt{1+y_w'^2}$ へ射影する。符号は
  `design/forge_design/metrics/sern_forces.py` の `twall_on_fluid` と同じ規約: node は False
  (twall = 壁ノードの CV に入る力、+x の付着流で twall_x>0) なので $\tau_{w,t}=+\,\mathbf{twall}\cdot\mathbf t$。
  $y_w(x)$ は `mesh/make_hill_mesh.py` の `wall_y` (メッシュ生成と同じ多項式) を中心差分する。
- $\rho_b,U_b$: 丘頂断面 (x=0 の節点列、y∈[h, 3.035h]) のバルク。各 z 列で y 方向に台形積分し
  $\rho_b=\int\rho\,dy/H_c$、$U_b=\int\rho U_x\,dy/\int\rho\,dy$ ([解釈] 質量流束で重み付けたバルク速度)、
  z 平均は一意 DOF 平均 (z 継ぎ目の重複 z=0/z=L_z は重み 1/2)。x=L_x の重複列は使わない。
- $C_f(x)$ の z 平均も同じく一意 DOF (z 継ぎ目は重み 1/2)。$x/h=0.5,2,6$ で壁ノード間の線形補間。
- $x_r$: 下壁 $C_f$ の負→正の最初のゼロ交差 ($x/h\in[1,8]$、線形補間)。交差なしは NaN
  (「未再付着」。§6 R1 は判定不能として R1 を落とさない)。
- 継ぎ目指標: z 継ぎ目の節点面 (z=0) と隣接内部の節点面 (z=Δz、[解釈] 第 1 内部層) を**同じ (x,y) 集合**で
  比べる。x 継ぎ目 (x=0, x=L_x) の節点は除く ([解釈] z 継ぎ目の効果だけを見るため)。壁ノードは含める。
  $|\nabla\mathbf u|$ = 速度勾配 9 成分の Frobenius ノルム、$|\nabla k|$・$|\nabla\omega|$ は `dKd*`・`dOmegad*`。
  比 = $\lVert g\rVert_2(\text{継ぎ目})/\lVert g\rVert_2(\text{隣接})$。
- $F_1$: **`sstF1` は出力変数に無い** (`variables.hpp` の `output_cellValNames` に含まれない) ので、
  `ransSource_d.cu` の `rans_sst_blend_f1_d` と同じ式・同じ定数を、同じ snapshot の出力 (`ro,k,omega,vis_lam,
  wall_dist,dKd*,dOmegad*`) から float64 で再計算する ([解釈] 出力時点の k/ω は勾配を作った後の更新値なので、
  カーネルが読んだ値と 1 step 分ずれる。定常到達後は同じ)。dF1_inf = 同じ (x,y) 集合での $\max|F_1^{seam}-F_1^{adj}|$。
"""
import argparse
import csv
import glob
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "mesh"))
from make_hill_mesh import H, LX, LY_TOP, LZ, wall_y  # noqa: E402

MU = 1.617e-4          # solverConfig の physProp.visc (Re の記録用)
TWALL_ON_FLUID = False  # node の規約 (sern_forces.py)
# ransSource_d.cu の定数
BETA_STAR, SIGMA_W2, K_SMALL = 0.09, 0.856, 1.0e-12
XS_CF = (0.5, 2.0, 6.0)


def steps_of(run):
    out = []
    for p in glob.glob(os.path.join(run, "res_*.h5")):
        m = re.match(r"res_(\d+)\.h5$", os.path.basename(p))
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def keyxy(x, y):
    return np.round(x / H, 7) * 1.0e8 + np.round(y / H, 7)


def dup_weight_z(z, tol):
    return np.where((np.abs(z) < tol) | (np.abs(z - LZ) < tol), 0.5, 1.0)


def wall_slope(x):
    d = 1.0e-7
    xp = np.clip(x + d, 0.0, LX)
    xm = np.clip(x - d, 0.0, LX)
    return np.array([(wall_y(a) - wall_y(b)) / (a - b) for a, b in zip(xp, xm)])


def bulk(f, tol):
    c = f["MESH/COORD"][:].astype(np.float64).reshape(-1, 3)
    ro = f["VALUE/ro"][:].astype(np.float64)
    ux = f["VALUE/Ux"][:].astype(np.float64)
    sel = np.abs(c[:, 0]) < tol
    cs, rs, us = c[sel], ro[sel], ux[sel]
    zs = np.round(cs[:, 2] / H, 6)
    num_r = num_m = wsum = 0.0
    hc = LY_TOP - H
    for zv in np.unique(zs):
        m = zs == zv
        o = np.argsort(cs[m, 1])
        y = cs[m, 1][o]
        if y[0] > H * (1 + 1e-4) or y[-1] < LY_TOP * (1 - 1e-4):
            sys.exit("丘頂列が y∈[h,3.035h] を覆っていない (z=%g h)" % zv)
        ir = np.trapezoid(rs[m][o], y)
        im = np.trapezoid(rs[m][o] * us[m][o], y)
        w = 0.5 if (abs(zv) < 1e-6 or abs(zv - LZ / H) < 1e-6) else 1.0
        num_r += w * ir
        num_m += w * im
        wsum += w
    rho_b = num_r / wsum / hc
    u_b = num_m / num_r
    return rho_b, u_b


def cf_line(fw, rho_b, u_b, tol):
    c = fw["MESH/COORD"][:].astype(np.float64).reshape(-1, 3)
    n = fw["VALUE/twall_x"].shape[0]
    if len(c) != n:
        sys.exit("境界出力の節点数と値の数が合わない (%d vs %d)" % (len(c), n))
    tx = fw["VALUE/twall_x"][:].astype(np.float64)
    ty = fw["VALUE/twall_y"][:].astype(np.float64)
    x, z = c[:, 0], c[:, 2]
    s = wall_slope(x)
    nt = np.sqrt(1.0 + s * s)
    sgn = -1.0 if TWALL_ON_FLUID else 1.0
    tau = sgn * (tx * 1.0 / nt + ty * s / nt)
    w = dup_weight_z(z, tol)
    xk = np.round(x / H, 7)
    xu = np.unique(xk)
    cf = np.empty(len(xu))
    for i, xv in enumerate(xu):
        m = xk == xv
        cf[i] = np.sum(w[m] * tau[m]) / np.sum(w[m])
    cf /= 0.5 * rho_b * u_b * u_b
    return xu, cf


def reattach(xh, cf):
    m = (xh >= 1.0) & (xh <= 8.0)
    xs, cs = xh[m], cf[m]
    for i in range(len(xs) - 1):
        if cs[i] < 0.0 <= cs[i + 1]:
            return xs[i] - cs[i] * (xs[i + 1] - xs[i]) / (cs[i + 1] - cs[i])
    return float("nan")


def f1_of(g):
    rho = np.maximum(g["ro"], K_SMALL)
    k = np.maximum(g["k"], 0.0)
    w = np.maximum(g["omega"], K_SMALL)
    y = np.maximum(g["wall_dist"], K_SMALL)
    gkw = g["dKdx"] * g["dOmegadx"] + g["dKdy"] * g["dOmegady"] + g["dKdz"] * g["dOmegadz"]
    cd = np.maximum(2.0 * rho * SIGMA_W2 / w * gkw, 1.0e-10)
    nu = g["vis_lam"] / rho
    a = np.sqrt(k) / (BETA_STAR * w * y)
    b = 500.0 * nu / (w * y * y)
    cc = 4.0 * rho * SIGMA_W2 * k / (cd * y * y)
    arg = np.minimum(np.maximum(a, b), cc)
    return np.tanh(arg ** 4)


GRADU = ["dU%sd%s" % (i, j) for i in "xyz" for j in "xyz"]
GRADK = ["dKdx", "dKdy", "dKdz"]
GRADW = ["dOmegadx", "dOmegady", "dOmegadz"]
F1IN = ["ro", "k", "omega", "vis_lam", "wall_dist"] + GRADK + GRADW


def seam(f, tol):
    c = f["MESH/COORD"][:].astype(np.float64).reshape(-1, 3)
    n = f["VALUE/ro"].shape[0]
    if len(c) != n:
        sys.exit("VALUE 長 %d が節点数 %d と違う (node の res でない?)" % (n, len(c)))
    zu = np.unique(np.round(c[:, 2] / H, 6))
    dz = (zu[1] - zu[0]) * H
    notx = (c[:, 0] > tol) & (c[:, 0] < LX - tol)
    s = notx & (np.abs(c[:, 2]) < tol)
    a = notx & (np.abs(c[:, 2] - dz) < tol)
    ks, ka = keyxy(c[s, 0], c[s, 1]), keyxy(c[a, 0], c[a, 1])
    os_, oa = np.argsort(ks), np.argsort(ka)
    if len(ks) != len(ka) or not np.array_equal(ks[os_], ka[oa]):
        sys.exit("継ぎ目面と隣接面の (x,y) 集合が一致しない")
    names = sorted(set(GRADU + GRADK + GRADW + F1IN))
    g = {nm: f["VALUE/" + nm][:].astype(np.float64) for nm in names}
    gs = {nm: v[s][os_] for nm, v in g.items()}
    ga = {nm: v[a][oa] for nm, v in g.items()}

    def ratio(cols):
        ns = np.sqrt(sum(gs[c] ** 2 for c in cols))
        na = np.sqrt(sum(ga[c] ** 2 for c in cols))
        return float(np.sqrt(np.sum(ns ** 2)) / np.sqrt(np.sum(na ** 2)))

    df1 = float(np.max(np.abs(f1_of(gs) - f1_of(ga))))
    return ratio(GRADU), ratio(GRADK), ratio(GRADW), df1, int(s.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default="r1_series.csv")
    ap.add_argument("--min-step", type=int, default=0)
    a = ap.parse_args()
    tol = 1.0e-6 * H
    rows, brows = [], []
    for st in steps_of(a.run):
        if st < a.min_step:
            continue
        pv = os.path.join(a.run, "res_%d.h5" % st)
        pw = os.path.join(a.run, "res_ylo_3_%d.h5" % st)
        if not os.path.exists(pw):
            print("skip step %d (下壁出力なし)" % st)
            continue
        with h5py.File(pv, "r") as f, h5py.File(pw, "r") as fw:
            rho_b, u_b = bulk(f, tol)
            xh, cf = cf_line(fw, rho_b, u_b, tol)
            cfs = [float(np.interp(x0, xh, cf)) for x0 in XS_CF]
            xr = reattach(xh, cf)
            ru, rk, rw, df1, nseam = seam(f, tol)
        rows.append([st] + cfs + [xr, ru, rk, rw, df1])
        brows.append([st, rho_b, u_b, rho_b * u_b * H / MU])
        print("step %6d  Cf(0.5,2,6)=%+.4e %+.4e %+.4e  x_r/h=%.3f  r(u,k,w)=%.4f %.4f %.4f  dF1=%.3e  (seam pts %d, U_b %.3f)"
              % (st, cfs[0], cfs[1], cfs[2], xr, ru, rk, rw, df1, nseam, u_b))
    out = os.path.join(a.run, a.out)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "Cf_x05", "Cf_x2", "Cf_x6", "xr_h", "r_gradu", "r_gradk", "r_gradw", "dF1_inf"])
        w.writerows(rows)
    with open(os.path.join(a.run, "r1_bulk.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "rho_b", "U_b", "Re_b"])
        w.writerows(brows)
    print("wrote", out)


if __name__ == "__main__":
    main()
