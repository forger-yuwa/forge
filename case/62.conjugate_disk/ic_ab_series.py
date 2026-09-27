#!/usr/bin/env python3
"""初期場 A/B (plan axisymmetric-graded-grid-static-gas §4.1) の比較量を毎 step の場から系列 CSV に集計する。

    python3 case/62.conjugate_disk/ic_ab_series.py <run>   -> <run>/ic_ab_series.csv

列: step, Uymax (max|U_y|), Uslip (slip 列 max|U|), Uint (内部 max|U|), dP_all (領域の Pmax−Pmin),
dP_mid (x=H/2 列の Pmax−Pmin), checker (共役壁 iface_q_eff の市松振幅 = 内部壁節点の値と左右隣接からの r 線形補間値の差の max)。
比較窓 step 250–500 の max の比 B/A が判定量 (登録 2026-09-27)。
"""
import glob, re, sys
import h5py, numpy as np
run = sys.argv[1].rstrip("/")
rows = []
for f in glob.glob(f"{run}/res_[0-9]*.h5"):
    st = int(re.search(r"res_(\d+)\.h5$", f).group(1))
    with h5py.File(f, "r") as h:
        c = h["MESH/COORD"][:].reshape(-1, 3); V = h["VALUE"]; Ux = V["Ux"][:]; Uy = V["Uy"][:]; P = V["P"][:]
    U = np.hypot(Ux, Uy); r = c[:, 1]
    slip = (np.abs(r - r.min()) < 1e-9) | (np.abs(r - r.max()) < 1e-9)
    xm = np.unique(np.round(c[:, 0], 9)); s = np.abs(c[:, 0] - xm[len(xm) // 2]) < 1e-9
    amp = np.nan
    try:
        with h5py.File(f"{run}/res_wall_cj_4_{st}.h5", "r") as w:
            rr = w["MESH/COORD"][:].reshape(-1, 3)[:, 1]; q = w["VALUE/iface_q_eff"][:]
        o = np.argsort(rr); rr, q = rr[o], q[o]
        amp = max(abs(q[i] - (q[i-1] + (q[i+1] - q[i-1]) * (rr[i] - rr[i-1]) / (rr[i+1] - rr[i-1]))) for i in range(1, len(rr) - 1))
    except OSError:
        pass
    rows.append((st, np.abs(Uy).max(), U[slip].max(), U[~slip].max(), P.max() - P.min(), P[s].max() - P[s].min(), amp))
rows.sort(); X = np.array(rows)
np.savetxt(f"{run}/ic_ab_series.csv", X, delimiter=",", header="step,Uymax,Uslip,Uint,dP_all,dP_mid,checker", comments="", fmt=["%d"] + ["%.10e"] * 6)
w = X[(X[:, 0] >= 250) & (X[:, 0] <= 500)]
print(f"{run}: {len(X)} steps, window 250-500 max: Uymax {w[:,1].max():.4e} checker {np.nanmax(w[:,6]):.4e}")
