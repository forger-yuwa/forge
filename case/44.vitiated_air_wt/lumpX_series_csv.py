"""run_0509–0511 用の報告量の時系列 (check_quasisteady.py --series-csv 用) と NaN/物理性の検査。
usage: python3 lumpX_series_csv.py RUN_DIR → RUN_DIR/lumpX_series.csv, 標準出力に最終値
量: mdot_in / mdot_out [kg/s] (列の ∫ρUx 2πr dr), exit_M_avg / exit_T_avg (出口 x_max−2 r_t 列の質量流束平均),
axis_M_exit, axis_dM_max (x_A..x_E−0.3 r_t の |M_axis − M_target| 最大), g_max (凝縮 run のみ)。"""
import sys, glob, os, json, h5py, numpy as np
trapz = getattr(np, "trapezoid", None) or np.trapz
run = sys.argv[1]; S = 0.205711
nc = h5py.File(os.path.join(run, "nozzle.h5"))["MESH/COORD"][:].reshape(-1, 3)
x, r = nc[:, 0], nc[:, 1]
cols = np.unique(np.round(x, 9))
def col(xv):
    i = np.where(np.abs(x - xv) < 1e-9)[0]; return i[np.argsort(r[i])]
cin, cex = col(cols[0]), col(cols[np.argmin(abs(cols - (x.max() - 2 * S)))])
ax = np.where(r < 1e-12)[0]; ax = ax[np.argsort(x[ax])]
tgt = np.loadtxt(os.path.join(run, "target_axis_M.csv"), delimiter=",", skiprows=1)
info = json.load(open(os.path.join(run, "prepare_info.json")))
xa, xe = info["x_A"] * S, (info["x_E"] - 0.3) * S
res = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
rows = []; bad = []
for p in res:
    v = h5py.File(p)["VALUE"]; step = int(os.path.basename(p)[4:-3])
    for k in v.keys():
        a = v[k][:]
        if not np.all(np.isfinite(a)): bad.append((step, k, int((~np.isfinite(a)).sum())))
    ro, ux, uy, T, P = (v[k][:] for k in ("ro", "Ux", "Uy", "T", "P"))
    M = np.hypot(ux, uy) / v["sonic"][:]
    def mdot(c): return trapz(ro[c] * ux[c] * 2 * np.pi * r[c], r[c])
    w = ro[cex] * ux[cex] * 2 * np.pi * r[cex]
    avg = lambda q: trapz(w * q[cex], r[cex]) / trapz(w, r[cex])
    m = (x[ax] >= xa) & (x[ax] <= xe)
    dM = np.abs(M[ax][m] - np.interp(x[ax][m], tgt[:, 0], tgt[:, 1])).max()
    row = dict(step=step, mdot_in=mdot(cin), mdot_out=mdot(cex), exit_M_avg=avg(M), exit_T_avg=avg(T),
               axis_M_exit=M[cex[0]], axis_dM_max=dM)
    if "g_0" in v: row["g_max"] = float(v["g_0"][:].max())
    row["ro_min"], row["T_min"], row["P_min"] = ro.min(), T.min(), P.min()
    rows.append(row)
keys = list(rows[0].keys())
with open(os.path.join(run, "lumpX_series.csv"), "w") as f:
    f.write(",".join(keys) + "\n")
    for rw in rows: f.write(",".join(f"{rw[k]:.10g}" for k in keys) + "\n")
print(os.path.basename(run), "NaN/Inf:", bad if bad else "none", "| snapshots", [rw["step"] for rw in rows])
print("  final:", {k: (f"{rows[-1][k]:.6g}") for k in keys})
