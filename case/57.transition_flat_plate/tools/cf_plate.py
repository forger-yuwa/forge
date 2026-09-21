#!/usr/bin/env python3
"""遷移平板の壁面摩擦係数 Cf(x) を forge / SU2 の場から**同じ式**で出し、実験・Blasius・乱流相関と重ねる。

  tau_w = mu * u(y1) / y1   (壁から第一内部節点。y1+ < 1 なので粘性底層の直線で十分)
  Cf    = tau_w / (0.5 rho_inf U_inf^2)   (基準は入口の一様流。plan turbulence-transition-lm2009 §6)
遷移開始位置 x_onset = Cf が層流側で最小になる x、遷移終了 x_end = その下流で Cf が最大になる x。

usage: cf_plate.py [--forge RUN[:label] ...] [--su2 DIR[:label] ...] [--exp ref/t3a_exp.dat] [--out cf.png] [--csv cf.csv]
"""
import argparse, csv, glob, os, sys
import numpy as np

RO, U, MU = 101325.0 / (287.0 * 300.0), 69.44, 2.269e-4
QINF = 0.5 * RO * U * U


def two_point(x, y, u, mu=None):
    mu = MU if mu is None else mu
    """構造格子前提: 壁節点 (y=0, x>=0) ごとに同じ x の最小 y>0 の節点を探す。"""
    wall = np.where((np.abs(y) < 1e-12) & (x >= -1e-12))[0]
    order = np.lexsort((y, np.round(x, 9)))
    xs, ys, us = np.round(x[order], 9), y[order], u[order]
    out = []
    for iw in wall:
        lo = np.searchsorted(xs, round(x[iw], 9), "left"); hi = np.searchsorted(xs, round(x[iw], 9), "right")
        yy, uu = ys[lo:hi], us[lo:hi]; m = yy > 1e-12
        if not m.any(): continue
        j = np.argmin(yy[m]); out.append((x[iw], mu * uu[m][j] / yy[m][j], yy[m][j]))
    a = np.array(sorted(out)); return a[:, 0], a[:, 1], a[:, 2]


def load_forge(run, step=None):
    import h5py
    fs = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
    f = os.path.join(run, f"res_{step}.h5") if step is not None else fs[-1]
    with h5py.File(f, "r") as h:
        c = h["MESH/COORD"][:].reshape(-1, 3); ux = h["VALUE/Ux"][:].astype(float)
    return two_point(c[:, 0].astype(float), c[:, 1].astype(float), ux), f


def load_su2(d):
    f = os.path.join(d, "restart_flow.csv")
    rows = list(csv.reader(open(f))); hd = [s.strip().strip('"') for s in rows[0]]; a = np.array(rows[1:], float)
    col = lambda n: a[:, hd.index(n)]
    return two_point(col("x"), col("y"), col("Momentum_x") / col("Density")), f


def onset(x, cf, xmin=0.01):
    """遷移開始 = 前縁から見て**最初の** Cf の極小、遷移終了 = その下流で最初の極大 (T3B のように遷移後の Cf が板端まで下がり続けても拾える)。
    極値が無ければ (完全乱流・完全層流) 板端を返す。"""
    m = x > xmin; xx, cc = x[m], cf[m]
    d = np.diff(cc)
    i0 = next((i for i in range(1, len(d)) if d[i - 1] < 0 <= d[i]), len(cc) - 1)
    i1 = next((i for i in range(i0 + 1, len(d)) if d[i - 1] > 0 >= d[i]), len(cc) - 1)
    return xx[i0], cc[i0], xx[i1], cc[i1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", nargs="*", default=[]); ap.add_argument("--su2", nargs="*", default=[])
    ap.add_argument("--exp", default=None); ap.add_argument("--out", default=None); ap.add_argument("--csv", default=None)
    ap.add_argument("--stations", default="0.3,0.9,1.3"); ap.add_argument("--mu", type=float, default=2.269e-4, help="粘性 [Pa s] (T3A 2.269e-4, T3B 1.3033e-4)")
    ap.add_argument("--series", nargs="*", default=None, help="時系列モード: 連続する run を順に並べる"); ap.add_argument("--series-out", default="cf_series.csv")
    a = ap.parse_args()
    globals()["MU"] = a.mu
    st = [float(s) for s in a.stations.split(",")]
    if a.series:
        # 連続する run の全スナップショットから報告量の時系列を書く (check_quasisteady.py --series-csv の入力)
        import h5py
        rows, off = [], 0
        for run in a.series:
            fs = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
            for f in fs:
                n = int(os.path.basename(f)[4:-3])
                if n == 0: continue
                (x, tw, y1), _ = load_forge(run, n); cf = tw / QINF; xo, cmin, xe, cmax = onset(x, cf)
                rows.append([off + n, xo, cmin, xe, cmax] + [float(np.interp(s_, x, cf)) for s_ in st])
            off += int(os.path.basename(fs[-1])[4:-3])
        with open(a.series_out, "w") as g:
            g.write("step,x_onset,cf_min,x_end,cf_max," + ",".join("cf_x%.1f" % s_ for s_ in st) + "\n")
            for r in rows: g.write(",".join(f"{v:.8g}" for v in r) + "\n")
        print("wrote", a.series_out, len(rows), "snapshots"); return
    curves = []
    for spec in a.forge:
        run, _, lab = spec.partition(":"); (x, tw, y1), f = load_forge(run); curves.append((lab or os.path.basename(run.rstrip("/")), x, tw / QINF, y1, f, "forge"))
    for spec in a.su2:
        d, _, lab = spec.partition(":"); (x, tw, y1), f = load_su2(d); curves.append((lab or os.path.basename(d.rstrip("/")), x, tw / QINF, y1, f, "su2"))
    print(f"{'curve':<28}{'x_onset[m]':>11}{'Cf_min':>10}{'x_end[m]':>10}{'Cf_max':>10}" + "".join(f"{'Cf@%.2f' % s:>11}" for s in st) + f"{'max y1+':>9}  file")
    for lab, x, cf, y1, f, kind in curves:
        xo, cmin, xe, cmax = onset(x, cf)
        utau = np.sqrt(np.maximum(cf, 0) * QINF / RO); yp = (y1 * utau * RO / MU)[x > 0.02].max()
        print(f"{lab:<28}{xo:11.4f}{cmin:10.5f}{xe:10.4f}{cmax:10.5f}" + "".join(f"{np.interp(s, x, cf):11.5f}" for s in st) + f"{yp:9.2f}  {f}")
    if a.csv:
        with open(a.csv, "w") as g:
            for lab, x, cf, *_ in curves:
                for xi, ci in zip(x, cf): g.write(f"{lab},{xi:.6e},{ci:.6e}\n")
    if a.out:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5.2))
        xx = np.linspace(0.005, 1.5, 400); rex = RO * U * xx / MU
        ax.plot(xx, 0.664 / np.sqrt(rex), "k:", lw=1, label="Blasius 0.664 Re_x^-1/2")
        ax.plot(xx, 0.0576 * rex ** -0.2, "k--", lw=1, label="turbulent 0.0576 Re_x^-1/5")
        if a.exp:
            e = np.genfromtxt(a.exp, comments="#"); ax.plot(e[:, 0] * 1e-3, e[:, 1], "ko", ms=5, mfc="none", label="ERCOFTAC T3A (exp.)")
        for lab, x, cf, y1, f, kind in curves:
            ax.plot(x, cf, "-" if kind == "forge" else "-.", lw=1.6, label=lab)
        ax.set_xlim(0, 1.5); ax.set_ylim(0, 0.01); ax.set_xlabel("x [m]  (Re/m = %.3g)" % (RO * U / MU)); ax.set_ylabel("Cf"); ax.grid(alpha=0.3)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False); fig.tight_layout(); fig.savefig(a.out, dpi=140); print("wrote", a.out)


if __name__ == "__main__":
    main()
