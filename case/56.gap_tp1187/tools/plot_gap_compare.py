#!/usr/bin/env python3
"""すきま内の鉛直速度 $U_y$ を forge と SU2 で比べる図。

**なぜこの量か**: 閉じた袋小路の定常解では、どの断面でも正味の質量流束が 0 でなければならない
(下向きがあれば戻りがある)。$U_y$ が幅全体で同符号なら $\\nabla\\cdot(\\rho u)\\neq0$ で、
その場は定常解ではなく**充填擬似過渡**である。深さ別の壁熱流束だけを見ていると、この区別がつかない。

規約は skill `forge-contour` に従う (カラーマップは turbo、カラーバーは図に重ねない)。

usage: plot_gap_compare.py FORGE_RES.h5 SU2_VOL.vtu --out fig.png
"""
import argparse, sys
from pathlib import Path
import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from su2_gap_eval import read_vtu_appended

W = 0.18e-2          # すきま幅 [m]
D = 6.35e-2          # すきま深さ [m]
XW = 0.5 * W         # 壁の x 位置


def load_forge(res):
    with h5py.File(res) as h:
        c = h["/MESH/COORD"][:].reshape(-1, 3)
        uy = h["/VALUE/Uy"][:].astype(float)
    m = (c[:, 1] < 0) & (np.abs(c[:, 0]) <= XW + 1e-9)
    return c[m, 0], c[m, 1], uy[m]


def load_su2(vtu):
    pts, pd = read_vtu_appended(vtu)
    c = pts[:, :2]
    V = np.asarray(pd["Velocity"])
    m = (c[:, 1] < 0) & (np.abs(c[:, 0]) <= XW + 1e-9)
    return c[m, 0], c[m, 1], V[m, 1].astype(float)


def band_mean(x, y, uy, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (-y > lo * W) & (-y <= hi * W)
        out.append(float(np.mean(uy[m])) if m.sum() else np.nan)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("forge"); ap.add_argument("su2"); ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in Path.home().glob(".fonts/NotoSansCJK*"):
        font_manager.fontManager.addfont(str(f))
    if any(f.name == "Noto Sans CJK JP" for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False

    fx, fy, fu = load_forge(a.forge)
    sx, sy, su = load_su2(a.su2)

    fig = plt.figure(figsize=(13.5, 7.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.5, 1.5], wspace=0.45)

    # --- (1)(2) すきま断面の Uy ---
    # **口元 (z/W<3) を外す**: そこは 178 m/s の剪断層で、深部 (mm/s 級) が色で潰れる。
    ZC = 3.0
    dm_f = -fy > ZC * W
    dm_s = -sy > ZC * W
    lim = float(np.percentile(np.abs(fu[dm_f]), 99.5))
    for k, (xx, yy, uu, lbl) in enumerate(((fx[dm_f], fy[dm_f], fu[dm_f], "forge"),
                                           (sx[dm_s], sy[dm_s], su[dm_s], "SU2"))):
        ax = fig.add_subplot(gs[0, k])
        sc = ax.scatter(xx * 1e3, yy * 1e2, c=uu * 1e3, s=3, cmap="turbo", vmin=-lim * 1e3, vmax=lim * 1e3)
        ax.set_xlabel("x [mm]"); ax.set_title(f"{lbl}\n$U_y$ [mm/s]  (z/W>3 のみ)", fontsize=10)
        ax.set_xlim(-1.0, 1.0); ax.set_ylim(-6.6, -ZC * W * 1e2)
        if k == 0:
            ax.set_ylabel("y [cm]  (0 = 板面、−6.35 = 床)")
        else:
            ax.set_yticklabels([])
        cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.09, fraction=0.05)
        cb.ax.tick_params(labelsize=7)

    # --- (3) 幅方向プロファイル ---
    ax = fig.add_subplot(gs[0, 2])
    for zW, ls in ((5, "-"), (15, "--"), (25, ":")):
        for xx, yy, uu, lbl, col in ((fx, fy, fu, "forge", "#c0504d"), (sx, sy, su, "SU2", "#1f4e79")):
            m = np.abs(-yy - zW * W) < 0.15e-3
            if m.sum() < 5:
                continue
            o = np.argsort(xx[m])
            ax.plot(xx[m][o] * 1e3, uu[m][o] * 1e3, ls, color=col, lw=1.4,
                    label=f"{lbl} z/W={zW}" if ls == "-" else None)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("$U_y$ [mm/s]")
    ax.set_title("幅方向プロファイル\n(実線 z/W=5、破線 15、点線 25)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="lower right")

    # --- (4) 帯平均 Uy vs 深さ ---
    ax = fig.add_subplot(gs[0, 3])
    edges = np.array([1, 3, 5, 10, 20, 35])
    mid = 0.5 * (edges[:-1] + edges[1:])
    for xx, yy, uu, lbl, col in ((fx, fy, fu, "forge", "#c0504d"), (sx, sy, su, "SU2", "#1f4e79")):
        bm = band_mean(xx, yy, uu, edges)
        ax.plot(np.abs(bm) * 1e3, mid, "o-", color=col, label=lbl)
    ax.set_xscale("log"); ax.invert_yaxis()
    ax.set_xlabel(r"帯平均 $|\overline{U_y}|$ [mm/s]"); ax.set_ylabel("z/W")
    ax.set_title("深さ別の正味の下向き流れ\n(定常解ならここは 0)", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    fig.suptitle("すきま内の鉛直速度 — forge は幅全体で同符号 (戻り流れ無し = 定常解でない)、SU2 は静止", fontsize=12)
    fig.savefig(a.out, dpi=120, bbox_inches="tight")
    print(f"-> {a.out}")
    print(f"   forge |Uy| 最大 {np.abs(fu).max()*1e3:.3f} mm/s / SU2 {np.abs(su).max()*1e3:.3e} mm/s")


if __name__ == "__main__":
    main()
