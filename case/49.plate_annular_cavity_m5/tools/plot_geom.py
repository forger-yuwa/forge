#!/usr/bin/env python3
"""case/49 の形状と計算領域を寸法つきで描く (cad/geom_used.json を読む)。

usage: python3 tools/plot_geom.py [--used cad/geom_used.json] [--out geom_layout.png]
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Rectangle

HERE = Path(__file__).resolve().parent
for cand in ("NotoSansCJKjp-Regular.otf", "NotoSansCJK-Regular.ttc", "NotoSansJP-Regular.otf"):
    for p in Path.home().joinpath(".fonts").glob(cand):
        font_manager.fontManager.addfont(str(p))
        matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(p)).get_name()
        break

C = dict(fluid="#dbeafe", solid="#9ca3af", wall_iso="#dc2626", wall_ad="#1d4ed8",
         slip="#9333ea", inlet="#059669", outlet="#ea580c")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--used", default=str(HERE.parent / "cad" / "geom_used.json"))
    ap.add_argument("--out", default=str(HERE.parent / "geom_layout.png"))
    a = ap.parse_args()
    P = json.load(open(a.used))
    Ro, Ri, off, d = P["Ro"], P["Ri"], P["x_off"], P["depth"]
    xi, xo, ym, zt = P["x_in"], P["x_out"], P["y_max"], P["z_top"]
    rp = P.get("r_patch", 0.0)

    fig = plt.figure(figsize=(15.5, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.28, wspace=0.18)

    # ---------- (a) 対称面 y=0 の全体断面 ----------
    ax = fig.add_subplot(gs[0, :])
    ax.add_patch(Rectangle((xi, 0), xo - xi, zt, fc=C["fluid"], ec="none"))
    ax.add_patch(Rectangle((xi, -d - 12), xo - xi, d + 12, fc=C["solid"], ec="none", alpha=.55))
    # すきま (流体) 2 枚
    for s in (-1, +1):
        x0 = off + s * Ri if s > 0 else -Ro
        w = (Ro - off - Ri) if s > 0 else (off - Ri + Ro)
        ax.add_patch(Rectangle((off + Ri, -d), Ro - off - Ri, d, fc=C["fluid"], ec="none"))
        ax.add_patch(Rectangle((-Ro, -d), off - Ri + Ro, d, fc=C["fluid"], ec="none"))
    ax.plot([xi, xo], [0, 0], color=C["wall_ad"], lw=2.5, zorder=5)
    for xs in ([-Ro, off - Ri], [off + Ri, Ro]):
        ax.plot(xs, [-d, -d], color=C["wall_iso"], lw=2.5, zorder=6)
    for xv in (-Ro, off - Ri, off + Ri, Ro):
        ax.plot([xv, xv], [-d, 0], color=C["wall_iso"], lw=2.5, zorder=6)
    ax.plot([off - Ri, off + Ri], [0, 0], color=C["wall_iso"], lw=2.5, zorder=7)
    ax.plot([xi, xo], [zt, zt], color=C["slip"], lw=2.2)
    ax.plot([xi, xi], [0, zt], color=C["inlet"], lw=3)
    ax.plot([xo, xo], [0, zt], color=C["outlet"], lw=3)
    ax.annotate("", xy=(xi + 8, zt * .72), xytext=(xi - 2, zt * .72),
                arrowprops=dict(arrowstyle="-|>", color=C["inlet"], lw=2))
    ax.text(xi + 11, zt * .72, "M∞ (BL 分布を inletProfile で与える)", color=C["inlet"],
            va="center", fontsize=10)
    ax.text(xi + 2, zt * .93, "inlet", color=C["inlet"], fontsize=9)
    ax.text(xo - 18, zt * .93, "outlet", color=C["outlet"], fontsize=9)
    ax.text((xi + xo) / 2, zt + 2, "top: slip", color=C["slip"], fontsize=9, ha="center")
    ax.text(xi + (xo - xi) * .18, 3, "plate: 断熱 no-slip", color=C["wall_ad"], fontsize=10)
    ax.text(40, -d / 2, "内部円柱 (固体)", color="#374151", fontsize=10)
    ax.set_xlim(xi - 6, xo + 6); ax.set_ylim(-d - 14, zt + 10)
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("z [mm]")
    ax.set_title("(a) 対称面 y=0 の計算領域 (流体=青, 固体=灰)", fontsize=12, loc="left")

    # ---------- (b) 平面図 z=0 ----------
    ax = fig.add_subplot(gs[1, 0])
    ax.add_patch(Rectangle((xi, 0), xo - xi, ym, fc=C["fluid"], ec="none"))
    if rp > Ro:
        ax.add_patch(Circle((0, 0), rp, fc="#bfdbfe", ec="#60a5fa", ls="--", lw=1.2))
    ax.add_patch(Circle((0, 0), Ro, fc=C["fluid"], ec=C["wall_iso"], lw=2))
    ax.add_patch(Circle((off, 0), Ri, fc=C["solid"], ec=C["wall_iso"], lw=2, alpha=.85))
    ax.plot([xi, xo], [0, 0], color="k", ls="-.", lw=1.2)
    ax.text(xo - 40, 1.5, "対称面 y=0 (slip)", fontsize=9)
    ax.text(xi + 4, ym - 6, "plate (断熱)", color=C["wall_ad"], fontsize=9)
    if rp > Ro:
        ax.text(rp * .62, rp * .62, "plate_in\n(細分パッチ)", color="#1e40af", fontsize=8.5)
    ax.plot([xi, xo], [ym, ym], color=C["slip"], lw=2.2)
    ax.text((xi + xo) / 2, ym + 1.5, "side: slip", color=C["slip"], fontsize=9, ha="center")
    ax.set_xlim(xi - 6, xo + 6); ax.set_ylim(-6, ym + 8)
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title("(b) 平面図 z=0 (半割: y≥0 を計算)", fontsize=12, loc="left")

    # ---------- (c) キャビティ拡大 ----------
    ax = fig.add_subplot(gs[1, 1])
    m = 12
    ax.add_patch(Rectangle((-Ro - m, 0), 2 * (Ro + m), m, fc=C["fluid"], ec="none"))
    ax.add_patch(Rectangle((-Ro - m, -d - 8), 2 * (Ro + m), d + 8, fc=C["solid"], ec="none", alpha=.55))
    ax.add_patch(Rectangle((off + Ri, -d), Ro - off - Ri, d, fc=C["fluid"], ec="none"))
    ax.add_patch(Rectangle((-Ro, -d), off - Ri + Ro, d, fc=C["fluid"], ec="none"))
    for xv in (-Ro, off - Ri, off + Ri, Ro):
        ax.plot([xv, xv], [-d, 0], color=C["wall_iso"], lw=2.5)
    for xs in ([-Ro, off - Ri], [off + Ri, Ro]):
        ax.plot(xs, [-d, -d], color=C["wall_iso"], lw=2.5)
    ax.plot([off - Ri, off + Ri], [0, 0], color=C["wall_iso"], lw=2.5)
    ax.plot([-Ro - m, -Ro], [0, 0], color=C["wall_ad"], lw=2.5)
    ax.plot([Ro, Ro + m], [0, 0], color=C["wall_ad"], lw=2.5)
    g_dn, g_up = Ro - off - Ri, Ro + off - Ri
    ax.annotate("", xy=(Ro, -d * .45), xytext=(off + Ri, -d * .45),
                arrowprops=dict(arrowstyle="<|-|>", color="k", lw=1.2))
    ax.text((Ro + off + Ri) / 2, -d * .41, "%.2f" % g_dn, ha="center", fontsize=9)
    ax.annotate("", xy=(-Ro, -d * .7), xytext=(off - Ri, -d * .7),
                arrowprops=dict(arrowstyle="<|-|>", color="k", lw=1.2))
    ax.text((-Ro + off - Ri) / 2, -d * .66, "%.2f" % g_up, ha="center", fontsize=9)
    ax.annotate("", xy=(Ro + 6, 0), xytext=(Ro + 6, -d),
                arrowprops=dict(arrowstyle="<|-|>", color="k", lw=1.2))
    ax.text(Ro + 7.5, -d / 2, "深さ %.0f" % d, rotation=90, va="center", fontsize=9)
    ax.annotate("", xy=(Ro, m * .7), xytext=(-Ro, m * .7),
                arrowprops=dict(arrowstyle="<|-|>", color="#374151", lw=1.2))
    ax.text(0, m * .78, "φ%.0f (キャビティ)" % (2 * Ro), ha="center", fontsize=9)
    ax.annotate("", xy=(off + Ri, -d - 5), xytext=(off - Ri, -d - 5),
                arrowprops=dict(arrowstyle="<|-|>", color="#374151", lw=1.2))
    ax.text(off, -d - 4.2, "φ%.0f (内部円柱)" % (2 * Ri), ha="center", fontsize=9)
    ax.text(0, -d * .2, "等温壁 %g K" % 500, color=C["wall_iso"], ha="center", fontsize=10)
    ax.set_xlim(-Ro - m, Ro + m + 10); ax.set_ylim(-d - 10, m)
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("z [mm]")
    ax.set_title("(c) キャビティ拡大 (x_off=%.2f → すきま %.2f〜%.2f mm)"
                 % (off, min(g_dn, g_up), max(g_dn, g_up)), fontsize=12, loc="left")

    fig.suptitle("case/49 超音速平板 環状深キャビティ — 形状と計算領域 [mm]", fontsize=14)
    fig.savefig(a.out, dpi=135, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
