#!/usr/bin/env python3
"""入口に与えている分布 (inlet_profile CSV) を描く。"""
import sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
for p in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(p))
    matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(p)).get_name()
csv = sys.argv[1]; out = sys.argv[2]
d = np.genfromtxt(csv, names=True)
z = d["z"] * 1e3
T = d["Ps"] / (287.0 * d["ro"])
fig, ax = plt.subplots(1, 5, figsize=(15, 4.2), sharey=True)
for a, (y, lab, c) in zip(ax, [(d["Ux"], "$U_x$ [m/s]", "#1d4ed8"), (T, "$T$ [K]", "#b91c1c"),
                               (d["ro"], r"$\rho$ [kg/m³]", "#047857"), (d["k"], "$k$ [m²/s²]", "#7c3aed"),
                               (d["omega"], r"$\omega$ [1/s]", "#c2410c")]):
    a.plot(y, z, lw=2, color=c); a.set_xlabel(lab); a.grid(alpha=.3); a.set_ylim(0, 12)
    if lab.startswith("$\\omega"): a.set_xscale("log")
ax[0].set_ylabel("壁からの距離 z [mm]")
ax[0].axhline(5.056, ls="--", c="k", lw=1); ax[0].text(50, 5.3, r"$\delta_{99}$=5.06 mm", fontsize=9)
fig.suptitle("3D 入口 (x=-80 mm) に与えている分布 — 2D 前駆計算 (M5 断熱平板, x=313 mm) から", fontsize=12)
fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); print("wrote", out)
