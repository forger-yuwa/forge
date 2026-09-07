#!/usr/bin/env python3
"""S4(b): NASA TM X-71972 系列 (run_0003–0007) の力係数を表と図にまとめる。"""
import json, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in Path.home().glob(".fonts/NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(f)); plt.rcParams["font.family"] = "Noto Sans CJK JP"
here = Path(__file__).parent
runs = [("run_0003_nasa_Lc20_euler", 2.0, 6), ("run_0004_nasa_Lc312_euler", 3.12, 6), ("run_0005_nasa_Lc45_euler", 4.5, 6),
        ("run_0006_nasa_cowl3deg_euler", 3.12, 3), ("run_0007_nasa_cowl12deg_euler", 3.12, 12)]
rows = []
for rd, Lc, ang in runs:
    mp = here / rd / "metrics.json"
    if not mp.exists():
        continue
    m = json.loads(mp.read_text()); mo = m["moc_forces"]
    rows.append((rd, Lc, ang, m["C_T"], m["C_L"], m["C_M"], mo["C_T"], mo["C_L"], mo["C_M"],
                 ",".join(v["verdict"][0] for v in m["steadiness"].values()), (m["convergence_verdict"] or ["?"])[-1][:20]))
print("| run | L_cowl/H | θ_c [°] | C_T forge | C_T MOC | C_L forge | C_L MOC | C_M forge | C_M MOC | steady | conv |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]:.4f} | {r[6]:.4f} | {r[4]:.4f} | {r[7]:.4f} | {r[5]:.3f} | {r[8]:.3f} | {r[9]} | {r[10]} |")
if len(rows) >= 2:
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    a = [r for r in rows if r[2] == 6]; b = [r for r in rows if r[1] == 3.12]
    for ax, data, xkey, xl in ((axs[0], a, 1, "L_cowl / H (θ_c = 6°)"), (axs[1], b, 2, "θ_c [°] (L_cowl = 3.12H)")):
        data = sorted(data, key=lambda r: r[xkey]); xs = [r[xkey] for r in data]
        ax.plot(xs, [r[3] for r in data], "ro-", label="C_T forge"); ax.plot(xs, [r[6] for r in data], "r--", label="C_T MOC")
        ax2 = ax.twinx(); ax2.plot(xs, [r[5] for r in data], "bs-", label="C_M forge"); ax2.plot(xs, [r[8] for r in data], "b--", label="C_M MOC")
        ax.set_xlabel(xl); ax.set_ylabel("C_T", color="r"); ax2.set_ylabel("C_M (頭上げ正)", color="b"); ax.grid(alpha=.3)
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="best")
    fig.suptitle("NASA TM X-71972 基準形 (平板ランプ 20°/18.54H, M∞10) の傾向: 非粘性 Euler vs MOC")
    fig.tight_layout(); fig.savefig(here / "nasa_trends.png", dpi=140); print("saved nasa_trends.png")
