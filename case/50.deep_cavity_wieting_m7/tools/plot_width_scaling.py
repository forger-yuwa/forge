#!/usr/bin/env python3
"""§4.4c の図 — 比で見たときと絶対量で見たときの幅依存性。"""
import json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in pathlib.Path.home().glob(".fonts/NotoSansCJK*"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

CASE = pathlib.Path(__file__).resolve().parents[1]
j = json.loads((CASE / "width_scaling.json").read_text(encoding="utf-8"))
r = j["rows"]
wd = np.array([x["wd"] for x in r]); w = np.array([x["w"] for x in r]) * 1e3
Qf = np.array([x["Qb"] for x in r]); Qm = np.array([x["Qm"] for x in r])
rl = np.array([x["ratio_lit"] for x in r]); rb = np.array([x["ratio_b"] for x in r])

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax[0].plot(wd, rl, "o-", color="#1f4e79", label="実測 [W70] Fig 12")
ax[0].plot(wd, rb, "s-", color="#c0504d", label="forge (前壁=Burggraf 置換)")
ax[0].set_xlabel("$w/d$"); ax[0].set_ylabel(r"$\bar q_c / q_{fp}$")
ax[0].set_title("比で見たとき — 狭いすきまほど食い違う", fontsize=11)
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3); ax[0].set_ylim(0, 1.2)
for x, a, b in zip(wd, rl, rb):
    ax[0].annotate(f"×{b/a:.2f}", (x, (a + b) / 2), fontsize=8, ha="center", color="#666")

ax[1].plot(Qf, Qm, "o", ms=9, color="#1f4e79", zorder=5, label="4 幅のデータ")
xs = np.linspace(0, 175, 50)
ax[1].plot(xs, xs, "k--", lw=1, label="$Q_{meas}=Q_{forge}$")
ax[1].plot(xs, xs + j["fit_offset"], "-", color="#c0504d",
           label=f"加算  $+{j['fit_offset']:.1f}$ W/m (RMS {j['rms_offset']:.1f})")
ax[1].plot(xs, j["fit_proportional"] * xs, ":", color="#7f7f7f",
           label=f"比例  ×{j['fit_proportional']:.3f} (RMS {j['rms_proportional']:.1f})")
for x, y, t in zip(Qf, Qm, wd):
    ax[1].annotate(f"$w/d$={t:.3f}", (x, y), textcoords="offset points",
                   xytext=(7, -11), fontsize=8, color="#444")
ax[1].set_xlabel("$Q_{forge}$ [W/m]"); ax[1].set_ylabel("$Q_{meas}$ [W/m]")
ax[1].set_title("絶対量で見たとき — 食い違いは幅に依らない加算量", fontsize=11)
ax[1].legend(fontsize=8.5, loc="upper left"); ax[1].grid(alpha=0.3)
ax[1].set_xlim(0, 175); ax[1].set_ylim(0, 185)

fig.suptitle("case/50 — すきま入熱の幅依存性 (前壁は W70 と同じく Burggraf 置換、4 run とも $Q_c$ STEADY)",
             fontsize=11.5)
fig.tight_layout()
out = CASE / "width_scaling.png"
fig.savefig(out, dpi=150)
print(f"→ {out}")
