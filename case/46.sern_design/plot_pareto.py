#!/usr/bin/env python3
"""MOO キャンペーンのパレート図: plot_pareto.py <campaign_dir> "<title>" """
import json, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in Path.home().glob(".fonts/NotoSansCJKjp-Regular.otf"): font_manager.fontManager.addfont(str(f)); plt.rcParams["font.family"]="Noto Sans CJK JP"
d = Path(sys.argv[1]); title = sys.argv[2] if len(sys.argv) > 2 else d.name
rows=[json.loads(l) for l in open(d/"ledger.jsonl")]; ok=[r for r in rows if r["status"]=="PASS"]; inf=[r for r in rows if r["status"]=="INFEASIBLE" and r.get("L_ramp")]
s=json.load(open(d/"pareto.json")); par=sorted(s["pareto"], key=lambda r: r["L_ramp"]); ops=list(ok[0]["ops"].keys())
fig,axs=plt.subplots(1,2,figsize=(12,4.3))
ax=axs[0]; sc=ax.scatter([r["L_ramp"] for r in ok],[r["C_T_w"] for r in ok], c=[r["C_M_w"] for r in ok], cmap="coolwarm", vmin=-6, vmax=0, s=45, edgecolor="k", lw=.4)
if inf: ax.scatter([r["L_ramp"] for r in inf],[r["C_T_w"] for r in inf], marker="x", c="gray", s=40, label="C_M 制約で除外")
ax.plot([r["L_ramp"] for r in par],[r["C_T_w"] for r in par],"k--",lw=1,label="パレート前線 (非劣解)")
for r in ok+inf: ax.annotate(r["tag"].replace("doe_","d").replace("inf_","i"), (r["L_ramp"], r["C_T_w"]), xytext=(4,3), textcoords="offset points", fontsize=7)
plt.colorbar(sc, ax=ax, label="C_M,w (x_ref = −20H, 頭上げ正)"); ax.set_xlabel("L_ramp / H"); ax.set_ylabel("C_T,w (重み付き)"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title(f"{title} ({s['n_eval']} 評価 / {s['n_pass']} PASS, HV {s['hv']:.3f})")
ax=axs[1]; ax.scatter([r["ops"][ops[0]]["C_T"] for r in ok],[r["ops"][ops[1]]["C_T"] for r in ok], c=[r["L_ramp"] for r in ok], cmap="viridis", s=45, edgecolor="k", lw=.4)
ax.set_xlabel(f"C_T {ops[0]}"); ax.set_ylabel(f"C_T {ops[1]}"); ax.grid(alpha=.3); ax.set_title("作動点ごとの推力 (色 = L_ramp/H)")
fig.tight_layout(); fig.savefig(d/"pareto.png", dpi=140); print("saved", d/"pareto.png")
