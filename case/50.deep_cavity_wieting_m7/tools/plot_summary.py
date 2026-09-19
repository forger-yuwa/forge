#!/usr/bin/env python3
"""幅依存のまとめ図: 開口面積平均比 q̄_c/q_fp vs w/d (forge / 理論 / 実験)。"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from burggraf import qs_over_qfp
from compare_ref import jp_font

CASE = Path(__file__).resolve().parents[1]
RUNS = {"0.063": "run_0006_T1_wd0063_long", "0.211": "run_0008_T1_wd0211_long", "0.383": "run_0009_T1_wd0383_long"}
EXP = {"0.063": 1.07, "0.211": 0.73, "0.383": 0.58, "0.524": 0.49}
d = 20.32e-3
W = json.loads((CASE / "geometry.json").read_text())["cavity"]["widths"]

def theory_avg(w):
    s = np.linspace(1e-6, 2 * d + w, 40001)
    return float(np.trapz(qs_over_qfp(s, w, d), s) / w)

jp_font()
import matplotlib.pyplot as plt
wd_t = np.linspace(0.05, 0.56, 60)
th = [theory_avg(x * d) for x in wd_t]
fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=140)
ax.plot(wd_t, th, "-.", lw=2, color="#4C7A34", label="Burggraf 理論 (周長積分)")
ax.plot(list(map(float, EXP)), list(EXP.values()), "o", ms=9, mfc="none", mew=1.8,
        color="#C4502A", label="実験 (W70 Fig 12)")
xs, ys = [], []
for k, r in RUNS.items():
    p = CASE / r / "cavity_eval.json"
    if p.exists():
        xs.append(float(k)); ys.append(json.loads(p.read_text())["qbar_over_lit"])
ax.plot(xs, ys, "s-", ms=8, lw=2, color="#2E6F9E", label="forge (層流 2D, 分母=文献 q_fp)")
ax.axhline(0.56, color="0.55", ls=":", lw=1.2)
ax.text(0.47, 0.575, "Chapman 0.56", fontsize=9, color="0.4")
ax.annotate("理論は適用外\n(粘性コア)", xy=(0.075, 1.02), xytext=(0.13, 1.12), fontsize=9.5,
            color="#C4502A", arrowprops=dict(arrowstyle="->", color="#C4502A", lw=1.2))
ax.set_xlabel("w / d  (幅 ÷ 深さ)"); ax.set_ylabel("$\\bar q_c / q_{fp}$  (開口面積平均)")
ax.set_xlim(0, 0.58); ax.set_ylim(0.3, 1.25); ax.grid(alpha=.25)
ax.legend(fontsize=9.5); ax.set_title("深キャビティの平均熱流束比 — 幅依存 (M 6.9, 層流, 冷壁)", fontsize=11)
fig.tight_layout(); fig.savefig(CASE / "width_sweep_summary.png")
print("→", CASE / "width_sweep_summary.png")
for k, y in zip(xs, ys):
    print(f"  w/d {k:.3f}: forge {y:.3f} / 理論 {theory_avg(k*d):.3f} / 実験 {EXP[f'{k:.3f}']:.2f}")
