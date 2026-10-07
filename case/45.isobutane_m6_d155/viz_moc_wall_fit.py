"""moc_wall_fit_probe の図: 補間壁 (MOC 1200 / 9600) と位置+壁角の同時当てはめ壁 (2400) の r″。
usage: design/.venv-opt/bin/python viz_moc_wall_fit.py → _band_ab/moc_wall_fit_r2.png
"""
import sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C))
from moc_wall_fit_probe import joint_fit, design_chain, load_problem  # noqa: E402 (probe を import すると測定も走る)

S = {}
for n in (1200, 9600, 2400):
    p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = n; d = design_chain(p); S[n] = d
fit, _ = joint_fit(S[2400]["wall_inv"], S[2400]["R"], lam=1e-7)
cur = [("補間 (生産, MOC 1200 → 334 点)", lambda x, k: S[1200]["wall"].r(x, k), "C0", 1.2),
       ("補間 (MOC 9600 → 2298 点)", lambda x, k: S[9600]["wall"].r(x, k), "C3", .8),
       ("位置+壁角の同時当てはめ (MOC 2400)", lambda x, k: fit(x, k), "k", 1.2)]
fig, ax = plt.subplots(3, 1, figsize=(11, 10))
for (xa, xb, ya, yb), a in zip(((0.5, 8, -0.02, 0.2), (20, 60, -0.008, 0), (60, 95, -0.0025, 0.0003)), ax):
    x = np.arange(xa, xb, 0.002)
    for lab, f, c, lw in cur:
        a.plot(x, f(x, 2), color=c, lw=lw, label=lab)
    a.set_xlim(xa, xb); a.set_ylim(ya, yb); a.grid(alpha=.3); a.set_ylabel("r″ [1/r_t]")
ax[0].set_title("壁の 2 階微分 r″ (形状のみ・境界層補正なし; problem_d155_ns_c2final)")
ax[2].set_xlabel("x / r_t")
fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(C / "_band_ab/moc_wall_fit_r2.png", dpi=120)
