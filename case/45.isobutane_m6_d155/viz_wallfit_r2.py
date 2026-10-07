"""補間壁 (生産) と位置+壁角の同時当てはめ壁の r″ を、縮流部 (上流 Hermite) から出口まで重ねて描く。Euler A/B (run_0053〜0064) と同じ壁 (MOC 2400 点)。
plan verification-m6-axis-wave-mesh-su2 §5.1 #13・#15。usage: design/.venv-opt/bin/python viz_wallfit_r2.py → _band_ab/wallfit_r2.png
"""
import copy, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from moc_wall_fit_ab import joint_fit  # noqa: E402

p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); d = design_chain(p)
WA = d["wall"]; WB = copy.copy(WA); WB._spl, _ = joint_fit(d["wall_inv"], d["R"], 1e-9)
tb = d["wall_inv"]
panels = [((-12.6, 96), (-0.3, 0.6), "全体 (x<0 は縮流部の 5 次 Hermite、x≥0 が MOC から作る区間)"),
          ((-0.4, 1.0), (0.3, 0.7), "スロート接合部の拡大 (x=0 で Hermite と接続)"),
          ((0.5, 8), (-0.02, 0.2), "スロート下流 x 0.5〜8"),
          ((20, 60), (-0.008, 0.0), "x 20〜60"),
          ((60, 95.2), (-0.0022, 0.0002), "x 60〜95 (出口側)")]
fig, ax = plt.subplots(len(panels), 1, figsize=(11, 16))
for a, ((xa, xb), yl, ttl) in zip(ax, panels):
    x = np.linspace(xa, xb, 20000)
    a.plot(x, WA.r(x, 2), "C0", lw=1.1, label="補間壁 (全点を通す; 現行の生産)")
    a.plot(x, WB.r(x, 2), "k", lw=1.1, label="当てはめ壁 (位置と壁角に近い滑らかな曲線)")
    m = (tb[:, 0] >= xa) & (tb[:, 0] <= xb)
    a.plot(tb[m, 0], np.zeros(m.sum()) + yl[0] + 0.03 * (yl[1] - yl[0]), "|", color="C3", ms=6, label="MOC の壁点の位置 (x)")
    a.axvline(0, color="gray", lw=.6, ls=":")
    a.set_xlim(xa, xb); a.set_ylim(*yl); a.grid(alpha=.3); a.set_title(ttl); a.set_ylabel("r″ [1/r_t]")
ax[-1].set_xlabel("x / r_t")
fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False)
fig.suptitle("壁の 2 階微分 r″ (設計壁 = 非粘性壁、problem_d155_euler_c2final_n2400)", y=0.995)
fig.tight_layout(rect=(0, 0.025, 1, 0.99)); fig.savefig(C / "_band_ab/wallfit_r2.png", dpi=110)
# 接合部の数値 (C² の確認)
for nm, W in (("interp", WA), ("fit", WB)):
    e = 1e-7
    print(nm, "x=0: r'' 左 %.6f 右 %.6f | r''' 左 %.4f 右 %.4f" % (W.r(np.r_[-e], 2)[0], W.r(np.r_[e], 2)[0], W.r(np.r_[-e], 3)[0], W.r(np.r_[e], 3)[0]),
          "| x_U=%.2f: r'' 左 %.2e 右 %.2e" % (-WA.up.L_U, W.r(np.r_[-WA.up.L_U - e], 2)[0], W.r(np.r_[-WA.up.L_U + e], 2)[0]))
