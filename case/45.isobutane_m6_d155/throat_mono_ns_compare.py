"""plan tooling-nozzle-throat-monotone-r2 §6 N・K の報告用: 旧壁と単調壁の NS / 凝縮 NS を同じ線上で比べる図 (AWS で回す)。
dry: 旧壁 run_0117 (res_40000〜60000) と単調壁 run_0147 + 延長 run_0149 (run_0147 res_60000 と run_0149 res_5000〜20000) の
     末尾 5 枚の平均で、M/M_d − 1 を r = 0 と r/r_w = 0.1 の線に沿って比べ、差 (単調 − 旧) も描く。
cond: 旧壁 run_0118 と単調壁 run_0148 の末尾 5 枚 (res_14000〜18000) の平均で、軸の凝縮質量分率 g_0 と過飽和度 S を比べる。
usage (AWS, case dir): python3 throat_mono_ns_compare.py [OUT_DIR (既定 _band_ab/throat_mono_compare)]
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402,F401  (load_field は日本語フォントも登録する)

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/throat_mono_compare"
OUT.mkdir(parents=True, exist_ok=True)
MD = 6.0
DRY = {"旧壁 run_0117": [("run_0117_ns_recal_final_ext", f"res_{s}.h5") for s in range(40000, 60001, 5000)],
       "単調壁 run_0147+0149": [("run_0147_ns_mono_final", "res_60000.h5")] + [("run_0149_ns_mono_final_ext", f"res_{s}.h5") for s in range(5000, 20001, 5000)]}
COND = {"旧壁 run_0118": [("run_0118_ns_recal_final_cond", f"res_{s}.h5") for s in range(14000, 18001, 1000)],
        "単調壁 run_0148": [("run_0148_ns_mono_final_cond", f"res_{s}.h5") for s in range(14000, 18001, 1000)]}


def mean_lines(items, keys, etas, xq):
    acc = {(k, e): [] for k in keys for e in etas}
    for run, res in items:
        F = load_field(C / run, res)
        for k in keys:
            for e in etas:
                acc[(k, e)].append(eta_line(F, k, e, xq))
    return {ke: np.mean(v, axis=0) for ke, v in acc.items()}, {ke: np.ptp(v, axis=0) for ke, v in acc.items()}


xq = np.linspace(0.0, 95.0, 1901)
D, Dr = {}, {}
for lab, items in DRY.items():
    D[lab], Dr[lab] = mean_lines(items, ["M"], [0.0, 0.1], xq)
(old, new) = DRY.keys()
fig, axs = plt.subplots(2, 1, figsize=(11, 7.2), sharex=True)
for lab, ls in ((old, "-"), (new, "--")):
    for e, c in ((0.0, "#1d2b3a"), (0.1, "#c2410c")):
        axs[0].plot(xq, 100 * (D[lab][("M", e)] / MD - 1), color=c, ls=ls, lw=1.1, label=f"{lab}  r/r_w = {e:g}")
axs[0].set_ylim(-0.4, 0.4); axs[0].set_ylabel("100 (M/M_d − 1)  [%]"); axs[0].grid(alpha=.3)
axs[0].set_title("試験部の M/M_d − 1 (末尾 5 枚の平均)。実線 = 旧壁、破線 = 単調壁")
for e, c in ((0.0, "#1d2b3a"), (0.1, "#c2410c")):
    d = 100 * (D[new][("M", e)] - D[old][("M", e)]) / MD
    noise = 100 * np.maximum(Dr[new][("M", e)], Dr[old][("M", e)]) / MD
    axs[1].plot(xq, d, color=c, lw=1.1, label=f"差 (単調 − 旧)  r/r_w = {e:g}")
    axs[1].fill_between(xq, -noise, noise, color=c, alpha=.12, label=f"末尾 5 枚の幅 (大きい方)  r/r_w = {e:g}")
axs[1].axhline(0, color="0.5", lw=.6); axs[1].set_ylabel("差  [%pt]"); axs[1].set_xlabel("x / r_t"); axs[1].grid(alpha=.3)
axs[1].set_ylim(-0.06, 0.06)
h0, l0 = axs[0].get_legend_handles_labels(); h1, l1 = axs[1].get_legend_handles_labels()
fig.legend(h0 + h1, l0 + l1, loc="lower center", ncol=3, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, .1, 1, 1)); fig.savefig(OUT / "fig_compare_dry_axis.png", dpi=150); plt.close(fig)
dsum = {}
for e in (0.0, 0.1):
    d = 100 * (D[new][("M", e)] - D[old][("M", e)]) / MD
    m = (xq >= 40) & (xq <= 94)
    dsum[f"eta{e:g}"] = {"max_abs_diff_pct_x0_95": float(np.abs(d).max()), "x_of_max": float(xq[np.argmax(np.abs(d))]),
                         "max_abs_diff_pct_test_section_40_94": float(np.abs(d[m]).max()),
                         "tail_range_pct_test_section_max": float(100 * max(Dr[new][("M", e)][m].max(), Dr[old][("M", e)][m].max()) / MD)}

C2, C2r = {}, {}
for lab, items in COND.items():
    C2[lab], C2r[lab] = mean_lines(items, ["g_0", "condS_0"], [0.0], xq)
(oc, nc) = COND.keys()
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
for lab, ls in ((oc, "-"), (nc, "--")):
    axs[0].plot(xq, C2[lab][("g_0", 0.0)], color="#0f766e", ls=ls, lw=1.2, label=lab)
    axs[1].plot(xq, C2[lab][("condS_0", 0.0)], color="#c2410c", ls=ls, lw=1.2, label=lab)
axs[0].set_title("軸の凝縮質量分率 g_0"); axs[1].set_title("軸の過飽和度 S")
axs[1].set_yscale("log")
for ax in axs:
    ax.set_xlabel("x / r_t"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_xlim(40, 95)
fig.tight_layout(); fig.savefig(OUT / "fig_compare_cond_axis.png", dpi=150); plt.close(fig)
m = (xq >= 40) & (xq <= 95)
csum = {"g0_max_abs_diff": float(np.abs(C2[nc][("g_0", 0.0)] - C2[oc][("g_0", 0.0)])[m].max()),
        "g0_exit_old": float(C2[oc][("g_0", 0.0)][-1]), "g0_exit_new": float(C2[nc][("g_0", 0.0)][-1])}
(OUT / "compare_summary.json").write_text(json.dumps({"dry": dsum, "cond": csum, "windows": {"dry": DRY, "cond": COND}}, indent=1, ensure_ascii=False))
print(json.dumps({"dry": dsum, "cond": csum}, indent=1, ensure_ascii=False))
