"""Euler 差分による排除厚さ抽出 (帯局所参照) の説明図。run_0022 の NS 場と run_0001 の Euler 場の実データで描く。

出力: report_axis_wave/band_fig*.png, band_numbers.json
usage: design/.venv-opt/bin/python viz_band_explain.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.metrics.deltastar import _load_structured, band_local_deficit  # noqa: E402

C = Path(__file__).resolve().parent
OUT = C / "report_axis_wave"
NS_RUN, E_RUN = "run_0022_ns_ib_pass0", "run_0001_euler_shortest_dry"
N = _load_structured(NS_RUN)
E = _load_structured(E_RUN)
xE = E["x"][:, 0]; rwE = E["r"][:, -1]; xs = N["x"][:, 0]
NUM = {}


def euler_at(x):
    k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
    w = float(np.clip((x - xE[k]) / (xE[k + 1] - xE[k]), 0, 1))
    f = lambda r: (1 - w) * np.interp(r, E["r"][k], E["q"][k]) + w * np.interp(r, E["r"][k + 1], E["q"][k + 1])
    return f, (1 - w) * rwE[k] + w * rwE[k + 1]


def column(xq):
    i = int(np.argmin(abs(xs - xq)))
    f, rw_e = euler_at(xs[i])
    r = N["r"][i]; q = N["q"][i]
    rr = np.linspace(0, r[-1], 20001)
    qn = np.interp(rr, r, q); qe = f(rr)
    return dict(x=float(xs[i]), r=r, q=q, f=f, rw_e=float(rw_e), rw=float(r[-1]), rr=rr, qn=qn, qe=qe,
                y=r[-1] - rr, d_in=float(r[-1] - rw_e))


def band_calc(c, yb):
    """band_local_deficit と同じ手順を図用に展開 (帯 [yb, 1.5 yb] の 1 次フィット → 壁へ延長 → 欠損 → δ_r)。"""
    y, qn, qe, rr = c["y"], c["qn"], c["qe"], c["rr"]
    mb = (y >= yb) & (y <= 1.5 * yb)
    ratio = qn / qe
    A = np.c_[np.ones(mb.sum()), y[mb] - yb]
    a, b = np.linalg.lstsq(A, ratio[mb], rcond=None)[0]
    ratio_fit = a + b * (y - yb)
    qref = ratio_fit * qe
    res = band_local_deficit(c["r"], c["q"], c["f"], c["rw_e"], y_b_fixed=yb)
    return dict(ratio_fit=ratio_fit, qref=qref, a=float(a), b=float(b), delta_r=res["delta_r"],
                D=res["mass_deficit"])


c21 = column(21.0)
c18 = column(18.0)
din = c21["d_in"]
yb_lo, yb_hi = 2.34 * din, 2.93 * din
B_lo = band_calc(c21, yb_lo); B_hi = band_calc(c21, yb_hi)
auto = band_local_deficit(c21["r"], c21["q"], c21["f"], c21["rw_e"])
NUM["x21"] = dict(x=c21["x"], delta_in=din, r_w=c21["rw"], yb_auto_over_din=auto["band_y_b"] / din,
                  delta_lo=B_lo["delta_r"], delta_hi=B_hi["delta_r"], ratio_lo_over_din=B_lo["delta_r"] / din,
                  ratio_hi_over_din=B_hi["delta_r"] / din, a_lo=B_lo["a"], a_hi=B_hi["a"])

# ---------------- Fig A: 断面の ρu (NS と Euler) ----------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
ax = axs[0]
ax.plot(c21["qe"], c21["rr"], color="k", lw=1.6, label="Euler (非粘性、設計壁の解)")
ax.plot(c21["qn"], c21["rr"], color="#d62728", lw=1.6, label="NS (粘性、物理壁 = 設計壁 + δ_in)")
ax.axhline(c21["rw_e"], color="k", ls=":", lw=0.9, label=f"設計壁 r = {c21['rw_e']:.3f} (Euler の壁)")
ax.axhline(c21["rw"], color="#d62728", ls=":", lw=0.9, label=f"物理壁 r = {c21['rw']:.3f} (= 設計壁 + δ_in)")
ax.set_xlabel("質量流束 ρu [kg/(m²·s)]"); ax.set_ylabel("r / r_t")
ax.set_title(f"(a) x = {c21['x']:.1f} r_t の断面全体 (軸 r=0 から壁まで)", fontsize=10); ax.set_ylim(0, 6.6)
ax.legend(fontsize=7.5, loc="center left", frameon=False)
ax = axs[1]
ym = 8 * din
m = c21["y"] <= ym
ax.plot(c21["qe"][m], c21["y"][m] / din, color="k", lw=1.6, label="Euler")
ax.plot(c21["qn"][m], c21["y"][m] / din, color="#d62728", lw=1.6, label="NS")
ax.fill_betweenx(c21["y"][m] / din, c21["qn"][m], c21["qe"][m], where=c21["qe"][m] > c21["qn"][m], color="#d62728", alpha=0.12)
ax.text(0.35, 0.18, "境界層で\n流れが足りない分", transform=ax.transAxes, fontsize=9, color="#a01818")
ax.set_ylim(0, 8); ax.set_xlabel("質量流束 ρu [kg/(m²·s)]"); ax.set_ylabel("壁からの距離 y / δ_in")
ax.set_title("(b) 壁の近く (y は NS の壁から測る。δ_in = 壁に足した排除厚)", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False)
fig.tight_layout(); fig.savefig(OUT / "band_figA_profiles.png", dpi=130); plt.close(fig)

# ---------------- Fig B: 比 NS/Euler と帯、延長線 ----------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
ratio = c21["qn"] / c21["qe"]
ax = axs[0]
ax.plot(ratio, c21["rr"], color="#1f5fa8", lw=1.4)
ax.axvline(1.0, color="0.6", lw=0.7)
ax.set_xlim(0.99, 1.01); ax.set_ylim(0, c21["rw"])
ax.set_xlabel("比 ρu_NS / ρu_Euler"); ax.set_ylabel("r / r_t")
ax.set_title("(a) 比を断面全体で見る: コアでも 1 ちょうどではない", fontsize=10)
core = c21["rr"] < 0.8 * c21["rw"]
NUM["x21"]["core_ratio_min"] = float(ratio[core].min()); NUM["x21"]["core_ratio_max"] = float(ratio[core].max())
ax = axs[1]
ym = 6 * din
m = c21["y"] <= ym
ax.plot(ratio[m], c21["y"][m] / din, color="#1f5fa8", lw=1.8, label="実際の比 NS/Euler")
for B, yb, col, lab in ((B_hi, yb_hi, "#2ca02c", "帯 A: y_b = 2.93 δ_in (x<20 で選ばれた)"),
                        (B_lo, yb_lo, "#ff7f0e", "帯 B: y_b = 2.34 δ_in (x≥20 で選ばれた)")):
    ax.axhspan(yb / din, 1.5 * yb / din, color=col, alpha=0.18)
    ax.plot(B["ratio_fit"][m], c21["y"][m] / din, color=col, lw=1.4, ls="--", label=lab + " の直線を壁まで延長")
ax.set_xlim(0.95, 1.01); ax.set_ylim(0, 6)
ax.set_xlabel("比 ρu_NS / ρu_Euler"); ax.set_ylabel("y / δ_in")
ax.set_title("(b) 壁の近く: 色の帯で比を直線フィット → 点線を壁まで延長", fontsize=10)
ax.legend(fontsize=7.8, loc="lower left", frameon=False)
fig.tight_layout(); fig.savefig(OUT / "band_figB_ratio.png", dpi=130); plt.close(fig)

# ---------------- Fig C: 基準 q_ref と欠損 (2 つの帯) ----------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
ym = 4 * din
m = c21["y"] <= ym
for ax, B, yb, col, lab in ((axs[0], B_hi, yb_hi, "#2ca02c", "帯 A (2.93 δ_in)"), (axs[1], B_lo, yb_lo, "#ff7f0e", "帯 B (2.34 δ_in)")):
    ax.plot(c21["qn"][m], c21["y"][m] / din, color="#d62728", lw=1.6, label="NS")
    ax.plot(B["qref"][m], c21["y"][m] / din, color=col, lw=1.6, ls="--", label="基準 = Euler × 延長した比")
    mz = c21["y"][m] <= yb
    ax.fill_betweenx(c21["y"][m][mz] / din, c21["qn"][m][mz], B["qref"][m][mz], color=col, alpha=0.25, label="欠損 (この面積から δ_r)")
    ax.axhspan(yb / din, 1.5 * yb / din, color=col, alpha=0.12)
    ax.set_ylim(0, 4); ax.set_xlabel("質量流束 ρu [kg/(m²·s)]")
    ax.set_title(f"{lab}: δ_r = {B['delta_r']:.5f} r_t (δ_in の {100 * B['delta_r'] / din:.1f} %)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left", frameon=False)
axs[0].set_ylabel("y / δ_in")
fig.tight_layout(); fig.savefig(OUT / "band_figC_deficit.png", dpi=130); plt.close(fig)

# ---------------- Fig D: x 方向 — 自動選択の帯と δ_r ----------------
xq = np.arange(10.0, 35.01, 0.25)
yb_auto, d_auto, d_fix = [], [], []
for xx in xq:
    c = column(xx)
    a = band_local_deficit(c["r"], c["q"], c["f"], c["rw_e"])
    yb_auto.append(a["band_y_b"] / c["d_in"]); d_auto.append(a["delta_r"] / c["d_in"])
    b = band_local_deficit(c["r"], c["q"], c["f"], c["rw_e"], y_b_fixed=3.66 * c["d_in"])
    d_fix.append(b["delta_r"] / c["d_in"])
fig, axs = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True)
axs[0].step(xq, yb_auto, where="mid", color="#1f5fa8", lw=1.6)
axs[0].axhline(2.93, color="#2ca02c", lw=0.8, ls="--"); axs[0].axhline(2.34, color="#ff7f0e", lw=0.8, ls="--")
axs[0].set_ylabel("自動で選ばれた y_b / δ_in"); axs[0].set_ylim(2.0, 3.3)
axs[0].set_title("(a) 断面ごとに自動で選ばれた帯の位置 (1.5 δ_in から 1.25 倍刻みで外へ)", fontsize=10)
axs[1].plot(xq, d_auto, color="#1f5fa8", lw=1.6, label="自動選択の帯で測った δ_r (実際に壁に使った)")
axs[1].plot(xq, d_fix, color="#2ca02c", lw=1.6, ls="--", label="帯を 3.66 δ_in に固定して測った δ_r")
axs[1].axhline(1.0, color="0.6", lw=0.6)
axs[1].set_ylabel("測った δ_r / δ_in"); axs[1].set_xlabel("x / r_t"); axs[1].set_ylim(0.975, 1.008)
axs[1].legend(fontsize=8.5, frameon=False, loc="lower right")
axs[1].set_title("(b) 測った排除厚さ (1 なら壁に足した値と同じ = 固定点)", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "band_figD_alongx.png", dpi=130); plt.close(fig)
NUM["along_x"] = {f"{x:.2f}": dict(yb=y, auto=a, fixed=f) for x, y, a, f in zip(xq, yb_auto, d_auto, d_fix)}

(OUT / "band_numbers.json").write_text(json.dumps(NUM, indent=1))
print(json.dumps(NUM["x21"], indent=1))
