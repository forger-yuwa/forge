"""CFD ピン (出口較正, run_0086) と V0 (run_0062) の比較図: コンタ・軸 M・壁 r″。plan tooling-nozzle-cfd-pinned-initial-line §9。
usage: design/.venv-opt/bin/python viz_cfdpin_compare.py [場のある case_dir] → _band_ab/cfdpin_fig_*.png"""
import sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from cfd_initial_line import pinned_factory  # noqa: E402
from moc_wall_fit_ab import joint_fit  # noqa: E402
D = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
OUT = C / "_band_ab"
RUNS = {"V0 (Hall 初期線, 現行の当てはめ壁)": "run_0062_euler_wallfit_fit_r1_ext6k", "CFD ピン (出口較正)": "run_0086_euler_wallfit_pincal_r1_ext6k"}
F = {k: load_field(D / v, "res_6000.h5") for k, v in RUNS.items()}
names = list(F)

# ---- 1. コンタ (3 量 × 2 腕 + 差) ----
def contour_set(fname, xlim, rmax_frac=1.0, title="", throat=False):
    qs = [("M", lambda V: V["M"], "M", None), ("dp", lambda V: V["dpdx_nd"], "(r_t/p) ∂p/∂x", None), ("sch", lambda V: V["schlieren"], "log10(|∇ρ| r_t/ρ)", None)] if throat else [("M", lambda V: 100 * (V["M"] / 6.0 - 1), "100(M/6−1) [%]", (-0.5, 0.5)),
          ("dp", lambda V: V["dpdx_nd"], "(r_t/p) ∂p/∂x", (-0.02, 0.02)),
          ("sch", lambda V: V["schlieren"], "log10(|∇ρ| r_t/ρ)", None)]
    fig, ax = plt.subplots(len(qs), 3, figsize=(18, 3.2 * len(qs)), constrained_layout=True)
    for i, (key, fn, lab, lim) in enumerate(qs):
        vals = [fn(F[n]["V"]) for n in names]
        X, R = F[names[0]]["X"], F[names[0]]["R"]
        m = (X[:, 0] >= xlim[0]) & (X[:, 0] <= xlim[1])
        if lim is None:
            v = np.concatenate([q[m].ravel() for q in vals]); lim = tuple(np.percentile(v[np.isfinite(v)], [0.2, 99.8]))
        for j, n in enumerate(names):
            Xn, Rn = F[n]["X"], F[n]["R"]; mm = (Xn[:, 0] >= xlim[0]) & (Xn[:, 0] <= xlim[1])
            pc = ax[i, j].pcolormesh(Xn[mm], Rn[mm], np.clip(vals[j][mm], *lim), cmap="turbo", vmin=lim[0], vmax=lim[1], shading="gouraud")
            ax[i, j].set_title(f"{lab} — {n}", fontsize=9); ax[i, j].set_aspect("auto")
            fig.colorbar(pc, ax=ax[i, j], location="bottom", shrink=0.8, pad=0.02)
        # 差: x 断面・η 位置が違うので、ピン側を V0 の格子へ (x, η) 線形補間
        Xa, Ra = F[names[0]]["X"], F[names[0]]["R"]; Xb, Rb = F[names[1]]["X"], F[names[1]]["R"]
        eb = Rb / Rb[:, -1:]; ea = Ra / Ra[:, -1:]
        vb_on_a = np.array([np.interp(Xa[:, 0], Xb[:, 0], vals[1][:, j]) for j in range(Xa.shape[1])]).T
        diff = vb_on_a - vals[0]; dl = float(np.percentile(np.abs(diff[m]), 99.5)) or 1e-9
        pc = ax[i, 2].pcolormesh(Xa[m], Ra[m], np.clip(diff[m], -dl, dl), cmap="turbo", vmin=-dl, vmax=dl, shading="gouraud")
        ax[i, 2].set_title(f"差 (ピン − V0)、±{dl:.3g} で頭打ち", fontsize=9); fig.colorbar(pc, ax=ax[i, 2], location="bottom", shrink=0.8, pad=0.02)
    for a in ax.ravel():
        a.set_xlim(*xlim); a.set_xlabel("x / r_t"); a.set_ylabel("r / r_t")
    fig.suptitle(title, fontsize=11); fig.savefig(OUT / fname, dpi=110); plt.close(fig)

contour_set("cfdpin_fig_contour_full.png", (-3, 95.2), title="コンタ (全体): V0 vs CFD ピン (Euler、res_6000 = 通算 18000 step)")
contour_set("cfdpin_fig_contour_throat.png", (-1, 12), title="コンタ (スロート〜膨張部 x −1〜12、色範囲は 0.2〜99.8 % で切る)", throat=True)
contour_set("cfdpin_fig_contour_test.png", (38, 95.2), title="コンタ (試験部 x 38〜95)")

# ---- 2. 軸 M ----
xq = np.arange(0.0, 95.0, 0.02)
fig, ax = plt.subplots(3, 1, figsize=(12, 11), constrained_layout=True)
for n, c in zip(names, ("C0", "C3")):
    for eta, ls in ((0.0, "-"), (0.1, "--")):
        d = 100 * (eta_line(F[n], "M", eta, xq) / 6.0 - 1)
        ax[1].plot(xq, d, c, ls=ls, lw=0.9, label=f"{n}, r/r_w={eta}")
        ax[2].plot(xq, 100 * (eta_line(F[n], "P", eta, xq) / 2237.0 - 1), c, ls=ls, lw=0.9, label=f"{n}, r/r_w={eta}")
    ax[0].plot(xq, eta_line(F[n], "M", 0.0, xq), c, lw=1, label=n)
ax[0].set_ylabel("軸の M"); ax[0].set_title("軸中心マッハ数 (全体)")
ax[1].set_xlim(30, 95); ax[1].set_ylim(-0.08, 0.08); ax[1].set_ylabel("100(M/6−1) [%]"); ax[1].set_title("試験部の M のずれ (軸と r/r_w=0.1)")
ax[2].set_xlim(38, 95); ax[2].set_ylim(-0.5, 0.5); ax[2].set_ylabel("100(P/P_ref−1) [%]"); ax[2].set_title("試験部の静圧のずれ (P_ref 2237 Pa)")
for a in ax:
    a.grid(alpha=.3); a.set_xlabel("x / r_t")
for a in ax[1:]:
    a.axvspan(41.82, 94.23, color="gray", alpha=0.08)
ax[0].legend(fontsize=8); ax[1].legend(fontsize=7, ncol=2, loc="lower left"); ax[2].legend(fontsize=7, ncol=2, loc="lower left")
fig.savefig(OUT / "cfdpin_fig_axis.png", dpi=110); plt.close(fig)

# ---- 3. 壁 r″ (設計壁 = 当てはめ) ----
W = {}
for n, prob, pin in ((names[0], "problem_d155_euler_c2final_n2400.yaml", False), (names[1], "problem_d155_euler_c2final_n2400_pincal.yaml", True)):
    H = RA.HallThroat
    if pin:
        RA.HallThroat = pinned_factory(D / RUNS[names[0]], ["res_6000.h5"])
    try:
        d = RA.design_chain(RA.load_problem(C / prob))
    finally:
        RA.HallThroat = H
    s, _ = joint_fit(d["wall_inv"], d["R"], 1e-9); w = d["wall"]; w._spl = s; W[n] = w
panels = [((-12.6, 95.3), (-0.3, 0.6), "全体"), ((-0.4, 1.0), (0.3, 0.6), "スロート付近"), ((0.5, 8), (-0.02, 0.2), "x 0.5〜8"),
          ((8, 45), (-0.012, 0.0), "x 8〜45"), ((45, 95.3), (-0.003, 0.0003), "x 45〜95")]
fig, ax = plt.subplots(len(panels) + 1, 1, figsize=(11, 18), constrained_layout=True)
for a, ((xa, xb), yl, t) in zip(ax, panels):
    x = np.linspace(xa, xb, 20000)
    for n, c in zip(names, ("C0", "C3")):
        a.plot(x, W[n].r(x, 2), c, lw=1.0, label=n)
    a.set_xlim(xa, xb); a.set_ylim(*yl); a.grid(alpha=.3); a.set_title(f"壁の r″ — {t}"); a.set_ylabel("r″ [1/r_t]")
x = np.linspace(0, 95.2, 40000)
ax[-1].plot(x, W[names[1]].r(x) - W[names[0]].r(x), "k", lw=1); ax[-1].set_title("壁の位置の差 r(ピン) − r(V0) [r_t] (m* の一様な膨らみと出口較正を含む)")
ax[-1].grid(alpha=.3); ax[-1].set_xlabel("x / r_t")
ax[0].legend(fontsize=8)
fig.savefig(OUT / "cfdpin_fig_wall_r2.png", dpi=110); plt.close(fig)
print("ok")
