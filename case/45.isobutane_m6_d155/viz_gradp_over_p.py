"""|∇p_s|/p_s [1/m] のコンタ (色範囲 0〜5, ユーザが以前使っていた見方)。V0 Euler・CFD ピン Euler・最終設計 NS (run_0051) を並べる。
usage: design/.venv-opt/bin/python viz_gradp_over_p.py [場のある case_dir] → _band_ab/gradp_over_p.png"""
import sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import load_field  # noqa: E402
D = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
RUNS = [("V0 (Hall 初期線, 当てはめ壁) — Euler", "run_0062_euler_wallfit_fit_r1_ext6k", "res_6000.h5"),
        ("CFD ピン (出口較正) — Euler", "run_0086_euler_wallfit_pincal_r1_ext6k", "res_6000.h5"),
        ("最終設計 run_0051 (補間壁, C2) — NS", "run_0051_ns_final_c2", None)]


def gradp_over_p(F):
    X, R, P = F["X"], F["R"], F["V"]["P"]
    xi, xj = np.gradient(X); ri, rj = np.gradient(R); J = xi * rj - xj * ri
    fi, fj = np.gradient(P)
    gx = (fi * rj - fj * ri) / J; gr = (-fi * xj + fj * xi) / J
    return np.hypot(gx, gr) / P / F["S"]          # X, R は r_t 単位 → 1/m


Fs = [(t, load_field(D / r, res)) for t, r, res in RUNS]
for xlim, name in (((-3, 96), "full"), ((30, 96), "test")):
    fig, ax = plt.subplots(len(Fs), 1, figsize=(15, 3.6 * len(Fs)), constrained_layout=True)
    for a, (t, F) in zip(ax, Fs):
        g = gradp_over_p(F); X, R = F["X"], F["R"]; m = (X[:, 0] >= xlim[0]) & (X[:, 0] <= xlim[1])
        pc = a.pcolormesh(X[m] * F["S"], R[m] * F["S"], g[m], cmap="turbo", vmin=0, vmax=5, shading="gouraud")
        a.set_title(f"|∇p_s|/p_s [1/m] (0〜5 で頭打ち) — {t}  (r_t {F['S']*1e3:.2f} mm, {F['res']})", fontsize=10)
        a.set_xlabel("x [m]"); a.set_ylabel("r [m]"); a.set_aspect("equal")
        fig.colorbar(pc, ax=a, location="right", shrink=0.9)
    fig.savefig(C / f"_band_ab/gradp_over_p_{name}.png", dpi=110); plt.close(fig)
print("ok")

# 試験部コア (x ∈ [41.8, 94.2] r_t, r/r_w ≤ 0.7) の |∇p|/p の最大・99 パーセンタイル [1/m] と、対数目盛の図
rows = []
fig, ax = plt.subplots(len(Fs), 1, figsize=(15, 3.6 * len(Fs)), constrained_layout=True)
for a, (t, F) in zip(ax, Fs):
    g = gradp_over_p(F); X, R = F["X"], F["R"]; eta = R / R[:, -1:]
    core = (X >= 41.82) & (X <= 94.23) & (eta <= 0.7)
    rows.append((t, float(g[core].max()), float(np.percentile(g[core], 99)), float(np.median(g[core]))))
    m = (X[:, 0] >= 30) & (X[:, 0] <= 96)
    pc = a.pcolormesh(X[m] * F["S"], R[m] * F["S"], np.log10(np.maximum(g[m], 1e-5)), cmap="turbo", vmin=-4, vmax=np.log10(5), shading="gouraud")
    a.set_title(f"log10(|∇p_s|/p_s [1/m]) (−4〜log10 5) — {t}", fontsize=10); a.set_xlabel("x [m]"); a.set_ylabel("r [m]"); a.set_aspect("equal")
    fig.colorbar(pc, ax=a, location="right", shrink=0.9)
fig.savefig(C / "_band_ab/gradp_over_p_test_log.png", dpi=110); plt.close(fig)
for r in rows:
    print(f"{r[0]}: core max {r[1]:.4f} /m, p99 {r[2]:.4f} /m, median {r[3]:.5f} /m")
