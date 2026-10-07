"""case/45 の生産の壁 (2026-10-07 採用: 上流 poly + MOC analytic/converge + 全域 1 本の B-spline) と旧生産の壁 (単調壁・ランプ・
legacy MOC、run_0147) を、スロート周りの r・r′・r″・r‴ で比べる図 (CFD 0 step)。

  新: _band_ab/prod_confirm/prep/wall_repr.json (生産の問題の準備、run_0167 の入力とビット同一) の物理壁 (1 本の B-spline) と設計壁
  旧: problem_d155_ns_finemesh_recal_final_mono.yaml を MOC のキー legacy・fixed2 を明示して作り直した設計壁と物理壁
      (作り直した物理壁と run_0147 の wall_physical.csv の差を記録する)

x と r は各壁の r_t で無次元化 (新 r_t 76.6654 mm・旧 76.6539 mm)。出力 _band_ab/prod_confirm/viz/*.png と viz_summary.json。
usage: design/.venv-opt/bin/python viz_prod_throat_r2.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
RUNS = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
OUT = HERE / "_band_ab" / "prod_confirm" / "viz"
OLD_PROB = "problem_d155_ns_finemesh_recal_final_mono.yaml"


def old_walls():
    from forge_design.evaluate.runner_axismach import build_physical_wall, design_chain, integral_delta_r, load_problem
    txt = (HERE / OLD_PROB).read_text().replace("initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k",
                                                f"initial_line_run: {RUNS / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    tmp = OUT / "old_problem.yaml"
    tmp.write_text(txt)
    p = load_problem(tmp)
    p.geometry["moc_axis_limit"], p.geometry["moc_corrector"] = "legacy", "fixed2"   # 2026-10-07 から既定が analytic・converge
    d = design_chain(p)
    _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    PW = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    return float(p.spec["r_throat"]), d["wall"], PW


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fp = Path.home() / ".fonts" / "NotoSansCJKjp-Regular.otf"
    if fp.is_file():
        font_manager.fontManager.addfont(str(fp))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    from forge_design.geometry.wall_axismach import load_wall_file
    OUT.mkdir(parents=True, exist_ok=True)
    W = load_wall_file(HERE / "_band_ab" / "prod_confirm" / "prep")
    rt_new = float(W["scale_m"])
    newP, newD = W["physical"], W["design"]
    rt_old, oldD, oldP = old_walls()
    summ = {"rt_new_m": rt_new, "rt_old_m": rt_old, "new_domain": list(map(float, W["domain"])), "new_throat": W["throat"],
            "old_throat": {"x": float(oldP.x_throat), "r": float(oldP.r_throat), "kappa": float(oldP.kappa_throat)}}
    # 旧の作り直しと run_0147 の保存物の一致
    csv = RUNS / "run_0147_ns_mono_final" / "wall_physical.csv"
    if csv.is_file():
        a = np.loadtxt(csv, delimiter=",", skiprows=1)
        summ["old_rebuild_vs_run_0147_wall_physical_max_abs_m"] = float(np.abs(oldP.r(a[:, 0] / rt_old) * rt_old - a[:, 1]).max())

    def ev(w, x, n):
        return np.asarray(w.r(np.asarray(x, float), n), float)

    C_NEW, C_OLD = "#0b6e4f", "#c0392b"
    # ---- 図 1: スロート周り (x ∈ [−2, 3]) の r・r′・r″・r‴ ----
    x = np.linspace(-2.0, 3.0, 4001)
    fig, ax = plt.subplots(4, 1, figsize=(9.5, 12.5), sharex=True, constrained_layout=True)
    labs = ["r / r_t", "dr/dx", "d²r/dx² [1/r_t]", "d³r/dx³ [1/r_t²]"]
    for n in range(4):
        a = ax[n]
        a.plot(x, ev(newP, x, n), color=C_NEW, lw=1.8, label="新 物理壁 (1 本の B-spline)")
        a.plot(x, ev(oldP, x, n), color=C_OLD, lw=1.4, label="旧 物理壁 (ランプ・legacy MOC)")
        a.plot(x, ev(newD, x, n), color=C_NEW, lw=1.0, ls="--", label="新 設計壁 (analytic + converge)")
        a.plot(x, ev(oldD, x, n), color=C_OLD, lw=1.0, ls="--", label="旧 設計壁 (legacy + fixed2)")
        a.axvline(0.0, color="0.6", lw=0.8)
        a.set_ylabel(labs[n])
        a.grid(alpha=0.3)
    ax[0].set_ylim(0.98, 1.25)
    ax[-1].set_xlabel("x / r_t (0 = 設計スロート)")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(OUT / "fig1_throat_r_derivs.png", dpi=150)
    plt.close(fig)
    # ---- 図 2: スロート直後の拡大 (x ∈ [−0.1, 0.45]): d²r/dx² と d³r/dx³ ----
    x2 = np.linspace(-0.1, 0.45, 5501)
    fig, ax = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True, constrained_layout=True)
    for n, a in ((2, ax[0]), (3, ax[1])):
        for w, c, ls, lab in ((newP, C_NEW, "-", "新 物理壁 (1 本の B-spline)"), (oldP, C_OLD, "-", "旧 物理壁 (ランプ・legacy MOC)"),
                              (newD, C_NEW, "--", "新 設計壁 (analytic + converge)"), (oldD, C_OLD, "--", "旧 設計壁 (legacy + fixed2)")):
            a.plot(x2, ev(w, x2, n), color=c, ls=ls, lw=1.8 if ls == "-" else 1.0, label=lab)
        a.axvline(0.0, color="0.6", lw=0.8)
        a.grid(alpha=0.3)
    ax[0].set_ylabel("d²r/dx² [1/r_t]"); ax[1].set_ylabel("d³r/dx³ [1/r_t²]"); ax[1].set_xlabel("x / r_t (0 = 設計スロート)")
    ax[0].set_ylim(0.40, 0.51)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(OUT / "fig2_throat_r2_zoom.png", dpi=150)
    plt.close(fig)
    # ---- 図 3: 入口から x = 12 までの r″ (上流の多項式 Q と旧ランプ) ----
    x3 = np.linspace(max(W["domain"][0], -12.5), 12.0, 8001)
    fig, ax = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True, constrained_layout=True)
    for w, c, ls, lab in ((newP, C_NEW, "-", "新 物理壁"), (oldP, C_OLD, "-", "旧 物理壁")):
        ax[0].plot(x3, ev(w, x3, 2), color=c, ls=ls, lw=1.4, label=lab)
    ax[0].axvspan(-11.0, -6.0, color="#f6e2df", zorder=0, label="旧 δ_r ランプ [−11, −6]")
    ax[0].set_ylabel("d²r/dx² [1/r_t]"); ax[0].grid(alpha=0.3); ax[0].legend(loc="upper left", frameon=False)
    dr_um = (ev(newP, x3, 0) * rt_new - ev(oldP, x3, 0) * rt_old) * 1e6
    ax[1].plot(x3, dr_um, color="#34495e", lw=1.4)
    ax[1].set_ylabel("新 − 旧 の物理壁の半径 [µm]")
    ax[1].axhline((rt_new - rt_old) * 1e6, color="0.5", lw=0.8, ls=":"); ax[1].set_xlabel("x / r_t"); ax[1].grid(alpha=0.3)
    fig.savefig(OUT / "fig3_upstream_r2_and_dr.png", dpi=150)
    plt.close(fig)
    # 数値の要約
    xs = np.linspace(0.0, 1.5, 30001)
    for lab, w in (("new_physical", newP), ("old_physical", oldP), ("new_design", newD), ("old_design", oldD)):
        r2 = ev(w, xs, 2)
        inc = float(np.max(np.maximum.accumulate(r2[::-1])[::-1] - r2)) if False else float(np.max(np.diff(r2).clip(min=0).cumsum()) if len(r2) else 0)
        summ[lab] = {"r2_at_0": float(ev(w, [0.0], 2)[0]), "r2_max_0_1p5": float(r2.max()), "x_r2_max": float(xs[np.argmax(r2)]),
                     "r2_total_rise_0_1p5": inc, "r3_max_0_1p5": float(ev(w, xs, 3).max())}
    xu = np.linspace(max(W["domain"][0], -12.5), -0.01, 20001)
    summ["upstream_r2_maxabs"] = {"new_physical": float(np.abs(ev(newP, xu, 2)).max()), "old_physical": float(np.abs(ev(oldP, xu, 2)).max())}
    summ["dr_um_range"] = [float(dr_um.min()), float(dr_um.max())]
    summ["rt_diff_um"] = (rt_new - rt_old) * 1e6
    xp = np.linspace(0.0, 0.3, 30001)
    for lab, w in (("new_physical", newP), ("old_physical", oldP)):
        r2 = ev(w, xp, 2)
        flat = xp[r2 >= r2.max() - 1e-4]
        summ[lab]["r2_plateau_within_1e-4_x"] = [float(flat.min()), float(flat.max())]
        summ[lab]["r3_min_0_0p3"] = float(ev(w, xp, 3).min())
        summ[lab]["x_r3_min"] = float(xp[np.argmin(ev(w, xp, 3))])
    (OUT / "viz_summary.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
    print(json.dumps(summ, indent=1, ensure_ascii=False, default=float))


if __name__ == "__main__" and (len(sys.argv) == 1 or sys.argv[1] != "cad"):
    main()


def fig_cad_overview():
    """生産の壁の全体像 (半断面の寸法つきと、3/4 を切り欠いた 3D の内面)。出力 _band_ab/prod_confirm/viz/fig4_cad_overview.png。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: F401
    fp = Path.home() / ".fonts" / "NotoSansCJKjp-Regular.otf"
    if fp.is_file():
        font_manager.fontManager.addfont(str(fp))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    from forge_design.geometry.wall_axismach import load_wall_file
    W = load_wall_file(HERE / "_band_ab" / "prod_confirm" / "prep")
    s = float(W["scale_m"]) * 1000.0
    x0, x1 = (float(v) for v in W["domain"])
    x = np.linspace(x0, x1, 3000)
    r = np.asarray(W["physical"].r(x)) * s
    xm = x * s
    th = W["throat"]
    # (a) 半断面 (縦横同じ縮尺)
    fig, a = plt.subplots(figsize=(12, 3.4), constrained_layout=True)
    a.fill_between(xm, r, r.max() * 1.1, color="#d8dfda")
    a.plot(xm, r, color="#0b6e4f", lw=1.6)
    a.axhline(0, color="0.4", lw=0.8, ls="-.")
    a.set_aspect("equal")
    a.set_xlabel("x [mm] (原点 = 設計スロート)"); a.set_ylabel("r [mm]")
    a.annotate(f"入口 r = {r[0]:.2f} mm", (xm[0], r[0]), (xm[0] + 300, r[0] + 120), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=10)
    a.annotate(f"スロート r = {th['r'] * s:.3f} mm (x = {th['x'] * s:.2f} mm)", (th["x"] * s, th["r"] * s), (350, 260),
               arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=10)
    a.annotate(f"出口 r = {r[-1]:.2f} mm", (xm[-1], r[-1]), (xm[-1] - 2000, r[-1] - 230), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=10)
    a.annotate("", (xm[0], -90), (xm[-1], -90), arrowprops=dict(arrowstyle="<->", lw=0.8))
    a.text(0.5 * (xm[0] + xm[-1]), -200, f"全長 {xm[-1] - xm[0]:.1f} mm (入口 → 出口)", ha="center", va="center", fontsize=10)
    a.set_ylim(-260, r.max() * 1.1)
    out = HERE / "_band_ab" / "prod_confirm" / "viz" / "fig4_section.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    # (b) 内面 (見る側の 1/4 を切り欠いた回転面)
    fig = plt.figure(figsize=(11, 5.2), constrained_layout=True)
    b = fig.add_subplot(projection="3d")
    xs = np.linspace(x0, x1, 260)
    rs = np.asarray(W["physical"].r(xs)) * s
    ph = np.linspace(np.pi, 2.5 * np.pi, 150)             # y < 0 かつ z > 0 の 1/4 (見る側) を除く
    X, PH = np.meshgrid(xs * s, ph, indexing="ij")
    R = np.repeat(rs[:, None], len(ph), axis=1)
    b.plot_surface(X, R * np.cos(PH), R * np.sin(PH), color="#7fb8a3", alpha=0.95, linewidth=0, antialiased=True, shade=True)
    b.set_box_aspect((xm[-1] - xm[0], 2 * r.max(), 2 * r.max()))
    b.set_xlabel("x [mm]", labelpad=10)
    b.set_yticks([-700, 0, 700]); b.set_zticks([-700, 0, 700])
    b.tick_params(axis="y", labelsize=8); b.tick_params(axis="z", labelsize=8)
    b.view_init(elev=24, azim=-58)
    out2 = HERE / "_band_ab" / "prod_confirm" / "viz" / "fig5_inner_surface.png"
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    print(out, out2)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "cad":
    fig_cad_overview()
