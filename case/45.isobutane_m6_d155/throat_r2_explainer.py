"""スロート近傍の壁の 2 階微分 r″ の山を説明する図と数値 (形状のみ, CFD 0 step)。
対象: 最終壁 problem_d155_ns_finemesh_recal_final.yaml (run_0117 と同じ設計チェーン) の設計壁と物理壁。
比較: 同じ MOC 点群・同じ始点条件 (r′(0)=0, r″(0)=1/R) で [0,1.5] の r″ を単調非増加に拘束した試し当てはめ C (生産には未使用) と、
      過去の対照 V4 (始点の r′・r″ 自由)。C は λ を 1e-10〜1e-8 で振って正則化依存でないことを確かめる。
usage: [CASE_RUNS=<run_* のある case dir>] python3 throat_r2_explainer.py [RUN_DIR (run_0117 の delta_r_initial.csv と nozzle.h5 を読む)] [OUT_DIR]
出力: OUT_DIR/throat_r2_explainer.json と fig1_overview / fig2_r2 / fig3_moc_curvature / fig4_variants / fig5_cost (*.png)
"""
import json, os, sys
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import BSpline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas, delta_r_from_table  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

RUN = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "run_0117_ns_recal_final_ext"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else C / "_band_ab/throat_r2_explainer"
OUT.mkdir(parents=True, exist_ok=True)
p = load_problem(C / "problem_d155_ns_finemesh_recal_final.yaml")
# 凍結源 run (initial_line_run) が別ツリーにあるときは CASE_RUNS で場所を渡す
p.geometry["initial_line_run"] = str((Path(os.environ.get("CASE_RUNS", C)) / p.geometry["initial_line_run"]).resolve())
d = design_chain(p)
W = d["wall"]; tb = d["wall_inv"]; R = float(d["R"]); rt_mm = float(p.spec["r_throat"]) * 1e3
spl = W._spl
x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]


def fit_variant(tb, start="V0", mono=None, lam=1e-9, k=5, sig_r=1e-6, sig_th=1e-4, h0=0.0125, h1=0.5, x_g=6.0):
    """joint_fit_wall と同じ目的関数・ノット・既定値で、始点条件と形状拘束だけを変えた当てはめ。
    start="V0": r(0)=r_0, r′(0)=0, r″(0)=1/R (= 生産, joint_fit_wall と同一) / "V4": r(0)=r_0 だけ (r′・r″ 自由, 過去の対照)。
    mono=(a, b): 台が [a, b] にかかる r‴ の B-spline 係数を ≤ 0 に拘束 (r″ が単調非増加の十分条件; 有効制約法の QP)。"""
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]; x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (x_g - x0), 1.0); xs.append(xs[-1] + h0 + (h1 - h0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * h1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1; E = np.eye(nc)
    D = lambda xq, dd: np.array([BSpline(t, E[i], k)(xq, dd) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1); w = np.gradient(x); w = w / w.mean()
    xq = np.linspace(x0, xe, 8000); B3 = D(xq, 3)
    A = (B0.T * (w / sig_r ** 2)) @ B0 + (B1.T * (w / sig_th ** 2)) @ B1 + lam * (B3.T * np.gradient(xq)) @ B3 / sig_r ** 2
    b = B0.T @ (w * r / sig_r ** 2) + B1.T @ (w * np.tan(th) / sig_th ** 2)
    sc = float(np.abs(np.diag(A)).max()); A = A / sc; b = b / sc
    rows = [D(np.r_[x0], 0)[0]]; dv = [r[0]]
    if start == "V0":
        rows += [D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0]]; dv += [x0 / R, 1.0 / R]
    rows += [D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]]; dv += [r[-1], np.tan(th[-1])]
    Ce = np.array(rows); de = np.array(dv); ne = len(de)
    if mono:
        # r‴ (2 次スプライン) の B-spline 係数 ≤ 0 は r‴ ≤ 0 の十分条件 (凸包性)。格子点で課すと制約が線形従属に近くなり有効制約法が循環する
        M3 = np.array([BSpline(t, E[i], k).derivative(3).c[:nc - 3] for i in range(nc)]).T
        t3 = t[3:len(t) - 3]
        G = M3[[j for j in range(nc - 3) if t3[j + 3] > mono[0] + 1e-12 and t3[j] < mono[1]]]
    else:
        G = np.zeros((0, nc))
    gs = np.abs(G).max(axis=1, keepdims=True) if len(G) else 1.0; G = G / gs
    act = []
    for it in range(3000):
        Cm = np.vstack([Ce, G[act]]) if act else Ce; dd = np.r_[de, np.zeros(len(act))]
        m = len(dd)
        sol = np.linalg.solve(np.block([[A, Cm.T], [Cm, np.zeros((m, m))]]), np.r_[b, dd])
        c, mu = sol[:nc], sol[nc + ne:]
        if len(act) and mu.min() < -1e-12 * max(1.0, np.abs(mu).max()):
            act.pop(int(np.argmin(mu))); continue
        if not len(G):
            break
        viol = G @ c
        if viol.max() <= 1e-10 * max(1.0, float(np.abs(G @ c).max())):
            break
        act.append(int(np.argmax(viol)))
    else:
        raise RuntimeError("有効制約法が 3000 回で収束しない")
    return BSpline(t, c, k), dict(n_active=len(act), iters=it + 1)


V0chk, _ = fit_variant(tb, "V0")
assert float(np.abs(V0chk(np.linspace(0, 3, 3001), 2) - W._spl(np.linspace(0, 3, 3001), 2)).max()) < 1e-6, "V0 が生産の joint 壁と一致しない"
MONO = (0.0, 1.5)
CV, c_info = fit_variant(tb, "V0", mono=MONO)
V4, _ = fit_variant(tb, "V4")
lam_sens = {}
for lam in (1e-10, 1e-9, 1e-8):
    s0, _ = fit_variant(tb, "V0", lam=lam); s1, _ = fit_variant(tb, "V0", mono=MONO, lam=lam)
    g = np.linspace(0, 0.3, 30001); nr = (x > 0) & (x < 0.3)
    lam_sens[f"{lam:g}"] = {nm: dict(r2_max=float(s(g, 2).max()), x_r2_max=float(g[np.argmax(s(g, 2))]),
                                     dth_max_lt0p3_deg=float(np.abs(np.degrees(np.arctan(s(x[nr], 1)) - th[nr])).max()),
                                     dr_max_lt0p3_um=float(np.abs(s(x[nr]) - r[nr]).max() * rt_mm * 1e3),
                                     dth_max_all_deg=float(np.abs(np.degrees(np.arctan(s(x[1:], 1)) - th[1:])).max()),
                                     dr_max_all_um=float(np.abs(s(x[1:]) - r[1:]).max() * rt_mm * 1e3))
                            for nm, s in (("V0", s0), ("C", s1))}
c0_deg = float(np.degrees(np.arctan((r[1] - r[0]) / (x[1] - x[0])) - 0.5 * (th[0] + th[1])))
# 物理壁 (run_0117 と同じ: 積分法初期 δ_r の平滑化済み表 + pw_ramp)
dr_tbl = np.loadtxt(RUN / "delta_r_initial.csv", delimiter=",", skiprows=1)
PW = PhysicalNozzleWall(W, tb, float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp,
                        offset="radial", delta_r_x=delta_r_from_table(dr_tbl[:, 0], dr_tbl[:, 1]),
                        ramp=tuple(float(v) for v in p.geometry["pw_ramp"]))
# 計算格子の壁節点 (run_0117 の nozzle.h5)
with h5py.File(RUN / "nozzle.h5") as f:
    xyz = np.asarray(f["MESH/COORD"]).reshape(-1, 3); wn = np.unique(np.asarray(f["BCONDS/3/vizBfaceNodes"]))
xw = np.sort(np.unique(np.round(xyz[wn, 0] / (float(p.spec["r_throat"])), 9)))

xg = np.linspace(-1.5, 3.0, 90001)
r2d = W.r(xg, 2); r3d = W.r(xg, 3)
mp = (xg >= 0) & (xg <= 0.3); i_max = np.argmax(np.where(mp, r2d, -1))
x_bump = float(xg[i_max]); r2_bump = float(r2d[i_max])
# 山の幅: r″ が 1/R を超えている区間
over = mp & (r2d > 1 / R); x_over = (float(xg[over].min()), float(xg[over].max())) if over.any() else (0, 0)
nodes_in_bump = int(((xw >= 0) & (xw <= 0.05)).sum())
dx_throat = float(np.diff(xw[(xw > -0.2) & (xw < 0.2)]).mean())
# MOC 点と当てはめの角度のずれ
near = (x > 0) & (x < 0.3)
dth_V0 = np.degrees(np.arctan(spl(x, 1)) - th); dth_C = np.degrees(np.arctan(CV(x, 1)) - th)
dr_V0 = spl(x) - r; dr_C = CV(x) - r
# 形状の差 (V0 − C)
xs = np.linspace(0, 1.0, 20001); dshape = (spl(xs) - CV(xs)) * rt_mm * 1e3   # µm
dslope = np.degrees(np.arctan(spl(xs, 1)) - np.arctan(CV(xs, 1)))
r2_C = CV(xs, 2)
# r″ の極値の数 [0,3]
xe3 = np.linspace(0, 3, 300001)


def n_extrema(s, tol=1e-6):
    """r″ の内部極値の数 = r‴ の符号反転の数 (|r‴| < tol の平坦区間は符号なしとして飛ばす)。"""
    g = s(xe3, 3); sg = np.sign(np.where(np.abs(g) < tol, 0.0, g)); sg = sg[sg != 0]
    return int(np.count_nonzero(np.diff(sg) != 0))


n_ext_V0 = n_extrema(spl); n_ext_C = n_extrema(CV)
# 物理壁の r″ (スロートは x_t にずれる)
r2p = PW.r(xg, 2)
ip = np.argmax(np.where(mp, r2p, -1))
out = dict(
    problem="problem_d155_ns_finemesh_recal_final.yaml", run=str(RUN.name), r_t_mm=rt_mm, R=R,
    design=dict(r2_at_0m=float(W.r(np.r_[-1e-9], 2)[0]), r2_at_0p=float(W.r(np.r_[1e-9], 2)[0]), r2_max=r2_bump, x_r2_max=x_bump,
                bump_height=r2_bump - 1 / R, x_r2_over_1overR=x_over, r3_max_0_0p05=float(r3d[(xg >= 0) & (xg < 0.05)].max()),
                r3_min_0_0p05=float(r3d[(xg >= 0) & (xg < 0.05)].min()), r2_extrema_0_3=n_ext_V0,
                r2_at=dict({f"{v}": float(W.r(np.r_[v], 2)[0]) for v in (-0.5, -0.25, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0)})),
    physical=dict(x_throat=float(PW.x_throat), r2_max_0_0p3=float(r2p[ip]), x_r2_max=float(xg[ip])),
    moc_points=dict(n_lt_0p3=int(near.sum()), x_first=[float(v) for v in x[:6]],
                    dth_V0_deg_lt0p3_absmax=float(np.abs(dth_V0[near]).max()), dth_C_deg_lt0p3_absmax=float(np.abs(dth_C[near]).max()),
                    dr_V0_lt0p3_absmax_um=float(np.abs(dr_V0[near]).max() * rt_mm * 1e3), dr_C_lt0p3_absmax_um=float(np.abs(dr_C[near]).max() * rt_mm * 1e3),
                    th_first_deg=[float(v) for v in np.degrees(th[:6])], th_parabola_deg=[float(v) for v in np.degrees(np.arctan(x[:6] / R))]),
    C_mono=dict(r2_at0=float(CV(np.r_[0.0], 2)[0]), r2_max_0_0p3=float(CV(np.linspace(0, .3, 3001), 2).max()), r2_extrema_0_3=n_ext_C,
                 shape_diff_um_absmax_0_1=float(np.abs(dshape).max()), slope_diff_deg_absmax_0_1=float(np.abs(dslope).max())),
    moc_segment_curvature_first4=[float(v) for v in (np.diff(np.tan(th)) / np.diff(x))[:4]],
    C_info=c_info, lambda_sensitivity=lam_sens, c0_first_interval_deg=c0_deg,
    V4=dict(r1_at0_deg=float(np.degrees(np.arctan(V4(np.r_[0.0], 1)[0]))), r2_at0=float(V4(np.r_[0.0], 2)[0])),
    mesh=dict(dx_wall_throat_rt=dx_throat, dx_wall_throat_mm=dx_throat * rt_mm, nodes_x_0_0p05=nodes_in_bump),
)
(OUT / "throat_r2_explainer.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))

INK, ACC, ALT, GR = "#1d2b3a", "#c2410c", "#0f766e", "#94a3b8"
# Fig1: 全体形状とスロート位置
fig, ax = plt.subplots(figsize=(10, 3.2))
xa = np.linspace(-12.5, float(x[-1]), 6000)
ax.plot(xa, W.r(xa), color=INK, lw=1.4, label="設計壁 (非粘性)"); ax.plot(xa, PW.r(xa), color=ACC, lw=1.0, ls="--", label="物理壁 (境界層ぶん外へ)")
ax.axvspan(-1.5, 3.0, color=ACC, alpha=.12); ax.annotate("この範囲を拡大 (x/r_t −1.5〜3)", (1.5, 1.2), (25, 2.2), color=ACC, arrowprops=dict(arrowstyle="->", color=ACC))
ax.set_xlabel("x / r_t  (r_t = %.2f mm、x=0 がスロート)" % rt_mm); ax.set_ylabel("r / r_t"); ax.grid(alpha=.3)
fig.legend(loc="lower center", ncol=2, frameon=False); fig.tight_layout(rect=(0, .08, 1, 1)); fig.savefig(OUT / "fig1_overview.png", dpi=160); plt.close(fig)

# Fig2: r″ (設計壁) — 全体と拡大
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1.4, 1]))
for ax, (a, b) in zip(axs, ((-1.5, 3.0), (-0.15, 0.25))):
    m = (xg >= a) & (xg <= b)
    up = m & (xg < 0); dn = m & (xg >= 0)
    ax.plot(xg[up], r2d[up], color=GR, lw=1.6, label="上流: 5 次 Hermite (収縮部)")
    ax.plot(xg[dn], r2d[dn], color=INK, lw=1.6, label="下流: MOC 点への当てはめ (5 次 B-spline)")
    ax.axhline(1 / R, color=ALT, lw=1, ls=":", label="1/R = 0.5 (スロート曲率の設計値)")
    ax.axvline(0, color=GR, lw=.6)
    ax.set_xlim(a, b); ax.set_xlabel("x / r_t"); ax.grid(alpha=.3)
axs[0].set_ylabel("r″  [1/r_t]"); axs[0].set_title("壁の 2 階微分 r″ (曲率にほぼ等しい)")
ax = axs[1]; ax.set_ylim(0.44, 0.52)
ax.plot([x_bump], [r2_bump], "o", color=ACC, ms=6)
ax.annotate(f"山 r″={r2_bump:.4f}\n(x={x_bump:.4f}, 1/R より +{r2_bump - 1/R:.4f})", (x_bump, r2_bump), (0.06, 0.508), color=ACC, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=ACC))
xm = xw[(xw >= -0.15) & (xw <= 0.25)]; ax.plot(xm, np.full_like(xm, 0.443), "|", color=ACC, ms=8, alpha=.7)
ax.text(-0.14, 0.4455, f"| = 計算格子の壁節点 (間隔 {dx_throat:.4f} r_t)", fontsize=8, color=ACC)
ax.set_title("拡大 (x −0.15〜0.25)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, .09, 1, 1)); fig.savefig(OUT / "fig2_r2.png", dpi=160); plt.close(fig)

# Fig3: MOC 点が区間ごとに要求する曲率 (隣り合う 2 点の壁角の差 / 間隔) と当てはめの r″
seg_x = 0.5 * (x[1:] + x[:-1]); seg_k = np.diff(np.tan(th)) / np.diff(x)
fig, ax = plt.subplots(figsize=(10, 4))
xx = np.linspace(0, 0.3, 3001)
for i in range(int((x < 0.3).sum())):
    ax.hlines(seg_k[i], x[i], x[i + 1], color=ACC, lw=2.4, label="MOC 点 2 つの壁角から決まる区間平均の曲率" if i == 0 else None)
ax.plot(xx, spl(xx, 2), color=INK, lw=1.5, label="現行の当てはめ r″ (始点 r″=1/R に固定)")
ax.axhline(1 / R, color=ALT, lw=1, ls=":", label="1/R = 0.5")
ax.plot(x[x < 0.3], np.full(int((x < 0.3).sum()), 0.405), "v", color=ACC, ms=6)
ax.text(0.003, 0.409, "▼ MOC 壁点の位置", color=ACC, fontsize=8)
ax.annotate(f"第 1 区間 {seg_k[0]:.3f}\n(1/R より +{seg_k[0] - 1/R:.3f})", (seg_x[0], seg_k[0]), (0.06, 0.515), color=ACC, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=ACC))
ax.set_ylim(0.40, 0.53); ax.set_xlim(-0.005, 0.3); ax.set_xlabel("x / r_t"); ax.set_ylabel("曲率 (≈ r″)"); ax.grid(alpha=.3)
ax.set_title("MOC の点は、スロート直後の最初の区間だけ 1/R より強い曲がりを要求している")
fig.legend(loc="lower center", ncol=3, frameon=False, fontsize=9); fig.tight_layout(rect=(0, .1, 1, 1)); fig.savefig(OUT / "fig3_moc_curvature.png", dpi=160); plt.close(fig)

# Fig4: 現行 V0 と試し C の r″・3 階微分 (λ 1e-9 = 生産値 と 1e-8)
fig, axs = plt.subplots(2, 1, figsize=(10, 6.6), sharex=True)
xx = np.linspace(0, 0.3, 6001)
for lam, ls in ((1e-9, "-"), (1e-8, "--")):
    s0, _ = fit_variant(tb, "V0", lam=lam); s1, _ = fit_variant(tb, "V0", mono=MONO, lam=lam)
    tag = "λ=1e-9 (生産値)" if lam == 1e-9 else "λ=1e-8 (平滑化 10 倍)"
    axs[0].plot(xx, s0(xx, 2), color=INK, ls=ls, lw=1.5, label=f"現行の当てはめ {tag}")
    axs[0].plot(xx, s1(xx, 2), color=ALT, ls=ls, lw=1.5, label=f"試し C (r″ 単調に拘束) {tag}")
    axs[1].plot(xx, s0(xx, 3), color=INK, ls=ls, lw=1.5); axs[1].plot(xx, s1(xx, 3), color=ALT, ls=ls, lw=1.5)
axs[0].axhline(1 / R, color=GR, lw=.8, ls=":"); axs[0].set_ylim(0.40, 0.52); axs[0].set_ylabel("r″")
axs[1].axhline(0, color=GR, lw=.6); axs[1].set_ylabel("d³r/dx³"); axs[1].set_xlabel("x / r_t")
for ax in axs:
    ax.grid(alpha=.3); ax.set_xlim(0, 0.3)
axs[0].set_title("r″ (上) と 3 階微分 (下): 山の高さは平滑化の強さ λ で変わる。試し C はどの λ でも山なし")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, .1, 1, 1)); fig.savefig(OUT / "fig4_variants.png", dpi=160); plt.close(fig)

# Fig5: 代償 — MOC 点からの壁角のずれ と 現行との形状差 (実寸)
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
mm = x < 0.3
ax = axs[0]
ax.plot(x[mm], dth_V0[mm], "o-", color=INK, ms=4, lw=1, label="現行")
ax.plot(x[mm], dth_C[mm], "s-", color=ALT, ms=4, lw=1, label="試し C")
ax.axhline(0, color=GR, lw=.6); ax.set_xlabel("x / r_t (MOC 壁点)"); ax.set_ylabel("当てはめの壁角 − MOC 点の壁角  [°]")
ax.set_title("MOC 点からの壁角のずれ"); ax.grid(alpha=.3); ax.legend(fontsize=9)
ax = axs[1]
ax.plot(xs, dshape, color=INK, lw=1.4); ax.axhline(0, color=GR, lw=.6)
ax.set_xlim(0, 0.5); ax.set_xlabel("x / r_t"); ax.set_ylabel("半径の差  [µm]"); ax.set_title(f"現行 − 試し C の半径差 (実寸, r_t = {rt_mm:.2f} mm)"); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(OUT / "fig5_cost.png", dpi=160); plt.close(fig)
print("figs ->", OUT)
