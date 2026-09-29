import sys, os, numpy as np, h5py
C = "/home/sano/work/forge-cht/case/63.graetz_cht"
sys.path.insert(0, C); os.chdir(C)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("/home/sano/.fonts/NotoSansCJKjp-Regular.otf", "/home/sano/.fonts/NotoSansCJKjp-Bold.otf"):
    fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Noto Sans CJK JP", "font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
import graetz_common as gc
OUT = sys.argv[1]
COL = {1: "#1f77b4", 2: "#9467bd", 3: "#7f7f7f", 4: "#d62728", 5: "#8c564b", 6: "#2ca02c"}
NAME = {1: "inlet (1) 放物速度・エントロピー固定", 2: "outlet (2) 静圧 101325 Pa", 3: "wall_up (3) 断熱 no-slip",
        4: "wall_heat (4) 共役壁 / 等温壁", 5: "wall_down (5) 断熱 no-slip", 6: "axis (6) 軸"}
run = "run_0016_g2_dT10_r32"

# ---------- Figure 1: 領域
m = h5py.File(f"{run}/mesh.h5", "r"); xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
sol = h5py.File(f"{run}/solid.h5", "r"); sc = np.asarray(sol["MESH/COORD"][:], float); tri = np.asarray(sol["MESH/TRIS"][:])
def patch(pid):
    b = m[f"BCONDS/{pid}"]; nodes = np.asarray(b["vizBfaceNodes"][:]); sz = np.asarray(b["vizBfaceSizes"][:])
    segs = []; k = 0
    for s in sz:
        segs.append(nodes[k:k+s]); k += s
    return segs
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.4), gridspec_kw={"width_ratios": [2.3, 1]})
for pid in range(1, 7):
    for sg in patch(pid):
        p = xyz[sg]; a1.plot(p[:, 0]*1e3, p[:, 1]*1e3, color=COL[pid], lw=3 if pid != 6 else 2.5, solid_capstyle="butt")
a1.fill_between([0, gc.L_HEAT*1e3], gc.R*1e3, gc.R*1e3+0.5, color="#d62728", alpha=0.15, lw=0)
a1.text(gc.L_HEAT*1e3/2, 1.28, "固体殻 (k_s 100 W/mK、厚さ 0.5 mm、外面 Robin h 1e8・T_c)", ha="center", fontsize=8.5, color="#8b1a1a")
a1.annotate("", xy=(40, 0.45), xytext=(5, 0.45), arrowprops=dict(arrowstyle="->", lw=1.4, color="#333"))
a1.text(8, 0.55, "流れ (M 0.05、Re 1000)", fontsize=8.5)
a1.set_xlim(-14, 187); a1.set_ylim(-0.1, 1.75); a1.set_xlabel("x [mm] (加熱開始点 x=0)"); a1.set_ylabel("r [mm]")
a1.set_title("(a) 全体 (縦を約 40 倍に拡大して表示)", fontsize=9.5, loc="left")
a1.add_patch(plt.Rectangle((-0.3, 0), 0.9, 1.5, fill=False, ec="#333", lw=0.8, ls="--")); a1.text(1.5, 1.52, "(b) の範囲", fontsize=8)
# (b) 実寸の拡大
X = xyz[:, 0]*1e3; R = xyz[:, 1]*1e3
xs = np.unique(np.round(X, 7)); rs = np.unique(np.round(R, 7))
for xv in xs[(xs > -0.35) & (xs < 0.65)]: a2.plot([xv, xv], [0, 1], color="#bbb", lw=0.4)
for rv in rs: a2.plot([-0.3, 0.6], [rv, rv], color="#bbb", lw=0.4)
a2.triplot(sc[:, 0]*1e3, sc[:, 1]*1e3, tri, color="#e8a0a0", lw=0.35)
for pid in (1, 3, 4, 6):
    for sg in patch(pid):
        p = xyz[sg]; a2.plot(p[:, 0]*1e3, p[:, 1]*1e3, color=COL[pid], lw=2.5)
a2.set_xlim(-0.3, 0.6); a2.set_ylim(0, 1.5); a2.set_aspect("equal"); a2.set_xlabel("x [mm]"); a2.set_ylabel("r [mm]")
a2.set_title("(b) x=0 付近 (実寸、格子 N_r=32)", fontsize=9.5, loc="left")
h = [plt.Line2D([], [], color=COL[p], lw=3, label=NAME[p]) for p in range(1, 7)]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.2), fontsize=8.5)
fig.savefig(f"{OUT}/fig1_domain.png"); plt.close(fig)

# ---------- Figure 2: 課した境界条件
fig, (b1, b2) = plt.subplots(1, 2, figsize=(10.5, 3.2))
prof = np.genfromtxt(f"{run}/inlet_profile_1.csv", names=True)
b1.plot(prof["Ux"], prof["y"]*1e3, color="#1f77b4", lw=2, label="指定 (inlet_profile_1.csv)")
V = h5py.File(f"{run}/res_600000.h5", "r")["VALUE"]; Ux = np.asarray(V["Ux"][:]); T = np.asarray(V["T"][:])
i = np.where(np.abs(X - X.min()) < 1e-6)[0]; o = np.argsort(R[i])
b1.plot(Ux[i][o], R[i][o], "o", ms=3, color="#333", label="実現 (入口面の節点、step 600000)")
b1.set_xlabel("u_x [m/s]"); b1.set_ylabel("r [mm]"); b1.set_title("(a) 入口で与えた速度分布 (十分発達した放物線)", fontsize=9.5, loc="left"); b1.legend(frameon=False, fontsize=8)
PE_M = 728.04
for rn, lab, c, tc in (("run_0016_g2_dT10_r32", "加熱 run (ΔT = 10 K、T_c = 310 K)", "#d62728", 310.0), ("run_0014_g2_dT0_r32", "対照 run (ΔT = 0、T_c = 300 K)", "#555", 300.0)):
    w = h5py.File(f"{rn}/res_wall_heat_4_600000.h5", "r"); cx = np.asarray(w["MESH/COORD"][:]).reshape(-1, 3)[:, 0]
    tw = np.asarray(w["VALUE/iface_Tw_bc"][:]); o = np.argsort(cx); cx, tw = cx[o], tw[o]; k = cx > 0
    b2.plot(cx[k] / (gc.D * PE_M), 1e3*(tw[k]-tc), color=c, lw=1.6, label=lab)
b2.set_xscale("log"); b2.set_xlabel("加熱開始点からの距離 x$^+$ = x/(D Pe) (対数軸)"); b2.set_ylabel("壁温 − 固体外面温度 T_c [mK]")
b2.set_title("(b) 管の内壁の温度は T_c からどれだけずれたか", fontsize=9.5, loc="left"); b2.legend(frameon=False, fontsize=8, loc="lower right")
fig.savefig(f"{OUT}/fig2_bc.png"); plt.close(fig)

# ---------- Figure 3: Nu
import graetz_ref
xs_, nu_ = graetz_ref.solve_march(400, 0.2, 3000)
fig, (c1, c2) = plt.subplots(2, 1, figsize=(7.6, 6.0), sharex=True, gridspec_kw={"height_ratios": [1.3, 1]})
c1.plot(xs_, nu_, color="#111", lw=1.4, label="古典 Graetz (基準、graetz_ref.py march)")
cols = {"16": "#f39c12", "32": "#2e86c1", "64": "#8e44ad"}
runs3 = {"16": "run_0020_g2_dT10_r16", "32": "run_0016_g2_dT10_r32", "64": "run_0024_g2_dT10_r64_ext"}
for k, rn in runs3.items():
    d = np.genfromtxt(f"{rn}/graetz_nu_600000.csv", delimiter=",", names=True); ok = d["xplus"] > 0
    c1.plot(d["xplus"][ok], d["Nu"][ok], "o", ms=2.2, color=cols[k], label=f"forge N_r={k} ({rn})")
    c2.plot(d["xplus"][ok], 100*d["err"][ok], "-", lw=1.1, color=cols[k], label=f"N_r={k}")
for ax in (c1, c2):
    ax.axvspan(3e-3, 0.1, color="#2e86c1", alpha=0.06, lw=0); ax.set_xscale("log")
c1.set_ylim(3, 16); c1.set_ylabel("局所 Nu"); c1.legend(frameon=False, fontsize=8); c1.set_title("窓 x$^+$ 3e-3〜0.1 (薄青) で判定、対照 (ΔT 0) を差し引いた Nu", fontsize=9.5, loc="left")
c2.axhline(0, color="#111", lw=0.8); c2.axhspan(-2, 2, color="#27ae60", alpha=0.05, lw=0)
for y in (-2, 2): c2.axhline(y, color="#27ae60", lw=0.8, ls="--")
c2.set_ylim(-2.5, 1.0); c2.set_xlim(1e-4, 0.13); c2.set_ylabel("Nu/Nu_ref − 1 [%]"); c2.set_xlabel("x$^+$ = x/(D Pe)、実測 Pe で算出 (対数軸)")
c2.text(1.2e-4, -2.25, "破線: 許容 ±2 % (窓全体)", fontsize=8, color="#1e8449"); c2.legend(frameon=False, fontsize=8, ncol=3, loc="lower right")
fig.savefig(f"{OUT}/fig3_nu.png"); plt.close(fig)

# ---------- Figure 4: 運動量収支
fig, ax = plt.subplots(figsize=(6.8, 3.2))
lab = []; acc = []; wsh = []; res = []; ramp = []; Adiff = []
for rn, t in (("run_0014_g2_dT0_r32", "対照 ΔT = 0"), ("run_0015_g2_dT5_r32", "加熱 ΔT = 5 K"), ("run_0016_g2_dT10_r32", "加熱 ΔT = 10 K")):
    S = np.genfromtxt(f"{rn}/momentum_balance_series.csv", delimiter=",", names=True); tail = S[len(S)//2:]   # 末尾窓 (後半) の平均
    dP = tail["dP"].mean(); lab.append(t)
    acc.append(100*tail["dM"].mean()/dP); wsh.append(100*(-tail["Fw"]-tail["A_pois"]).mean()/dP)
    res.append(100*tail["resid_o2"].mean()/dP); ramp.append(100*np.ptp(tail["resid_o2"])/2/dP); Adiff.append(100*(tail["dP"]-tail["A_pois"]).mean()/dP)
y = np.arange(3)
ax.barh(y, acc, color="#2e86c1", label="運動量流束の増分 (加熱・膨張による加速)")
ax.barh(y, wsh, left=acc, color="#e67e22", label="壁せん断の Poiseuille からの超過")
ax.plot(Adiff, y, "k|", ms=18, mew=2, label="測った差 (区間の圧力力 − Poiseuille の予測)")
for yy, r_, m_, a_ in zip(y, res, ramp, Adiff): ax.text(a_ + 0.05, yy - 0.05, f"収支の不一致 {r_:+.4f} ± {m_:.1e} %", fontsize=7.5, color="#333", va="center")
ax.set_yticks(y, lab); ax.set_xlabel("区間の圧力力に対する割合 [%] (区間 x 43〜130 mm、N_r=32、末尾窓の平均)"); ax.set_xlim(0, 3.1)
ax.legend(frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.45, -0.28), ncol=1)
fig.savefig(f"{OUT}/fig4_momentum.png"); plt.close(fig)

# ---------- Figure 5: 温度再現
from temp_reproduce import Grid, source, solve_T
rn = "run_0014_g2_dT0_r32"; G = Grid(xyz)
V = h5py.File(f"{rn}/res_600000.h5", "r")["VALUE"]; ro, u, v, p, T = (np.asarray(V[k][:], float) for k in ("ro", "Ux", "Uy", "P", "T"))
phi, work = source(G, ro, u, v, p); Ts = T[G.idx]; heat = (G.xs >= -1e-7) & (G.xs <= gc.L_HEAT + 1e-7)
wall = np.where(heat, Ts[:, -1], np.nan)
TA = solve_T(G, ro, u, v, np.zeros_like(phi), Ts[0], wall); TB = solve_T(G, ro, u, v, phi + work, Ts[0], wall)
i0 = int(np.argmin(np.abs(G.xs))); rr = G.r*1e3
fig, ax = plt.subplots(figsize=(6.0, 3.4))
ax.plot(Ts[i0] - 300, rr, "o", ms=3.5, color="#111", label="forge の保存場 (run_0014、対照 ΔT 0)")
ax.plot(TA[i0] - 300, rr, "-", color="#7f8c8d", lw=1.6, label="A: 散逸・圧力仕事なし")
ax.plot(TB[i0] - 300, rr, "--", color="#c0392b", lw=1.6, label="B: 散逸・圧力仕事あり")
ax.set_xlabel("T − 300 K [K] (x=0 断面)"); ax.set_ylabel("r [mm]"); ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("正式判定: 判定不能 (x=0 の後処理不確かさ 2.8e-3 K > 上限 1.7e-3 K)", fontsize=9, loc="left", color="#8a4b00")
fig.savefig(f"{OUT}/fig5_temp.png"); plt.close(fig)
print("ok", np.ptp(Ts[i0, :-1]), np.ptp(TA[i0, :-1]), np.ptp(TB[i0, :-1]))

# ---------- Figure: 問題の模式図
fig, ax = plt.subplots(figsize=(10.5, 3.3))
ax.set_xlim(-2.2, 11.4); ax.set_ylim(-1.35, 2.0); ax.axis("off")
ax.fill_between([0, 9], 1.0, 1.35, color="#f4c7c3"); ax.text(4.5, 1.47, "固体殻の外面を T_c に保つ (T_c = T_in + ΔT)", ha="center", fontsize=9, color="#8b1a1a")
ax.plot([-2, 10.2], [1, 1], color="#555", lw=2); ax.plot([-2, 10.2], [-1, -1], color="#555", lw=2)
ax.plot([0, 9], [1, 1], color="#d62728", lw=3); ax.plot([0, 9], [-1, -1], color="#d62728", lw=3)
ax.plot([-2, 10.2], [0, 0], color="#999", lw=0.8, ls="-.")
rr = np.linspace(-1, 1, 41)
for x0 in (-1.6,):
    ax.plot(x0 + 0.9*(1-rr**2), rr, color="#1f77b4", lw=1.5)
    for r0 in (-0.6, 0, 0.6): ax.annotate("", xy=(x0 + 0.9*(1-r0**2), r0), xytext=(x0, r0), arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1))
ax.text(-2.1, 1.47, "入口: 放物速度 (十分発達)\n一様温度 T_in", fontsize=8.5, color="#1f77b4", va="bottom")
# 熱境界層の発達
xx = np.linspace(0, 9, 200); dl = 1 - np.minimum(1, 0.5*np.cbrt(xx/9*8))
ax.fill_between(xx, dl, 1, color="#d62728", alpha=0.12); ax.fill_between(xx, -1, -dl, color="#d62728", alpha=0.12)
ax.text(2.2, 0.62, "温められた層 (下流ほど厚くなる)", fontsize=8.5, color="#8b1a1a")
for xq in (0.6, 2.5, 5.5):
    ax.annotate("", xy=(xq, 0.82), xytext=(xq, 1.18), arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.2))
ax.text(0.75, 1.08, "q_w", fontsize=9, color="#c0392b")
ax.text(0, -1.25, "x = 0 (加熱開始)", fontsize=8.5, ha="center"); ax.text(-1.5, -1.25, "断熱壁", fontsize=8.5, ha="center", color="#555"); ax.text(9, -1.25, "x$^+$ = 0.12", fontsize=8.5, ha="center")
ax.text(9.35, 0.25, "混合平均温度 T_b(x) は\n下流ほど T_w に近づく", fontsize=8.5)
fig.savefig(f"{OUT}/figS_problem.png"); plt.close(fig)

# ---------- Figure: 差し引きの中身
h = np.genfromtxt("run_0016_g2_dT10_r32/graetz_nu_600000.csv", delimiter=",", names=True); k = (h["xplus"] > 0) & (h["xplus"] < 0.115)
fig, (d1, d2) = plt.subplots(1, 2, figsize=(10.5, 3.4))
d1.plot(h["xplus"][k], -h["q"][k], color="#d62728", lw=1.6, label="加熱 run (ΔT = 10 K)")
d1.plot(h["xplus"][k], -h["q0"][k], color="#555", lw=1.6, label="対照 run (ΔT = 0)")
d1.plot(h["xplus"][k], -(h["q"]-h["q0"])[k], "--", color="#2e86c1", lw=1.4, label="差 (加熱 − 対照)")
d1.set_xscale("log"); d1.set_yscale("symlog", linthresh=10); d1.set_xlabel("x$^+$"); d1.set_ylabel("壁からガスへの熱流束 [W/m²]")
d1.set_title("(a) 壁熱流束: 対照 run でも負の値 (ガス → 壁) が出る", fontsize=9.5, loc="left"); d1.legend(frameon=False, fontsize=8)
d2.plot(h["xplus"][k], (h["Tw"]-h["Tb"])[k], color="#d62728", lw=1.6, label="加熱 run")
d2.plot(h["xplus"][k], (h["Tw0"]-h["Tb0"])[k], color="#555", lw=1.6, label="対照 run")
d2.plot(h["xplus"][k], ((h["Tw"]-h["Tb"])-(h["Tw0"]-h["Tb0"]))[k], "--", color="#2e86c1", lw=1.4, label="差 (加熱 − 対照)")
d2.set_xscale("log"); d2.set_xlabel("x$^+$"); d2.set_ylabel("壁温 − 混合平均温度 [K]")
d2.set_title("(b) 壁温と混合平均温度の差: 対照 run でも最大 0.17 K", fontsize=9.5, loc="left"); d2.legend(frameon=False, fontsize=8)
fig.savefig(f"{OUT}/figD_subtract.png"); plt.close(fig)
