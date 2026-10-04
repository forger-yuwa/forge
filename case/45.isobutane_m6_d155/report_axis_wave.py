"""軸 M の山の報告用の図と数値 (plan verification-m6-axis-wave-mesh-su2)。

出力: case/45.isobutane_m6_d155/report_axis_wave/fig*.png と numbers.json
usage: design/.venv-opt/bin/python report_axis_wave.py
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402

fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

C = Path(__file__).resolve().parent
OUT = C / "report_axis_wave"
OUT.mkdir(exist_ok=True)
NUM = {}

# 壁の系列: (run, problem, δ 入力 CSV, 列, ラベル, 色)
SERIES = [
    ("run_0022_ns_ib_pass0", "problem_d155_ns.yaml", "run_0022_ns_ib_pass0/delta_r_initial.csv", 1,
     "pass 0: 積分法 (CONTUR)", "#1f77b4"),
    ("run_0023_ns_ib_pass1", "problem_d155_ns.yaml", "run_0022_ns_ib_pass0/delta_r_next.csv", 1,
     "pass 1: run_0022 から抽出 (旧 3 次平滑化)", "#d62728"),
    ("run_0025_ns_ib_pass2_q", "problem_d155_ns.yaml", "run_0023_ns_ib_pass1/delta_r_next.csv", 1,
     "pass 2: run_0023 から抽出 (5 次 P-spline λ=1)", "#ff7f0e"),
    ("run_0038_ns_final_rt77p02", "problem_d155_ns_rt77p02.yaml", "run_0025_ns_ib_pass2_q/delta_r_next.csv", 1,
     "最終: run_0025 から抽出 (P-spline λ=100)", "#2ca02c"),
]
EULER_OF = {"run_0022_ns_ib_pass0": "run_0001_euler_shortest_dry", "run_0023_ns_ib_pass1": "run_0001_euler_shortest_dry",
            "run_0025_ns_ib_pass2_q": "run_0001_euler_shortest_dry", "run_0038_ns_final_rt77p02": "run_0037_euler_rt77p02"}
FINAL_RES = {"run_0022_ns_ib_pass0": "res_24000.h5", "run_0023_ns_ib_pass1": "res_24000.h5",
             "run_0025_ns_ib_pass2_q": "res_24000.h5", "run_0038_ns_final_rt77p02": "res_12000.h5",
             "run_0001_euler_shortest_dry": None, "run_0037_euler_rt77p02": None}


def last_res(run):
    import re
    fs = sorted((C / run).glob("res_[0-9]*.h5"), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    return fs[-1].name


def load_field(run, res=None):
    rd = C / run
    info = json.loads((rd / "prepare_info.json").read_text())
    S = info["scale_m"]; ni = info["mesh"]["ni"]
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    nj = nc.shape[0] // ni
    res = res or last_res(run)
    with h5py.File(rd / res) as f:
        V = {k: f["/VALUE/" + k][:].reshape(ni, nj).astype(float) for k in ("Ux", "Uy", "sonic", "P")}
        if "dPdx" in f["/VALUE"]:
            V["dPdx"] = f["/VALUE/dPdx"][:].reshape(ni, nj).astype(float)
    X = (nc[:, 0] / S).reshape(ni, nj); R = (nc[:, 1] / S).reshape(ni, nj)
    V["M"] = np.hypot(V["Ux"], V["Uy"]) / V["sonic"]
    V["theta"] = np.arctan2(V["Uy"], V["Ux"])
    return X, R, V, S, res


def eta_prof(X, R, F, xq, eta):
    xa = X[:, 0]
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2)
    w = (xq - xa[i]) / (xa[i + 1] - xa[i])
    out = np.empty(len(xq))
    for k in range(len(xq)):
        v = [np.interp(eta, R[ii] / R[ii, -1], F[ii]) for ii in (i[k], i[k] + 1)]
        out[k] = (1 - w[k]) * v[0] + w[k] * v[1]
    return out


XQ = np.linspace(-2.0, 94.0, 961)


def bump(xq, d, lo, hi, blo, bhi):
    s = (xq >= lo) & (xq <= hi); b = (xq >= blo) & (xq <= bhi)
    k = np.argmax(d[s]); return float(d[s][k] - np.median(d[b])), float(xq[s][k])


# ---------------------------------------------------------------- 壁の再構成 (prepare_ns と同じ経路)
_designs = {}
walls = {}
for run, prob, csv, col, lab, colr in SERIES:
    if prob not in _designs:
        p = load_problem(C / prob); _designs[prob] = (p, design_chain(p))
    p, d = _designs[prob]
    tbl = np.loadtxt(C / csv, delimiter=",", skiprows=1)
    fdel = (lambda x, _t=tbl, _c=col: np.interp(x, _t[:, 0], _t[:, _c]))
    w = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]),
                           float(p.spec["Tt"]), _gam_or_gas(p), p.cp, offset="radial", delta_r_x=fdel)
    ref = np.loadtxt(C / run / "wall_physical.csv", delimiter=",", skiprows=1)
    S = float(p.spec["r_throat"])
    err = float(np.abs(w.r(ref[:, 0] / S) * S - ref[:, 1]).max())
    walls[run] = dict(w=w, inv=d["wall"], tbl=tbl, col=col, lab=lab, c=colr, S=S, recon_err_m=err)
NUM["wall_reconstruction_max_err_m"] = {r: v["recon_err_m"] for r, v in walls.items()}

# ---------------------------------------------------------------- Fig 1: 解析領域と境界
X, R, V, S38, _ = load_field("run_0038_ns_final_rt77p02")
with h5py.File(C / "run_0038_ns_final_rt77p02/nozzle.h5") as f:
    bnodes = {}
    for pid in ("1", "2", "3", "4"):
        g = f["BCONDS"][pid]
        key = "vizBfaceNodes" if "vizBfaceNodes" in g else list(g.keys())[0]
        bnodes[pid] = np.unique(g[key][:].ravel())
nc = np.c_[X.ravel(), R.ravel()]
BC = {"1": ("inlet (physID 1): 全圧入口 Pt 5.5 MPa, Tt 1600 K", "#1f77b4"),
      "2": ("outlet (physID 2): 静圧出口 2237 Pa (超音速流出では外挿)", "#9467bd"),
      "3": ("wall (physID 3): 断熱 no-slip (低 Re SST)", "#d62728"),
      "4": ("axis (physID 4): 軸対称の軸", "#2ca02c")}
fig, axs = plt.subplots(2, 1, figsize=(11, 6.2), gridspec_kw=dict(height_ratios=[1.5, 1]))
for ax, xl, yl, si, sj in ((axs[0], (-13, 96), (-0.5, 14), 10, 4), (axs[1], (-3, 6), (-0.1, 2.0), 1, 2)):
    m = (X[:, 0] >= xl[0] - 1) & (X[:, 0] <= xl[1] + 1)
    for i in np.where(m)[0][::si]:
        ax.plot(X[i], R[i], color="0.82", lw=0.3)
    for j in range(0, X.shape[1], sj):
        ax.plot(X[m, j], R[m, j], color="0.82", lw=0.3)
    for pid, (lab, colr) in BC.items():
        pts = nc[bnodes[pid]] if bnodes[pid].max() < len(nc) else None
        if pts is None:
            continue
        o = np.argsort(pts[:, 0] if pid in ("3", "4") else pts[:, 1])
        ax.plot(pts[o, 0], pts[o, 1], color=colr, lw=2.2, label=lab)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_xlabel("x / r_t"); ax.set_ylabel("r / r_t")
axs[0].set_aspect("equal"); axs[1].set_aspect("equal")
axs[0].annotate("", xy=(8, 6.5), xytext=(-8, 6.5), arrowprops=dict(arrowstyle="->", lw=1.5))
axs[0].text(0, 7.2, "主流", ha="center")
axs[0].set_title(f"全体 (格子線は i 方向 10 本に 1 本・j 方向 4 本に 1 本を表示; 実格子 {X.shape[0]}×{X.shape[1]})", fontsize=9)
axs[1].set_title("スロート近傍の拡大 (x = −3〜6 r_t; i 方向は全格子線)", fontsize=9)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=2, fontsize=8.5, frameon=False)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(OUT / "fig1_domain.png", dpi=130); plt.close(fig)
NUM["mesh_run_0038"] = dict(ni=int(X.shape[0]), nj=int(X.shape[1]), x_min=float(X.min()), x_max=float(X.max()),
                             r_max=float(R.max()), axis_gap_frac=float(R[0, 1] / R[0, -1]))

# ---------------------------------------------------------------- δ_r: 各手法の分布
xs = np.linspace(-1.0, 94.0, 3801)
fig, axs = plt.subplots(3, 1, figsize=(10, 10.5), sharex=False)
NUM["delta_r"] = {}
for run, prob, csv, col, lab, colr in SERIES:
    W = walls[run]
    dr = np.interp(xs, W["tbl"][:, 0], W["tbl"][:, col])
    axs[0].plot(xs, dr, color=colr, lw=1.4, label=lab + f"  → {run} の壁")
    NUM["delta_r"][run] = {f"x{int(x)}": float(np.interp(x, xs, dr)) for x in (0, 5, 10, 15, 20, 25, 30, 40, 60, 80, 94)}
# 各 run から抽出した生値 (点) — 次の壁の材料
for run, colr in (("run_0022_ns_ib_pass0", "#d62728"), ("run_0023_ns_ib_pass1", "#ff7f0e"), ("run_0025_ns_ib_pass2_q", "#2ca02c")):
    t = np.loadtxt(C / run / "delta_r_equiv.csv", delimiter=",", skiprows=1)
    axs[0].plot(t[::3, 0], t[::3, 1], ".", ms=2.0, color=colr, alpha=0.5)
axs[0].set_ylabel("δ_r / r_t"); axs[0].set_xlim(-1, 95)
axs[0].legend(fontsize=8, loc="upper left", frameon=False)
axs[0].set_title("(a) 壁に載せた半径方向排除厚 δ_r(x)。点 = 各 run の場から Euler 差で抽出した生値 (同色の線がそれを平滑化して次の壁に使った値)", fontsize=8.5)
base = np.interp(xs, walls["run_0022_ns_ib_pass0"]["tbl"][:, 0], walls["run_0022_ns_ib_pass0"]["tbl"][:, 1])
for run, prob, csv, col, lab, colr in SERIES[1:]:
    W = walls[run]
    dr = np.interp(xs, W["tbl"][:, 0], W["tbl"][:, col])
    axs[1].plot(xs, 100 * (dr / np.maximum(base, 1e-6) - 1), color=colr, lw=1.4, label=lab)
axs[1].axhline(0, color="0.5", lw=0.6)
axs[1].set_xlim(-1, 95); axs[1].set_ylim(-40, 40); axs[1].set_ylabel("δ_r / δ_r(積分法) − 1 [%]")
axs[1].legend(fontsize=8, frameon=False, loc="upper right")
axs[1].set_title("(b) 積分法 (pass 0) に対する比", fontsize=9)
for run, prob, csv, col, lab, colr in SERIES:
    W = walls[run]
    dr = np.interp(xs, W["tbl"][:, 0], W["tbl"][:, col])
    axs[2].plot(xs, dr, color=colr, lw=1.4)
axs[2].set_xlim(-1, 30); axs[2].set_ylim(0, 0.35)
axs[2].set_xlabel("x / r_t"); axs[2].set_ylabel("δ_r / r_t")
axs[2].set_title("(c) (a) の x = −1〜30 拡大", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "fig3_delta_r.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 壁の 2 階微分 (実際のスプライン)
xw = np.linspace(0.05, 94.0, 18801)
fig, axs = plt.subplots(3, 1, figsize=(10, 9.5))
NUM["wall_curvature"] = {}
for run, prob, csv, col, lab, colr in SERIES:
    W = walls[run]
    r2p = W["w"].r(xw, 2); r2i = W["inv"].r(xw, 2)
    dd = r2p - r2i                                   # = 壁に乗った δ_r'' (スプライン経由の実効値)
    axs[0].plot(xw, r2p, color=colr, lw=0.9, label=lab)
    axs[1].plot(xw, dd, color=colr, lw=0.9, label=lab)
    # 高周波成分: 幅 2 r_t の移動平均からの偏差
    k = int(round(2.0 / (xw[1] - xw[0])))
    sm = np.convolve(dd, np.ones(k) / k, mode="same")
    hf = dd - sm
    m = (xw > 3) & (xw < 90)
    axs[2].plot(xw, hf, color=colr, lw=0.6, label=lab)
    NUM["wall_curvature"][run] = dict(dd_rms_3_90=float(np.sqrt(np.mean(dd[m] ** 2))), hf_rms_3_90=float(np.sqrt(np.mean(hf[m] ** 2))),
                                      hf_max_3_90=float(np.abs(hf[m]).max()),
                                      r2_inv_absmax_3_90=float(np.abs(r2i[m]).max()))
r2i = walls["run_0022_ns_ib_pass0"]["inv"].r(xw, 2)
axs[0].plot(xw, r2i, color="k", lw=0.9, ls="--", label="非粘性設計壁 r_inv''")
axs[0].set_xlim(0, 94); axs[0].set_ylim(-0.02, 0.06); axs[0].set_ylabel("r_w'' [1/r_t]")
axs[0].legend(fontsize=7.5, frameon=False, ncol=2)
axs[0].set_title("(a) 物理壁の 2 階微分 r_w''(x) (メッシュ生成が使う 5 次スプラインをそのまま微分)", fontsize=9)
axs[1].axhline(0, color="0.5", lw=0.5)
axs[1].set_xlim(0, 94); axs[1].set_ylim(-0.006, 0.006); axs[1].set_ylabel("r_w'' − r_inv'' [1/r_t]")
axs[1].set_title("(b) 物理壁と設計壁の 2 階微分の差 (= 壁に乗った δ_r'')", fontsize=9)
axs[2].axhline(0, color="0.5", lw=0.5)
axs[2].set_xlim(0, 94); axs[2].set_ylim(-0.004, 0.004); axs[2].set_xlabel("x / r_t"); axs[2].set_ylabel("高周波成分 [1/r_t]")
axs[2].set_title("(c) (b) から幅 2 r_t の移動平均を引いた高周波成分", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "fig4_wall_curvature.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 軸 M
cache = {}


def prof(run, eta, res=None):
    k = (run, eta, res)
    if k not in cache:
        Xr, Rr, Vr, _, _ = load_field(run, res)
        cache[k] = eta_prof(Xr, Rr, Vr["M"], XQ, eta)
    return cache[k]


fig, axs = plt.subplots(3, 1, figsize=(10, 10))
for e in ("run_0037_euler_rt77p02",):
    axs[0].plot(XQ, 100 * (prof(e, 0.0) / 6 - 1), color="k", lw=1.2, label="Euler run_0037 (= run_0001)")
NUM["axis"] = {}
for run, prob, csv, col, lab, colr in SERIES:
    m0 = prof(run, 0.0); e0 = prof(EULER_OF[run], 0.0)
    m1 = prof(run, 0.10); e1 = prof(EULER_OF[run], 0.10)
    axs[0].plot(XQ, 100 * (m0 / 6 - 1), color=colr, lw=1.1, label=lab)
    d0 = 100 * (m0 / e0 - 1); d1 = 100 * (m1 / e1 - 1)
    axs[1].plot(XQ, d0, color=colr, lw=1.1, label=lab)
    axs[2].plot(XQ, d1, color=colr, lw=1.1, label=lab)
    w = (XQ >= 42) & (XQ <= 94)
    NUM["axis"][run] = dict(b70_eta0=bump(XQ, d0, 60, 80, 45, 90), b35_eta0=bump(XQ, d0, 25, 45, 20, 50),
                            b70_eta01=bump(XQ, d1, 60, 80, 45, 90), pp_eta0=float(d0[w].max() - d0[w].min()),
                            pp_eta01=float(d1[w].max() - d1[w].min()))
ew = (XQ >= 42) & (XQ <= 94)
me = 100 * (prof("run_0037_euler_rt77p02", 0.0) / 6 - 1)
NUM["axis"]["euler_pp_42_94"] = float(me[ew].max() - me[ew].min())
axs[0].set_xlim(10, 94); axs[0].set_ylim(-1.0, 0.8); axs[0].axhline(0, color="0.6", lw=0.5)
axs[0].set_ylabel("M_axis / 6 − 1 [%]"); axs[0].legend(fontsize=7.5, frameon=False, ncol=2, loc="lower right")
axs[0].set_title("(a) 軸上 (r=0) の M と設計値 6 の差", fontsize=9)
for ax, t in ((axs[1], "(b) 軸上 (r=0): 各 NS と同じ壁設計の Euler との比"), (axs[2], "(c) r/r_w = 0.1: 同上")):
    ax.set_xlim(10, 94); ax.set_ylim(-0.4, 0.6); ax.axhline(0, color="0.6", lw=0.5)
    ax.set_ylabel("M_NS / M_Euler − 1 [%]"); ax.set_title(t, fontsize=9)
axs[2].set_xlabel("x / r_t")
fig.tight_layout(); fig.savefig(OUT / "fig5_axis_mach.png", dpi=130); plt.close(fig)

# 山の時系列 (各 run の全スナップショット)
NUM["bump_series"] = {}
for run, *_ in SERIES:
    import re
    fs = sorted((C / run).glob("res_[0-9]*.h5"), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    ser = []
    for f in fs:
        st = int(re.findall(r"\d+", f.name)[0])
        if st == 0:
            continue
        d0 = 100 * (prof(run, 0.0, f.name) / prof(EULER_OF[run], 0.0) - 1)
        ser.append((st, *bump(XQ, d0, 60, 80, 45, 90)))
    NUM["bump_series"][run] = ser

# ---------------------------------------------------------------- 特性線 (C−) の追跡


def trace_cminus(run, x_starts, eta0=0.985):
    """壁近傍 (η=eta0, 境界層の外側) から下流・軸方向の C− 特性線 dr/dx = tan(θ−μ) を RK2 で追う。"""
    Xr, Rr, Vr, _, _ = load_field(run)
    xa = Xr[:, 0]; eta = Rr[0] / Rr[0, -1]
    if not np.allclose(Rr / Rr[:, -1:], eta[None, :], atol=1e-6):
        raise SystemExit("η が x で一定でない")
    mu = np.arcsin(np.clip(1.0 / np.maximum(Vr["M"], 1.0001), -1, 1))
    slope = np.tan(Vr["theta"] - mu)
    f_s = RegularGridInterpolator((xa, eta), slope, bounds_error=False, fill_value=None)
    rw = lambda x: np.interp(x, xa, Rr[:, -1])
    lines = []
    for x0 in x_starts:
        x, r = x0, eta0 * rw(x0); pts = [(x, r)]
        h = 0.02
        while r > 0 and x < xa[-1] - 0.1:
            s1 = float(f_s([[x, r / rw(x)]])[0])
            xm, rm = x + 0.5 * h, r + 0.5 * h * s1
            s2 = float(f_s([[xm, max(rm, 0) / rw(xm)]])[0])
            x, r = x + h, r + h * s2
            pts.append((x, r))
        lines.append(np.array(pts))
    return lines


x_starts = np.arange(4.0, 41.0, 2.0)
lines22 = trace_cminus("run_0022_ns_ib_pass0", x_starts)
NUM["cminus_wall_to_axis"] = {f"{x0:.0f}": float(L[-1, 0]) for x0, L in zip(x_starts, lines22)}

# ---------------------------------------------------------------- コンタ
X22, R22, V22, _, _ = load_field("run_0022_ns_ib_pass0")
X23, R23, V23, _, _ = load_field("run_0023_ns_ib_pass1")
X38, R38, V38, _, _ = load_field("run_0038_ns_final_rt77p02")


def pctl(a, lo=0.2, hi=99.8):
    return np.percentile(a, lo), np.percentile(a, hi)


# Fig 6: マッハ数 (最終 run_0038)
fig, axs = plt.subplots(2, 1, figsize=(11, 6.4), gridspec_kw=dict(height_ratios=[1.3, 1]))
vmin, vmax = 0.0, 6.2
for ax, xl, yl in ((axs[0], (-13, 95), (0, 10.5)), (axs[1], (20, 95), (0, 10.5))):
    cs = ax.pcolormesh(X38, R38, V38["M"], cmap="turbo", vmin=vmin, vmax=vmax, shading="gouraud", rasterized=True)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
axs[1].set_xlabel("x / r_t")
axs[0].set_title(f"マッハ数 (run_0038 res_12000, 最大 {V38['M'].max():.3f})", fontsize=9)
axs[1].set_title("x = 20〜95 r_t (試験部側)", fontsize=9)
fig.colorbar(cs, ax=axs, orientation="horizontal", fraction=0.05, pad=0.1, label="M")
fig.savefig(OUT / "fig6_mach.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# Fig 7: 試験部の M の微細構造 (M − 6 を ±1 % で)、3 run
fig, axs = plt.subplots(3, 1, figsize=(11, 8.6), constrained_layout=True)
for ax, (run, Xr, Rr, Vr) in zip(axs, (("run_0022_ns_ib_pass0 (積分法壁)", X22, R22, V22),
                                         ("run_0023_ns_ib_pass1 (抽出壁 pass 1)", X23, R23, V23),
                                         ("run_0038_ns_final_rt77p02 (最終壁)", X38, R38, V38))):
    dm = 100 * (Vr["M"] / 6 - 1)
    cs = ax.pcolormesh(Xr, Rr, dm, cmap="turbo", vmin=-0.8, vmax=0.8, shading="gouraud", rasterized=True)
    ax.set_xlim(40, 95); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(f"{run}: M/6 − 1 [%] (±0.8 % で頭打ち。左上の暗部は膨張途中 (M<5.95) と境界層)", fontsize=9)
for L in lines22:
    for ax in axs:
        ax.plot(L[:, 0], L[:, 1], color="w", lw=0.5, alpha=0.6)
axs[-1].set_xlabel("x / r_t")
fig.colorbar(cs, ax=axs, orientation="horizontal", fraction=0.04, pad=0.02, label="M/6 − 1 [%]")
fig.savefig(OUT / "fig7_mach_finestructure.png", dpi=130); plt.close(fig)

# Fig 8: 壁を替えた応答 ΔM (積分法壁 run_0022 基準、同一トポロジ (i,j) で対応; r_t 単位の座標は 0.1 % 以内で一致)
X25, R25, V25, _, _ = load_field("run_0025_ns_ib_pass2_q")
NUM["dM_vs_0022"] = {}
pairs = (("run_0023 (pass 1, 旧 3 次平滑化) − run_0022", V23, X23, R23),
         ("run_0025 (pass 2, P-spline) − run_0022", V25, X25, R25),
         ("run_0038 (最終, P-spline) − run_0022", V38, X38, R38))
fig, axs = plt.subplots(4, 1, figsize=(11, 11.5), constrained_layout=True)
for ax, (lab, Vb, Xb, Rb) in zip(axs[:3], pairs):
    dM = 100 * (Vb["M"] / V22["M"] - 1)
    m = (Xb[:, :1] > 2.0) & (np.arange(dM.shape[1])[None, :] < int(0.9 * dM.shape[1]))
    NUM["dM_vs_0022"][lab] = float(np.abs(dM[m]).max())
    cs = ax.pcolormesh(Xb, Rb, dM, cmap="turbo", vmin=-0.5, vmax=0.5, shading="gouraud", rasterized=True)
    for L in lines22:
        ax.plot(L[:, 0], L[:, 1], color="k", lw=0.4, alpha=0.45)
    ax.set_xlim(-1, 95); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(f"{lab}: 100·(M/M_0022 − 1) [%] (±0.5 で頭打ち)", fontsize=9)
dM = 100 * (V23["M"] / V22["M"] - 1)
cs = axs[3].pcolormesh(X23, R23, dM, cmap="turbo", vmin=-0.5, vmax=0.5, shading="gouraud", rasterized=True)
for L in lines22:
    axs[3].plot(L[:, 0], L[:, 1], color="k", lw=0.5, alpha=0.5)
axs[3].set_xlim(0, 40); axs[3].set_ylim(0, 9); axs[3].set_aspect("equal"); axs[3].set_ylabel("r / r_t")
axs[3].set_xlabel("x / r_t")
axs[3].set_title("run_0023 − run_0022 の x = 0〜40 拡大。黒線 = run_0022 の場で追った C− 特性線 (壁 x = 4, 6, …, 40 発)", fontsize=9)
fig.colorbar(cs, ax=axs, orientation="horizontal", fraction=0.03, pad=0.02, label="ΔM [%]")
fig.savefig(OUT / "fig8_dM_wallchange.png", dpi=130); plt.close(fig)

# Fig 9: 無次元圧力勾配 (r_t/p) ∂p/∂x
fig, axs = plt.subplots(3, 1, figsize=(11, 8.6))
for ax, (lab, Xr, Rr, Vr) in zip(axs, (("run_0022 (積分法壁)", X22, R22, V22), ("run_0023 (抽出壁 pass 1)", X23, R23, V23),
                                        ("差 run_0023 − run_0022", X23, R23, None))):
    if Vr is None:
        g = V23["dPdx"] / V23["P"] - V22["dPdx"] / V22["P"]
        lim = 0.004
    else:
        g = Vr["dPdx"] / Vr["P"]
        lim = 0.02
    g = g * walls["run_0022_ns_ib_pass0"]["S"]           # (r_t/p) ∂p/∂x (dPdx は [Pa/m])
    cs = ax.pcolormesh(Xr, Rr, g, cmap="turbo", vmin=-lim, vmax=lim, shading="gouraud", rasterized=True)
    ax.set_xlim(20, 95); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(f"{lab}: (r_t/p) ∂p/∂x  (±{lim} で頭打ち)", fontsize=9)
    fig.colorbar(cs, ax=ax, orientation="vertical", fraction=0.02, pad=0.01)
axs[-1].set_xlabel("x / r_t")
fig.tight_layout(); fig.savefig(OUT / "fig9_dpdx.png", dpi=130); plt.close(fig)

(OUT / "numbers.json").write_text(json.dumps(NUM, indent=1, default=float))
print(json.dumps(NUM, indent=1, default=float)[:6000])
