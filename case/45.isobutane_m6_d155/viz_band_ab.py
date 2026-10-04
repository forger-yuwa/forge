"""帯修正の CFD A/B (plan verification-m6-axis-wave-mesh-su2 §4.4) の図と ParaView 用ファイル。

出力:
  _band_ab/viz/fig_*.png            軸 M・壁 (r, r', r'')・圧力勾配と M のコンタ
  _band_ab/paraview/<tag>.h5/.xmf   構造格子 (2DSMesh) に派生量 (Mach, M/6−1, M/M_Euler−1, (r_t/p)∂p/∂x) を載せたもの
  _band_ab/paraview/B1_minus_B0.*   同じトポロジの差 (B1 − B0)
usage: design/.venv-opt/bin/python viz_band_ab.py
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
VIZ = C / "_band_ab" / "viz"; VIZ.mkdir(parents=True, exist_ok=True)
PV = C / "_band_ab" / "paraview"; PV.mkdir(parents=True, exist_ok=True)
RUNS = {  # tag: (run, 壁の δ_r CSV, 線色, ラベル)
    "Euler": ("run_0047_euler_rt77p02_newbin", None, "k", "Euler (設計壁、新バイナリ) run_0047"),
    "A0p": ("run_0042_ns_restart_ctrl", "run_0025_ns_ib_pass2_q/delta_r_next.csv", "#7f7f7f", "旧最終壁 (adaptive 3 回目) run_0042"),
    "B0": ("run_0048_ns_band_adaptive_ext", "_band_ab/adaptive/delta_r_next.csv", "#d62728", "B0: adaptive で作り直した壁 run_0045→0048"),
    "B1": ("run_0046_ns_band_edge", "_band_ab/edge_tests_v2/T5_F3_c1.25/delta_r_next.csv", "#1f77b4", "B1: 方式 E で作り直した壁 run_0046"),
}
S = 0.07702
NUM = {}


def last_res(run):
    import re
    fs = sorted((C / run).glob("res_[0-9]*.h5"), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    return fs[-1]


def field(run):
    info = json.loads((C / run / "prepare_info.json").read_text()); ni = info["mesh"]["ni"]
    with h5py.File(C / run / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    nj = nc.shape[0] // ni
    res = last_res(run)
    with h5py.File(res) as f:
        V = {k: f["/VALUE/" + k][:].astype(float).reshape(ni, nj) for k in ("ro", "Ux", "Uy", "P", "T", "sonic")}
    V["M"] = np.hypot(V["Ux"], V["Uy"]) / V["sonic"]
    V["theta"] = np.arctan2(V["Uy"], V["Ux"])
    X = (nc[:, 0] / S).reshape(ni, nj); R = (nc[:, 1] / S).reshape(ni, nj)
    # ∂p/∂x (x は r_t 単位) を構造格子の座標変換で: p_x = (p_ξ r_η − p_η r_ξ) / (x_ξ r_η − x_η r_ξ)
    px_i, px_j = np.gradient(V["P"]); x_i, x_j = np.gradient(X); r_i, r_j = np.gradient(R)
    V["dPdx_rt"] = (px_i * r_j - px_j * r_i) / (x_i * r_j - x_j * r_i)
    return X, R, V, res.name


F = {k: field(v[0]) for k, v in RUNS.items()}
# Euler の M を (x, η) で補間する関数 (Euler 格子も η が全断面共通)
XE, RE, VE, _ = F["Euler"]
etaE = RE[0] / RE[0, -1]
assert np.allclose(RE / RE[:, -1:], etaE[None, :], atol=1e-6)
ME_interp = RegularGridInterpolator((XE[:, 0], etaE), VE["M"], bounds_error=False, fill_value=None)


def m_vs_euler(X, R):
    eta = R / R[:, -1:]
    return ME_interp(np.c_[X.ravel(), eta.ravel()]).reshape(X.shape)


# ---------------------------------------------------------------- 壁の再構成 (prepare_ns と同じ経路)
p = load_problem(C / "problem_d155_ns_rt77p02.yaml"); d = design_chain(p)
walls = {}
for k, (run, csv, col, lab) in RUNS.items():
    if csv is None:
        continue
    t = np.loadtxt(C / csv, delimiter=",", skiprows=1)
    w = PhysicalNozzleWall(d["wall"], d["wall_inv"], S, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp,
                           offset="radial", delta_r_x=lambda x, _t=t: np.interp(x, _t[:, 0], _t[:, 1]))
    ref = np.loadtxt(C / run / "wall_physical.csv", delimiter=",", skiprows=1)
    walls[k] = w
    NUM.setdefault("wall_recon_err_m", {})[k] = float(np.abs(w.r(ref[:, 0] / S) * S - ref[:, 1]).max())

# ---------------------------------------------------------------- 図 1: 軸 M
xq = np.linspace(-2, 95, 1941)


def eta_line(X, R, Mf, eta):
    xa = X[:, 0]; out = np.empty(len(xq))
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2); wgt = (xq - xa[i]) / (xa[i + 1] - xa[i])
    for k in range(len(xq)):
        v = [np.interp(eta, R[ii] / R[ii, -1], Mf[ii]) for ii in (i[k], i[k] + 1)]
        out[k] = (1 - wgt[k]) * v[0] + wgt[k] * v[1]
    return out


fig, axs = plt.subplots(3, 1, figsize=(10.5, 10.5))
ME0 = eta_line(XE, RE, VE["M"], 0.0); ME1 = eta_line(XE, RE, VE["M"], 0.1)
for k, (run, csv, col, lab) in RUNS.items():
    X, R, V, _ = F[k]
    m0 = eta_line(X, R, V["M"], 0.0); m1 = eta_line(X, R, V["M"], 0.1)
    axs[0].plot(xq, m0, color=col, lw=1.2, label=lab)
    if k != "Euler":
        axs[1].plot(xq, 100 * (m0 / ME0 - 1), color=col, lw=1.2, label=lab)
        axs[2].plot(xq, 100 * (m1 / ME1 - 1), color=col, lw=1.2, label=lab)
axs[0].set_xlim(-2, 95); axs[0].set_ylabel("軸上の M"); axs[0].legend(fontsize=8, frameon=False, loc="lower right")
axs[0].set_title("(a) 軸上 (r = 0) のマッハ数", fontsize=10)
for ax, t in ((axs[1], "(b) 軸上: M_NS / M_Euler − 1 [%] (Euler = run_0047、同じバイナリ)"),
              (axs[2], "(c) r / r_w = 0.1: 同上 (軸ノードの値の影響を避けた線)")):
    ax.set_xlim(10, 95); ax.set_ylim(-0.6, 0.6); ax.axhline(0, color="0.6", lw=0.5); ax.axvspan(42, 94, color="0.92", zorder=0)
    ax.set_ylabel("[%]"); ax.set_title(t, fontsize=10)
axs[2].set_xlabel("x / r_t"); axs[1].legend(fontsize=8, frameon=False, loc="lower left")
axs[1].text(43, 0.5, "灰色 = 判定区間 x = 42〜94", fontsize=8, color="0.4")
fig.tight_layout(); fig.savefig(VIZ / "fig_axis_mach.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 図 2: 壁 r, r', r''
xw = np.linspace(-12.4, 95.0, 21481)
fig, axs = plt.subplots(3, 2, figsize=(13, 10.5))
r_inv = [d["wall"].r(xw, k) for k in (0, 1, 2)]
for row, (k_der, nm) in enumerate(((0, "r_w / r_t"), (1, "dr_w/dx"), (2, "d²r_w/dx² [1/r_t]"))):
    axs[row, 0].plot(xw, r_inv[k_der], color="k", lw=1.0, ls="--", label="設計壁 (非粘性)")
    for k, w in walls.items():
        v = w.r(xw, k_der)
        axs[row, 0].plot(xw, v, color=RUNS[k][2], lw=0.9, label=RUNS[k][3])
        axs[row, 1].plot(xw, v - r_inv[k_der], color=RUNS[k][2], lw=0.9, label=RUNS[k][3])
    axs[row, 0].set_ylabel(nm); axs[row, 1].set_ylabel(nm.split(" [")[0] + " − 設計壁")
    axs[row, 1].axhline(0, color="0.6", lw=0.5)
for a in axs[:, 0]:
    a.set_xlim(-12.5, 95)
for a in axs[:, 1]:
    a.set_xlim(-12.5, 95)
axs[0, 0].set_title("(a) 壁の形 (左: そのもの / 右: 設計壁との差)", fontsize=10)
axs[1, 0].set_ylim(-0.4, 0.35); axs[1, 1].set_ylim(-0.004, 0.004)
axs[2, 0].set_ylim(-0.03, 0.08); axs[2, 1].set_ylim(-0.003, 0.003)
axs[0, 1].set_ylim(-0.03, 0.75)
axs[1, 0].set_title("(b) 1 階微分 (壁の傾き)", fontsize=10); axs[2, 0].set_title("(c) 2 階微分 (壁の曲率)", fontsize=10)
axs[0, 0].legend(fontsize=7.5, frameon=False, loc="upper right")
for a in axs[2]:
    a.set_xlabel("x / r_t")
fig.tight_layout(); fig.savefig(VIZ / "fig_wall_derivatives.png", dpi=130); plt.close(fig)
# 拡大: x = -1〜60 の差分
fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
for row, k_der in enumerate((0, 1, 2)):
    for k, w in walls.items():
        axs[row].plot(xw, w.r(xw, k_der) - r_inv[k_der], color=RUNS[k][2], lw=1.0, label=RUNS[k][3])
    axs[row].axhline(0, color="0.6", lw=0.5)
axs[0].set_ylabel("r_w − 設計壁"); axs[1].set_ylabel("r_w′ − 設計壁"); axs[2].set_ylabel("r_w″ − 設計壁")
axs[0].set_ylim(-0.02, 0.45); axs[1].set_ylim(-0.0005, 0.009); axs[2].set_ylim(-0.0015, 0.0015)
axs[2].set_xlim(-1, 60); axs[2].set_xlabel("x / r_t"); axs[0].legend(fontsize=8, frameon=False)
axs[0].set_title("設計壁との差 (x = −1〜60 の拡大)。B0 は x ≈ 35〜50、旧最終壁は x ≈ 17〜29 にへこみ", fontsize=10)
fig.tight_layout(); fig.savefig(VIZ / "fig_wall_derivatives_zoom.png", dpi=130); plt.close(fig)
for k, w in walls.items():
    m = (xw > 3) & (xw < 90); dd = w.r(xw[m], 2) - r_inv[2][m]
    kk = int(2.0 / (xw[1] - xw[0])); hf = dd - np.convolve(dd, np.ones(kk) / kk, "same")
    NUM.setdefault("wall_d2_hf_rms", {})[k] = float(np.sqrt(np.mean(hf[kk:-kk] ** 2)))


# ---------------------------------------------------------------- C− 特性線 (B1 の場)
def trace_cminus(X, R, V, x_starts, eta0=0.985):
    xa = X[:, 0]; eta = R[0] / R[0, -1]
    mu = np.arcsin(np.clip(1.0 / np.maximum(V["M"], 1.0001), -1, 1))
    f_s = RegularGridInterpolator((xa, eta), np.tan(V["theta"] - mu), bounds_error=False, fill_value=None)
    rw = lambda x: np.interp(x, xa, R[:, -1])
    lines = []
    for x0 in x_starts:
        x, r = x0, eta0 * rw(x0); pts = [(x, r)]; h = 0.02
        while r > 0 and x < xa[-1] - 0.1:
            s1 = float(f_s([[x, r / rw(x)]])[0]); xm, rm = x + 0.5 * h, r + 0.5 * h * s1
            s2 = float(f_s([[xm, max(rm, 0) / rw(xm)]])[0]); x, r = x + h, r + h * s2; pts.append((x, r))
        lines.append(np.array(pts))
    return lines


XB1, RB1, VB1, _ = F["B1"]; XB0, RB0, VB0, _ = F["B0"]; XA, RA, VA, _ = F["A0p"]
lines = trace_cminus(XB1, RB1, VB1, np.arange(4.0, 41.0, 3.0))
NUM["cminus_B1"] = {f"{L[0,0]:.0f}": float(L[-1, 0]) for L in lines}

# ---------------------------------------------------------------- 図 3: 圧力勾配のコンタ
def nd_dpdx(V):
    """(r_t/p) ∂p/∂x。∂p/∂x は格子から計算 (x が r_t 単位なのでそのまま r_t を掛けた形)。"""
    return V["dPdx_rt"] / V["P"]


fig, axs = plt.subplots(4, 1, figsize=(11.5, 13), constrained_layout=True)
panels = (("B1 (方式 E の壁): (r_t/p) ∂p/∂x", XB1, RB1, nd_dpdx(VB1), 0.02),
          ("B0 (adaptive の壁): (r_t/p) ∂p/∂x", XB0, RB0, nd_dpdx(VB0), 0.02),
          ("差 B0 − B1 (同じ格子トポロジの節点どうし)", XB1, RB1, nd_dpdx(VB0) - nd_dpdx(VB1), 0.003),
          ("差 旧最終壁 A0′ − B1", XB1, RB1, nd_dpdx(VA) - nd_dpdx(VB1), 0.003))
for ax, (t, X, R, Z, lim) in zip(axs, panels):
    cs = ax.pcolormesh(X, R, Z, cmap="turbo", vmin=-lim, vmax=lim, shading="gouraud", rasterized=True)
    for L in lines:
        ax.plot(L[:, 0], L[:, 1], color="k", lw=0.4, alpha=0.45)
    ax.set_xlim(0, 95); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(t + f"  (±{lim} で頭打ち)", fontsize=9.5)
    fig.colorbar(cs, ax=ax, orientation="vertical", fraction=0.02, pad=0.01)
axs[-1].set_xlabel("x / r_t")
fig.savefig(VIZ / "fig_dpdx_contours.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 図 4: M のコンタ (Euler 比)
fig, axs = plt.subplots(3, 1, figsize=(11.5, 10), constrained_layout=True)
for ax, (t, X, R, V) in zip(axs, (("B1 (方式 E の壁)", XB1, RB1, VB1), ("B0 (adaptive の壁)", XB0, RB0, VB0), ("旧最終壁 A0′", XA, RA, VA))):
    Z = 100 * (V["M"] / m_vs_euler(X, R) - 1)
    cs = ax.pcolormesh(X, R, Z, cmap="turbo", vmin=-0.5, vmax=0.5, shading="gouraud", rasterized=True)
    for L in lines:
        ax.plot(L[:, 0], L[:, 1], color="k", lw=0.4, alpha=0.45)
    ax.set_xlim(0, 95); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(f"{t}: M / M_Euler − 1 [%] (±0.5 で頭打ち、境界層は下限に張り付く)", fontsize=9.5)
axs[-1].set_xlabel("x / r_t")
fig.colorbar(cs, ax=axs, orientation="horizontal", fraction=0.03, pad=0.02, label="[%]")
fig.savefig(VIZ / "fig_mach_vs_euler_contours.png", dpi=130); plt.close(fig)


# ---------------------------------------------------------------- ParaView 用 (2DSMesh)
def write_pv(tag, X, R, fields, note):
    h5p = PV / f"{tag}.h5"
    ni, nj = X.shape
    with h5py.File(h5p, "w") as f:
        f["x"] = X.astype("f8"); f["r"] = R.astype("f8")
        for k, v in fields.items():
            f[k] = v.astype("f8")
        f.attrs["note"] = note
    att = "\n".join(f"""      <Attribute Name="{k}" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{ni} {nj}" NumberType="Float" Precision="8" Format="HDF">{h5p.name}:/{k}</DataItem>
      </Attribute>""" for k in fields)
    (PV / f"{tag}.xmf").write_text(f"""<?xml version="1.0" ?>
<!-- {note} -->
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="{tag}" GridType="Uniform">
      <Topology TopologyType="2DSMesh" Dimensions="{ni} {nj}"/>
      <Geometry GeometryType="X_Y">
        <DataItem Dimensions="{ni} {nj}" NumberType="Float" Precision="8" Format="HDF">{h5p.name}:/x</DataItem>
        <DataItem Dimensions="{ni} {nj}" NumberType="Float" Precision="8" Format="HDF">{h5p.name}:/r</DataItem>
      </Geometry>
{att}
    </Grid>
  </Domain>
</Xdmf>
""")


for k, (run, csv, col, lab) in RUNS.items():
    X, R, V, resname = F[k]
    flds = dict(Mach=V["M"], M_over_6_minus1_pct=100 * (V["M"] / 6 - 1), dpdx_nd=nd_dpdx(V), P=V["P"], T=V["T"], rho=V["ro"],
                Ux=V["Ux"], Uy=V["Uy"])
    if k != "Euler":
        flds["M_vs_Euler_pct"] = 100 * (V["M"] / m_vs_euler(X, R) - 1)
    write_pv(k, X, R, flds, f"{lab}: {run}/{resname}。座標は r_t 単位 (r_t = {S} m)")
write_pv("B1_minus_B0", XB1, RB1,
         dict(dM_pct=100 * (VB1["M"] / VB0["M"] - 1), d_dpdx_nd=nd_dpdx(VB1) - nd_dpdx(VB0)),
         "B1 (run_0046) − B0 (run_0048) を同じ (i, j) 節点で。座標は B1")
(VIZ / "numbers.json").write_text(json.dumps(NUM, indent=1))
print(json.dumps(NUM, indent=1))
