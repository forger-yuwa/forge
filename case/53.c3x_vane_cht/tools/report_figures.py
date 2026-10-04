#!/usr/bin/env python3
r"""報告 (Artifact「Cooled Vane CHT Validation」) の**線グラフを全部ここから作る**。

2026-09-21 に `/tmp` の作業領域が再起動で消え、図ごとに散らばっていた作図スクリプトを失った。
以後は本ファイルが唯一の作図元。**どの図がどの run から出るかは下の `RUNS` が正本**で、
報告の run 一覧・キャプションはこれと一致させる。

図 (出力ファイル名):
  h_c3x.png / h_mk.png         実測 $T_w$ を課した $h$ (灰 = 全壁節点, 黒 = 熱電対位置で読んだ同じ解)
  press_c3x.png / press_mk.png 表面静圧
  su2.png                      C3X 一様壁: forge と SU2 (同一メッシュ・同一節点)
  turb.png                     C3X 入口乱流スイープ
  lam_c3x.png / lam_mk.png     層流対照 (層流解が定常な区間だけ描く)
  trans_plate.png              遷移平板 T3A の $C_f$ (forge 3 格子・SU2 LM・実験・遷移なし SST・層流)
  trans_c3x.png / trans_mk.png 遷移モデル ON/OFF と入口粘性比の感度 (実測 $T_w$、同一メッシュ)
  trans_su2.png                C3X 一様壁: forge と SU2 の遷移モデルどうし (同一メッシュ)
  trans_sweep.png              入口乱流粘性比と Re_θt 下限の掃引 (領域別 h 偏差。compare_h.py と同じ規約)

規約: 熱流束は全図 `iface_q_eff`。実測点には報告の表 V/VI の不確かさを付ける。
凡例はデータに重ねない (skill `forge-contour` と同じ)。コンタ図は
`solver_density_cuda/tools/plot_field_contours.py`、形状重ね合わせは `check_geometry.py`、
連成壁温は `plot_tw.py`。

usage: python3 case/53.c3x_vane_cht/tools/report_figures.py --out <dir> [--only h_c3x su2 ...]
"""
from __future__ import annotations

import argparse
import glob
import re
import struct
import sys
from pathlib import Path

import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from compare_h import arc_map                       # noqa: E402
from ref_data import PRESS                          # noqa: E402
from run_data import TABLES                         # noqa: E402
from uncertainty import h_uncertainty_pct           # noqa: E402

H0 = 1135.0
FLUX = "iface_q_eff"
C3 = ROOT / "case/53.c3x_vane_cht"
MK = ROOT / "case/54.markii_vane_cht"

# ---- どの図がどの run か (報告の run 一覧と一致させること) ----
# 2026-09-21: node 内部双対面の fx=0.5 固定 (plan discretization-node-face-weight-midpoint) 後の run に差し替え。
# `run_0125` / `run_0024` は環境変数 FORGE_NODE_FX_HALF=1 の試行バイナリ (恒久実装と同じスキーム。run_0127 で照合)。
RUNS = {
    "c3x": dict(base=C3, key="run108", Tg=786.0, Pt=319500.0, table="VI",
                prod="run_0125_prod_fxhalf", laminar="run_0134_laminar_fx05",
                prod_note="cs400 mesh, reconT=1, SU2-matched turbulence, exit 188.9 kPa",
                name="C3X run 108"),
    "markii": dict(base=MK, key="run42", Tg=788.0, Pt=337100.0, table="V",
                   prod="run_0024_tecut_fxhalf", laminar="run_0026_laminar_fx05",
                   prod_note="cs400 mesh with the base cut, reconT=1, exit 169.9 kPa",
                   name="Mark II run 42"),
}
SU2_CTRL = dict(forge="run_0128_cf0_fx05", su2="su2_smooth/vol_solution.vtu", Tw=566.0)
TURB = [("run_0133_tu3_fx05", "Tu 3 %", "#2a9d8f", "--"),
        ("run_0128_cf0_fx05", r"Tu 6.5 %, $\mu_t/\mu$=10  (used here)", "k", "-"),
        ("run_0132_tu15_fx05", "Tu 15 %", "#e9a23b", "--"),
        ("run_0130_mut1_fx05", r"$\mu_t/\mu$ = 1", "tab:blue", "-"),
        ("run_0131_mut100_fx05", r"$\mu_t/\mu$ = 100", "tab:red", "-")]


def last_wall(run_dir):
    fs = sorted(glob.glob(str(run_dir / "res_wall_5_*.h5")),
                key=lambda p: int(re.findall(r"(\d+)\.h5", p)[0]))
    if not fs:
        sys.exit(f"[report_figures] 壁ダンプが無い: {run_dir}")
    return fs[-1]


def wall(run_dir):
    f = last_wall(run_dir)
    w = h5py.File(f, "r")
    C = np.array(w["MESH/COORD"]).reshape(-1, 3)[:, :2]
    s, ss = arc_map(C)
    step = int(re.findall(r"(\d+)\.h5", f)[0])
    return s, ss, w["VALUE"], C, step


def measured(key):
    rows = [r for r in TABLES[key]["rows"] if r[3] is not None]
    sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows])
    return sd, hd, int(np.argmin(sd))


def plot_measured(ax, vane, key, table):
    sd, hd, i = measured(key)
    for sgn, sl, lab, col in ((+1, slice(i, None), "SS", "tab:red"), (-1, slice(0, i + 1), "PS", "tab:blue")):
        eu = h_uncertainty_pct(vane, lab == "SS", np.abs(sd[sl])) / 100.0 * hd[sl]
        ax.errorbar(sgn * sd[sl], hd[sl], yerr=eu, fmt="o", ms=4.6, color=col, elinewidth=1.0,
                    capsize=2.5, label=f"measured {lab}  (± report Table {table})", zorder=6)
    return sd, hd, i


def sides(ax, s, ss, y, mask=None, **kw):
    lab = kw.pop("label", None)
    for side, sgn in (("SS", +1), ("PS", -1)):
        m = ss if side == "SS" else ~ss
        if mask is not None:
            m = m & mask(side)
        o = np.argsort(s[m])
        ax.plot(sgn * s[m][o], y[m][o], label=(lab if side == "SS" else None), **kw)


def finish(ax, fig, out, title, legend_below=True):
    ax.set_xlabel("$-s/S$  (pressure side)   |   $+s/S$  (suction side)")
    ax.grid(alpha=.3)
    ax.set_title(title, fontsize=10)
    if legend_below:            # 凡例はデータに重ねない
        ax.legend(fontsize=8.6, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=115, bbox_inches="tight")
    print(f"[report_figures] -> {out}")


def fig_h(vane, out, plt):
    R = RUNS[vane]
    s, ss, V, _, step = wall(R["base"] / R["prod"])
    Ts = np.array(V["Ts"]); h = np.array(V[FLUX]) / (R["Tg"] - Ts) / H0
    fig, a = plt.subplots(figsize=(9.2, 5.6))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    sd, hd, i = plot_measured(a, vane, R["key"], R["table"])
    sides(a, s, ss, h, color="#b8b8b8", lw=1.0, zorder=2, label=f"forge, every wall node ({len(s)})")
    for side, sgn, sl in (("SS", +1, slice(i, None)), ("PS", -1, slice(0, i + 1))):
        m = ss if side == "SS" else ~ss; o = np.argsort(s[m])
        x = np.abs(sd[sl]); k = np.argsort(x)
        a.plot(sgn * x[k], np.interp(x, s[m][o], h[m][o])[k], "-", color="k", lw=1.9, zorder=4,
               label=("same solution, read only where the thermocouples are" if side == "SS" else None))
    axt = a.twinx()
    sides(axt, s, ss, Ts, color="#c9b6e4", lw=2.0, alpha=.6, zorder=1, label="imposed $T_w$ (right axis)")
    axt.set_ylabel("imposed $T_w$  [K]", color="#7f62a8"); axt.tick_params(axis="y", colors="#7f62a8")
    a.set_zorder(axt.get_zorder() + 1); a.patch.set_visible(False)
    a.set_ylabel("$h/h_0$")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = axt.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=8.4, loc="upper center", bbox_to_anchor=(0.5, -0.13),
             ncol=2, frameon=False)
    finish(a, fig, out, f"{R['name']} — heat transfer coefficient, measured $T_w$ imposed\n"
           f"{R['prod']}, step {step}  ({R['prod_note']})", legend_below=False)


def fig_press(vane, out, plt):
    R = RUNS[vane]
    s, ss, V, _, step = wall(R["base"] / R["prod"])
    P = np.array(V["Ps"]) / R["Pt"]
    fig, a = plt.subplots(figsize=(8.6, 5.3))
    for side, sgn, col in (("SS", +1, "tab:red"), ("PS", -1, "tab:blue")):
        d = np.array(PRESS[R["key"]][side])
        a.plot(sgn * d[:, 0], d[:, 2], "o", ms=5.5, color=col, label="measured " + side, zorder=5)
        m = ss if side == "SS" else ~ss; o = np.argsort(s[m])
        e = np.interp(d[:, 0], s[m][o], P[m][o]) - d[:, 2]
        print(f"   {vane} {side} n={len(d)} bias {e.mean():+.4f} rms {np.sqrt((e**2).mean()):.4f}")
    sides(a, s, ss, P, color="k", lw=1.7, label="forge")
    a.set_ylabel("$P_s/P_{T1}$")
    finish(a, fig, out, f"{R['name']} — surface static pressure\n{R['prod']}, step {step}")


def read_su2(path):
    DT = {"Float32": ("f", 4), "Float64": ("d", 8)}
    raw = open(path, "rb").read()
    hdr = raw[:raw.index(b"<AppendedData")].decode("utf8", "replace")
    A = re.findall(r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents=\s*"(\d+)" offset="(\d+)"', hdr)
    st = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    def g(name):
        for t, n, c, o in A:
            if n == name:
                off = st + int(o); nb = struct.unpack_from("<Q", raw, off)[0]; f, sz = DT[t]
                arr = np.frombuffer(raw, dtype=np.dtype("<" + f), count=nb // sz, offset=off + 8)
                return arr.reshape(-1, int(c)) if int(c) > 1 else arr
        sys.exit(f"[report_figures] SU2 出力に {name!r} が無い")
    return g("")[:, :2], np.abs(g("Heat_Flux"))


def fig_su2(out, plt):
    from scipy.spatial import cKDTree
    s, ss, V, C, step = wall(C3 / SU2_CTRL["forge"])
    XY, QS = read_su2(C3 / SU2_CTRL["su2"])
    d, idx = cKDTree(XY).query(C)
    assert d.max() == 0.0, f"forge と SU2 の壁節点が一致しない (max {d.max()})"
    sc = 1.0 / ((RUNS["c3x"]["Tg"] - SU2_CTRL["Tw"]) * H0)
    fig, a = plt.subplots(figsize=(9.2, 5.6))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    plot_measured(a, "c3x", "run108", "VI")
    sides(a, s, ss, np.array(V[FLUX]) * sc, color="k", lw=2.0, zorder=4, label="forge")
    sides(a, s, ss, QS[idx] * sc, color="tab:green", lw=1.7, zorder=3, label="SU2 8.5, same mesh")
    a.set_ylabel("$h/h_0$")
    finish(a, fig, out, "C3X run 108 — forge and SU2 on the identical mesh, uniform $T_w$ = 566 K\n"
           f"{SU2_CTRL['forge']}, step {step} . su2_smooth  (wall nodes coincide exactly)")


# ---- 遷移モデル (plan turbulence-transition-lm2009)。run は case README の「遷移モデル」節と一致させる ----
T57 = ROOT / "case/57.transition_flat_plate"
TRANS = {
    "c3x": [("run_0135_hwall1um_fx05", "no transition model, inlet $\\mu_t/\\mu$ = 10", "#1f77b4", "-"),
            ("run_0149_sst_1um_mur100_ctrl", "no transition model, inlet $\\mu_t/\\mu$ = 100", "#1f77b4", "--"),
            ("run_0148_lm_1um_mur1", "transition model, $\\mu_t/\\mu$ = 1", "#f2b134", "-"),
            ("run_0146_lm_1um_cont", "transition model, $\\mu_t/\\mu$ = 10", "#e8590c", "-"),
            ("run_0147_lm_1um_mur100", "transition model, $\\mu_t/\\mu$ = 100", "#862e9c", "-")],
    "markii": [("run_0036_sst_1um", "no transition model, inlet $\\mu_t/\\mu$ = 10", "#1f77b4", "-"),
               ("run_0037_lm_1um", "transition model, $\\mu_t/\\mu$ = 10", "#f2b134", "-"),
               ("run_0039_lm_1um_mur20", "transition model, $\\mu_t/\\mu$ = 20", "#e8590c", "-"),
               ("run_0040_lm_1um_mur40", "transition model, $\\mu_t/\\mu$ = 40", "#c92a2a", "-"),
               ("run_0038_lm_1um_mur100", "transition model, $\\mu_t/\\mu$ = 100 (pressure side still drifting)", "#862e9c", "-")],
}
TRANS_SU2 = dict(forge="run_0153_lm_cf0_su2ctrl", su2="su2_smooth_lm/vol_solution.vtu", forge_off="run_0128_cf0_fx05", Tw=566.0)


def regional_bias(run_dir, vane, key="run108"):
    """領域別の $h$ 平均偏差 [%] を `compare_h.py` と**同じ規約**で出す (最後の壁ダンプ、iface_ok、s/S<=0.87)。
    戻り値 {PS, SS_lam, SS_post, all}。plan §6.3/§6.5 とレポートの表はこの値。"""
    from compare_h import TG
    s, ss, V, C, step = wall(Path(run_dir))
    h = np.array(V[FLUX]) / (TG[key] - np.array(V["Ts"]))
    ok = np.array(V["iface_ok"]) > 0.5
    coarse = {i for i, d in getattr(__import__("run_data"), "RUN42_PARTIAL", {}).items()
              if key == "run42" and d["tol"] > 0.01}
    rows = [r for i, r in enumerate(TABLES[key]["rows"]) if r[3] is not None and i not in coarse]
    sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows]) * H0; i0 = int(np.argmin(sd))
    exp = {"PS": (sd[:i0 + 1][::-1], hd[:i0 + 1][::-1]), "SS": (sd[i0:], hd[i0:])}
    he = np.array([np.interp(s[k], *exp["SS" if ss[k] else "PS"]) for k in range(len(s))])
    out = {}
    for name, m in (("PS", ok & ~ss & (s <= 0.87)), ("SS_lam", ok & ss & (s < 0.25)),
                    ("SS_post", ok & ss & (s >= 0.25) & (s <= 0.87)), ("all", ok & (s <= 0.87))):
        out[name] = 100.0 * float(((h[m] - he[m]) / he[m]).mean())
    return out


# 掃引の点 (x 値, run)。数字は図の中で regional_bias が読み直すので、ここは run の対応表だけ。
SWEEP = {
    "c3x_ratio": ("C3X", "run108", C3, "inlet eddy viscosity ratio", "log",
                  [(1, "run_0148_lm_1um_mur1"), (10, "run_0146_lm_1um_cont"), (30, "run_0152_lm_1um_mur30_cont"), (100, "run_0147_lm_1um_mur100")],
                  [(10, "run_0135_hwall1um_fx05"), (100, "run_0149_sst_1um_mur100_ctrl")]),
    "mk_ratio": ("Mark II", "run42", MK, "inlet eddy viscosity ratio", "log",
                 [(10, "run_0037_lm_1um"), (20, "run_0039_lm_1um_mur20"), (40, "run_0040_lm_1um_mur40"), (100, "run_0038_lm_1um_mur100")],
                 [(10, "run_0036_sst_1um")]),
    "mk_reth": ("Mark II", "run42", MK, "lower bound on $\\widetilde{Re}_{\\theta t}$", "linear",
                [(20, "run_0037_lm_1um"), (130, "run_0042_lm_rt130"), (200, "run_0043_lm_rt200")], []),
}
REGIONS = [("PS", "pressure side", "#1f77b4", "o"), ("SS_lam", "suction, $s/S<0.25$", "#c92a2a", "s"),
           ("SS_post", "suction, after transition", "#2b8a3e", "^")]


def fig_trans_sweep(out, plt):
    fig, axs = plt.subplots(1, 3, figsize=(12.4, 4.6), sharey=True)
    for ax, keyname in zip(axs, ("c3x_ratio", "mk_ratio", "mk_reth")):
        name, key, base, xlab, xscale, on, off = SWEEP[keyname]
        ax.axhline(0, color="k", lw=1.0, zorder=1)
        ax.axhspan(-10, 10, color="0.85", alpha=.45, zorder=0)
        for rk, rlab, col, mk in REGIONS:
            xs = [x for x, r in on if (base / r).exists()]
            ys = [regional_bias(base / r, name, key)[rk] for x, r in on if (base / r).exists()]
            ax.plot(xs, ys, "-", marker=mk, ms=6, color=col, lw=1.8, zorder=4,
                    label=(rlab if keyname == "c3x_ratio" else None))
            xo = [x for x, r in off if (base / r).exists()]
            yo = [regional_bias(base / r, name, key)[rk] for x, r in off if (base / r).exists()]
            if xo:
                ax.plot(xo, yo, linestyle="none", marker=mk, ms=7, mfc="none", mec=col, mew=1.4, zorder=3,
                        label=("same, no transition model" if (keyname == "c3x_ratio" and rk == "PS") else None))
        ax.set_xscale(xscale); ax.grid(alpha=.3)
        ax.set_xlabel(xlab)
        ax.set_title(f"{name} — {'inlet turbulence decay' if xscale == 'log' else 'the model constant'}", fontsize=10)
        if xscale == "log":
            ax.set_xticks([x for x, _ in on]); ax.set_xticklabels([str(x) for x, _ in on])
        else:
            ax.set_xticks([x for x, _ in on])
    axs[0].set_ylabel("mean departure of $h$ from the measurement  [%]")
    axs[0].set_ylim(-45, 90)
    h1, l1 = axs[0].get_legend_handles_labels()
    fig.legend(h1, l1, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=115, bbox_inches="tight"); print(f"[report_figures] -> {out}")


def fig_trans_plate(out, plt):
    sys.path.insert(0, str(T57 / "tools")); import cf_plate as cp
    fig, a = plt.subplots(figsize=(9.2, 5.4))
    xx = np.linspace(0.004, 1.5, 400); rex = cp.RO * cp.U * xx / cp.MU
    a.plot(xx, 0.664 / np.sqrt(rex) * 1e3, ":", color="0.35", lw=1.1, label="Blasius, laminar")
    a.plot(xx, 0.0576 * rex ** -0.2 * 1e3, "--", color="0.35", lw=1.1, label="$0.0576\\,Re_x^{-1/5}$, turbulent")
    e = np.genfromtxt(T57 / "ref/t3a_exp.dat", comments="#")
    a.plot(e[:, 0] * 1e-3, e[:, 1] * 1e3, "o", ms=5.5, mfc="w", mec="k", zorder=6, label="ERCOFTAC T3A, measured")
    for run, lab, col, ls, lw in (("run_0003_t3a_sst", "forge, no transition model", "#1f77b4", "-", 1.5), ("run_0010_t3a_lam", "forge, laminar", "#2b8a3e", "-", 1.5),
                                  ("run_0007_t3a_lm_coarse", "forge + transition model, coarse grid", "#ffa94d", "-", 1.3), ("run_0009_t3a_lm_fine", "forge + transition model, fine grid", "#c92a2a", "-", 1.3),
                                  ("run_0005_t3a_lm_cont", "forge + transition model, base grid", "k", "-", 2.0)):
        (x, tw, _), _f = cp.load_forge(str(T57 / run)); a.plot(x, tw / cp.QINF * 1e3, ls, color=col, lw=lw, label=lab, zorder=(5 if col == "k" else 3))
    (x, tw, _), _f = cp.load_su2(str(T57 / "su2_t3a_lm")); a.plot(x, tw / cp.QINF * 1e3, "-.", color="tab:green", lw=1.8, label="SU2 8.5 + same model, base grid", zorder=4)
    a.set_xlim(0, 1.5); a.set_ylim(0, 8); a.set_xlabel("distance from the leading edge  $x$  [m]   ($Re_x$ = 3.6e5 · $x$)"); a.set_ylabel("$C_f \\times 10^3$"); a.grid(alpha=.3)
    a.set_title("Flat plate, free-stream turbulence 3.3 % (ERCOFTAC T3A) — skin friction\n"
                "run_0005 / 0007 / 0009 (transition model) . run_0003 (none) . run_0010 (laminar) . su2_t3a_lm", fontsize=10)
    a.legend(fontsize=8.6, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=115, bbox_inches="tight"); print(f"[report_figures] -> {out}")


def fig_trans(vane, out, plt):
    R = RUNS[vane]
    fig, a = plt.subplots(figsize=(9.6, 5.8))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    plot_measured(a, vane, R["key"], R["table"])
    names = []
    for run, lab, col, ls in TRANS[vane]:
        if not (R["base"] / run).exists(): print(f"[report_figures] skip (no run yet): {run}"); continue
        s, ss, V, _, step = wall(R["base"] / run)
        h = np.array(V[FLUX]) / (R["Tg"] - np.array(V["Ts"])) / H0
        sides(a, s, ss, h, color=col, ls=ls, lw=1.6, alpha=.95, zorder=4, label=lab); names.append(run)
    a.set_ylabel("$h/h_0$"); a.set_ylim(0, 1.15 if vane == "c3x" else 1.6)
    finish(a, fig, out, f"{R['name']} — with and without the transition model, measured $T_w$ imposed, same mesh\n" + " . ".join(names))


def fig_trans_su2(out, plt):
    from scipy.spatial import cKDTree
    s, ss, V, C, step = wall(C3 / TRANS_SU2["forge"])
    XY, QS = read_su2(C3 / TRANS_SU2["su2"])
    d, idx = cKDTree(XY).query(C)
    assert d.max() == 0.0, f"forge と SU2 の壁節点が一致しない (max {d.max()})"
    sc = 1.0 / ((RUNS["c3x"]["Tg"] - TRANS_SU2["Tw"]) * H0)
    fig, a = plt.subplots(figsize=(9.2, 5.6))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    plot_measured(a, "c3x", "run108", "VI")
    s0, ss0, V0, _, _ = wall(C3 / TRANS_SU2["forge_off"])
    sides(a, s0, ss0, np.array(V0[FLUX]) * sc, color="#9aa0a6", lw=1.3, zorder=2, label="forge, no transition model")
    sides(a, s, ss, np.array(V[FLUX]) * sc, color="k", lw=2.0, zorder=4, label="forge + transition model")
    sides(a, s, ss, QS[idx] * sc, color="tab:green", lw=1.7, zorder=3, label="SU2 8.5 + same model, same mesh")
    a.set_ylabel("$h/h_0$"); a.set_ylim(0, 1.15)
    finish(a, fig, out, "C3X run 108 — the transition model in forge and in SU2 on the identical mesh, uniform $T_w$ = 566 K\n"
           f"{TRANS_SU2['forge']}, step {step} . su2_smooth_lm  (wall nodes coincide exactly)")


def fig_turb(out, plt):
    fig, a = plt.subplots(figsize=(9.6, 5.8))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    plot_measured(a, "c3x", "run108", "VI")
    for run, lab, col, ls in TURB:
        s, ss, V, _, _ = wall(C3 / run)
        h = np.array(V[FLUX]) / (786.0 - SU2_CTRL["Tw"]) / H0
        sides(a, s, ss, h, color=col, ls=ls, lw=(2.0 if col == "k" else 1.5), alpha=.92, zorder=4, label=lab)
    a.set_ylabel("$h/h_0$")
    finish(a, fig, out, "C3X run 108 — what the inlet turbulence can and cannot move\n"
           "uniform $T_w$=566 K, cs2000 mesh; only the inlet $k$ and $\\omega$ differ")


def fig_lam(vane, out, plt):
    R = RUNS[vane]
    fig, a = plt.subplots(figsize=(9.4, 5.6))
    a.axvspan(0, 0.25, color="#f0d9a8", alpha=.35, zorder=0)
    sd, hd, _ = plot_measured(a, vane, R["key"], R["table"])
    for run, lab, col in ((R["prod"], "turbulent (SST)", "k"), (R["laminar"], "laminar (model off)", "tab:purple")):
        s, ss, V, _, _ = wall(R["base"] / run)
        h = np.array(V[FLUX]) / (R["Tg"] - np.array(V["Ts"])) / H0
        # 層流解は遷移後の負圧面で剥離して非定常になる。**その区間と後縁近傍は描かない**
        mask = (lambda side: s < (0.25 if side == "SS" else 0.90)) if col != "k" else None
        sides(a, s, ss, h, mask=mask, color=col, lw=1.9, alpha=.92, zorder=4, label=lab)
    a.set_ylim(0, max(1.15, hd.max() * 1.35)); a.set_ylabel("$h/h_0$")
    finish(a, fig, out, f"{R['name']} — the measurement sits between the two limits\n"
           f"{R['prod']}  vs  {R['laminar']}")


FIGS = {"h_c3x": lambda o, p: fig_h("c3x", o, p), "h_mk": lambda o, p: fig_h("markii", o, p),
        "press_c3x": lambda o, p: fig_press("c3x", o, p), "press_mk": lambda o, p: fig_press("markii", o, p),
        "su2": fig_su2, "turb": fig_turb,
        "lam_c3x": lambda o, p: fig_lam("c3x", o, p), "lam_mk": lambda o, p: fig_lam("markii", o, p),
        "trans_plate": fig_trans_plate, "trans_c3x": lambda o, p: fig_trans("c3x", o, p),
        "trans_mk": lambda o, p: fig_trans("markii", o, p), "trans_su2": fig_trans_su2,
        "trans_sweep": fig_trans_sweep}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", nargs="*", default=None, choices=list(FIGS))
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for name in (a.only or FIGS):
        FIGS[name](out / f"{name}.png", plt)


if __name__ == "__main__":
    main()
