"""ノズル設計の標準出力 (図・数値・条件表) を 1 つの NS run から作る。

規約の正本は procedures/nozzle-design-outputs.md (2026-10-05 ユーザ指示でルール化)。
出力 (既定 <run>/report/):
  fig_*.png        解析領域・コンタ・線グラフ・壁形状
  report.json      条件 (境界条件・数値設定・物性・ゲート) と評価量。build_pptx.py が読む
  <run>_report.pptx  build_pptx.py が `.venv-pptx` の python-pptx で作る (無ければ図と json だけ)

usage:
  design/.venv-opt/bin/python -m forge_design.report.nozzle_report RUN_DIR [--euler EULER_RUN] [--out DIR] [--no-pptx]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
from scipy.interpolate import RegularGridInterpolator, make_interp_spline  # noqa: E402

_FONT = "/home/sano/.fonts/NotoSansCJKjp-Regular.otf"
if os.path.exists(_FONT):
    fm.fontManager.addfont(_FONT)
    matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path(__file__).resolve().parents[3]
TOOLS = REPO / "solver_density_cuda" / "tools"
PPTX_PY = REPO / ".venv-pptx" / "bin" / "python"
CMAP = "turbo"          # skill forge-contour: 全量 turbo


# ------------------------------------------------------------------ 入力
def _res_files(run):
    fs = [f for f in os.listdir(run) if re.fullmatch(r"res_\d+\.h5", f) and f != "res_0.h5"]
    return sorted(fs, key=lambda f: int(re.findall(r"\d+", f)[0]))


def _yaml_flow(path):
    """forge の config は flow style の YAML。そのまま読む。"""
    return yaml.safe_load(Path(path).read_text())


def load_field(run, res=None):
    run = Path(run)
    info = json.loads((run / "prepare_info.json").read_text())
    S = float(info["scale_m"]); ni = int(info["mesh"]["ni"])
    with h5py.File(run / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    nj = nc.shape[0] // ni
    res = res or _res_files(run)[-1]
    V = {}
    with h5py.File(run / res) as f:
        for k in f["/VALUE"].keys():
            ds = f["/VALUE/" + k]
            if ds.shape and ds.shape[0] == ni * nj and ds.ndim == 1:
                V[k] = ds[:].astype(float).reshape(ni, nj)
    X = (nc[:, 0] / S).reshape(ni, nj); R = (nc[:, 1] / S).reshape(ni, nj)
    V["M"] = np.hypot(V["Ux"], V["Uy"]) / V["sonic"]
    V["U2"] = V["Ux"] ** 2 + V["Uy"] ** 2
    V["q_dyn"] = 0.5 * V["ro"] * V["U2"]
    V["ke"] = 0.5 * V["U2"]
    # 勾配: 構造格子の座標変換 (x, r は r_t 単位)
    xi, xj = np.gradient(X); ri, rj = np.gradient(R); J = xi * rj - xj * ri

    def ddx(F):
        fi, fj = np.gradient(F); return (fi * rj - fj * ri) / J

    def ddr(F):
        fi, fj = np.gradient(F); return (-fi * xj + fj * xi) / J

    V["dpdx_nd"] = ddx(V["P"]) / V["P"]
    gx, gr = ddx(V["ro"]), ddr(V["ro"])
    V["schlieren"] = np.log10(np.maximum(np.hypot(gx, gr) / V["ro"], 1e-8))
    if "g_0" in V and "condS_0" not in V:
        try:
            _saturation(run, V)
        except BaseException as e:  # noqa: BLE001
            print(f"[nozzle_report] 飽和量を作れない: {e}")
    if "condTsat_0" in V:
        V["dT_sub"] = np.where(V["condTsat_0"] > 0, V["condTsat_0"] - V["T"], np.nan)
    return dict(X=X, R=R, V=V, S=S, ni=ni, nj=nj, info=info, res=res)


# ---- 飽和量 (凝縮 ON で condS_0 / condTsat_0 が出力に無いとき)。solver_density_cuda/tools/paraview/forge_filters.py の
#      _h2o_psat_liquid / _tsat_newton / H2O_RV の写し (あちらは ParaView を import するのでここでは読めない)。
#      forge の condensationProperties_d.cuh (Murphy & Koop 2005 過冷却液) と同じ式。
H2O_RV = 461.5


def _h2o_psat_liquid(T):
    Tc = np.maximum(T, 120.0)
    lnp = (54.842763 - 6763.22 / Tc - 4.210 * np.log(Tc) + 0.000367 * Tc
           + np.tanh(0.0415 * (Tc - 218.8)) * (53.878 - 1331.22 / Tc - 9.44523 * np.log(Tc) + 0.014025 * Tc))
    return np.exp(lnp)


def _tsat_newton(psat_fn, pv, T_guess):
    pv = np.asarray(pv, dtype=np.float64)
    valid = pv > 1.0e-6
    lnpv = np.log(np.where(valid, pv, 1.0))
    T = np.where((T_guess > 50.0) & (T_guess < 1000.0), T_guess, 250.0).astype(np.float64)
    h = 0.01
    for _ in range(25):
        f = np.log(np.maximum(psat_fn(T), 1.0e-300)) - lnpv
        dfdT = (np.log(np.maximum(psat_fn(T + h), 1.0e-300)) - np.log(np.maximum(psat_fn(T - h), 1.0e-300))) / (2 * h)
        T = np.clip(T - np.clip(f / np.maximum(dfdT, 1.0e-6), -0.3 * T, 0.3 * T), 50.0, 1000.0)
    return np.where(valid, T, 0.0)


def _saturation(run, V):
    """H2O 凝縮の S と T_sat を T・ρ・Y・g から作る (forge cond_vapor_state と同じ p_v = ρ (Y_v − g) R_v T)。"""
    sys.path.insert(0, str(TOOLS))
    from forge_species import species_info
    yv_name = species_info(str(run))["vapor_array"]
    pv = V["ro"] * np.maximum(V[yv_name] - V.get("g_0", 0.0), 0.0) * H2O_RV * V["T"]
    V["condS_0"] = pv / np.maximum(_h2o_psat_liquid(V["T"]), 1.0e-300)
    V["condTsat_0"] = _tsat_newton(_h2o_psat_liquid, pv, V["T"])
    V["_saturation_source"] = "post"


def eta_line(F, key, eta, xq):
    X, R, Q = F["X"], F["R"], F["V"][key]
    xa = X[:, 0]
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2); w = (xq - xa[i]) / (xa[i + 1] - xa[i])
    out = np.empty(len(xq))
    for k in range(len(xq)):
        v = [np.interp(eta, R[ii] / R[ii, -1], Q[ii]) for ii in (i[k], i[k] + 1)]
        out[k] = (1 - w[k]) * v[0] + w[k] * v[1]
    return out


def pspline(x, v, knot=10.0, lam=1.0):
    from forge_design.metrics.deltastar import pspline_uniform
    return pspline_uniform(x, v, knot=knot, lam=lam)


# ------------------------------------------------------------------ 条件 (境界条件・数値設定・物性・ゲート)
def conditions(run, F):
    run = Path(run)
    cfg = _yaml_flow(run / "solverConfig.yaml"); bc = _yaml_flow(run / "bcondConfig.yaml")
    info = F["info"]
    pp = cfg.get("physProp", {}); tm = cfg.get("time", {}); dT = tm.get("deltaT", {})
    turb = cfg.get("turbulence", {}); sp = cfg.get("space", {}); mesh = cfg.get("mesh", {})
    species = pp.get("species", [])
    sp_txt = ", ".join(s if isinstance(s, str) else s.get("name", "?") for s in species) or "—"
    gas = {0: "熱量的完全気体 (CPG)", 2: "熱的完全気体 (NASA-9、比熱は温度依存)"}.get(pp.get("thermalMethod"), f"thermalMethod {pp.get('thermalMethod')}")
    bcrows = []
    for name, b in bc.items():
        fl = b.get("floats") or {}
        val = ", ".join(f"{k} {v:g}" if isinstance(v, (int, float)) else f"{k} {v}" for k, v in fl.items()) or "—"
        bcrows.append([str(b.get("physID")), name, str(b.get("kind")), val])
    num = [
        ["離散化", f"{mesh.get('discretization', '?')} 中心、軸対称 {'あり' if mesh.get('isAxisymmetric') else 'なし'}"],
        ["対流スキーム", f"{cfg.get('solver', '?')}、convMethod {sp.get('convMethod')} (2 次)、limiter {sp.get('limiter')}" + (f"、limiterScaled {sp['limiterScaled']}" if 'limiterScaled' in sp else "")],
        ["時間積分", f"定常 (擬似時間)、timeIntegration {tm.get('timeIntegration')}、blockDPLUR {dT.get('blockDPLUR')}、CFL {dT.get('cfl')}、implicitRelax {dT.get('implicitRelax', '—')}"],
        ["step 数", f"{tm.get('last', {}).get('nStepOuter')} (出力間隔 {tm.get('outStepInterval')})"],
        ["乱流モデル", f"{turb.get('model', 'なし')}" + (f" (dilatationCorrection {turb.get('dilatationCorrection')}, KL {turb.get('katoLaunder')}, 壁処理 {'低 Re' if turb.get('wallTreatmentSST', 0) == 0 else '壁関数'}, Pr_t {turb.get('turbulentPrandtl')})" if turb.get('model') not in (None, 'none') else "")],
        ["格子", f"{F['ni']} × {F['nj']} (x × r)、第 1 セル {info['mesh'].get('wall_first_frac')} r_w"],
        ["初期場", str(info.get("ic_from") or "—")],
        ["壁の δ_r", str(info.get("dstar_source", "—"))[:120]],
    ]
    if "condensation" in cfg:
        c = cfg["condensation"]
        num.append(["凝縮", ", ".join(f"{k} {v}" for k, v in c.items())])
    phys = [
        ["気体", f"{gas}; 化学種 {sp_txt}"],
        ["粘性", {0: "定数", 1: "Sutherland"}.get(pp.get("viscMethod"), str(pp.get("viscMethod")))],
        ["熱伝導", f"Pr {pp.get('prandtlLam', '—')} 一定 (thermCondMethod {pp.get('thermCondMethod', '—')})"],
        ["スロート半径 r_t", f"{F['S'] * 1e3:.2f} mm"],
        ["設計マッハ数", f"{info.get('Md', '—')}"],
    ]
    gates = []
    cv = run / "CONVERGENCE_VERDICT.txt"
    if cv.exists():
        t = cv.read_text().splitlines()
        overall = next((l.strip() for l in t if "OVERALL" in l), "?")
        worst = [l.strip() for l in t if l.strip().startswith("rms_") and "<--" in l][:3]
        gates.append(["収束 (check_convergence)", overall + (" — " + "; ".join(worst) if worst else "")])
    mq = run / "MESH_QUALITY.txt"
    if mq.exists():
        gates.append(["メッシュ品質 (check_mesh_quality)", next((l.strip() for l in mq.read_text().splitlines()[::-1] if "VERDICT" in l), "?")])
    nonfin = int(sum((~np.isfinite(F["V"][k])).sum() for k in ("ro", "P", "T", "Ux", "Uy")))
    gates.append(["NaN / Inf (最終場)", str(nonfin)])
    return dict(bc=bcrows, numerics=num, physics=phys, gates=gates,
                run=str(run), res=F["res"], case=run.parent.name, run_name=run.name)


# ------------------------------------------------------------------ 評価量
def metrics(run, F, euler=None):
    info = F["info"]; Md = float(info.get("Md", 6.0))
    xE = float(info.get("x_E", 40.0)); xF = float(F["X"][-1, 0])
    xq = np.linspace(float(F["X"][0, 0]), xF, 2401)
    out = dict(Md=Md, x_E=xE, x_F=xF)
    w_test = (xq >= xE + 2) & (xq <= xF - 1); w_ov = (xq >= xE - 15) & (xq <= xF)
    for eta in (0.0, 0.1):
        d = 100 * (eta_line(F, "M", eta, xq) / Md - 1)
        xx, v = xq[w_test], d[w_test]
        cf = np.polyfit(xx, v, 1)
        out[f"eta{eta}"] = dict(wave_pct=float(np.abs(v - pspline(xx, v)).max()), overshoot_pct=float(d[w_ov].max()),
                                x_overshoot=float(xq[w_ov][np.argmax(d[w_ov])]), slope_pct=float(cf[0] * (xx[-1] - xx[0])),
                                range_pct=float(v.max() - v.min()), test_window=[float(xx[0]), float(xx[-1])])
    eta_last = F["R"][-1] / F["R"][-1, -1]; core = (eta_last >= 0.05) & (eta_last <= 0.7)
    Mx = F["V"]["M"][-1][core]
    out["exit_core"] = dict(M_mean=float(Mx.mean()), M_min=float(Mx.min()), M_max=float(Mx.max()), eta_range=[0.05, 0.7])
    if euler is not None:
        try:
            from forge_design.metrics.deltastar import massflow_ratio
            # 無次元 (r_t 単位) の流量 2π∫ρU_x r dr を x∈(−2.5, x_E) の中央値で比べる (r_t が違う run 同士でも比べられる)。
            # metrics.deltastar.massflow_ratio は実寸化で r_t を 3 乗で掛けており (正しくは 2 乗)、r_t が同じ run 同士でしか比が正しくない。
            from forge_design.metrics.deltastar import _load_structured
            md = []
            for rd in (run, euler):
                d_ = _load_structured(rd); xs_ = d_["x"][:, 0]
                m_ = 2 * np.pi * np.trapezoid(d_["q"] * d_["r"], d_["r"], axis=1)
                md.append(float(np.median(m_[(xs_ > -2.5) & (xs_ < d_["info"]["x_E"])])))
            out["mdot_ratio_vs_euler"] = md[0] / md[1]
        except Exception as e:  # noqa: BLE001
            out["mdot_ratio_vs_euler"] = f"未計算 ({e})"
    if "g_0" in F["V"]:
        g0 = eta_line(F, "g_0", 0.0, xq)
        on = np.where(g0 > 1e-4)[0]
        out["condensation"] = dict(onset_x_axis=float(xq[on[0]]) if len(on) else None,
                                   exit_g_axis=float(g0[-1]), onset_threshold_g=1e-4, exit_g_core_mean=float(F["V"]["g_0"][-1][core].mean()),
                                   S_max=float(np.nanmax(F["V"]["condS_0"])) if "condS_0" in F["V"] else None)
    # 準定常: 末尾 5 スナップショットの変動
    ser = []
    for r in _res_files(run)[-5:]:
        G = load_field(run, r)
        d = 100 * (eta_line(G, "M", 0.1, xq) / Md - 1)
        xx, v = xq[w_test], d[w_test]
        ser.append((d[w_ov].max(), np.abs(v - pspline(xx, v)).max(), G["V"]["M"][-1][core].mean()))
    a = np.array(ser)
    out["tail5_range"] = dict(overshoot_pct=float(np.ptp(a[:, 0])), wave_pct=float(np.ptp(a[:, 1])), exit_core_M=float(np.ptp(a[:, 2])), n=len(ser))
    return out


# ------------------------------------------------------------------ 図
def fig_domain(run, F, path):
    run = Path(run); X, R = F["X"], F["R"]
    bc = _yaml_flow(run / "bcondConfig.yaml")
    nc = np.c_[X.ravel(), R.ravel()]
    cols = ["#1f77b4", "#9467bd", "#d62728", "#2ca02c", "#ff7f0e", "#8c564b"]
    fig, axs = plt.subplots(2, 1, figsize=(12, 6.2), gridspec_kw=dict(height_ratios=[1.4, 1]))
    with h5py.File(run / "nozzle.h5") as f:
        for ax, xl, yl, si, sj in ((axs[0], (X.min() - 0.5, X.max() + 0.5), (-0.4, R.max() * 1.05), 10, 4), (axs[1], (-3, 6), (-0.1, 2.0), 1, 2)):
            m = (X[:, 0] >= xl[0] - 1) & (X[:, 0] <= xl[1] + 1)
            for i in np.where(m)[0][::si]:
                ax.plot(X[i], R[i], color="0.82", lw=0.3)
            for j in range(0, X.shape[1], sj):
                ax.plot(X[m, j], R[m, j], color="0.82", lw=0.3)
            for k, (name, b) in enumerate(bc.items()):
                pid = str(b.get("physID"))
                if pid not in f["BCONDS"]:
                    continue
                nodes = np.unique(f["BCONDS"][pid]["vizBfaceNodes"][:].ravel())
                pts = nc[nodes]; o = np.argsort(pts[:, 0] if b.get("kind") in ("wall", "axis") else pts[:, 1])
                fl = b.get("floats") or {}
                lab = f"{name} (physID {pid}): {b.get('kind')}" + (" — " + ", ".join(f"{a} {v:g}" for a, v in list(fl.items())[:3]) if fl else "")
                ax.plot(pts[o, 0], pts[o, 1], color=cols[k % len(cols)], lw=2.2, label=lab if ax is axs[0] else None)
            ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.set_xlabel("x / r_t"); ax.set_ylabel("r / r_t")
    axs[0].set_title(f"解析領域 (軸対称の上半分; 格子線は 10×4 本に 1 本表示、実格子 {F['ni']}×{F['nj']})", fontsize=10)
    axs[1].set_title("スロート近傍の拡大 (x = −3〜6 r_t)", fontsize=10)
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(path, dpi=140); plt.close(fig)


def _contour(ax, F, key, title, vmin, vmax, xlim, unit=""):
    cs = ax.pcolormesh(F["X"], F["R"], F["V"][key] if isinstance(key, str) else key, cmap=CMAP, vmin=vmin, vmax=vmax,
                       shading="gouraud", rasterized=True)
    ax.set_xlim(*xlim); ax.set_ylim(0, F["R"].max() * 1.02); ax.set_aspect("equal"); ax.set_ylabel("r / r_t")
    ax.set_title(title, fontsize=9.5)
    cb = plt.colorbar(cs, ax=ax, orientation="vertical", fraction=0.025, pad=0.01); cb.set_label(unit, fontsize=8)


def _pct(a, lo=0.2, hi=99.8):
    a = a[np.isfinite(a)]
    return float(np.percentile(a, lo)), float(np.percentile(a, hi))


def fig_contours(F, Md, out_dir):
    xl = (float(F["X"].min()), float(F["X"].max())); xt = (max(0.0, xl[0]), xl[1])
    V = F["V"]; figs = []
    groups = [
        ("fig_contour_mach.png", [("M", f"マッハ数 (最大 {np.nanmax(V['M']):.4f})", 0.0, float(np.nanmax(V['M'])), xl, "M"),
                                   (100 * (V["M"] / Md - 1), f"M / M_d − 1 [%] (M_d = {Md}; ±0.5 % で頭打ち。境界層と膨張部は下限に張り付く)", -0.5, 0.5, xt, "%")]),
        ("fig_contour_gradients.png", [("dpdx_nd", "無次元の軸方向圧力勾配 (r_t/p) ∂p/∂x (±0.02 で頭打ち)", -0.02, 0.02, xt, "—"),
                                        ("schlieren", "数値シュリーレン log10(|∇ρ| r_t / ρ)", *_pct(V["schlieren"], 5, 99.5), xl, "log10")]),
        ("fig_contour_PT.png", [(V["P"] / 1e3, "静圧 [kPa] (対数目盛ではない; 0.2–99.8 % で頭打ち)", *_pct(V["P"] / 1e3), xl, "kPa"),
                                (V["T"], "静温 [K]", *_pct(V["T"]), xl, "K")]),
    ]
    if "g_0" in V:
        groups.append(("fig_contour_condensation.png", [("g_0", "液滴の質量分率 g", 0.0, _pct(V["g_0"], 0, 99.9)[1], xt, "—"),
                                                        ("dT_sub", "過冷却度 T_sat − T [K] (正 = 過冷却。過冷却側に合わせて頭打ち)",
                                                         -0.25 * max(_pct(V["dT_sub"], 0.5, 99.5)[1], 1.0), max(_pct(V["dT_sub"], 0.5, 99.5)[1], 1.0), xt, "K"),
                                                        ("condS_0", "過飽和度 S", *_pct(V["condS_0"], 0.5, 99.5), xt, "—")]))
    for name, panels in groups:
        fig, axs = plt.subplots(len(panels), 1, figsize=(12, 3.3 * len(panels)), constrained_layout=True)
        for ax, (key, title, vmin, vmax, xlim, unit) in zip(np.atleast_1d(axs), panels):
            _contour(ax, F, key, title, vmin, vmax, xlim, unit)
        np.atleast_1d(axs)[-1].set_xlabel("x / r_t")
        fig.savefig(out_dir / name, dpi=130); plt.close(fig); figs.append(name)
    return figs


def _line_quantities(V):
    q = [("M", "M", 1.0), ("P", "静圧 [kPa]", 1e-3), ("T", "静温 [K]", 1.0), ("ro", "密度 [kg/m³]", 1.0),
         ("q_dyn", "動圧 0.5ρv² [kPa]", 1e-3), ("ke", "0.5 v² [kJ/kg]", 1e-3), ("P0", "全圧 [kPa]", 1e-3)]
    if "g_0" in V:
        q += [("g_0", "液滴の質量分率 g", 1.0), ("dT_sub", "過冷却度 T_sat − T [K]", 1.0), ("condS_0", "過飽和度 S", 1.0)]
    return [x for x in q if x[0] in V]


def fig_axis_lines(F, Md, path, path_dev, euler=None):
    xq = np.linspace(float(F["X"][0, 0]), float(F["X"][-1, 0]), 2401)
    qs = _line_quantities(F["V"]); n = len(qs); nc = 2; nr = (n + 1) // 2
    fig, axs = plt.subplots(nr, nc, figsize=(13, 2.6 * nr), sharex=True)
    for ax, (k, lab, sc) in zip(axs.ravel(), qs):
        for eta, ls in ((0.0, "-"), (0.1, "--")):
            ax.plot(xq, sc * eta_line(F, k, eta, xq), ls=ls, lw=1.1, label=f"r/r_w = {eta}")
        ax.set_ylabel(lab, fontsize=9)
        if k in ("P", "ro", "q_dyn"):
            ax.set_yscale("log")
        if k == "P0":
            ax.set_title("軸ノード (r=0) は既知の軸の値の癖 (静圧が 0.5〜0.7 % 低い) で全圧が高めに出る", fontsize=8)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    axs.ravel()[0].legend(fontsize=8, frameon=False)
    for ax in axs[-1]:
        ax.set_xlabel("x / r_t")
    fig.suptitle("軸中心 (r = 0、実線) と r/r_w = 0.1 (破線) の分布。r/r_w = 0.1 は軸ノード固有の値の影響を避けた線", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    # 試験部の M/M_d − 1 の拡大 (評価量の図)
    fig, ax = plt.subplots(figsize=(12, 4.2))
    for eta, ls in ((0.0, "-"), (0.1, "--")):
        ax.plot(xq, 100 * (eta_line(F, "M", eta, xq) / Md - 1), ls=ls, lw=1.2, label=f"NS r/r_w = {eta}")
    if euler is not None:
        for eta, ls in ((0.0, "-"), (0.1, "--")):
            ax.plot(xq, 100 * (eta_line(euler, "M", eta, xq) / Md - 1), ls=ls, color="k", lw=0.9, alpha=0.6, label=f"Euler (設計壁) r/r_w = {eta}")
    ax.axhline(0, color="0.6", lw=0.5); ax.set_ylim(-0.6, 0.6); ax.set_xlim(max(10, xq[0]), xq[-1])
    ax.set_xlabel("x / r_t"); ax.set_ylabel("M / M_d − 1 [%]"); ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.set_title("試験部の M の設計値からのずれ (圧力波・オーバーシュート・傾きを見る図)", fontsize=10)
    fig.tight_layout(); fig.savefig(path_dev, dpi=130); plt.close(fig)


def fig_exit_lines(F, path):
    eta = F["R"][-1] / F["R"][-1, -1]; V = F["V"]
    qs = _line_quantities(V) + [("flow_angle", "流れ角 atan(v/u) [deg]", 1.0)]
    V["flow_angle"] = np.degrees(np.arctan2(V["Uy"], V["Ux"]))
    n = len(qs); nr = (n + 2) // 3
    fig, axs = plt.subplots(nr, 3, figsize=(13, 2.8 * nr))
    for ax, (k, lab, sc) in zip(axs.ravel(), qs):
        ax.plot(eta, sc * V[k][-1], lw=1.2); ax.set_ylabel(lab, fontsize=9); ax.set_xlim(0, 1)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    for ax in axs[-1]:
        ax.set_xlabel("r / r_w (出口断面)")
    fig.suptitle(f"出口断面 (x = {F['X'][-1, 0]:.2f} r_t) の半径方向分布", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def fig_wall_lines(F, path):
    """壁面: 壁節点 (j=-1) と第 1 内部節点から。y1+ は第 1 内部節点までの壁法線距離と接線壁応力 (AGENTS.md: ソルバの ypls は使わない)。"""
    X, R, V = F["X"], F["R"], F["V"]; S = F["S"]
    x = X[:, 0]
    th = np.arctan2(np.gradient(R[:, -1]), np.gradient(X[:, -1]))
    ut = V["Ux"] * np.cos(th)[:, None] + V["Uy"] * np.sin(th)[:, None]
    yn = (R[:, -1] - R[:, -2]) * np.cos(th) * S
    mu = V.get("vis_lam", np.full_like(V["P"], np.nan))
    tw = mu[:, -1] * ut[:, -2] / np.maximum(yn, 1e-30)
    je = np.argmin(np.abs(R / R[:, -1:] - 0.8), axis=1)
    qe = np.array([V["q_dyn"][i, j] for i, j in enumerate(je)])
    Cf = tw / qe
    y1p = yn * np.sqrt(V["ro"][:, -1] * np.abs(tw)) / mu[:, -1]
    # この y₁⁺ は速度差からの近似で**図示だけ**に使う。判定・報告の数値は正式ツール check_wall_resolution.py
    # (wall_resolution(); plan tooling-nozzle-cfd-pinned-initial-line §5.1 #10 — codex result M4: 近似は x>0 の節点割合で最大 8.7/18 % と出ていた)
    rows = [(V["P"][:, -1] / 1e3, "壁圧 [kPa]"), (V["T"][:, -1], "壁温 [K]"), (Cf * 1e3, "C_f × 10³ (接線壁応力 / 動圧@r/r_w=0.8)"),
            (y1p, "y₁⁺ 近似 (図示のみ; 判定は check_wall_resolution.py)")]
    fig, axs = plt.subplots(2, 2, figsize=(13, 6), sharex=True)
    for ax, (v, lab) in zip(axs.ravel(), rows):
        ax.plot(x, v, lw=1.1); ax.set_ylabel(lab, fontsize=9)
    axs[1, 1].axhline(1.0, color="0.5", lw=0.6, ls="--")
    mC = x > 0.5
    axs[1, 0].set_ylim(0, 1.3 * float(np.nanmax(Cf[mC] * 1e3)))
    axs[0, 0].set_yscale("log")
    for ax in axs[-1]:
        ax.set_xlabel("x / r_t")
    fig.suptitle("壁面の分布", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    m = x > 0
    return dict(y1p_approx_note="速度差からの近似 (x>0 の節点、図示用)。判定には metrics.wall_resolution (正式ツール) を使う",
                y1p_approx_max=float(np.nanmax(y1p[m])), y1p_approx_x_max=float(x[m][np.nanargmax(y1p[m])]))


def wall_resolution(run, F, over_frac=None):
    """壁解像の正式値: `solver_density_cuda/tools/check_wall_resolution.py` を no-slip 壁全部 (bcondConfig の kind が wall*) で
    実行し、出力 (壁ごとの y₁⁺ 平均・p99・最大・評価面積・目標超過面積・最大の位置・VERDICT) を解析する。
    最大の位置は壁ダンプ `res_<群>_<physID>_<step>.h5` の座標 (index) から x/r_t, r/r_t に直す。
    over_frac: 目標超過を許す面積割合 [%] (None = ツール既定)。"""
    run = Path(run)
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text()) or {}
    walls = [k for k, v in bc.items() if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")]
    if not walls:
        return {"verdict": "INDETERMINATE", "reason": "no-slip 壁が無い"}
    cmd = [sys.executable, str(TOOLS / "check_wall_resolution.py"), str(run), "--groups", ",".join(walls)]
    if over_frac is not None:
        cmd += ["--over-frac", str(float(over_frac))]
    r = subprocess.run(cmd, capture_output=True, text=True)
    txt = r.stdout + r.stderr
    out = {"tool": "solver_density_cuda/tools/check_wall_resolution.py", "cmd": " ".join(cmd[1:]), "returncode": r.returncode,
           "groups": walls, "per_wall": {}}
    m = re.search(r"^VERDICT: (\S+)(.*)$", txt, re.M)
    out["verdict"] = m.group(1) if m else "INDETERMINATE"
    out["verdict_line"] = m.group(0).strip() if m else None
    m = re.search(r"目標 ([0-9.eE+-]+) を超える面積 最大 ([0-9.]+) % \(許容 ([0-9.eE+-]+) %\)", txt)
    if m:
        out.update(target=float(m.group(1)), over_area_pct=float(m.group(2)), over_area_allow_pct=float(m.group(3)))
    m = re.search(r"最大 y1\+ = ([0-9.]+)", txt)
    if m:
        out["y1p_max"] = float(m.group(1))
    S = F["S"]
    for mm in re.finditer(r"^  (\S+)\s+step\s+(\d+)\s+y1 = ([0-9.e+-]+) m\s+y1\+ 平均\s+([0-9.]+) / p99\s+([0-9.]+) / 最大\s+([0-9.]+)\n"
                          r"\s+評価できた面積割合 ([0-9.]+) % ; y1\+ > \S+ が ([0-9.]+) % ; 最大の位置 index (\d+)", txt, re.M):
        name, step = mm.group(1), int(mm.group(2))
        w = dict(step=step, y1_median_m=float(mm.group(3)), y1p_mean=float(mm.group(4)), y1p_p99=float(mm.group(5)),
                 y1p_max=float(mm.group(6)), evaluated_area_pct=float(mm.group(7)), over_area_pct=float(mm.group(8)),
                 index_max=int(mm.group(9)))
        pid = int(bc[name]["physID"])
        dump = run / f"res_{name}_{pid}_{step}.h5"
        if dump.exists():
            with h5py.File(dump) as f:
                xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
            if w["index_max"] < len(xyz):
                w["x_max_rt"], w["r_max_rt"] = float(xyz[w["index_max"], 0] / S), float(xyz[w["index_max"], 1] / S)
        out["per_wall"][name] = w
    if not out["per_wall"]:
        out["raw"] = txt[-3000:]
    return out


def fig_wall_shape(run, F, path):
    """物理壁 r_w (wall_physical.csv) と設計壁 r_inv (wall_design.csv の逆 MOC 点列) の形と 1・2 階微分。
    2 階微分の高周波 (0.5 r_t 移動平均の残差) を壁全体と δ_r 部分 (r_w − r_inv) に分けて出す。"""
    run = Path(run); S = F["S"]
    w = np.loadtxt(run / "wall_physical.csv", delimiter=",", skiprows=1)
    x, r = w[:, 0] / S, w[:, 1] / S
    sp = make_interp_spline(x, r, k=5)
    xs = np.linspace(x[0], x[-1], 20001)
    k = int(round(0.5 / (xs[1] - xs[0])))
    hfun = lambda v: v - np.convolve(v, np.ones(k) / k, "same")
    spi = None
    wd = run / "wall_design.csv"
    if wd.exists():
        t = np.loadtxt(wd, delimiter=",", skiprows=1); xd, rd = t[:, 0] / S, t[:, 1] / S
        o = np.argsort(xd); xd, rd = xd[o], rd[o]
        keep = np.r_[True, np.diff(xd) > 1e-9]
        spi = make_interp_spline(xd[keep], rd[keep], k=5)
    fig, axs = plt.subplots(3, 2, figsize=(13, 8.5), sharex=True)
    for row, der, lab in ((0, 0, "r / r_t"), (1, 1, "dr/dx"), (2, 2, "d²r/dx² [1/r_t]")):
        axs[row, 0].plot(xs, sp(xs, der), lw=1.0, label="物理壁 r_w")
        if spi is not None:
            axs[row, 0].plot(xs, spi(np.clip(xs, xd[0], xd[-1]), der), lw=0.9, ls="--", color="k", label="設計壁 r_inv (逆 MOC)")
            axs[row, 1].plot(xs, sp(xs, der) - spi(np.clip(xs, xd[0], xd[-1]), der), lw=1.0)
        axs[row, 0].set_ylabel(lab); axs[row, 1].set_ylabel(lab.split(" [")[0] + " (r_w − r_inv)")
    lo, hi = np.percentile(sp(xs[xs > 2], 2), [0.5, 99.5]); axs[2, 0].set_ylim(lo - 0.2 * abs(lo), hi + 0.2 * abs(hi))
    axs[0, 0].legend(fontsize=8, frameon=False); axs[0, 0].set_title("壁の形と微分 (左: そのもの / 右: 設計壁との差 = 境界層の補正分)", fontsize=10)
    for a_ in axs[-1]:
        a_.set_xlabel("x / r_t")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    m = (xs > 2) & (xs < xs[-1] - 0.5)
    out = dict(r2_highfreq_max_x_gt2=float(np.abs(hfun(sp(xs, 2))[m][k:-k]).max()), exit_radius_m=float(w[-1, 1]), x_F_rt=float(x[-1]))
    if spi is not None:
        dd = sp(xs, 2) - spi(np.clip(xs, xd[0], xd[-1]), 2)
        out["r2_highfreq_max_x_gt2_delta_part"] = float(np.abs(hfun(dd)[m][k:-k]).max())
        out["r2_highfreq_max_x_gt2_design_wall"] = float(np.abs(hfun(spi(np.clip(xs, xd[0], xd[-1]), 2))[m][k:-k]).max())
    return out


# ------------------------------------------------------------------ 本体
def make_report(run, euler=None, out=None, pptx=True, wall_over_frac=None):
    run = Path(run).resolve(); out = Path(out) if out else run / "report"; out.mkdir(parents=True, exist_ok=True)
    F = load_field(run)
    try:
        sys.path.insert(0, str(TOOLS))
        from total_quantities import total_state
        ts = total_state(str(run), str(run / F["res"]))
        F["V"]["P0"] = np.asarray(ts["P0"], float).reshape(F["ni"], F["nj"])
    except BaseException as e:  # noqa: BLE001  (total_quantities は h0 が無いと SystemExit で抜ける)
        print(f"[nozzle_report] 全圧は未計算 (h0 の無い古い res など): {e}")
    E = load_field(euler) if euler else None
    Md = float(F["info"].get("Md", 6.0))
    rep = dict(conditions=conditions(run, F), metrics=metrics(run, F, euler), figures={})
    fig_domain(run, F, out / "fig_domain.png"); rep["figures"]["domain"] = "fig_domain.png"
    rep["figures"]["contours"] = fig_contours(F, Md, out)
    fig_axis_lines(F, Md, out / "fig_axis_lines.png", out / "fig_axis_deviation.png", E)
    rep["figures"].update(axis="fig_axis_lines.png", axis_dev="fig_axis_deviation.png")
    fig_exit_lines(F, out / "fig_exit_lines.png"); rep["figures"]["exit"] = "fig_exit_lines.png"
    rep["metrics"]["wall"] = fig_wall_lines(F, out / "fig_wall_lines.png"); rep["figures"]["wall"] = "fig_wall_lines.png"
    rep["metrics"]["wall_resolution"] = wall_resolution(run, F, over_frac=wall_over_frac)
    rep["metrics"]["wall_shape"] = fig_wall_shape(run, F, out / "fig_wall_shape.png"); rep["figures"]["wall_shape"] = "fig_wall_shape.png"
    rep["euler_ref"] = str(euler) if euler else None
    (out / "report.json").write_text(json.dumps(rep, indent=1, ensure_ascii=False, default=float))
    if pptx:
        if PPTX_PY.exists():
            r = subprocess.run([str(PPTX_PY), str(Path(__file__).with_name("build_pptx.py")), str(out)], capture_output=True, text=True)
            print(r.stdout.strip() or r.stderr.strip()[-2000:])
        else:
            print(f"[nozzle_report] {PPTX_PY} が無いので pptx は作らない (図と report.json のみ)")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run"); ap.add_argument("--euler", default=None); ap.add_argument("--out", default=None)
    ap.add_argument("--no-pptx", action="store_true")
    ap.add_argument("--wall-over-frac", type=float, default=None,
                    help="壁解像: y1+ > 1 を許す面積割合 [%%] (check_wall_resolution.py --over-frac; 既定はツール既定)")
    a = ap.parse_args(argv)
    print(make_report(a.run, a.euler, a.out, pptx=not a.no_pptx, wall_over_frac=a.wall_over_frac))


if __name__ == "__main__":
    main()
