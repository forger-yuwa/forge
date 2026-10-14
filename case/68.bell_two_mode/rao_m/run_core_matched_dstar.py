#!/usr/bin/env python3
r"""③ベル (M) L100 (CFD ピン版) — コアで Euler を NS に合わせて基準を固定し、壁から不足流量を積分する δ_r の抽出 (ユーザ指定の方法)。

2026-10-12、plan tooling-nozzle-bell-two-mode §5.1 #7 (7b) の比較。NS = run_0019 (res_60000)、Euler = run_0018 (res_30000、同じ物理壁)。

方法 (ユーザ指定):
  q = ρ·Ux (軸方向速度)。NS は一定 x の断面 (`_load_structured(remap_constant_x=True)`)、Euler は同じ x・同じ r に線形補間。
  a(x) = 〈q_NS〉core / 〈q_E〉core、コア = 0 ≤ r/R ≤ f (f = 0.30 標準、0.20・0.40 は感度)、平均は 2πr 重みの面積平均。
  q_ref = a·q_E を断面全体に固定。Δm(Y) = 2π∫_{R−Y}^{R}(q_ref − q_NS) r dr (符号付き、Y は壁からの半径方向距離)。
  δ_r: 2π∫_{R−δ_r}^{R} q_ref r dr = Δm_adopt (負の Δm は符号付きの負の厚さ)。

主セッションが置いた判定 (仮、結果と一緒に明示する):
  平らな所の入口 Y_p = 「[Y, 2Y] での Δm の変化 (最大 − 最小) が |Δm(Y)| の 1 % 以内」を最初に満たす Y (2Y ≤ (1 − f)·R の範囲、
  合わせ込みのコアに入る手前まで)。採用値 = [Y_p, 2Y_p] での Δm の平均。入口が無ければ「平らにならない」と記録する。
  積分は両方を NS 節点の 4 分割の共通 r 格子に載せて台形則。

出力 (`case/68.bell_two_mode/rao_m/`): core_matched_dstar_l100.{json,npz}、core_matched_stations.png、core_matched_dstar_x.png。
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.metrics import deltastar as ds  # noqa: E402

NS_RUN = HERE.parent / "run_0019_tp_m_ns_l100"
EU_RUN = HERE.parent / "run_0018_tp_m_euler_l100_pin"
CORE_FRACS = (0.20, 0.30, 0.40)
STATIONS = (0.70, 3.99, 10.89, 21.51)
PLATEAU_TOL = 0.01
PLATEAU_SPAN = 2.0
REFINE = 4


def section(ns, eu, i):
    """断面 i の (r, qN, qE, R) を細かい共通の r 格子で返す (r は軸 0 → 壁 R)。"""
    xi = float(ns["x"][i, 0])
    rn, qn = ns["r"][i, :], ns["q"][i, :]
    xe = eu["x"][:, 0]
    j = int(np.clip(np.searchsorted(xe, xi), 1, len(xe) - 1))
    w = (xi - xe[j - 1]) / (xe[j] - xe[j - 1])
    re = (1 - w) * eu["r"][j - 1, :] + w * eu["r"][j, :]
    qe = (1 - w) * eu["q"][j - 1, :] + w * eu["q"][j, :]
    t = np.linspace(0.0, 1.0, REFINE, endpoint=False)
    rr = np.concatenate([(rn[:-1, None] + np.diff(rn)[:, None] * t[None, :]).ravel(), [rn[-1]]])
    qN = np.interp(rr, rn, qn)
    qE = np.interp(rr, re, qe)
    return xi, rr, qN, qE, float(rn[-1]), float(re[-1])


def cumtrapz0(f, r):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])


def core_matched(rr, qN, qE, R, f):
    m = rr <= f * R
    a = float(np.trapezoid(qN[m] * rr[m], rr[m]) / np.trapezoid(qE[m] * rr[m], rr[m]))
    qref = a * qE
    # 壁からの累積: Y = R − r。rr は昇順 (軸 → 壁) なので逆順に積む
    rd = rr[::-1]
    Y = R - rd
    dm = 2 * np.pi * cumtrapz0(((qref - qN) * rr)[::-1], rd) * -1.0   # rd は降順なので符号を戻す
    F = 2 * np.pi * cumtrapz0((qref * rr)[::-1], rd) * -1.0
    # 平らな所の入口
    Ymax = (1.0 - f) * R
    Yp = np.nan; adopt = np.nan; var = np.nan
    for k in range(1, len(Y)):
        y0 = Y[k]
        if PLATEAU_SPAN * y0 > Ymax:
            break
        sel = (Y >= y0) & (Y <= PLATEAU_SPAN * y0)
        if sel.sum() < 4 or dm[k] == 0:
            continue
        span = float(dm[sel].max() - dm[sel].min())
        if span <= PLATEAU_TOL * abs(dm[k]):
            Yp = float(y0); adopt = float(dm[sel].mean()); var = span / abs(dm[k])
            break
    if np.isfinite(adopt):
        if adopt >= 0:
            dr = float(np.interp(adopt, F, Y))
        else:   # 符号付きの負の厚さ: 壁の外へ q_ref(R) の延長で近似
            dr = float(adopt / (2 * np.pi * R * qref[-1]))
    else:
        dr = np.nan
    return dict(a=a, Y=Y, dm=dm, F=F, Yp=Yp, adopt=adopt, var=var, dr=dr, qref=qref, total=float(dm[-1]))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fpath in font_manager.findSystemFonts([os.path.expanduser("~/.fonts")]):
        font_manager.fontManager.addfont(fpath)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False

    ns = ds._load_structured(NS_RUN, remap_constant_x=True)
    eu = ds._load_structured(EU_RUN)
    x = ns["x"][:, 0]
    ii = np.where(x >= 0.0)[0]
    res = {f: dict(a=[], dr=[], Yp=[], adopt=[], var=[], total=[]) for f in CORE_FRACS}
    for i in ii:
        xi, rr, qN, qE, R, Re = section(ns, eu, i)
        for f in CORE_FRACS:
            o = core_matched(rr, qN, qE, R, f)
            for k in res[f]:
                res[f][k].append(o[k])
    xs = x[ii]
    out = {"method": __doc__.split("\n\n")[1], "plateau_tol": PLATEAU_TOL, "plateau_span": PLATEAU_SPAN, "refine": REFINE,
           "ns_run": str(NS_RUN.relative_to(ROOT)), "eu_run": str(EU_RUN.relative_to(ROOT)), "ns_res": ns["res"], "eu_res": eu["res"],
           "core_fracs": CORE_FRACS, "n": int(len(ii))}
    arr = {"x": xs}
    for f in CORE_FRACS:
        for k, v in res[f].items():
            arr[f"{k}_{int(f*100)}"] = np.asarray(v, dtype=float)
    # 既存法
    ed = np.load(HERE / "dstar_l100_ns_run0019_constx_edge.npz", allow_pickle=True)
    ad = np.load(HERE / "dstar_l100_ns_run0019_constx.npz", allow_pickle=True)
    co = np.load(HERE / "contur_l100_euler_run0018.npz", allow_pickle=True)
    d_ed = np.interp(xs, ed["x"], ed["delta_r_use"]); d_ad = np.interp(xs, ad["x"], ad["delta_r_use"]); d_co = np.interp(xs, co["x"], co["delta_r_ann"])
    d30 = arr["dr_30"]
    sens = np.nanmax(np.abs(np.c_[arr["dr_20"], arr["dr_40"]] / d30[:, None] - 1.0), axis=1)
    noplat = {int(f*100): int(np.sum(~np.isfinite(arr[f"dr_{int(f*100)}"]))) for f in CORE_FRACS}
    regions = [(0, 1), (1, 3), (3, 8), (8, 16), (16, 23)]
    summary = []
    for lo, hi in regions:
        m = (xs >= lo) & (xs < hi) & np.isfinite(d30)
        summary.append({"x": [lo, hi], "n": int(((xs >= lo) & (xs < hi)).sum()), "n_plateau_30": int(m.sum()),
                        "core_sens_max": float(np.nanmax(sens[m])) if m.any() else None,
                        "vs_edge_median": float(np.nanmedian(d30[m] / d_ed[m] - 1)) if m.any() else None,
                        "vs_edge_maxabs": float(np.nanmax(np.abs(d30[m] / d_ed[m] - 1))) if m.any() else None,
                        "vs_contur_median": float(np.nanmedian(d30[m] / d_co[m])) if m.any() else None})
    out.update(no_plateau=noplat, regions=summary)
    st = {}
    for s_ in STATIONS:
        i = int(ii[np.argmin(abs(xs - s_))])
        xi, rr, qN, qE, R, Re = section(ns, eu, i)
        row = {"x": xi, "R": R, "R_euler": Re}
        for f in CORE_FRACS:
            o = core_matched(rr, qN, qE, R, f)
            row[f"core{int(f*100)}"] = {k: o[k] for k in ("a", "Yp", "adopt", "var", "dr", "total")}
        k = int(np.argmin(abs(xs - xi)))
        row["edge_dr"] = float(d_ed[k]); row["adaptive_dr"] = float(d_ad[k]); row["contur_dr"] = float(d_co[k])
        st[f"{xi:.2f}"] = row
    out["stations"] = st
    (HERE / "core_matched_dstar_l100.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    np.savez(HERE / "core_matched_dstar_l100.npz", **arr, d_edge=d_ed, d_adaptive=d_ad, d_contur=d_co, core_sens=sens)

    # ---- 図 1: 4 断面
    col = {0.20: "#7a5195", 0.30: "#d9822b", 0.40: "#2f6db5"}
    fig, ax = plt.subplots(3, 4, figsize=(16, 11))
    for c, s_ in enumerate(STATIONS):
        i = int(ii[np.argmin(abs(xs - s_))])
        xi, rr, qN, qE, R, Re = section(ns, eu, i)
        o3 = core_matched(rr, qN, qE, R, 0.30)
        a0 = ax[0, c]; rat = qN / qE
        a0.plot(rr / R, rat, color="#222", lw=1.4, label="q_NS / q_Euler")
        for f in CORE_FRACS:
            o = core_matched(rr, qN, qE, R, f)
            a0.axhline(o["a"], color=col[f], lw=1.0, ls="--", label=f"a (コア {int(f*100)} %) = {o['a']:.4f}")
            a0.axvline(f, color=col[f], lw=0.8, ls=":")
        a0.axvspan(0, 0.30, color=col[0.30], alpha=0.07, lw=0)
        rr_core = rat[rr <= 0.95 * R]
        lo = max(0.0, float(np.nanmin(rr_core)) - 0.01); hi = float(np.nanmax(rr_core)) + 0.01
        a0.set_ylim(lo, hi); a0.set_title(f"x = {xi:.2f} r_t  (R = {R:.3f})", fontsize=11); a0.grid(alpha=0.3)
        a0.set_xlabel("r / R (軸 0、壁 1)")
        if c == 0: a0.set_ylabel("q_NS / q_Euler")
        a0.legend(fontsize=7.5, loc="lower left")
        a1 = ax[1, c]; y = R - rr
        Ymax = max(4 * (o3["Yp"] if np.isfinite(o3["Yp"]) else 0.05 * R), 0.06 * R)
        sel = y <= Ymax
        a1.plot(y[sel], qN[sel], color="#222", lw=1.5, label="q_NS")
        a1.plot(y[sel], o3["qref"][sel], color=col[0.30], lw=1.3, ls="--", label="q_ref = a·q_Euler (コア 30 %)")
        a1.set_xlabel("壁からの距離 Y [r_t]"); a1.grid(alpha=0.3); a1.legend(fontsize=8)
        if c == 0: a1.set_ylabel("ρu [kg/(m²·s)]")
        a2 = ax[2, c]
        for f in CORE_FRACS:
            o = core_matched(rr, qN, qE, R, f)
            a2.plot(o["Y"], o["dm"], color=col[f], lw=1.3, label=f"コア {int(f*100)} %")
            if np.isfinite(o["Yp"]):
                a2.axvspan(o["Yp"], PLATEAU_SPAN * o["Yp"], color=col[f], alpha=0.10, lw=0)
                a2.axhline(o["adopt"], color=col[f], lw=0.7, ls=":")
        a2.axvline((1 - 0.30) * R, color="0.5", lw=0.8, ls="--")
        a2.set_xscale("log"); a2.set_xlim(max(1e-4, y[y > 0].min()), R)
        tot = max(abs(o3["total"]), abs(np.nanmax(o3["dm"])))
        a2.set_ylim(min(0, np.nanmin(o3["dm"])) - 0.1 * tot, 1.3 * tot)
        a2.set_xlabel("壁からの距離 Y [r_t] (対数)"); a2.grid(alpha=0.3, which="both")
        if c == 0: a2.set_ylabel("Δm(Y) [kg/s]")
        a2.legend(fontsize=8, loc="upper left")
        txt = ["δ_r [r_t]:"]
        for f in CORE_FRACS:
            o = core_matched(rr, qN, qE, R, f)
            txt.append(f"コア {int(f*100)} %: {o['dr']:.5f}" if np.isfinite(o["dr"]) else f"コア {int(f*100)} %: 平らにならない")
        k = int(np.argmin(abs(xs - xi)))
        txt += [f"edge: {d_ed[k]:.5f}", f"adaptive: {d_ad[k]:.5f}"]
        a2.text(0.98, 0.04, "\n".join(txt), transform=a2.transAxes, ha="right", va="bottom", fontsize=8.5,
                bbox=dict(fc="white", ec="0.7", alpha=0.9))
    fig.suptitle("コアで合わせた基準 (q_ref = a·q_Euler) による不足流量。上: 比と倍率 a、中: 基準と NS、下: 壁から積んだ Δm(Y) (影 = 平らな所の判定区間 [Y_p, 2Y_p]、破線 = コア 30 % の端)", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(HERE / "core_matched_stations.png", dpi=120); plt.close(fig)

    # ---- 図 2: 全断面
    fig, ax = plt.subplots(3, 1, figsize=(10, 10.5), sharex=True, gridspec_kw={"height_ratios": [1.3, 1, 1]})
    a = ax[0]
    for f in CORE_FRACS:
        a.plot(xs, arr[f"dr_{int(f*100)}"], color=col[f], lw=1.4 if f == 0.30 else 0.9, label=f"コアで合わせる方法 (コア {int(f*100)} %)")
    a.plot(xs, d_ed, color="#555", lw=1.0, ls="--", label="既存: edge")
    a.plot(xs, d_ad, color="#999", lw=1.0, ls=":", label="既存: adaptive")
    a.plot(xs, d_co, color="#3a8f5c", lw=1.0, ls="-.", label="理論式 contur_v2 (較正なし)")
    a.set_ylabel("δ_r [r_t]"); a.grid(alpha=0.3); a.legend(fontsize=8.5, loc="upper left")
    a.set_title("L100 (CFD ピン版) の δ_r(x): コアで合わせる方法と既存の帯の方法", fontsize=11)
    b = ax[1]
    for f in CORE_FRACS:
        b.plot(xs, arr[f"dr_{int(f*100)}"] / d_ed, color=col[f], lw=1.2, label=f"コア {int(f*100)} % / edge")
    b.plot(xs, d_co / d_ed, color="#3a8f5c", lw=1.0, ls="-.", label="理論式 / edge")
    b.axhline(1, color="0.4", lw=0.8); b.set_ylim(0.7, 1.3); b.set_ylabel("edge との比"); b.grid(alpha=0.3); b.legend(fontsize=8.5)
    c = ax[2]
    for f in CORE_FRACS:
        c.plot(xs, arr[f"a_{int(f*100)}"], color=col[f], lw=1.2, label=f"倍率 a (コア {int(f*100)} %)")
    c.set_ylabel("倍率 a(x)"); c.grid(alpha=0.3)
    c2 = c.twinx(); c2.semilogy(xs, np.maximum(sens, 1e-5), color="#c0392b", lw=0.9, alpha=0.8, label="コア範囲への感度 (20・40 % と 30 % の δ_r の差の最大)")
    c2.axhline(0.05, color="#c0392b", lw=0.8, ls="--"); c2.set_ylim(1e-4, 1); c2.set_ylabel("感度 (相対)")
    nan30 = ~np.isfinite(d30)
    c.fill_between(xs, 0, 1, where=nan30, transform=c.get_xaxis_transform(), color="#c0392b", alpha=0.12, lw=0)
    h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels(); c.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    c.set_xlabel("x [r_t] (スロート = 0、出口 = 22.66)。薄い赤 = コア 30 % で Δm が平らにならない断面")
    fig.tight_layout(); fig.savefig(HERE / "core_matched_dstar_x.png", dpi=130); plt.close(fig)
    print(json.dumps({"no_plateau": noplat, "regions": summary, "stations": {k: {kk: (vv if not isinstance(vv, dict) else {q: (round(w, 6) if isinstance(w, float) else w) for q, w in vv.items()}) for kk, vv in v.items()} for k, v in st.items()}}, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
