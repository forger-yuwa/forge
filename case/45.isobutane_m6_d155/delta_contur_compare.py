"""生産の NS の実際の排除厚さ δ と、積分法 (CONTUR) の式の δ・C_f を比べる (plan verification-m6-axis-wave-mesh-su2 §5.1 #8d の前段:
4 係数化が要るかを、今の生産で試験部に沿った δ_E/δ_C の傾きを測って決める。2026-10-08 ユーザ「測定して、プロットしてほしい。
実態としての δ と、CONTUR の式とかを比べてほしい」)。CFD 0 step。

  extract (AWS、生産の run と同じ環境):
    NS の場 = 生産の dry (N2) の判定窓の 5 枚: run_0167_ns_n012_N2/res_80000 と run_0179_ns_n012_N2_ext/res_5000〜20000 (通算 85000〜100000)。
    Euler 参照 = run_0174_euler_v5d_M_r1 (N2 の δ_E の参照、E5 で参照の時刻への感度は合格)。
    - 実際の δ: `deltastar_from_core_matched_euler(..., band_select="edge")` の半径方向の等価排除厚さ δ_r (生のもの・平滑したもの・採用したもの)
      と、生産の判定の経路 (`extract_and_merge` → delta_r_next.csv の delta_E)。
    - 壁の摩擦: res_wall_3_<step>.h5 の壁面せん断 → C_f = |τ_w| / (½ ρ_e u_e²) (ρ_e・u_e は CONTUR と同じ非粘性の縁の状態)。
    - CONTUR: 生産の問題 (problem_d155_ns_prod.yaml) の設計壁で integral_bl を k_f = 1 (較正なし) と生産の k_f で。δ_r・C_f・H・N・θ。
    → _band_ab/delta_contur/extract.npz と extract.json
  plot (手元): 図 3 枚 → _band_ab/delta_contur/fig_*.png と summary.json (試験部 [40, 94] の δ_E/δ_C の 1 次の傾き・端から端の振れ)。
  tw (手元): CONTUR に NS の壁温 (ns_wall_T.csv) を与えた場合と断熱の式の場合を、k_f = 1 と出口合わせの k_f で比べる
    → fig4_wall_temperature.png・fig5_wall_T.png・tw_experiment.json。
  taw (手元): T_aw の式 (局所 γ_e の現行 / 全温基準) の違いが δ_r を動かす量を、断熱と等温 1000/600/300 K で → taw_sensitivity.json。
  knobs (手元): 較正の係数 (k_f・k_N・a・Δm) と第 1 層 (エンタルピー形の熱閉包 + 混合気の μ) の効き方を、出口で合わせ直して比べる
    → knobs.json・knobs_profiles.json (説明ページ用)。
  hform (手元): 熱閉包の A/B (温度形 / エンタルピー形、isothermal plan §5.1 #10a の事前登録) → hform_ab.json。
  predict (手元): V-c45 の CONTUR の予測 (300 K/断熱 の δ_r の比、k_f 1 と 1.0541、両腕の精度の検査) → cooling_ratio_predictions.json。

usage: python3 delta_contur_compare.py extract   (AWS の case dir)
       /home/sano/work/forge/design/.venv-opt/bin/python delta_contur_compare.py {plot|tw|taw|knobs|hform|predict}   (手元)
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
OUT = HERE / "_band_ab" / "delta_contur"
PROB = "problem_d155_ns_prod.yaml"
EU = "run_0174_euler_v5d_M_r1"
SNAPS = [("run_0167_ns_n012_N2", "res_80000.h5", 80000)] + [("run_0179_ns_n012_N2_ext", f"res_{s}.h5", 80000 + s) for s in (5000, 10000, 15000, 20000)]
TEST = (40.0, 94.0)


def _ns_copy(td: Path, run: Path, res: str) -> Path:
    dd = td / "ns"
    dd.mkdir()
    for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
        shutil.copy(run / f, dd / f)
    os.symlink((run / "nozzle.h5").resolve(), dd / "nozzle.h5")
    os.symlink((run / res).resolve(), dd / res)
    return dd


def _wall_cf(run: Path, step_local: int, edge, scale: float):
    """壁面せん断から C_f。res_wall_3_<step>.h5 のデータセット名は実物で確かめて使う (見つからなければ None と名前の一覧)。"""
    import h5py
    f = run / f"res_wall_3_{step_local}.h5"
    if not f.is_file():
        return None, f"{f.name} が無い"
    with h5py.File(f, "r") as h:
        names = []
        h.visit(names.append)
        V = h["VALUE"] if "VALUE" in h else h
        keys = list(V.keys())
        tx = next((k for k in ("twall_x", "tauw_x", "twall") if k in keys), None)
        if tx is None:
            return None, f"壁面せん断の名前が無い: {keys}"
        txv = np.asarray(V[tx][()], float)
        tyv = np.asarray(V["twall_y"][()], float) if "twall_y" in keys else np.zeros_like(txv)
        coords = None
        for cand in ("MESH/COORD", "COORD", "centCoords", "MESH/centCoords"):
            if cand in h:
                coords = np.asarray(h[cand][()], float).reshape(-1, 3)
                break
        if coords is None or len(coords) != len(txv):
            return None, f"座標が見つからない・長さが違う ({names[:20]})"
    x = coords[:, 0] / scale
    o = np.argsort(x)
    x, tau = x[o], np.hypot(txv[o], tyv[o])
    rho_e = np.asarray(edge["rho_e"](x)) if callable(edge.get("rho_e")) else np.interp(x, edge["x"], edge["rho_e"])
    ue = np.asarray(edge["ue"](x)) if callable(edge.get("ue")) else np.interp(x, edge["x"], edge["ue"])
    return {"x": x, "tau": tau, "cf": tau / (0.5 * rho_e * ue ** 2), "field": tx}, None


def extract():
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    from forge_design.metrics.deltastar import deltastar_from_core_matched_euler, smooth_delta_quintic
    OUT.mkdir(parents=True, exist_ok=True)
    p = load_problem(HERE / PROB)
    d = design_chain(p)
    rt = float(p.spec["r_throat"])
    kf = float(p.raw["deltastar_initializer"]["cf_scale"])
    xF = float(d["wall_inv"][-1, 0])
    xq = np.linspace(-12.0, xF, 2401)
    cont = {}
    from forge_design.evaluate.runner_axismach import integral_delta_r
    init = dict(p.raw["deltastar_initializer"])
    for lab, k in (("k1", 1.0), ("kprod", kf)):
        # 生産と同じ経路 (integral_delta_r: integral_bl → 5 次 P-spline 平滑化 → 壁に渡す δ_r 関数)。k1 は cf_scale だけ 1 にしたもの
        r, drx, _ = integral_delta_r(p, d, {**init, "cf_scale": k})
        xr = np.asarray(r["x"])
        cont[lab] = {"x": xr, "delta_r": np.asarray(r["delta_r_raw_integral"]), "delta_r_smooth": np.asarray(drx(xr)),
                     "Cf": np.asarray(r["Cf"]), "H": np.asarray(r["H"]), "N": np.asarray(r["N"]), "theta_rt": np.asarray(r["theta"]), "k_f": k}  # integral_bl の theta は既に θ/r_t (2026-10-08 codex 指摘で修正: 以前は r_t でもう一度割っていた。それで作った extract.npz の cont_*_theta_rt は 1/r_t = 13.04 倍)
    # 縁の状態 (CONTUR と同じ非粘性の縁): integral_bl の EdgeConditions
    from forge_design.feedback.deltastar_integral import EdgeConditions
    ec = EdgeConditions(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt)
    ed = {"x": xq, "rho_e": np.array([ec.at(float(x))["rho_e"] for x in xq]), "ue": np.array([ec.at(float(x))["ue"] for x in xq])}
    snaps = []
    for run_name, res, step in SNAPS:
        run = HERE / run_name
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            dd = _ns_copy(td, run, res)
            x = deltastar_from_core_matched_euler(dd, HERE / EU, band_select="edge")
            extract_and_merge(dd, HERE / EU, band_select="edge")
            nx = read_delta_r_next(dd / "delta_r_next.csv")
        cf, why = _wall_cf(run, int(res.split("_")[1].split(".")[0]), ed, rt)
        snaps.append({"run": run_name, "res": res, "step_total": step,
                      "x": np.asarray(x["x"]), "draw": np.asarray(x["delta_r_raw"]), "dsm": np.asarray(x["delta_r_smooth"]),
                      "duse": np.asarray(x["delta_r_use"]), "ok": np.asarray(x["ok"]),
                      "merged_x": np.asarray(nx["x_rt"]), "merged_dE": np.asarray(nx["delta_E"]),
                      "cf": cf, "cf_note": why})
        print(f"[extract] {run_name}/{res}: δ 列 {len(x['x'])}、C_f {'ok' if cf else why}", flush=True)
    np.savez(OUT / "extract.npz",
             **{f"cont_{lab}_{k}": v for lab, c in cont.items() for k, v in c.items() if isinstance(v, np.ndarray)},
             **{f"snap{i}_{k}": v for i, s in enumerate(snaps) for k, v in s.items() if isinstance(v, np.ndarray)},
             **{f"snap{i}_cf_{k}": v for i, s in enumerate(snaps) if s["cf"] for k, v in s["cf"].items() if isinstance(v, np.ndarray)},
             edge_x=ed["x"], edge_rho_e=ed["rho_e"], edge_ue=ed["ue"])
    meta = {"problem": PROB, "euler": EU, "r_t_m": rt, "k_f_prod": kf, "x_F": xF, "test": TEST,
            "snaps": [{k: s[k] for k in ("run", "res", "step_total", "cf_note")} | {"cf_field": (s["cf"] or {}).get("field")} for s in snaps],
            "environment": {"numpy": np.__version__}}
    (OUT / "extract.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    print(json.dumps(meta, indent=1, ensure_ascii=False))


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fp = Path.home() / ".fonts" / "NotoSansCJKjp-Regular.otf"
    if fp.is_file():
        font_manager.fontManager.addfont(str(fp))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    Z = np.load(OUT / "extract.npz")
    meta = json.loads((OUT / "extract.json").read_text())
    rt_mm = meta["r_t_m"] * 1000.0
    n = len(meta["snaps"])
    xC = Z["cont_kprod_x"]
    # 窓の 5 枚の δ (採用値) を共通の x に
    xs = Z["snap0_x"]
    D = np.array([np.interp(xs, Z[f"snap{i}_x"], Z[f"snap{i}_duse"]) for i in range(n)])
    Dm, Dlo, Dhi = D.mean(0), D.min(0), D.max(0)
    M = np.array([np.interp(xs, Z[f"snap{i}_merged_x"], Z[f"snap{i}_merged_dE"]) for i in range(n)]).mean(0)
    c1 = np.interp(xs, xC, Z["cont_k1_delta_r"]); cp_ = np.interp(xs, xC, Z["cont_kprod_delta_r"]); cps = np.interp(xs, xC, Z["cont_kprod_delta_r_smooth"])
    C_NS, C_K1, C_KP = "#1f3b73", "#c0392b", "#0b6e4f"
    sel = (xs >= 2.0)
    # 図 1: δ_r [mm]
    fig, a = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    a.fill_between(xs[sel], Dlo[sel] * rt_mm, Dhi[sel] * rt_mm, color=C_NS, alpha=0.25, lw=0, label="NS 実測 δ_E (窓 5 枚の幅)")
    a.plot(xs[sel], Dm[sel] * rt_mm, color=C_NS, lw=1.8, label="NS 実測 δ_E (窓 5 枚の平均、帯 E の採用値)")
    a.plot(xs[sel], c1[sel] * rt_mm, color=C_K1, lw=1.3, ls="--", label="CONTUR 式 (較正なし、k_f = 1)")
    a.plot(xs[sel], cps[sel] * rt_mm, color=C_KP, lw=1.6, label=f"CONTUR 式 (生産の較正、k_f = {meta['k_f_prod']:.4f}、平滑後 = 壁に使う値)")
    a.axvspan(*meta["test"], color="0.93", zorder=0, label="試験部 [40, 94] r_t")
    a.set_xlabel("x / r_t (0 = 設計スロート)"); a.set_ylabel("半径方向の排除厚さ δ_r [mm]"); a.grid(alpha=0.3)
    a.legend(loc="upper left", frameon=False, fontsize=9)
    fig.savefig(OUT / "fig1_delta.png", dpi=150); plt.close(fig)
    # 図 2: 比 δ_E/δ_C と試験部の 1 次の傾き
    ratio = Dm / cps
    rlo, rhi = Dlo / cps, Dhi / cps
    t = (xs >= meta["test"][0]) & (xs <= meta["test"][1]) & np.isfinite(ratio)
    co = np.polyfit(xs[t], ratio[t], 1)
    swing = abs(co[0]) * (meta["test"][1] - meta["test"][0])
    fig, a = plt.subplots(figsize=(10, 4.6), constrained_layout=True)
    a.fill_between(xs[sel], rlo[sel], rhi[sel], color=C_NS, alpha=0.25, lw=0, label="窓 5 枚の幅")
    a.plot(xs[sel], ratio[sel], color=C_NS, lw=1.6, label="δ_E (NS) / δ_C (CONTUR、生産の較正)")
    a.plot(xs[sel], (Dm / c1)[sel], color=C_K1, lw=1.0, ls="--", label="δ_E / δ_C (較正なし)")
    xx = np.array(meta["test"])
    a.plot(xx, np.polyval(co, xx), color="#e67e22", lw=2.0, label=f"試験部の 1 次の当てはめ: 傾き {co[0]*100:+.4f} %/r_t、端から端 {swing*100:.2f} %")
    a.axhline(1.0, color="0.5", lw=0.8)
    a.axvspan(*meta["test"], color="0.93", zorder=0)
    a.set_ylim(0.9, 1.15)
    a.set_xlabel("x / r_t"); a.set_ylabel("δ_E / δ_C"); a.grid(alpha=0.3); a.legend(loc="lower right", frameon=False, fontsize=9)
    fig.savefig(OUT / "fig2_ratio.png", dpi=150); plt.close(fig)
    # 図 3: C_f
    fig, a = plt.subplots(figsize=(10, 4.6), constrained_layout=True)
    have = [i for i in range(n) if f"snap{i}_cf_x" in Z]
    if have:
        xw = Z[f"snap{have[-1]}_cf_x"]
        CF = np.array([np.interp(xw, Z[f"snap{i}_cf_x"], Z[f"snap{i}_cf_cf"]) for i in have])
        sw = xw >= 0.5
        a.plot(xw[sw], CF.mean(0)[sw] * 1e3, color=C_NS, lw=1.4, label="NS 壁面せん断から (窓の平均)")
    a.plot(xC[xC >= 0.5], Z["cont_k1_Cf"][xC >= 0.5] * 1e3, color=C_K1, lw=1.2, ls="--", label="CONTUR の C_f 式 (k_f = 1)")
    a.plot(xC[xC >= 0.5], Z["cont_kprod_Cf"][xC >= 0.5] * 1e3, color=C_KP, lw=1.4, label=f"CONTUR の C_f 式 × k_f ({meta['k_f_prod']:.4f})")
    a.axvspan(*meta["test"], color="0.93", zorder=0)
    a.set_xlabel("x / r_t"); a.set_ylabel("C_f × 10³"); a.grid(alpha=0.3); a.legend(loc="upper right", frameon=False, fontsize=9)
    a.set_ylim(0, None)
    fig.savefig(OUT / "fig3_cf.png", dpi=150); plt.close(fig)
    xF = meta["x_F"]
    summ = {"test": meta["test"], "slope_per_rt": float(co[0]), "swing_end_to_end": float(swing), "criterion_swing": 0.003,
            "within_0p3pct": bool(swing <= 0.003),
            "ratio_mean_test": float(ratio[t].mean()), "ratio_range_test": [float(ratio[t].min()), float(ratio[t].max())],
            "ratio_at_xF_window_mean": float(np.interp(xF, xs, ratio)),
            "merged_dE_over_dC_at_xF": float(np.interp(xF, xs, M) / np.interp(xF, xs, cps)),
            "time_width_ratio_test_max": float((rhi - rlo)[t].max()),
            "ratio_k1_mean_test": float((Dm / c1)[t].mean()),
            "cf_available": bool(have)}
    if have:
        tt = (xw >= meta["test"][0]) & (xw <= meta["test"][1])
        cfk = np.interp(xw, xC, Z["cont_kprod_Cf"]); cf1 = np.interp(xw, xC, Z["cont_k1_Cf"])
        summ["cf_ns_over_contur_kprod_test_mean"] = float((CF.mean(0) / cfk)[tt].mean())
        summ["cf_ns_over_contur_k1_test_mean"] = float((CF.mean(0) / cf1)[tt].mean())
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False))
    print(json.dumps(summ, indent=1, ensure_ascii=False))


def tw_experiment():
    """CONTUR に NS の実際の壁温 (ns_wall_T.csv、断熱壁の NS の窓の平均) を与えると δ のずれと傾きがどう変わるか (手元、CFD 0 step)。
    断熱壁の式の T_aw = T_e(1 + Pr^(1/3)(γ_e − 1)/2 M²) は局所の γ_e を使うため、燃焼ガスでは全温を超える (試験部で約 +190 K)。"""
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback.deltastar_integral import integral_bl
    from forge_design.metrics.deltastar import smooth_delta_quintic
    runs = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
    txt = (HERE / PROB).read_text().replace("initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k",
                                            f"initial_line_run: {runs / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    tmp = OUT / "prod_local.yaml"; tmp.write_text(txt)
    p = load_problem(tmp); d = design_chain(p); rt = float(p.spec["r_throat"])
    tw = np.loadtxt(OUT / "ns_wall_T.csv", delimiter=",", skiprows=1)
    Z = np.load(OUT / "extract.npz"); n = len(json.loads((OUT / "extract.json").read_text())["snaps"])
    xs = Z["snap0_x"]
    Dm = np.array([np.interp(xs, Z[f"snap{i}_x"], Z[f"snap{i}_duse"]) for i in range(n)]).mean(0)
    xF = float(d["wall_inv"][-1, 0])
    res = {}
    def run(lab, k, tbc):
        r = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt, thermal_bc=tbc, cf_scale=k)
        f_s, _ = smooth_delta_quintic(r["x"], r["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)
        dc = np.asarray(f_s(xs))
        ratio = Dm / dc
        t = (xs >= TEST[0]) & (xs <= TEST[1])
        co = np.polyfit(xs[t], ratio[t], 1)
        res[lab] = {"k_f": k, "ratio_test_mean": float(ratio[t].mean()), "slope_per_rt": float(co[0]),
                    "swing": float(abs(co[0]) * (TEST[1] - TEST[0])), "ratio_xF": float(np.interp(xF, xs, ratio)),
                    "Tw_test": [float(np.interp(TEST[0], r["x"], r["Tw"])), float(np.interp(TEST[1], r["x"], r["Tw"]))],
                    "Cf_test_mean": float(np.interp(np.linspace(*TEST, 50), r["x"], r["Cf"]).mean()), "_ratio": ratio}
        return res[lab]
    ad = {"mode": "adiabatic"}
    nsT = {"mode": "prescribed_temperature", "Tw_table": tw.tolist()}
    run("adiabatic_k1", 1.0, ad)
    run("nsTw_k1", 1.0, nsT)
    # 出口で δ を合わせる k_f (C2 と同じ考え): ratio_xF ∝ 1/δ_C(k)。2 点の割線で解く
    for lab, tbc in (("adiabatic_kfit", ad), ("nsTw_kfit", nsT)):
        k0, k1 = 1.0, 1.06
        f0 = run(lab, k0, tbc)["ratio_xF"] - 1.0
        f1 = run(lab, k1, tbc)["ratio_xF"] - 1.0
        for _ in range(4):
            k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
            k0, f0, k1 = k1, f1, k2
            f1 = run(lab, k1, tbc)["ratio_xF"] - 1.0
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fp = Path.home() / ".fonts" / "NotoSansCJKjp-Regular.otf"
    if fp.is_file():
        font_manager.fontManager.addfont(str(fp)); plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    fig, a = plt.subplots(figsize=(10, 4.8), constrained_layout=True)
    sel = xs >= 2.0
    sty = {"adiabatic_k1": ("#c0392b", "--", "断熱の式の壁温、k_f = 1"), "adiabatic_kfit": ("#0b6e4f", "-", "断熱の式の壁温、k_f を出口に合わせる (= 生産)"),
           "nsTw_k1": ("#8e44ad", "--", "NS の壁温を与える、k_f = 1"), "nsTw_kfit": ("#1f6fb2", "-", "NS の壁温を与える、k_f を出口に合わせる")}
    for lab, (c, ls, t_) in sty.items():
        a.plot(xs[sel], res[lab]["_ratio"][sel], color=c, ls=ls, lw=1.5, label=f"{t_} (k_f {res[lab]['k_f']:.4f}、試験部の振れ {res[lab]['swing']*100:.2f} %)")
    a.axhline(1.0, color="0.5", lw=0.8); a.axvspan(*TEST, color="0.93", zorder=0)
    a.set_ylim(0.9, 1.12); a.set_xlabel("x / r_t"); a.set_ylabel("δ_E (NS) / δ_C (CONTUR)"); a.grid(alpha=0.3)
    a.legend(loc="upper right", frameon=False, fontsize=8.5)
    fig.savefig(OUT / "fig4_wall_temperature.png", dpi=150); plt.close(fig)
    fig, a = plt.subplots(figsize=(10, 3.8), constrained_layout=True)
    from forge_design.feedback.deltastar_integral import EdgeConditions
    ec = EdgeConditions(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt)
    xx = np.linspace(0.5, xF, 400)
    taw = [ec.at(float(x))["Te"] * (1 + ec.Pr ** (1 / 3) * 0.5 * (ec.at(float(x))["gam"] - 1) * ec.at(float(x))["M"] ** 2) for x in xx]
    a.plot(xx, taw, color="#c0392b", lw=1.5, label="CONTUR の断熱壁温 T_aw (一定 γ の式に局所 γ_e)")
    a.plot(tw[:, 0][tw[:, 0] >= 0.5], tw[:, 1][tw[:, 0] >= 0.5], color="#1f3b73", lw=1.5, label="NS の壁温 (断熱壁、窓の平均)")
    a.axhline(float(p.spec["Tt"]), color="0.4", lw=0.9, ls=":", label=f"全温 T_t = {float(p.spec['Tt']):.0f} K")
    a.set_xlabel("x / r_t"); a.set_ylabel("壁温 [K]"); a.grid(alpha=0.3); a.legend(loc="lower right", frameon=False, fontsize=9)
    fig.savefig(OUT / "fig5_wall_T.png", dpi=150); plt.close(fig)
    out = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in res.items()}
    (OUT / "tw_experiment.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


def taw_sensitivity():
    """T_aw の式の違いが CONTUR の δ_r をどれだけ動かすか (断熱と等温 1000/600/300 K、生産の k_f、手元、CFD 0 step)。
    現行: T_aw = T_e(1 + Pr^(1/3)(γ_e − 1)/2 M²) (局所 γ_e)。比較: T_aw = T_e + Pr^(1/3)(T_t − T_e) (全温基準)。
    等温壁でも T_aw は Crocco の温度分布の (T_aw − T_w) の項に入るので、壁温を与えても式の差は残る。
    比較の式は EdgeConditions.at の戻り値を差し替えて試すだけ (コードの既定は変えない)。"""
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback import deltastar_integral as DI
    p = load_problem(OUT / "prod_local.yaml"); d = design_chain(p); rt = float(p.spec["r_throat"])
    Tt = float(p.spec["Tt"]); k = float(p.raw["deltastar_initializer"]["cf_scale"])
    orig_at = DI.EdgeConditions.at

    def at_tt(self, x):
        e = orig_at(self, x)
        e["Taw"] = e["Te"] + self.Pr ** (1 / 3) * (Tt - e["Te"])
        return e
    xF = float(d["wall_inv"][-1, 0]); xq = np.array([TEST[0], 70.0, xF])
    out = {"x_query_rt": xq.tolist(), "k_f": k}
    cases = (("adiabatic", {"mode": "adiabatic"}), ("Tw1000", {"mode": "prescribed_temperature", "Tw": 1000.0}),
             ("Tw600", {"mode": "prescribed_temperature", "Tw": 600.0}), ("Tw300", {"mode": "prescribed_temperature", "Tw": 300.0}))
    try:
        for lab, tbc in cases:
            row = {}
            for form in ("current", "tt"):
                DI.EdgeConditions.at = orig_at if form == "current" else at_tt
                r = DI.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), Tt, rt, thermal_bc=tbc, cf_scale=k)
                row[form] = {"delta_r_mm": [float(np.interp(x, r["x"], r["delta_r"])) * rt * 1e3 for x in xq],
                             "Taw_test": [float(np.interp(x, r["x"], r["Taw"])) for x in xq[:2]],
                             "Tw_test": [float(np.interp(x, r["x"], r["Tw"])) for x in xq[:2]]}
            row["rel_change_tt_vs_current"] = [b / a - 1 for a, b in zip(row["current"]["delta_r_mm"], row["tt"]["delta_r_mm"])]
            out[lab] = row
    finally:
        DI.EdgeConditions.at = orig_at
    (OUT / "taw_sensitivity.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


# 混合気 (case/45 の組成) の μ_mix/μ_Sutherland(空気)。solver_density_cuda/tests/unit/transport_reference.py (NS と同じ CEA 輸送物性の
# 独立参照実装) で 2026-10-08 に計算した値。knobs の「第 1 層の一部」の試算にだけ使う
MU_RATIO_T = [250, 300, 350, 400, 600, 800, 1000, 1200, 1400, 1470, 1600]
MU_RATIO = [0.956, 0.965, 0.973, 0.980, 1.007, 1.030, 1.050, 1.068, 1.084, 1.089, 1.099]


def knobs():
    """較正の係数ごとの効き方 (説明用、手元、CFD 0 step)。生産の CONTUR (k_f 1.0541, k_N 1, a 1) を基準に、係数を 1 つ動かしたときの
    δ_r の変化を x に沿って出す。k_f 以外は出口の δ_r を基準と同じに戻すよう k_f を解き直す (生産の C2 と同じく出口で合わせる)。
    Δm (C_f の Re 指数ずらし、#8d の案) と第 1 層 (エンタルピー形の熱閉包 = hform の B 腕 + 混合気の μ) は実装に無いので、ここで差し替えて試す。"""
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback import deltastar_integral as DI
    from forge_design.metrics.deltastar import smooth_delta_quintic
    p = load_problem(OUT / "prod_local.yaml"); d = design_chain(p); rt = float(p.spec["r_throat"]); Tt = float(p.spec["Tt"])
    k0 = float(p.raw["deltastar_initializer"]["cf_scale"])
    orig_closure, orig_at, orig_mu, orig_prof = DI.closure_contur, DI.EdgeConditions.at, DI._sutherland, DI._profile_integrals
    st = {"dm": 0.0, "re_ref": 1.0e4}

    def closure_dm(theta_m, e, Tw, a=1.0, cf_scale=1.0, n_scale=1.0):
        c = orig_closure(theta_m, e, Tw, a=a, cf_scale=cf_scale, n_scale=n_scale)
        if st["dm"]:
            rei = max(c["F_Rdelta"] * c["Re_theta_c"], 300.0)
            c["Cf"] *= (rei / st["re_ref"]) ** (-st["dm"])
        return c
    mu_mix = lambda T: orig_mu(T) * np.interp(T, MU_RATIO_T, MU_RATIO)
    xF = float(d["wall_inv"][-1, 0])

    def run(kf, kN=1.0, a=1.0, dm=0.0, layer1=False):
        st["dm"] = dm
        DI.closure_contur = closure_dm
        if layer1:
            DI.EdgeConditions.at, DI._profile_integrals = hform_patches(DI, _gam_or_gas(p), float(p.cp), Tt)
            DI._sutherland = mu_mix
        try:
            return DI.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), Tt, rt,
                                  thermal_bc={"mode": "adiabatic"}, a_crocco=a, cf_scale=kf, n_scale=kN)
        finally:
            DI.closure_contur, DI.EdgeConditions.at, DI._sutherland, DI._profile_integrals = orig_closure, orig_at, orig_mu, orig_prof
    base = run(k0); xs = base["x"]; dF0 = float(np.interp(xF, xs, base["delta_r"]))

    def pinned(**kw):
        # 出口の δ_r を基準に戻す k_f を、挟み込み (brentq) で解く。δ_r は k_f に単調増加
        from scipy.optimize import brentq
        g = lambda k: float(np.interp(xF, xs, run(k, **kw)["delta_r"])) / dF0 - 1.0
        lo, hi = 0.8 * k0, 1.25 * k0
        k = brentq(g, lo, hi, xtol=1e-7)
        return k, run(k, **kw)
    variants = {"kf_plus5": (k0 * 1.05, run(k0 * 1.05))}
    for lab, kw in (("kN_plus10", dict(kN=1.1)), ("a_0p5", dict(a=0.5)), ("dm_plus0p05", dict(dm=0.05)), ("dm_minus0p05", dict(dm=-0.05)),
                    ("layer1", dict(layer1=True))):
        variants[lab] = pinned(**kw)
    Z = np.load(OUT / "extract.npz"); n = len(json.loads((OUT / "extract.json").read_text())["snaps"])
    xe = Z["snap0_x"]
    Dm = np.array([np.interp(xe, Z[f"snap{i}_x"], Z[f"snap{i}_duse"]) for i in range(n)]).mean(0)
    sm = lambda r: np.asarray(smooth_delta_quintic(r["x"], r["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)[0](xe))
    dC0 = sm(base); ratio0 = Dm / dC0
    t = (xe >= TEST[0]) & (xe <= TEST[1])
    out = {"k_f_base": k0, "x_F": xF, "test": list(TEST), "variants": {}}
    sel = (xs >= 2.0)
    prof = {"x": xs[sel].tolist(), "base": {k: base[k][sel].tolist() for k in ("delta_r", "H", "N", "Cf", "Re_theta_c")}}
    for lab, (kf, r) in variants.items():
        rel = r["delta_r"] / base["delta_r"] - 1.0
        ratio = Dm / sm(r); co = np.polyfit(xe[t], ratio[t], 1)
        out["variants"][lab] = {"k_f": kf, "rel_at": {str(xq): float(np.interp(xq, xs, rel)) for xq in (5.0, 20.0, 40.0, 70.0, xF)},
                                "ratio_swing_test": float(abs(co[0]) * (TEST[1] - TEST[0])), "ratio_slope_per_rt": float(co[0]),
                                "ratio_xF": float(np.interp(xF, xe, ratio))}
        prof[lab] = {"rel": rel[sel].tolist(), "ratio": ratio.tolist()}
    co0 = np.polyfit(xe[t], ratio0[t], 1)
    out["base"] = {"ratio_swing_test": float(abs(co0[0]) * (TEST[1] - TEST[0])), "ratio_slope_per_rt": float(co0[0])}
    prof["xe"] = xe.tolist(); prof["ratio0"] = ratio0.tolist()
    (OUT / "knobs.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    (OUT / "knobs_profiles.json").write_text(json.dumps(prof))
    print(json.dumps(out, indent=1, ensure_ascii=False))


def hform_patches(DI, gas, cp_const: float, Tt: float, n_table: int = 8000, T_range=(100.0, 1800.0)):
    """熱閉包のエンタルピー形 (isothermal plan §5.1 #10a の B 腕) に差し替える (EdgeConditions.at, _profile_integrals) の組を返す。
    h(T) は気体の c_p(T) (CPG なら一定 c_p) の積分。h_aw = h_e + r(h_0 − h_e)、h(v) = h_w + a(h_aw − h_w)v + [h_e − a(h_aw − h_w) − h_w]v²。"""
    orig_at = DI.EdgeConditions.at
    Tg = np.linspace(float(T_range[0]), float(T_range[1]), int(n_table))
    cp = np.asarray(gas.cp_mass(Tg)) if hasattr(gas, "cp_mass") else np.full_like(Tg, cp_const)
    hg = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(Tg))])
    def H(T):
        T = np.asarray(T, dtype=float)
        if np.any(T < Tg[0]) or np.any(T > Tg[-1]):
            raise ValueError(f"hform: 温度 {float(np.min(T)):.1f}〜{float(np.max(T)):.1f} K が表の範囲 {Tg[0]}〜{Tg[-1]} K の外 (端値で外挿しない)")
        return np.interp(T, Tg, hg)

    def Hinv(h):
        h = np.asarray(h, dtype=float)
        if np.any(h < hg[0]) or np.any(h > hg[-1]):
            raise ValueError("hform: エンタルピーが表の範囲の外 (端値で外挿しない)")
        return np.interp(h, hg, Tg)

    def at_h(self, x):
        e = orig_at(self, x); he = H(e["Te"]); h0 = H(Tt)
        e["Taw"] = float(Hinv(he + self.Pr ** (1.0 / 3.0) * (h0 - he)))
        return e

    def prof_h(delta, N, Tw, Taw, Te, rw, cos_phi, a):
        u = DI._GL_U; w = DI._GL_W
        hw, haw, he = H(Tw), H(Taw), H(Te)
        T = Hinv(hw + a * (haw - hw) * u + (he - a * (haw - hw) - hw) * u ** 2)
        rho_rel = Te / np.maximum(T, 1e-30)
        z = delta * u ** N
        dz = N * delta * u ** (N - 1.0)
        curv = 1.0 - z * cos_phi / rw
        theta = float(np.sum(w * curv * rho_rel * u * (1.0 - u) * dz))
        dstar = float(np.sum(w * curv * (1.0 - rho_rel * u) * dz))
        theta_c = float(np.sum(w * rho_rel * u * (1.0 - u) * dz))
        Fc = float(np.sum(w * np.sqrt(rho_rel))) ** -2
        return theta, dstar, theta_c, Fc
    return at_h, prof_h


def hform():
    """熱閉包の A/B (plan tooling-nozzle-isothermal-wall-chain §5.1 #10a、2026-10-08 事前登録、codex diagnose の判別 A/B)。CFD 0 step。
    A = 今の温度形 (Eq. 69 を T で、T_aw = T_e(1 + r(γ_e−1)/2 M²))。
    B = エンタルピー形: h_aw = h_e + r(h_0 − h_e)、h(v) = h_w + a(h_aw − h_w)v + [h_e − a(h_aw − h_w) − h_w]v²、T(v) = h⁻¹(h(v))。
    変える因子は熱閉包の表現だけ (r = 0.72^(1/3)、a = 1、k_f = k_N = 1、今の μ・入口 θ₀・縁条件・運動量式・求積は固定)。
    判定: 試験部 [40, 94] で |(δ_r(300 K)/δ_r(断熱))_B / (同)_A − 1| の最大 ≥ 1 % なら第 1 仮説を支持、< 1 % なら棄却。
    前提: CPG 極限で A と B の差 < 0.1 %、積分精度 (rtol 1e-6 と 1e-8) の差 < 0.1 %。"""
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback import deltastar_integral as DI
    from forge_design.metrics.deltastar import smooth_delta_quintic
    p = load_problem(OUT / "prod_local.yaml"); d = design_chain(p); rt = float(p.spec["r_throat"]); Tt = float(p.spec["Tt"])
    orig_prof, orig_at = DI._profile_integrals, DI.EdgeConditions.at
    gas_sp = _gam_or_gas(p)

    def run(arm, tbc, gas=gas_sp, rtol=1e-6):
        if arm == "B":
            DI.EdgeConditions.at, DI._profile_integrals = hform_patches(DI, gas, float(p.cp), Tt)
        try:
            return DI.integral_bl(d["wall"], d["wall_inv"], gas, p.cp, float(p.spec["Pt"]), Tt, rt, thermal_bc=tbc, rtol=rtol)
        finally:
            DI.EdgeConditions.at, DI._profile_integrals = orig_at, orig_prof
    walls = {"adiabatic": {"mode": "adiabatic"}, "Tw1000": {"mode": "prescribed_temperature", "Tw": 1000.0},
             "Tw600": {"mode": "prescribed_temperature", "Tw": 600.0}, "Tw300": {"mode": "prescribed_temperature", "Tw": 300.0}}
    R = {arm: {k: run(arm, tbc) for k, tbc in walls.items()} for arm in ("A", "B")}
    xs = R["A"]["adiabatic"]["x"]; t = (xs >= TEST[0]) & (xs <= TEST[1])
    out = {"test": list(TEST), "r": 0.72 ** (1 / 3), "k_f": 1.0, "k_N": 1.0, "a": 1.0}
    # 前提 1: CPG 極限 (気体を設計点の γ・c_p 一定にすると、温度形とエンタルピー形は式として一致する)
    pre = {}
    for k in ("adiabatic", "Tw300"):
        a_, b_ = run("A", walls[k], gas=float(p.gamma)), run("B", walls[k], gas=float(p.gamma))
        pre[f"cpg_{k}_max_rel"] = float(np.max(np.abs(b_["delta_r"] / a_["delta_r"] - 1.0)[xs >= 0.5]))
    # 前提 2: 積分精度
    for k in ("adiabatic", "Tw300"):
        fine = run("A", walls[k], rtol=1e-8)
        pre[f"rtol_{k}_max_rel_test"] = float(np.max(np.abs(fine["delta_r"] / R["A"][k]["delta_r"] - 1.0)[t]))
    pre["ok"] = bool(max(pre.values()) < 1e-3)
    out["preconditions"] = pre
    # 判定
    ch = {}
    for k in ("Tw1000", "Tw600", "Tw300"):
        ra = R["A"][k]["delta_r"] / R["A"]["adiabatic"]["delta_r"]
        rb = R["B"][k]["delta_r"] / R["B"]["adiabatic"]["delta_r"]
        c = rb / ra - 1.0
        ch[k] = {"max_abs_change_test": float(np.max(np.abs(c[t]))), "x_at_max": float(xs[t][np.argmax(np.abs(c[t]))]),
                 "ratio_A_at": {str(xq): float(np.interp(xq, xs, ra)) for xq in (40.0, 70.0, 94.0)},
                 "ratio_B_at": {str(xq): float(np.interp(xq, xs, rb)) for xq in (40.0, 70.0, 94.0)}}
    out["cooling_ratio_change"] = ch
    out["verdict"] = ("第 1 仮説を支持 (≥ 1 %)" if ch["Tw300"]["max_abs_change_test"] >= 0.01 else "第 1 仮説を棄却 (< 1 %)") if pre["ok"] else "判定不能 (前提を満たさない)"
    # 記録のみ: 試験部の平均の量と、NS (断熱) との比
    stat = {}
    for arm in ("A", "B"):
        for k in walls:
            r = R[arm][k]
            stat[f"{arm}_{k}"] = {q: float(np.mean(np.asarray(r[q])[t])) for q in ("delta_r", "theta", "H", "Cf", "Tw", "Taw")}
    out["test_means"] = stat
    Z = np.load(OUT / "extract.npz"); n = len(json.loads((OUT / "extract.json").read_text())["snaps"])
    xe = Z["snap0_x"]; te = (xe >= TEST[0]) & (xe <= TEST[1])
    Dm = np.array([np.interp(xe, Z[f"snap{i}_x"], Z[f"snap{i}_duse"]) for i in range(n)]).mean(0)
    xF = float(d["wall_inv"][-1, 0]); ns = {}
    for arm in ("A", "B"):
        r = R[arm]["adiabatic"]
        dc = np.asarray(smooth_delta_quintic(r["x"], r["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)[0](xe))
        ratio = Dm / dc; co = np.polyfit(xe[te], ratio[te], 1)
        ns[arm] = {"ratio_test_mean": float(ratio[te].mean()), "swing_test": float(abs(co[0]) * (TEST[1] - TEST[0])),
                   "ratio_xF": float(np.interp(xF, xe, ratio))}
    out["ns_adiabatic_k1_record_only"] = ns
    prof = {"x": xs.tolist()}
    for arm in ("A", "B"):
        for k in walls:
            prof[f"{arm}_{k}"] = {q: np.asarray(R[arm][k][q]).tolist() for q in ("delta_r", "H", "Cf", "Taw", "Tw")}
    (OUT / "hform_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    (OUT / "hform_ab_profiles.json").write_text(json.dumps(prof))
    print(json.dumps(out, indent=1, ensure_ascii=False))


def predict():
    """V-c45 (plan tooling-nozzle-isothermal-wall-chain §6) の CONTUR の予測: 冷却 300 K / 断熱 の δ_r の比 R_A (温度形)・R_B (エンタルピー形)。
    熱閉包だけを替える (hform と同じ経路、同じ壁)。k_f は診断用の固定値として 1 と生産の 1.0541 を両腕に同じく使う (§5.1 #19)。
    両腕の数値精度 (codex plan m8): RK45 の rtol 1e-6 → 1e-8、B 腕のエンタルピー表 8000 → 32000 点で、試験部の R の変化 < 0.1 % を確かめる。
    出力: _band_ab/delta_contur/cooling_ratio_predictions.json (x、各 k_f・各腕の δ_r と R、精度の検査)。"""
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    from forge_design.feedback import deltastar_integral as DI
    p = load_problem(OUT / "prod_local.yaml"); d = design_chain(p); rt = float(p.spec["r_throat"]); Tt = float(p.spec["Tt"])
    kprod = float(p.raw["deltastar_initializer"]["cf_scale"])
    orig_at, orig_prof = DI.EdgeConditions.at, DI._profile_integrals
    gas = _gam_or_gas(p)
    walls = {"adiabatic": {"mode": "adiabatic"}, "Tw300": {"mode": "prescribed_temperature", "Tw": 300.0}}

    def run(arm, tbc, kf, rtol=1e-6, n_table=8000):
        if arm == "B":
            DI.EdgeConditions.at, DI._profile_integrals = hform_patches(DI, gas, float(p.cp), Tt, n_table=n_table)
        try:
            return DI.integral_bl(d["wall"], d["wall_inv"], gas, p.cp, float(p.spec["Pt"]), Tt, rt, thermal_bc=tbc, cf_scale=kf, rtol=rtol)
        finally:
            DI.EdgeConditions.at, DI._profile_integrals = orig_at, orig_prof
    out = {"k_f": [1.0, kprod], "test": list(TEST), "arms": {}, "precision": {}}
    xs = None
    for kf in (1.0, kprod):
        for arm in ("A", "B"):
            R = {w: run(arm, tbc, kf) for w, tbc in walls.items()}
            xs = R["adiabatic"]["x"]
            ratio = R["Tw300"]["delta_r"] / R["adiabatic"]["delta_r"]
            out["arms"][f"{arm}_kf{kf:.6f}"] = {"k_f": kf, "arm": arm, "delta_r_ad": R["adiabatic"]["delta_r"].tolist(),
                                                 "delta_r_tw300": R["Tw300"]["delta_r"].tolist(), "R": ratio.tolist()}
            # 精度: rtol と (B は) 表の細分
            t = (xs >= TEST[0]) & (xs <= TEST[1])
            Rf = {w: run(arm, tbc, kf, rtol=1e-8) for w, tbc in walls.items()}
            ch = {"rtol_1e-8": float(np.max(np.abs((Rf["Tw300"]["delta_r"] / Rf["adiabatic"]["delta_r"]) / ratio - 1.0)[t]))}
            if arm == "B":
                Rt = {w: run(arm, tbc, kf, n_table=32000) for w, tbc in walls.items()}
                ch["table_32000"] = float(np.max(np.abs((Rt["Tw300"]["delta_r"] / Rt["adiabatic"]["delta_r"]) / ratio - 1.0)[t]))
            ch["ok"] = bool(max(ch.values()) < 1e-3)
            out["precision"][f"{arm}_kf{kf:.6f}"] = ch
    out["x"] = xs.tolist()
    out["precision_ok"] = bool(all(v["ok"] for v in out["precision"].values()))
    (OUT / "cooling_ratio_predictions.json").write_text(json.dumps(out, ensure_ascii=False))
    print(json.dumps({"precision": out["precision"], "precision_ok": out["precision_ok"]}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    {"extract": extract, "plot": plot, "tw": tw_experiment, "taw": taw_sensitivity, "knobs": knobs, "hform": hform,
     "predict": predict}[sys.argv[1]]()
