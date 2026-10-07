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

usage: python3 delta_contur_compare.py extract   (AWS の case dir)
       /home/sano/work/forge/design/.venv-opt/bin/python delta_contur_compare.py plot   (手元)
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
                     "Cf": np.asarray(r["Cf"]), "H": np.asarray(r["H"]), "N": np.asarray(r["N"]), "theta_rt": np.asarray(r["theta"]) / rt, "k_f": k}
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
    a.set_yscale("log")
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


if __name__ == "__main__":
    {"extract": extract, "plot": plot}[sys.argv[1]]()
