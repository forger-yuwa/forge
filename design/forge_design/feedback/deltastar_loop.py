r"""排除厚さの固定点反復ドライバ (固定 Euler 基準・コア整合抽出)。

計画: plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5–4.6。

1 pass = 前 pass の NS run から $\delta_r(x)$ を抽出 (`metrics.deltastar.deltastar_from_core_matched_euler`)
→ 緩和して次の入力 $\delta_{in}^{k+1} = (1-\omega)\delta_{in}^k + \omega\,\delta_{use}^k$ (hard 不合格の断面は前回値保持)
→ 半径方向オフセットの物理壁で `prepare_ns` → 段階起動 NS → `collect` (質量流量帳簿込み)
→ 新 run からも抽出して固定点差を記録 (`fixed_point.json`)。

使い方 (リポジトリルートで):
  design/.venv-opt/bin/python -m forge_design.feedback.deltastar_loop \
      --problem case/45.isobutane_m6_d155/problem_d155_ns.yaml \
      --euler-ref case/45.isobutane_m6_d155/run_0001_euler_shortest_dry \
      --prev case/45.isobutane_m6_d155/run_0004_ns_v3 \
      --run-dir case/45.isobutane_m6_d155/run_0019_ns_cm_pass1 [--omega 0.5] [--prepare-only]

推奨 (2026-09-04, case/45 run_0026 で検証): 収束済み NS 場からの warm start では
  `--stages none --cfl 5 --implicit-relax 0.7 --steps 12000`
(soft/mid 6000 step は不要、cfl 5 + relax 0.7 で 5000〜9000 step で cfl1/24000 step と同じ残差水準・出口 M は 8000 step で凍結。
 NS 1 本 ≈ 95 s)。cold start (Euler → 中継 → 初回 NS) では `--stages full` (既定) か `--stages ramp --ramp 1,2,3.5`。
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..metrics.deltastar import deltastar_from_core_matched_euler, massflow_ratio


# delta_r_next.csv の列 (2026-10-05, plan tooling-nozzle-cfd-pinned-initial-line §5.1 #8 — codex result M1):
# 第 2 列は**緩和・再平滑化後の次 pass の壁入力** δ_target、第 4 列が**未緩和の平滑化抽出値** δ_E。
# 2026-10-05 以前のファイルは同じ並びで名前だけ旧称 (delta_r / delta_r_use_smoothed)。列の位置は変えていないので、
# 第 2 列を壁入力として読む既存の読み手 (prepare_ns の delta_r_csv 等) の挙動は変わらない。
DELTA_R_NEXT_COLUMNS = ("x_rt", "delta_target", "delta_in_prev", "delta_E", "held")
_DELTA_R_NEXT_OLD = ("x_rt", "delta_r", "delta_in_prev", "delta_r_use_smoothed", "held")


def read_delta_r_next(path) -> dict:
    """`delta_r_next.csv` を列名で読む。戻り: x_rt / delta_target (緩和後の壁入力) / delta_in_prev /
    delta_E (未緩和の平滑化抽出) / held と header_kind ("named" か "legacy")。並びが想定と違えば拒否する。"""
    path = Path(path)
    with open(path) as f:
        head = f.readline().strip()
    names = tuple(c.strip().split(" ")[0] for c in head.split(","))
    if names[:5] == DELTA_R_NEXT_COLUMNS:
        kind = "named"
    elif names[:5] == _DELTA_R_NEXT_OLD:
        kind = "legacy"
    else:
        raise ValueError(f"{path}: delta_r_next.csv の列が想定外 ({names})")
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    out = {k: a[:, i] for i, k in enumerate(DELTA_R_NEXT_COLUMNS)}
    out["header_kind"] = kind
    return out


def extract_and_merge(prev_run, euler_run, omega: float = 0.5, smooth_lam: float = 1.0,
                      knot_spacing: float = 2.0, out_dir=None, max_lam_factor: float | None = None, **kw) -> dict:
    """前 pass の NS run から抽出し、次 pass の入力 δ_r(x) を作る。
    出力 (out_dir、既定 prev_run): delta_r_equiv.csv / delta_r_equiv_diag.json / delta_r_next.csv / delta_r_extract_summary.json。
    out_dir は同じ場から別設定 (例 band_select) で抽出するとき互いに上書きしないために使う。"""
    prev_run = Path(prev_run)
    od = Path(out_dir) if out_dir is not None else prev_run
    od.mkdir(parents=True, exist_ok=True)
    d = deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=od, **kw)
    d_in = d["delta_in"]
    d_use = d["delta_r_use"]
    held = ~np.isfinite(d_use)
    # 平滑化 (2026-09-04 ユーザ指摘: 3 次平滑化では壁曲率が凸凹): 抽出 raw を 5 次 P-spline (3 階差分ペナルティ,
    # ノット 2 r_t, λ=1) で「曲率が滑らか」な δ_ext(x) にしてから緩和する。hard 不合格の点は重み 0。
    from ..metrics.deltastar import smooth_delta_quintic
    # 符号付き δ_r (2026-09-12): 等温壁 (prepare_info.json の wall_thermal) では負値を許す。断熱は従来どおり ≥0 に丸める
    try:
        _wt = json.loads((prev_run / "prepare_info.json").read_text()).get("wall_thermal", {"mode": "adiabatic"})
    except Exception:
        _wt = {"mode": "adiabatic"}
    positive = (str(_wt.get("mode", "adiabatic")) == "adiabatic")
    x = d["x"]; r_inv = d["r_wall_euler"]
    raw = d["delta_r_raw"]; wts = d["hard_ok"].astype(float) & np.isfinite(raw) if False else (d["hard_ok"] & np.isfinite(raw)).astype(float)
    lam_used = None; diag = None
    for lam in (smooth_lam, smooth_lam * 10, smooth_lam * 100, smooth_lam * 1000):
        f_s, diag = smooth_delta_quintic(x, raw, weights=wts, knot_spacing=knot_spacing, lam=lam, positive=positive)
        d_ext = f_s(x)
        d_next = (1.0 - omega) * d_in + omega * d_ext
        if omega < 1.0:
            # 前回入力 δ_in (旧方式の壁など) の凸凹を引き継がないよう、緩和後も同じ P-spline を通す
            f_b, _ = smooth_delta_quintic(x, d_next, knot_spacing=knot_spacing, lam=lam, positive=positive)
            d_next = f_b(x)
        # 単調性ガード: 新しい物理壁 r_inv + δ_next がスロート下流で非単調なら λ を 10 倍して再平滑化
        r_new = r_inv + d_next
        m = x >= 0.5
        if np.all(np.diff(r_new[m]) > -1e-9):
            lam_used = lam; break
    else:
        lam_used = f"{lam} (still non-monotone)"
    # 単調性ガードの λ 昇格に上限 (2026-10-04, diagnostician): 昇格しすぎた平滑化は抽出値を追わない壁になるので pass を失敗にする
    if max_lam_factor is not None and (not isinstance(lam_used, float) and not isinstance(lam_used, int)
                                       or lam_used > smooth_lam * max_lam_factor):
        raise RuntimeError(f"extract_and_merge: 単調性ガードで λ が {lam_used} まで上がった (上限 smooth_lam×{max_lam_factor})")
    d_use = d_ext                                   # 以降の帳簿は平滑化後の抽出値で取る
    held = ~np.isfinite(d["delta_r_use"])
    np.savetxt(od / "delta_r_next.csv", np.c_[d["x"], d_next, d_in, d_use, held.astype(int)],
               delimiter=",", comments="",
               header=",".join(DELTA_R_NEXT_COLUMNS) + f" (delta_target = relaxed wall input; delta_E = unrelaxed smoothed extraction; omega={omega}; quintic P-spline knot={knot_spacing} lam={lam_used}; resid_rel_rms={diag['resid_rel_rms']:.4f})")
    fin = np.isfinite(d_use) & (d_in > 1e-4)
    ratio = d_use[fin] / d_in[fin]
    summary = {"omega": omega, "smooth": {"kind": "quintic_pspline", "knot_spacing": knot_spacing, "lam": str(lam_used), **diag},
               "n_stations": int(len(d["x"])), "n_ok": int(d["ok"].sum()),
               "n_hard_ok": int(d["hard_ok"].sum()), "n_held": int(held.sum()),
               "use_over_in_median": float(np.median(ratio)) if fin.any() else None,
               "use_over_in_p10_p90": [float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))] if fin.any() else None,
               "max_abs_change": float(np.nanmax(np.abs(d_use - d_in))) if fin.any() else None,
               "massflow": d["massflow"],
               "delta_r_throat_use": float(np.interp(0.0, d["x"], np.where(np.isfinite(d_use), d_use, d_in))),
               "delta_r_throat_in": float(np.interp(0.0, d["x"], d_in))}
    (od / "delta_r_extract_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def _sizing_delta_supplier(p, init_cfg: dict, rtol: float | None = None):
    """寸法の逆算の CFD 前の δ_r の供給: `prepare_ns` と同じ経路 (`integral_delta_r`: k_f の cf_scale・熱条件・5 次 P-spline の平滑化・
    `delta_r_from_table`) で、反復のたびに r_t を変えて作り直す (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)。
    rtol: 診断用 (同 plan §6 U2c)。`integral_delta_r(rtol=...)` へそのまま渡す (None = `integral_bl` の既定)。
    戻り: supplier(d, rt) → (delta_r_x, 記録)。記録の `integral_rtol` は solve_ivp に実際に渡った値。"""
    from ..evaluate.runner_axismach import integral_delta_r

    def supplier(d, rt):
        res, drx, info = integral_delta_r(p, d, init_cfg, scale=rt, rtol=rtol)
        return drx, {"kind": "integral_delta_r (prepare_ns と同じ経路)", "settings": {k: info.get(k) for k in
                     ("model", "thermal_bc", "cf_scale", "n_scale", "a_crocco", "closure", "theta0_m", "x_virtual_m")},
                     "smooth": info.get("smooth", {}).get("kind"),
                     "integral_rtol": float(res["solve_ivp"]["rtol"]), "integral_rtol_injected": rtol is not None}
    return supplier


SIZING_TOL_M = 1e-9      # 寸法の逆算の数値解法の許容差 [m] (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2・U2b の合格条件と同じ値)。
                         # 既定で使うのは NS 後の `solve_rt_throat` (δ_E の全分布) と試験用の `_delta_supplier`。CFD 前の既定は下の 2 つ
SIZING_MAX_ITER = 30
# CFD 前の寸法 (積分法の δ_r) は**初期見積もり**: その既定の許容差 [m] (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d、
# 諮問 notes/reviews/2026-10-07-upoly-sizing-after-u2c-diagnose.md の P1)。CFD 前の δ_r は r_t に対して滑らかでなく、局所の散らばり
# (rtol 1e-6: スロート半径 7.41e-9 m・出口半径 2.64e-6 m、U2c) が SIZING_TOL_M より大きくて届かないため。真値の誤差の保証ではなく、
# 精密な寸法・感度の評価には未検証 (同 plan §6 既知の制約)。SIZING_TOL_M (共通) は変えない
SIZING_TOL_PRE_CFD_THROAT_M = 1e-7
SIZING_TOL_PRE_CFD_EXIT_M = 1e-5


def _resolve_sizing_tol(tol_R_m, target: str, cfd_before: bool) -> tuple[float, str]:
    """寸法の逆算の実効の許容差 [m] とその出典。明示 (tol_R_m) が優先。既定は CFD 前なら初期見積もりの定数 (target 別)、
    それ以外 (NS 後・試験用の供給) は SIZING_TOL_M。"""
    if tol_R_m is not None:
        if isinstance(tol_R_m, bool) or not np.isfinite(float(tol_R_m)) or float(tol_R_m) <= 0.0:
            raise ValueError(f"tol_R_m は正の有限の数値 [m] (受け取った値: {tol_R_m!r})")
        return float(tol_R_m), "explicit"
    if cfd_before:
        if target == "exit":
            return SIZING_TOL_PRE_CFD_EXIT_M, "default: SIZING_TOL_PRE_CFD_EXIT_M (CFD 前の初期見積もり)"
        return SIZING_TOL_PRE_CFD_THROAT_M, "default: SIZING_TOL_PRE_CFD_THROAT_M (CFD 前の初期見積もり)"
    return SIZING_TOL_M, "default: SIZING_TOL_M"


class SizingNotConverged(RuntimeError):
    """寸法の逆算が反復の上限で許容差に収まらなかった (不合格)。`history` に各反復の r_t・残差、`tol_R_m` に実効の許容差を持つ。"""

    def __init__(self, msg: str, history: list, tol_R_m: float | None = None) -> None:
        super().__init__(msg)
        self.history = history
        self.tol_R_m = tol_R_m


def _fixed_point_rt(p, d, R_target_m: float, target: str, supplier, rt0: float, max_iter: int, tol_R_m: float) -> dict:
    """r_t を固定点反復で解く: 反復ごとに δ_r を作り直して物理壁を組み (`build_physical_wall`、prepare_ns と同じ)、
    target = 'exit' なら r_t·r_W(x_e)、'throat' なら r_t·min r_W (物理壁の大域最小) を目標に合わせる。
    各反復で今の r_t の残差 rt·r(rt) − R を評価し、|残差| ≤ tol_R_m で止める (返す r_t の残差は最後に評価した値そのもの)。
    max_iter 回で収まらなければ例外 (plan §4.2「反復の上限に達したら不合格」)。"""
    from ..evaluate.runner_axismach import build_physical_wall
    rt = float(rt0)
    hist = []
    for k in range(int(max_iter)):
        drx, src = supplier(d, rt)
        W = build_physical_wall(p, d, rt, delta_r_x=drx, offset="radial")
        x_e = float(W.x_e)
        r_e = float(W.r(np.array([x_e]))[0])
        r_val = r_e if target == "exit" else float(W.r_throat)
        res = rt * r_val - R_target_m
        hist.append({"iter": k, "r_t_m": rt, "r_rt": r_val, "residual_m": res, "x_throat_rt": float(W.x_throat),
                     "r_throat_rt": float(W.r_throat), "exit_r_rt": r_e, "delta_r_exit_rt": float(drx(np.array([x_e]))[0]),
                     "delta_r_x0_rt": float(drx(np.array([0.0]))[0])})
        if abs(res) <= tol_R_m:
            return {"r_t_m": rt, "residual_m": res, "wall": W, "delta_r_x": drx, "delta_r_source": src, "iters": hist,
                    "n_iter": k + 1, "tol_R_m": float(tol_R_m)}
        rt = R_target_m / r_val
    res = [h["residual_m"] for h in hist]
    raise SizingNotConverged(f"寸法の逆算 ({target}): {max_iter} 回で |r_t·r − R| ≤ {tol_R_m:g} m に収まらない "
                             f"(最後の残差 {res[-1]:.3e} m、後半の |残差| の範囲 {min(map(abs, res[len(res) // 2:])):.2e}〜"
                             f"{max(map(abs, res[len(res) // 2:])):.2e} m) — 不合格 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)",
                             hist, tol_R_m=float(tol_R_m))


def _sizing_result(fp: dict, p, R_target_m: float, target: str, tol_source: str) -> dict:
    """`solve_rt`・`solve_rt_throat` の戻りの共通部分: 解いた r_t、最終の寸法での残差、実効の許容差とその出典、物理スロートの位置と半径、
    出口半径、δ_r の出典。"""
    W, rt = fp["wall"], fp["r_t_m"]
    x_e = float(W.x_e)
    return {"r_t_m": rt, "sizing": target, "target_m": R_target_m, "residual_m": fp["residual_m"], "n_iter": fp["n_iter"],
            "tol_R_m": fp["tol_R_m"], "tol_R_m_source": tol_source,
            "pw_upstream": W.pw_upstream, "physical_wall_repr": type(W).__name__,
            "physical_throat": {"x_rt": float(W.x_throat), "r_rt": float(W.r_throat), "kappa_rt": float(W.kappa_throat),
                                "x_m": float(W.x_throat) * rt, "r_m": float(W.r_throat) * rt},
            "exit_radius_m": float(W.r(np.array([x_e]))[0]) * rt, "x_e_rt": x_e, "length_m": x_e * rt,
            "delta_r_source": fp["delta_r_source"], "iter_history": fp["iters"]}


def solve_rt(problem, R_exit_m: float, prev_run=None, euler_run=None, n_iter: int = 6,
             max_iter: int = SIZING_MAX_ITER, tol_R_m: float | None = None, integral_rtol: float | None = None) -> dict:
    r"""出口の物理半径 $R$ を仕様に合わせる **スロート半径 $r_t$ の 1 変数解**。

    設計は $r_t$ 無次元で不変なので $R = r_t\,r_W(x_e; r_t)$ の $r_t$ だけを解く。
    - prev_run なし (CFD 前): `prepare_ns` と同じ δ_r の経路 (`integral_delta_r`: problem の `deltastar_initializer` の k_f・熱条件・
      平滑化) と同じ壁の構築 (`build_physical_wall`、pw_upstream に従う) で、反復のたびに δ_r と物理壁を作り直し、
      $r_t\,r_W(x_e) = R$ を |残差| ≤ tol_R_m まで解く (2026-10-07 変更、plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2:
      以前は未較正・未平滑化の `integral_bl` を直接呼んでいて、生産の壁と δ_r が違った)。max_iter 回で収まらなければ例外
      (`SizingNotConverged`)。tol_R_m の既定 (None) は `SIZING_TOL_PRE_CFD_EXIT_M` (1e-5 m)、明示で渡せる。
      **CFD 前の寸法は初期見積もり。精密な寸法・感度の評価には未検証** (同 plan §6 U2d・既知の制約): `integral_bl` (RK45 rtol 1e-6) の
      δ_r は r_t に対して滑らかでなく (r_t を 1e-12 m 変えただけで δ_r(x_e) が ~1e-5 r_t 揺れる)、r_t の ±1e-10 m の 11 点での出口半径の
      局所の散らばりは 2.64e-6 m (1 次の傾向を除いた max − min、U2c `_band_ab/upoly/U2c.json`。rtol 1e-10 でも 5.09e-8 m)。
      `SIZING_TOL_M` (1e-9 m) には届かない (U2b)。原因は未確定、真値の誤差は未評価。
    - prev_run あり (NS 後の最終補正、変えていない): その run の抽出 δ_r(x_F) を使い、$r_t$ 依存は $Re^{-0.2}$ で補正。n_iter 回の反復。
      **この経路は固定回数 (n_iter) の反復で、tol_R_m も max_iter も判定しない** (渡しても使わない。残差も返さない)。NS 後の 1e-9 m の
      保証は `solve_rt_throat` の δ_E の全分布の経路だけ (同 plan §6「NS 後の経路の保証の範囲の訂正」)。
    - integral_rtol: **診断用** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c)。CFD 前の経路の `integral_bl` の RK45 の
      相対許容差を明示の引数で注入する (None = 既定 1e-6、今の振る舞いのまま)。積分を使わない NS 後の経路で指定したら例外
      (渡っていない口を作らない)。実効値は戻りの `delta_r_source.integral_rtol` (solve_ivp に渡った値)。
    戻り: dict(r_t_m, delta_exit_rt, source, delta_column [prev_run ありで読んだ列], iters, …)。CFD 前は残差・実効の許容差
    (`tol_R_m`・`tol_R_m_source`)・物理スロート・出口半径・δ_r の出典も (`_sizing_result`)。"""
    from ..probdef import load_problem
    from ..evaluate.runner_axismach import design_chain
    if integral_rtol is not None and prev_run is not None:
        raise ValueError("solve_rt: integral_rtol は CFD 前の経路 (prev_run なし) 専用 — NS 後の経路は積分法を使わない")
    p = load_problem(problem); d = design_chain(p)
    S0 = float(p.spec["r_throat"])
    rF = float(d["wall_inv"][-1, 1]); xF = float(d["wall_inv"][-1, 0])
    hist = []
    if prev_run is None:
        init = resolve_integral_initializer(problem)
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "exit", cfd_before=True)
        fp = _fixed_point_rt(p, d, float(R_exit_m), "exit", _sizing_delta_supplier(p, init, rtol=integral_rtol), S0, max_iter, tol)
        out = _sizing_result(fp, p, float(R_exit_m), "exit", tol_src)
        last = fp["iters"][-1]
        out.update(r_t_prev_m=S0, delta_exit_rt=last["delta_r_exit_rt"], r_F_rt=rF, x_F_rt=xF, R_exit_m=R_exit_m,
                   source="integral_delta_r (prepare_ns と同じ経路: deltastar_initializer の k_f・熱条件・5 次 P-spline 平滑化) + "
                          "build_physical_wall", delta_column=None, delta_measured_xF_rt=None,
                   iters=[(h["r_t_m"], h["delta_r_exit_rt"], h["residual_m"]) for h in fp["iters"]],
                   initializer=init)
        return out
    else:
        prev_run = Path(prev_run)
        if not (prev_run / "delta_r_equiv.csv").exists():
            if euler_run is None: raise ValueError("prev_run に抽出結果が無く euler_run も未指定")
            deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=prev_run)
        # 未緩和の平滑化抽出値 δ_E (delta_r_next.csv の delta_E 列) の x_F 端を使う。2026-10-05 訂正 (codex result M1):
        # 以前は第 2 列 (ω 緩和後の次 pass 壁入力 δ_target) を読んでいて、出口の δ を約半分しか動かしていなかった。
        # delta_r_next.csv が無ければ生抽出の x_F−0.3
        if (prev_run / "delta_r_next.csv").exists():
            nx = read_delta_r_next(prev_run / "delta_r_next.csv")
            d_meas = float(np.interp(xF, nx["x_rt"], nx["delta_E"]))
            d_col = f"delta_E (unrelaxed smoothed extraction; delta_r_next.csv {nx['header_kind']} header, col 4)"
        else:
            e = np.genfromtxt(prev_run / "delta_r_equiv.csv", delimiter=",", names=True)
            d_meas = float(np.interp(xF - 0.3, e["x_rt"], e["delta_r_raw"]))
            d_col = "delta_r_raw at x_F-0.3 (delta_r_equiv.csv)"
        S_prev = float(json.loads((prev_run / "prepare_info.json").read_text())["scale_m"])
        rt = S_prev
        for k in range(n_iter):
            de = d_meas * (rt / S_prev) ** -0.2
            rt_new = R_exit_m / (rF + de); hist.append((rt, de, rt_new)); rt = rt_new
        src = f"measured delta_r from {prev_run.name} (Re^-0.2 scaling)"
    return dict(r_t_m=float(rt), r_t_prev_m=S0, delta_exit_rt=float(hist[-1][1]), r_F_rt=rF, x_F_rt=xF,
                R_exit_m=R_exit_m, source=src, delta_column=d_col, delta_measured_xF_rt=d_meas,
                iters=hist, length_m=float(xF * rt))


def solve_rt_throat(problem, R_throat_m: float, prev_run=None, delta_r_out=None, delta_next=None,
                    max_iter: int = SIZING_MAX_ITER, tol_R_m: float | None = None, _delta_supplier=None) -> dict:
    r"""**物理スロート径から寸法を決める**: $r_t\cdot\min_x r_W(x; r_t) = R_{throat}$ を $r_t$ について解く
    (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2、ユーザ決定 2026-10-07)。$r_W$ は物理壁 (r_t 単位、pw_upstream に従う)、
    最小は物理壁の大域最小 (`wall_global_min`)。スロート径を固定するときは出口径を同時に固定条件にしない (出口半径は戻り値に記録するだけ)。

    - prev_run なし (CFD 前): 反復のたびに `prepare_ns` と同じ経路 (`integral_delta_r`: k_f の cf_scale・熱条件・5 次 P-spline の
      平滑化・`delta_r_from_table`) で δ_r を作り、同じ壁の構築 (`build_physical_wall`) で最小半径を求める。k_f を較正し直したら
      寸法も解き直すこと (k_f はスロートの δ_r も変える)。
    - prev_run あり (NS 後): その run の `delta_r_next.csv` (delta_next で別の抽出 [例 `_extract_edge/delta_r_next.csv`] を指定可。
      S_prev は prev_run の prepare_info.json の scale_m) の **δ_E の全分布** (未緩和の平滑化抽出) から補正関数
      $\delta(x; r_t) = \delta_E(x)\,(r_t/S_{prev})^{-0.2}$ を作り、同じ関数で寸法と壁を作る ($Re^{-0.2}$ の換算は予測の近似)。
      解いた r_t での補正関数の表を delta_r_out (必須) に書く — 次の壁はこの表を `prepare_ns(delta_r_csv=...)` に渡して作る
      (寸法と壁で同じ関数)。delta_r_next.csv が無ければ例外。
    - 反復: 各反復で今の r_t の残差 r_t·min r_W − R を評価し、|残差| ≤ tol_R_m で止める (返す残差 = 最終の寸法で評価し直した値)。
      max_iter 回で収まらなければ例外 (`SizingNotConverged`、不合格)。tol_R_m の既定 (None) は経路で違う: CFD 前は
      `SIZING_TOL_PRE_CFD_THROAT_M` (1e-7 m)、NS 後と `_delta_supplier` は `SIZING_TOL_M` (1e-9 m、今までどおり)。明示で渡せる。
      **CFD 前の寸法は初期見積もり。精密な寸法・感度の評価には未検証** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d・既知の
      制約): CFD 前の δ_r (`integral_bl`、RK45 rtol 1e-6) は r_t を 1e-12 m 変えただけで出口で ~1e-5 r_t・スロートで ~1e-7 r_t 揺れ、
      r_t の ±1e-10 m の 11 点での物理スロート半径の局所の散らばりは 7.41e-9 m (U2c)。1e-9 m には届かない。原因は未確定、真値の誤差は未評価。
    - `_delta_supplier` (試験用): supplier(d, rt) → (delta_r_x, 記録) で δ_r の供給経路だけを差し替える (判別 A/B の腕 A)。
    戻り: r_t_m・残差・実効の許容差 (`tol_R_m`・`tol_R_m_source`)・反復の履歴・δ_r の出典と設定・物理スロートの位置と半径・
    そのとき決まる出口半径 ほか。"""
    from ..probdef import load_problem
    from ..evaluate.runner_axismach import delta_r_from_table, design_chain
    p = load_problem(problem); d = design_chain(p)
    S0 = float(p.spec["r_throat"])
    R = float(R_throat_m)
    extra = {}
    if _delta_supplier is not None:
        supplier, rt0 = _delta_supplier, S0
        extra["delta_r_route"] = "test supplier (_delta_supplier)"
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=False)
    elif prev_run is None:
        init = resolve_integral_initializer(problem)
        supplier, rt0 = _sizing_delta_supplier(p, init), S0
        extra.update(delta_r_route="cfd_before: integral_delta_r (prepare_ns と同じ経路)", initializer=init)
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=True)
    else:
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=False)
        prev_run = Path(prev_run)
        if delta_r_out is None:
            raise ValueError("solve_rt_throat (NS 後): delta_r_out (次の壁に渡す補正関数の表の出力先) が必要 — 寸法と壁を同じ関数で作るため")
        dn = Path(delta_next) if delta_next is not None else prev_run / "delta_r_next.csv"
        if not dn.exists():
            raise ValueError(f"solve_rt_throat (NS 後): {dn} が無い (extract_and_merge で作る — δ_E の全分布が要る)")
        nx = read_delta_r_next(dn)
        xE, dE = np.asarray(nx["x_rt"], dtype=float), np.asarray(nx["delta_E"], dtype=float)
        if not np.all(np.isfinite(dE)):
            raise ValueError(f"{dn} の delta_E に非有限値がある")
        S_prev = float(json.loads((prev_run / "prepare_info.json").read_text())["scale_m"])
        out_csv = Path(delta_r_out)

        def supplier(d_, rt, _x=xE, _dE=dE, _S=S_prev, _out=out_csv):
            # 表を書いて読み直し、prepare_ns(delta_r_csv=...) と同じ読み方・同じ関数で壁を作る (寸法と壁で同じ関数)
            fac = (rt / _S) ** -0.2
            np.savetxt(_out, np.c_[_x, _dE * fac], delimiter=",", comments="",
                       header=f"x_rt,delta_r (= delta_E x (r_t/S_prev)^-0.2; r_t={rt!r} m, S_prev={_S!r} m, factor={fac!r}; "
                              f"from {dn} col delta_E; solve_rt_throat)")
            tbl = np.loadtxt(_out, delimiter=",", skiprows=1)
            return delta_r_from_table(tbl[:, 0], tbl[:, 1]), {
                "kind": "delta_E 全分布 × (r_t/S_prev)^-0.2 (Re^-0.2 換算は予測の近似)", "prev_run": str(prev_run), "delta_next": str(dn),
                "delta_column": f"delta_E ({nx['header_kind']} header)", "S_prev_m": _S, "factor": fac, "table": str(_out)}
        rt0 = S_prev
        extra.update(delta_r_route="ns_after: delta_E 全分布", delta_r_out=str(out_csv))
    fp = _fixed_point_rt(p, d, R, "throat", supplier, rt0, max_iter, tol)
    out = _sizing_result(fp, p, R, "throat", tol_src)
    out.update(r_t_problem_m=S0, R_throat_m=R, **extra)
    return out


def run_pass(problem, euler_ref, prev_run, run_dir, omega: float = 0.5, ic_from=None,
             prepare_only: bool = False, nsteps=None, cfl_main=None, implicit_relax=None,
             stages: str = "full", ramp=None, ramp_steps: int = 1000) -> dict:
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    prev_run = Path(prev_run); run_dir = Path(run_dir)
    summ = extract_and_merge(prev_run, euler_ref, omega=omega)
    print("extract(prev):", json.dumps({k: v for k, v in summ.items() if k != "massflow"}), flush=True)
    print("massflow(prev):", json.dumps(summ["massflow"]), flush=True)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from or prev_run,
                      delta_r_csv=prev_run / "delta_r_next.csv", offset="radial",
                      euler_ref=euler_ref, omega=omega, prev_run=prev_run,
                      cfl_main=cfl_main, implicit_relax=implicit_relax)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "mesh", "nStepOuter", "cfl_main")}, indent=1), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    # 固定点差: 新 run から抽出し、入力 (= 新 run の壁 − Euler 壁) と比較
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("fixed_point(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    _design_report(run_dir, euler_ref)
    return m


def resolve_integral_initializer(problem) -> dict:
    """`--init-integral` の積分法初期化: problem YAML の `deltastar_initializer` があればそれ、無ければ断熱 contur (cf_scale 1)。
    2026-10-07 までは YAML を読まずに常に {"model": "contur"} を渡していた (help・docstring と食い違い)。
    C2 で較正した k_f (cf_scale) を YAML に書いた生産問題で、prep_c2pin.py (k_f を明示) と同じ物理壁にするため
    (plan tooling-nozzle-throat-monotone-r2 §9 2026-10-07: 読まないと出口半径が 2.4 mm ずれた)。"""
    from ..evaluate.runner_axismach import load_problem
    cfg = load_problem(problem).raw.get("deltastar_initializer")
    if cfg is None:
        return {"model": "contur"}
    if not isinstance(cfg, dict) or "model" not in cfg:
        raise ValueError(f"deltastar_initializer は model を持つ mapping であること ({cfg!r})")
    return dict(cfg)


def run_pass0_integral(problem, euler_ref, run_dir, ic_from, initializer=None, prepare_only: bool = False,
                       nsteps=None, omega: float = 0.5, cfl_main=None, implicit_relax=None,
                       stages: str = "full", ramp=None, ramp_steps: int = 1000) -> dict:
    """pass 0: 積分法 (CONTUR) 初期壁で NS を立て、終了後に抽出して次 pass 用 delta_r_next.csv まで作る。
    initializer=None なら problem YAML の `deltastar_initializer` (無ければ断熱 contur)。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    run_dir = Path(run_dir)
    # thermal_bc は prepare_ns 側で spec.wall_thermal から決める (ここでは model のみ)
    init = initializer if initializer is not None else resolve_integral_initializer(problem)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from, initializer=init,
                      euler_ref=euler_ref, omega=omega, cfl_main=cfl_main, implicit_relax=implicit_relax)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "initializer", "mesh", "nStepOuter", "cfl_main")},
                     indent=1, default=str), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("extract(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    _design_report(run_dir, euler_ref)
    return m


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="δ_r 固定点反復 (固定 Euler 基準)")
    ap.add_argument("--problem", required=True)
    ap.add_argument("--euler-ref", required=True)
    ap.add_argument("--prev", default=None, help="前 pass の NS run (抽出元・IC 既定)。--init-integral では不要")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--omega", type=float, default=0.5)
    ap.add_argument("--ic-from", default=None)
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--extract-only", action="store_true", help="抽出と delta_r_next.csv だけ作る")
    ap.add_argument("--init-integral", action="store_true",
                    help="pass 0: 積分法 (CONTUR) 初期壁で NS を立てる (YAML の deltastar_initializer を使う)")
    ap.add_argument("--init-thermal", default=None, help="廃止 (2026-09-12): 壁温は problem YAML の spec.wall_thermal が単一ソース")
    ap.add_argument("--cfl", type=float, default=None, help="本段 cfl (YAML evaluate.cfl_main を上書き)")
    ap.add_argument("--implicit-relax", type=float, default=None, help="implicitRelax (YAML evaluate.implicit_relax を上書き)")
    ap.add_argument("--stages", default="full", choices=("full", "none", "ramp"), help="起動: full=soft/mid/本段, none=本段のみ, ramp=cfl を段階的に上げて本段")
    ap.add_argument("--ramp", default=None, help="--stages ramp の cfl 列 (例 '1,2,3.5')")
    ap.add_argument("--ramp-steps", type=int, default=1000)
    ap.add_argument("--solve-rt", type=float, default=None, metavar="R_EXIT_M",
                    help="出口物理半径 [m] を与えて r_t を解くだけ (--prev があれば実測 δ_r で、無ければ prepare_ns と同じ積分法の経路で)")
    ap.add_argument("--solve-rt-throat", type=float, default=None, metavar="R_THROAT_M",
                    help="物理スロート半径 [m] を与えて r_t を解くだけ (--prev があれば δ_E の全分布から作る補正関数で。"
                         "そのとき --delta-r-out に次の壁に渡す表を書く)")
    ap.add_argument("--delta-r-out", default=None, help="--solve-rt-throat --prev のとき、解いた r_t での補正関数の表の出力先")
    ap.add_argument("--delta-next", default=None, help="--solve-rt-throat --prev のとき、δ_E を読む delta_r_next.csv (既定 PREV/delta_r_next.csv)")
    a = ap.parse_args(argv)
    if a.solve_rt is not None:
        r = solve_rt(a.problem, a.solve_rt, prev_run=a.prev, euler_run=a.euler_ref)
        print(json.dumps(r, indent=1, default=str))
        return 0
    if a.solve_rt_throat is not None:
        r = solve_rt_throat(a.problem, a.solve_rt_throat, prev_run=a.prev, delta_r_out=a.delta_r_out, delta_next=a.delta_next)
        print(json.dumps(r, indent=1, default=str))
        return 0
    if a.extract_only:
        s = extract_and_merge(a.prev, a.euler_ref, omega=a.omega)
        print(json.dumps(s, indent=1))
        return 0
    if a.init_integral:
        init = None
        if a.init_thermal:
            ap.error("--init-thermal は廃止。spec.wall_thermal (problem YAML) で指定する")
        ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
        run_pass0_integral(a.problem, a.euler_ref, a.run_dir, a.ic_from, initializer=init,
                           prepare_only=a.prepare_only, nsteps=a.steps, omega=a.omega, cfl_main=a.cfl,
                           implicit_relax=a.implicit_relax, stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps)
        return 0
    if not a.prev:
        ap.error("--prev が必要 (--init-integral でなければ)")
    ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
    run_pass(a.problem, a.euler_ref, a.prev, a.run_dir, omega=a.omega, ic_from=a.ic_from,
             prepare_only=a.prepare_only, nsteps=a.steps, cfl_main=a.cfl, implicit_relax=a.implicit_relax,
             stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps)
    return 0



def _design_report(run_dir, euler_ref):
    """ノズル設計の標準出力 (procedures/nozzle-design-outputs.md) を NS pass の後に自動で作る。失敗しても chain は止めない。"""
    if os.environ.get("FORGE_NO_DESIGN_REPORT"):
        return
    try:
        from ..report.nozzle_report import make_report
        print("design report:", make_report(run_dir, euler_ref), flush=True)
    except BaseException as e:  # noqa: BLE001
        print(f"design report: 失敗 ({type(e).__name__}: {e}) — chain は続行", flush=True)


if __name__ == "__main__":
    sys.exit(main())
