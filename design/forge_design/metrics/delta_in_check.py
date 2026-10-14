"""抽出器の δ_in と NS の壁の入力 δ_r の照合 (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36b、
codex diagnose `notes/reviews/2026-10-14-half-inlet-370k-pass-inputs-diagnose.md`)。

δ* ループの抽出器 (`deltastar.deltastar_from_core_matched_euler`) は NS の各列で
δ_in = (NS の壁の節点の半径) − (Euler の参照の壁の半径をその列の x へ列間で線形補間した値) と定義する
(`deltastar.py` の `euler_q_at` と `band_local_deficit` の `d_in`)。固定の Euler の参照 (補正 δ_r ≡ 0 の物理壁、
`runner_axismach.euler_reference_wall`) が NS の物理壁と同じ設計壁 S の上に作られていれば、x ≥ 0 では δ_in は NS の壁に入れた
δ_r(x) と一致するはず。違いになりうるのは Euler の壁の列間の線形補間・座標の丸め・1 本の B-spline の表現の差だけ。
この照合を格子の座標だけで行う (場は要らないので CFD 0 step でよい):

- 被覆: NS のすべての列が Euler の補間できる範囲の内にあること (外の列は抽出器が抽出しない)。範囲は抽出器・受け渡しの抽出と同じ
  共通の関数 `deltastar.euler_cover_range` (§5.1 #36g、remap_edge="missing" の経路と同じ edge_tol=True): Euler の端の外側だけに
  幾何の許容 (端の列の間隔 × 1e−4) を足し、許容の内側で端を越える列は端の列の値で評価する (同一の端点の丸めを吸収する許容 =
  微小区間の定数の延長。2026-10-14 の B の dry-run で出口の列が 1.68e−9 r_t 越えて被覆が落ちた件: 原因は Euler の格子の msh の 10 桁の丸め)。
- x ≥ x_min (既定 0) の各列で |δ_in − δ_r(x)| ≤ tol。tol の既定 `DELTA_IN_TOL_MM` = 0.001 mm は、測る前に決めた値:
  codex diagnose 2 が参照の取り替えによる δ_E の差と δ_E・δ_target の時間の変動に置いた暫定の予算 0.01 mm の 1/10
  (参照そのものの幾何の食い違いが、予算の 1 割を超えて使わないこと)。
- x < 0 は判定しない (pw_upstream poly の物理壁は上流に δ_r を足さず、縮流部は 5 次多項式 Q。δ_in は NS と Euler の Q の差
  = x = 0 で δ_r の (値, 1 階, 2 階) に C² でつながる 5 次 Hermite になる)。記録だけ出す。

`compare_delta_in` は配列 (r_t 単位) の核、`check_runs` は prepare 済みの run (nozzle.h5 の /MESH/COORD と prepare_info.json) を読む。
`python -m forge_design.metrics.delta_in_check NS_RUN EULER_RUN [--delta-r-csv CSV] [--out JSON]` で回す (合格 0・不合格 1)。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

DELTA_IN_TOL_MM = 1.0e-3      # [mm] 予算 0.01 mm の 1/10 (上の docstring)
REGIONS = {"throat_0_1": (0.0, 1.0), "x_1_40": (1.0, 40.0), "test_40_94": (40.0, 94.0), "end_94_exit": (94.0, np.inf)}


def euler_wall_at(xE, rwE, x: float, rng: dict | None = None) -> float:
    """抽出器 (`deltastar_from_core_matched_euler` の `euler_q_at`) と同じ式で Euler の壁の半径を x へ補間する
    (列 k = clip(searchsorted(xE, x) − 1)、重み w を [0, 1] に丸めた線形補間)。範囲外 (抽出器が列を外す x) は NaN。
    範囲 rng は `deltastar.euler_cover_range(xE, edge_tol=True)` (受け渡しの抽出の remap_edge="missing" と同じ; 省けばここで作る)。
    許容の内側で端を越える x は w の丸めで端の列の値になる。"""
    from .deltastar import euler_cover_range, euler_covers
    x = float(x)
    if rng is None:
        rng = euler_cover_range(xE, edge_tol=True)
    if not euler_covers(x, rng):
        return float("nan")
    k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
    w = (x - xE[k]) / max(xE[k + 1] - xE[k], 1e-30)
    w = float(np.clip(w, 0.0, 1.0))
    return float((1 - w) * rwE[k] + w * rwE[k + 1])


def _max_rec(err, x, mask, scale_m) -> dict:
    if not np.any(mask):
        return {"n": 0}
    e = np.abs(err[mask])
    i = int(np.argmax(e))
    return {"n": int(mask.sum()), "max_abs_rt": float(e[i]), "max_abs_mm": float(e[i] * scale_m * 1e3), "x_at_max": float(x[mask][i]),
            "rms_rt": float(np.sqrt(np.mean(e ** 2)))}


def compare_delta_in(x_ns, rw_ns, x_eu, rw_eu, delta_r, scale_m: float, x_min: float = 0.0, tol_mm: float = DELTA_IN_TOL_MM,
                     r_ns_exact=None, r_eu_exact=None, return_arrays: bool = False) -> dict:
    """抽出器の δ_in と NS の壁の入力 δ_r を照合する (r_t 単位の配列)。

    x_ns・rw_ns: NS の列の x (軸の節点 = 壁の station) と壁の節点の半径。x_eu・rw_eu: Euler の列の x と壁の節点の半径。
    delta_r: callable(x) → NS の壁に入れた δ_r(x) [r_t] (prepare_ns の `delta_r_from_table(CSV)`)。scale_m: r_t [m] (mm 換算だけ)。
    r_ns_exact・r_eu_exact: 任意の callable(x) (NS・Euler の物理壁の r(x))。あれば内訳 (Euler の壁の列間の補間の誤差・NS の壁の
    節点と S + δ_r の差・真の壁の間隔と δ_r の差) も出す。return_arrays: 列ごとの x・δ_in・δ_r を "arrays" に入れる (試験・作図用)。
    戻り: {"pass", "covered" (全列が Euler の範囲の内), "judged" (x ≥ x_min の列の max|δ_in − δ_r| ほか), "by_region",
    "upstream_x_lt_min" (記録だけ), "decomposition" (exact があれば), "tol_mm", ...}。"""
    x_ns = np.asarray(x_ns, dtype=float); rw_ns = np.asarray(rw_ns, dtype=float)
    x_eu = np.asarray(x_eu, dtype=float); rw_eu = np.asarray(rw_eu, dtype=float)
    if x_ns.shape != rw_ns.shape or x_eu.shape != rw_eu.shape or x_ns.ndim != 1 or x_eu.ndim != 1 or len(x_eu) < 2:
        raise ValueError("compare_delta_in: x・r は同じ長さの 1 次元の配列 (Euler は 2 列以上)")
    if not np.all(np.diff(x_eu) > 0):
        raise ValueError("compare_delta_in: Euler の列の x が単調増加でない")
    from .deltastar import euler_cover_range
    rng = euler_cover_range(x_eu, edge_tol=True)
    rwE_at = np.array([euler_wall_at(x_eu, rw_eu, x, rng) for x in x_ns])
    d_in = rw_ns - rwE_at                                    # 抽出器の δ_in (同じ式・同じ順の演算)
    dr = np.asarray(delta_r(x_ns), dtype=float)
    err = d_in - dr
    unc = ~np.isfinite(rwE_at)
    et = ~unc & ((x_ns < rng["x0"]) | (x_ns > rng["x1"]))     # 許容の内側で Euler の端を越え、端の列の値で評価した列 (記録)
    tol_rt = float(tol_mm) * 1e-3 / float(scale_m)
    sel = (x_ns >= x_min) & ~unc
    judged = _max_rec(err, x_ns, sel, scale_m)
    ok_cov = not bool(np.any(unc))
    ok_err = bool(judged.get("n", 0) > 0 and judged["max_abs_rt"] <= tol_rt)
    out = {"pass": bool(ok_cov and ok_err), "covered": ok_cov, "n_cols_ns": int(len(x_ns)), "n_cols_euler": int(len(x_eu)),
           "n_uncovered": int(unc.sum()), "uncovered_x": [float(v) for v in x_ns[unc][:5]],
           "range_ns": [float(x_ns[0]), float(x_ns[-1])], "range_euler": [float(x_eu[0]), float(x_eu[-1])],
           "euler_cover": rng, "n_edge_tol": int(et.sum()), "edge_tol_x": [float(v) for v in x_ns[et][:5]],
           "edge_tol_beyond": [float(max(rng["x0"] - v, v - rng["x1"])) for v in x_ns[et][:5]],
           "x_min": float(x_min), "tol_mm": float(tol_mm), "tol_rt": tol_rt, "scale_m": float(scale_m), "judged": judged,
           "by_region": {k: _max_rec(err, x_ns, sel & (x_ns >= lo) & (x_ns < hi), scale_m) for k, (lo, hi) in REGIONS.items()},
           "upstream_x_lt_min": {**_max_rec(err, x_ns, (x_ns < x_min) & ~unc, scale_m),
                                 "note": "判定しない (poly の物理壁は上流に δ_r を足さない; δ_in は NS と Euler の縮流部の多項式の差)"}}
    if r_ns_exact is not None and r_eu_exact is not None:
        rN = np.asarray(r_ns_exact(x_ns), dtype=float)
        rE = np.asarray(r_eu_exact(x_ns), dtype=float)
        rEn = np.asarray(r_eu_exact(x_eu), dtype=float)
        gap = rN - rE                                        # 真の壁の間隔 (補間なし)
        out["decomposition"] = {
            "euler_interp_minus_exact": _max_rec(rwE_at - rE, x_ns, ~unc, scale_m),
            "euler_nodes_minus_exact": _max_rec(rw_eu - rEn, x_eu, np.ones(len(x_eu), bool), scale_m),
            "ns_nodes_minus_exact": _max_rec(rw_ns - rN, x_ns, np.ones(len(x_ns), bool), scale_m),
            "exact_gap_minus_delta_r_x_ge_min": _max_rec(gap - dr, x_ns, x_ns >= x_min, scale_m),
            "delta_in_minus_exact_gap_all": _max_rec(d_in - gap, x_ns, ~unc, scale_m),
            "note": "δ_in − δ_r = (Euler の補間 − 真の Euler の壁) の符号反転 + (NS の節点 − 真の NS の壁) + (真の壁の間隔 − δ_r)"}
    if return_arrays:
        out["arrays"] = {"x": x_ns, "delta_in": d_in, "delta_r": dr, "r_wall_euler": rwE_at}
    return out


def wall_columns_from_run(run_dir) -> dict:
    """prepare 済みの node 構造格子の run (index = i*nj + j) の列の x と壁の節点の半径 [r_t] を読む。場 (res_*.h5) は要らない。
    読み方は抽出器の `_load_structured` と同じ (prepare_info.json の mesh.ni・nj・scale_m、nozzle.h5 の /MESH/COORD、列 x = 軸の節点)。"""
    import h5py
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5", "r") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
        dtype = str(f["/MESH/COORD"].dtype)
    if len(nc) != ni * nj:
        raise ValueError(f"{run_dir}: node 構造格子 run でない (COORD {len(nc)}, ni*nj {ni * nj})")
    x = nc[:, 0].reshape(ni, nj) / S
    r = nc[:, 1].reshape(ni, nj) / S
    return {"x": x[:, 0].copy(), "rw": r[:, -1].copy(), "scale_m": S, "info": info, "coord_dtype": dtype}


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_runs(ns_run, euler_run, delta_r_csv=None, tol_mm: float = DELTA_IN_TOL_MM, x_min: float = 0.0) -> dict:
    """prepare 済み (CFD 0 step でよい) の NS run と Euler の参照の run を照合する。delta_r_csv を省けば NS の prepare_info.json の
    `delta_r_csv`。δ_r は prepare_ns と同じ読み方 (`runner_axismach.delta_r_from_table`、CSV の先頭 2 列)。
    NS と Euler の r_t (scale_m) が違えば記録する (抽出は r_t 単位なので δ_in には効かないが、流量の比 `massflow_ratio` は実寸で
    取るので比が (S_NS/S_E)² だけずれる)。"""
    from ..evaluate.runner_axismach import delta_r_from_table
    N = wall_columns_from_run(ns_run)
    E = wall_columns_from_run(euler_run)
    csv = delta_r_csv if delta_r_csv is not None else N["info"].get("delta_r_csv")
    if not csv:
        raise ValueError(f"{ns_run}: delta_r_csv が無い (prepare_info.json にも無い) — 照合する入力が決まらない")
    tbl = np.loadtxt(csv, delimiter=",", skiprows=1)
    dr = delta_r_from_table(tbl[:, 0], tbl[:, 1])
    out = compare_delta_in(N["x"], N["rw"], E["x"], E["rw"], lambda x: dr(x), N["scale_m"], x_min=x_min, tol_mm=tol_mm)
    eu_info = E["info"]
    out.update({"ns_run": str(ns_run), "euler_run": str(euler_run), "delta_r_csv": str(csv), "delta_r_csv_sha256": _sha256(csv),
                "scale_m_ns": N["scale_m"], "scale_m_euler": E["scale_m"], "scale_same": bool(N["scale_m"] == E["scale_m"]),
                "coord_dtype": {"ns": N["coord_dtype"], "euler": E["coord_dtype"]},
                "euler_wall": {"physical_upstream": eu_info.get("physical_upstream"),
                               "kind": (eu_info.get("euler_wall") or {}).get("kind")}})
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="抽出器の δ_in と NS の壁の入力 δ_r の照合 (CFD 0 step、plan tooling-nozzle-isothermal-wall-chain §5.1 #36b)")
    ap.add_argument("ns_run")
    ap.add_argument("euler_run")
    ap.add_argument("--delta-r-csv", default=None)
    ap.add_argument("--tol-mm", type=float, default=DELTA_IN_TOL_MM)
    ap.add_argument("--x-min", type=float, default=0.0)
    ap.add_argument("--out", default=None, help="結果の JSON を書くファイル")
    a = ap.parse_args(argv)
    try:
        r = check_runs(a.ns_run, a.euler_run, a.delta_r_csv, tol_mm=a.tol_mm, x_min=a.x_min)
    except Exception as e:  # noqa: BLE001 — 入力の欠けは判定不能 (合格にしない)
        print(f"VERDICT: 判定不能 ({type(e).__name__}: {e})")
        return 2
    if a.out:
        Path(a.out).write_text(json.dumps(r, indent=1, ensure_ascii=False))
    j = r["judged"]
    print(json.dumps({k: r[k] for k in ("covered", "n_uncovered", "judged", "by_region", "scale_same")}, indent=1, ensure_ascii=False))
    cov = "全列" if r["covered"] else "外れ %d 列" % r["n_uncovered"]
    mx = j.get("max_abs_mm", float("nan"))
    cmp = "≤" if mx <= r["tol_mm"] else ">"
    print(f"VERDICT: {'PASS' if r['pass'] else 'FAIL'} (x ≥ {r['x_min']:g}: max|δ_in − δ_r| = {mx:.3g} mm {cmp} {r['tol_mm']:g} mm"
          f" (x = {j.get('x_at_max', float('nan')):.4g})、被覆 {cov})")
    return 0 if r["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
