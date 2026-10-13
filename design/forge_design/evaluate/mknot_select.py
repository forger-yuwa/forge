r"""軸 Mach 則 (knot) の knot のマッハ数 M_knot を、長さを固定したまま選ぶ (最大壁角の最小化、規則 `min_theta_max_v1`)。

plan: plans/active/tooling-nozzle-mknot-auto.md §4.1・§4.2 (§5.1 #2・#2b)。設計計算だけで CFD は回さない。

評価 (`evaluate_mknot`) は case/70.m5_air_wt_demo/mknot_sweep.py の `one()` と同じ組み合わせ:
design_chain + `inverse_design` の捕捉 + `characteristic_topology_diagnostics` (seam の除外は
case/42.isobutane_wt/compare_axislaw_D.py と同じ) + `wall_margin_diagnostics` + `wall_curvature_diagnostics`
+ prepare の入口の MOC のゲート (`runner_axismach.require_moc_gate` をそのまま呼ぶ)。新しい診断は足さない。

規則 (§4.1・§4.2):
  - 目的: 最大壁角 θ_max (`wall_qa` の `theta_max_deg`) の最小化 (M6 と同じ)。
  - 制約 (全部を数える。例外が出なかっただけでは成立にしない): 設計チェーンが通る (壁 QA・設計壁の検査・軸則のゲート)、
    網の内部の向きの反転 0・特性線の交差 0 (標本)・点の重なり 0、網の最小 sin 角 ≥ 0.02、壁での余裕 min(μ_w − θ_w) ≥ 1°、
    軸則の `mpp_quality_ok` (M″ の符号反転 ≤ 1)、prepare の入口と同じ MOC のゲート (`require_moc_gate` と同じ判定;
    `moc_corrector: fixed2` は `pass = None` のまま通す)。
  - 探索: 格子の原点を共通にした粗い格子 M_A + 0.1 + 0.1 k (M_d − 0.1 まで) を全部評価し、成立点の上の θ_max の局所極小を
    全部 (成立の区間の端を含む) 拾って、それぞれ ±0.1 を同じ原点の 0.01 の格子で詰める。
  - 境目: 評価した点を M_K の順に並べ、成立と不成立が隣り合う区間を二分法で幅 ≤ 0.001 まで詰める。境目の位置は
    区間の成立側の端 (保守的な端)。候補から境目までの距離は全部の境目の成立側の端までの距離の最小。
  - 選ぶ順序 (2026-10-14 確定、諮問 mknot-result-v6 の Major): (1) 制約が全部成立し、境目から ΔM_safe (0.05) 以上の候補に
    絞る (各境目 + ΔM_safe の点を評価して候補に足す) → (2) その中の θ_max の最小 → (3) 同点の規則: (1) の候補で (2) の最小との
    差 ≤ 1e-3° の集合から境目から最も遠い点。距離で除外される点の θ_max は選択を左右しない (記録の theta_min は距離の前の素の最小)。
  - 選択値は評価した点の値のまま (格子の原点 M_A + 0.1 からの端数を丸めない; 表示だけ丸める)。
  - 成立点が無ければ例外 (表つき)。黙って 0.45 × M_d には戻さない。
  - 解像度のゲート (`resolution_gate`): n_axis_inv を 2 倍にして選んだ点の ±0.1 を 0.01 で詰め直し、境目も求め直す。
    N と 2N の境目の区間の保守的な側を共通の境目 (`conservative_boundaries`) とし、その境目で N と 2N の選択をそれぞれ
    選び直して比べる: 選択値の差 ≤ 0.01、両方の解像度で制約と境目からの距離が成立、最小が詰めの端に無いこと。
    満たさなければ「選択は未確定」。

候補 (選ばれうる点) は粗い格子・詰めの格子・境目 + ΔM_safe の点 (と解像度のゲート・V6 の格子の点)。二分法の点は境目を
決めるためだけに使い、候補にしない。
"""
from __future__ import annotations

import copy
import functools
import math
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml

RULE_NAME = "min_theta_max_v1"
SIN_MIN = 0.02            # 網の最小 sin 角の床 (sweep_axislaw_A_MK.py の TOPO_FLOOR)
MARGIN_MIN_DEG = 1.0      # 壁での余裕 min(μ_w − θ_w) の床 [°] (同 MARGIN_FLOOR_DEG)
CANDIDATE_ROLES = ("coarse", "fine", "safe", "grid", "window")
_EPS = 1e-9


class MknotSelectionError(ValueError):
    """選べない (成立点が無い・境目からの距離を満たす点が無い)。`table` に全評価点。"""

    def __init__(self, msg: str, table: list | None = None):
        self.table = list(table or [])
        super().__init__(msg + ("\n" + format_table(self.table) if self.table else ""))


class MknotUndeterminedError(MknotSelectionError):
    """解像度のゲートが不合格 (選択は未確定)。`result` にゲートの結果。"""

    def __init__(self, msg: str, result: dict):
        self.result = result
        super().__init__(msg, (result.get("rows_N") or []) + (result.get("rows_2N") or []))


# --- 評価 -----------------------------------------------------------------------------------------------------------
def feasibility_reasons(row: dict) -> list[str]:
    """§4.1 の制約を全部数えて、成立しない理由を返す (空なら成立)。値が無い・有限でない制約は「判定不能」として不成立。"""
    rs = []
    if row.get("error"):
        rs.append(f"設計チェーンが例外: {row['error']}")
        return rs
    if row.get("qa_violations"):
        rs.append("壁 QA の違反: " + "; ".join(map(str, row["qa_violations"])))
    if row.get("law_violations"):
        rs.append("軸則のゲートの違反: " + "; ".join(map(str, row["law_violations"])))

    def _count(key, label):
        v = row.get(key)
        if v is None or (isinstance(v, (int, float)) and v < 0):
            rs.append(f"{label}が判定不能 ({key} = {v!r})")
        elif v != 0:
            rs.append(f"{label} {v}")

    _count("n_flip_interior", "網の内部の向きの反転")
    _count("n_crossings_sampled", "特性線の交差 (標本)")
    _count("pileup_pairs_interior", "点の重なり")
    s = row.get("min_sin_angle")
    if s is None or not np.isfinite(s):
        rs.append(f"網の最小 sin 角が判定不能 ({s!r})")
    elif s < SIN_MIN:
        rs.append(f"網の最小 sin 角 {s:.4f} < {SIN_MIN}")
    m = row.get("margin_deg")
    if m is None or not np.isfinite(m):
        rs.append(f"壁での余裕が判定不能 ({m!r})")
    elif m < MARGIN_MIN_DEG:
        rs.append(f"壁での余裕 {m:.3f}° < {MARGIN_MIN_DEG}°")
    q = row.get("mpp_quality_ok")
    if q is not True:
        rs.append(f"軸則の M″ の符号反転 {row.get('mpp_sign_changes')!r} > 1 (mpp_quality_ok = {q!r})")
    g = row.get("moc_gate_ok")                      # moc_gate_result (require_moc_gate そのもの) の判定
    if g is None:
        rs.append("prepare の入口の MOC のゲートが判定不能 (moc_gate_ok が無い)")
    elif g is not True:
        rs.append(f"prepare の入口の MOC のゲート (require_moc_gate) が不成立: {row.get('moc_gate_reason')}")
    t = row.get("theta_max_deg")
    if t is None or not np.isfinite(t):
        rs.append(f"最大壁角が有限でない ({t!r})")
    return rs


def moc_gate_result(p, d: dict) -> dict:
    """prepare の入口の MOC のゲートを、計算準備と同じ関数 `runner_axismach.require_moc_gate(p, d)` で判定する
    (判定の条件を写して二重に持たない)。`moc_corrector: converge` で診断が無い・`pass` が True でなければ不成立、
    `fixed2` (`pass = None`) は計算準備と同じに通す。p は `geometry` を持つもの (Problem)、d は design_chain の戻り。
    戻り: {"moc_gate_ok": bool, "moc_gate_reason": 不成立の理由 (成立なら None)}。"""
    import forge_design.evaluate.runner_axismach as RA
    try:
        RA.require_moc_gate(p, d)
    except ValueError as e:
        return {"moc_gate_ok": False, "moc_gate_reason": str(e)}
    return {"moc_gate_ok": True, "moc_gate_reason": None}


def evaluate_mknot(problem_dict: dict, M_K: float, n_axis_inv: int | None = None, *,
                   base_dir: str | Path | None = None, n_crossing_sample: int = 20) -> dict:
    """問題 (dict) の `geometry.M_knot` を M_K に (と `geometry.n_axis_inv` を上書き) して design_chain を回し、
    制約の値と成立の判定を返す。例外は行に記録する (成立しない理由)。

    base_dir: 問題の相対パス (initial_line_run など) の基準。一時の YAML はシステムの一時ディレクトリ ($TMPDIR) に書き、
    読んだ後の `Problem.path` を base_dir の中の仮の名前にする (case の中にファイルを作らない)。
    n_crossing_sample: 特性線の交差の標本の数 (mknot_sweep.py の one() と同じ 20)。"""
    import forge_design.evaluate.runner_axismach as RA
    from forge_design.geometry.moc_diagnostics import (characteristic_topology_diagnostics, wall_curvature_diagnostics,
                                                        wall_margin_diagnostics)
    from forge_design.probdef import load_problem

    t0 = time.time()
    raw = copy.deepcopy(problem_dict)
    raw.setdefault("geometry", {})["M_knot"] = float(M_K)
    if n_axis_inv:
        raw["geometry"]["n_axis_inv"] = int(n_axis_inv)
    raw.pop("_card", None)
    row = {"M_K": float(M_K), "n_axis": None}
    caps = []
    orig = RA.inverse_design

    def wrapped(*a, **k):
        r = orig(*a, **k)
        caps.append((r, k.get("n_axis")))
        return r

    RA.inverse_design = wrapped
    try:
        with tempfile.TemporaryDirectory(prefix="forge_mknot_") as td:
            tmp = Path(td) / "problem.yaml"
            tmp.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
            p = load_problem(str(tmp))
            if base_dir is not None:
                p.path = str(Path(base_dir).resolve() / "_mknot_select_eval.yaml")
            d = RA.design_chain(p)
        w = d["wall_inv"]
        x_A = float(d["x_A"])
        # from_length は逆 MOC を何回も回す。design_chain が採用した壁 (d["wall_inv"]) を返した呼び出しの網を見る
        hit = [(r, n) for r, n in caps if r.get("wall") is w]
        if not hit:
            raise RuntimeError("design_chain の壁に対応する逆 MOC の網が捕捉できない (網の診断の対象が不明)")
        res, n_cap = hit[-1]
        n_axis = int(n_cap or p.geometry.get("n_axis_inv", 500))
        topo = characteristic_topology_diagnostics(res["levels"], n_axis - 1, n_crossing_sample=n_crossing_sample,
                                                   seam_i_pad=90, seam_k_max=110, seam_x_max=x_A + 1.0, wall=w[:, :2])
        mg = wall_margin_diagnostics(w)
        wc = wall_curvature_diagnostics(d["wall"], 0.0, float(d["wall"].x_e), x_split=0.25)
        gates = d.get("gates") or {}
        moc_gate = (d.get("moc") or {}).get("gate") or {}
        row.update({
            "L_c": float(d["L_c"]), "x_K": float(d["x_K"]), "x_F": float(d["qa"]["x_F"]),
            "theta_max_deg": float(d["qa"]["theta_max_deg"]),
            "margin_deg": float(mg["min_mu_minus_theta_deg"]), "x_at_min_margin": float(mg["x_at_min_margin"]),
            "min_sin_angle": float(topo.get("min_sin_angle_interior", float("nan"))),
            "n_flip_interior": int(topo.get("n_orientation_flip_interior", -1)),
            "n_crossings_sampled": int(topo.get("n_crossings_sampled", -1)),
            "pileup_pairs_interior": int(topo.get("pileup_pairs_interior", -1)),
            "mpp_sign_changes": (None if gates.get("mpp_sign_changes") is None else int(gates["mpp_sign_changes"])),
            "mpp_quality_ok": (None if gates.get("mpp_quality_ok") is None else bool(gates["mpp_quality_ok"])),
            "qa_violations": list(d["qa"].get("violations") or []),
            "law_violations": list(gates.get("violations") or []),
            "n_axis": n_axis,
            # 格子の原点と上端 (select_mknot の M_A・M_d): 軸のアンカーの M と、軸則・MOC に渡した M_d (出口較正の Md_moc_offset 込み)
            "M_A": float(d["anchor"][0]), "M_d_moc": float(d["Md_moc"]),
            # 制約: prepare の入口の MOC のゲート (require_moc_gate と同じ判定; 2026-10-14 諮問 mknot-result-v6 で追加)。
            # moc_gate_pass は診断の値そのもの (記録; fixed2 では None)
            "moc_corrector": RA._moc_keys(p.geometry)[1], "moc_gate_pass": moc_gate.get("pass"),
            **moc_gate_result(p, d),
            # 記録だけ (制約にしない): 下流の曲率
            "max_kappa_down": float(wc.get("downstream_max_abs_kappa", float("nan"))),
            "J_dkappa_down": float(wc.get("downstream_J_dkappa_ds2", float("nan"))),
        })
    except Exception as e:  # noqa: BLE001 — 失敗も表の行として残す
        row["error"] = f"{type(e).__name__}: {str(e)[:300]}"
    finally:
        RA.inverse_design = orig
    row["reasons"] = feasibility_reasons(row)
    row["ok"] = not row["reasons"]
    row["time_s"] = time.time() - t0
    return row


def make_eval_fn(problem_dict: dict, n_axis_inv: int | None = None, base_dir: str | Path | None = None):
    """`select_mknot` に渡す評価 (M_K → 行)。ProcessPoolExecutor に渡せるよう partial で作る。"""
    return functools.partial(evaluate_mknot, problem_dict, n_axis_inv=n_axis_inv, base_dir=base_dir)


def _feasible(row: dict) -> bool:
    t = row.get("theta_max_deg")
    return row.get("ok") is True and not row.get("reasons") and t is not None and bool(np.isfinite(t))


def _key(M: float) -> float:
    return round(float(M), 9)


def _call(args):
    fn, M = args
    try:
        r = dict(fn(M))
    except Exception as e:  # noqa: BLE001 — 評価の関数そのものの例外も不成立の行にする
        r = {"error": f"{type(e).__name__}: {str(e)[:300]}"}
        r["reasons"] = [f"評価が例外: {r['error']}"]
        r["ok"] = False
    r["M_K"] = float(M)
    return r


class Evaluator:
    """M_K → 行 の評価の窓口 (同じ M_K は 1 回だけ評価する)。workers > 1 なら ProcessPoolExecutor で並列。"""

    def __init__(self, fn, workers: int = 1, log=None):
        self.fn, self.workers, self.log = fn, int(max(1, workers)), log
        self.rows: dict[float, dict] = {}
        self.n_eval = 0
        self._pool = None

    def __enter__(self):
        if self.workers > 1:
            self._pool = ProcessPoolExecutor(max_workers=self.workers)
        return self

    def __exit__(self, *exc):
        if self._pool is not None:
            self._pool.shutdown()
            self._pool = None
        return False

    def seed(self, rows) -> None:
        """既に評価した行を入れる (同じ評価の関数・同じ解像度のものだけ)。"""
        for r in rows:
            self.rows.setdefault(_key(r["M_K"]), dict(r))

    def many(self, Ms, role: str) -> list[dict]:
        keys = [_key(m) for m in Ms]
        todo = sorted({k for k in keys if k not in self.rows})
        if todo:
            t0 = time.time()
            jobs = [(self.fn, k) for k in todo]
            if self._pool is not None and len(todo) > 1:
                out = list(self._pool.map(_call, jobs))
            else:
                out = [_call(j) for j in jobs]
            for k, r in zip(todo, out):
                r["role"] = role
                self.rows[k] = r
            self.n_eval += len(todo)
            if self.log:
                self.log(f"  [{role}] {len(todo)} 点を評価 ({time.time() - t0:.1f} s): "
                         + ", ".join(f"{k:.6g}{'' if _feasible(self.rows[k]) else '×'}" for k in todo))
        return [self.rows[k] for k in keys]

    def table(self) -> list[dict]:
        return [self.rows[k] for k in sorted(self.rows)]


# --- 規則の部品 -------------------------------------------------------------------------------------------------------
def coarse_local_minima(rows: list[dict]) -> list[int]:
    """M_K の順の行 (粗い格子) の上で、成立点の θ_max の局所極小の添字を全部返す。隣が不成立・格子の端の側は比べない
    (成立の区間の端を含む)。等しい値の並びは全部を極小に数える (≤)。"""
    out = []
    for i, r in enumerate(rows):
        if not _feasible(r):
            continue
        t = r["theta_max_deg"]
        nb = [rows[j] for j in (i - 1, i + 1) if 0 <= j < len(rows) and _feasible(rows[j])]
        if all(t <= q["theta_max_deg"] for q in nb):
            out.append(i)
    return out


def find_boundaries(ev: Evaluator, bisect_tol: float, role: str = "bisect") -> list[dict]:
    """評価した点を M_K の順に並べ、成立と不成立が隣り合う区間を全部、二分法で幅 ≤ bisect_tol まで詰める。

    戻り: [{lo, hi, width, kind, feasible_end, reasons_infeasible}]。kind = 'lower' (下が不成立・上が成立 =
    成立の範囲の下の境目) / 'upper'。feasible_end は区間の成立側の端 (保守的な境目)。"""
    while True:
        pts = ev.table()
        brk = [(a, b) for a, b in zip(pts[:-1], pts[1:]) if _feasible(a) != _feasible(b)]
        open_ = [(a, b) for a, b in brk if b["M_K"] - a["M_K"] > bisect_tol + _EPS]
        if not open_:
            break
        mids = [0.5 * (a["M_K"] + b["M_K"]) for a, b in open_]
        ev.many(mids, role)
    out = []
    for a, b in brk:
        lower = _feasible(b)
        inf = a if lower else b
        out.append({"lo": float(a["M_K"]), "hi": float(b["M_K"]), "width": float(b["M_K"] - a["M_K"]),
                    "kind": "lower" if lower else "upper",
                    "feasible_end": float(b["M_K"] if lower else a["M_K"]),
                    "reasons_infeasible": list(inf.get("reasons") or [])})
    return out


def _same_boundaries(a: list[dict], b: list[dict]) -> bool:
    return [(x["kind"], x["lo"], x["hi"]) for x in a] == [(x["kind"], x["lo"], x["hi"]) for x in b]


def boundary_distance(M: float, boundaries: list[dict]) -> float:
    """候補 M から境目 (全部の境目の成立側の端) までの距離の最小。境目が無ければ inf。"""
    if not boundaries:
        return float("inf")
    return float(min(abs(float(M) - b["feasible_end"]) for b in boundaries))


def _pick(cands: list[dict], boundaries: list[dict], tie_deg: float) -> tuple[dict, dict, list[dict]]:
    """候補の中の θ_max の最小と、その最小との差 ≤ tie_deg の集合 (同点の集合) から境目から最も遠い点
    (同じ距離なら θ_max の小さい方・M_K の小さい方)。戻り: (選んだ点, 最小の点, 同点の集合)。"""
    tmin_row = min(cands, key=lambda c: (c["theta_max_deg"], c["M_K"]))
    ties = [c for c in cands if c["theta_max_deg"] <= tmin_row["theta_max_deg"] + tie_deg]
    best = min(ties, key=lambda c: (-boundary_distance(c["M_K"], boundaries), c["theta_max_deg"], c["M_K"]))
    return best, tmin_row, ties


def safe_points_for(boundaries: list[dict], dM_safe: float, safe_range=None) -> list[float]:
    """各境目の成立側の端から成立の側へ ΔM_safe の点 (境目 + ΔM_safe の点; 下の境目なら +、上の境目なら −)。
    safe_range = (lo, hi) の外の点は除く (探索の範囲・解像度のゲートの窓の外を評価しない)。ΔM_safe ≤ 0 なら空
    (距離の制約が無いので足す点も無い。二分法の点を候補にしないため)。"""
    if dM_safe <= 0.0:
        return []
    lo, hi = safe_range if safe_range is not None else (-math.inf, math.inf)
    Ms = sorted({_key(b["feasible_end"] + (dM_safe if b["kind"] == "lower" else -dM_safe)) for b in boundaries})
    return [m for m in Ms if lo - _EPS <= m <= hi + _EPS]


def apply_rule(rows: list[dict], boundaries: list[dict], *, tie_deg: float, dM_safe: float, safe_eval=None,
               candidate_roles=CANDIDATE_ROLES, safe_range=None) -> dict:
    """§4.1・§4.2 の選び方を評価済みの行に当てる (2026-10-14 確定の順序、諮問 mknot-result-v6 の Major):
      (1) 制約が全部成立し、境目 (各境目の成立側の端) から ΔM_safe 以上の候補に絞る。各境目 + ΔM_safe の点を
          safe_eval で評価し、成立すれば候補に足す (足すが、必ずそれを選ぶわけではない)。
      (2) (1) の候補の中の θ_max の最小。
      (3) 同点の規則: (1) の候補のうち (2) の最小との差 ≤ tie_deg の集合から境目から最も遠い点。
    距離で除外される点の θ_max は選択を左右しない (先に全体で同点の集合を作ると、除外される点の値で選択が変わる反例がある)。
    記録の theta_min は距離の前の素の最小 (参考)。constrained は素の最小が距離で除外されたか。

    safe_eval: M_K の列 → 行の列 (境目 + ΔM_safe の点の評価)。None なら評価せず、評価済みの候補だけで選ぶ。
    safe_range: 境目 + ΔM_safe の点を評価する範囲 (lo, hi)。None なら全部の境目。
    候補が無ければ MknotSelectionError。選択値は評価した点の M_K のまま (丸めない)。"""
    cands = [r for r in rows if r.get("role") in candidate_roles and _feasible(r)]
    if not cands:
        raise MknotSelectionError("M_knot を選べない: 成立する候補が無い (§4.1 の制約を満たす点が 0)",
                                  sorted(rows, key=lambda r: r["M_K"]))
    raw_min = min(cands, key=lambda c: (c["theta_max_deg"], c["M_K"]))       # 素の最小 (距離の前; 記録だけ)
    safe_M = safe_points_for(boundaries, dM_safe, safe_range)
    safe_rows = list(safe_eval(safe_M)) if (safe_eval and safe_M) else []
    pool = {_key(c["M_K"]): c for c in cands}
    for r in safe_rows:
        if _feasible(r):
            pool.setdefault(_key(r["M_K"]), r)

    def dist(c):
        return boundary_distance(c["M_K"], boundaries)

    adm = sorted((c for c in pool.values() if dist(c) >= dM_safe - _EPS), key=lambda c: c["M_K"])       # (1)
    excl = sorted((c for c in pool.values() if dist(c) < dM_safe - _EPS), key=lambda c: c["M_K"])
    raw_excluded = dist(raw_min) < dM_safe - _EPS
    out = {"theta_min": {"M_K": raw_min["M_K"], "theta_max_deg": raw_min["theta_max_deg"], "distance": dist(raw_min),
                         "excluded_by_distance": raw_excluded},
           "safe_points": [{"M_K": r["M_K"], "ok": _feasible(r), "theta_max_deg": r.get("theta_max_deg"),
                            "distance": boundary_distance(r["M_K"], boundaries), "reasons": r.get("reasons")} for r in safe_rows],
           "excluded_by_distance": [c["M_K"] for c in excl], "n_admissible": len(adm)}
    if not adm:
        raise MknotSelectionError(
            f"M_knot を選べない: 境目から ΔM_safe = {dM_safe:g} 以上離れた成立点が無い (成立の候補 {len(pool)} 点は全部が境目から "
            f"ΔM_safe 未満; 素の最小 {raw_min['M_K']:.6g} は境目から {dist(raw_min):.4g}; 境目 + ΔM_safe の点 "
            f"{', '.join(f'{m:.6g}' for m in safe_M) or '—'} の成立 {[p['ok'] for p in out['safe_points']]})",
            sorted(rows + safe_rows, key=lambda r: r["M_K"]))
    best, amin, ties = _pick(adm, boundaries, tie_deg)                                                     # (2)・(3)
    d_best = dist(best)
    why = ("θ_max の最小" if best is amin else
           f"同点の規則 (最小 {amin['M_K']:.6g} との差 ≤ {tie_deg:g}° の中で境目から最も遠い点)")
    scope = f"境目から ΔM_safe {dM_safe:g} 以上の候補 {len(adm)} 点の中"
    if excl:
        scope += (f"; ΔM_safe 未満の {len(excl)} 点を除外"
                  + (f" (素の最小 {raw_min['M_K']:.6g} も境目から {dist(raw_min):.4g} で除外)" if raw_excluded else ""))
    if safe_M:
        scope += f"; 境目 + ΔM_safe の点 {', '.join(f'{m:.6g}' for m in safe_M)} を評価"
    out.update(chosen=best["M_K"], chosen_row=best, distance=d_best, constrained=raw_excluded,
               theta_min_admissible={"M_K": amin["M_K"], "theta_max_deg": amin["theta_max_deg"], "distance": dist(amin)},
               tie_set=[c["M_K"] for c in ties], reason=f"{why} ({scope})")
    return out


def format_table(rows: list[dict]) -> str:
    """評価点の表 (例外のメッセージ・ログ用)。"""
    def f(v, fmt):
        return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, fmt)
    lines = ["  M_K        成立  θ_max[°]   余裕[°]  min sin  M″反転  役割    理由"]
    for r in sorted(rows, key=lambda r: r["M_K"]):
        lines.append(f"  {r['M_K']:<10.6f} {'○' if _feasible(r) else '×':<4} {f(r.get('theta_max_deg'), '9.5f'):>9}  "
                     f"{f(r.get('margin_deg'), '7.3f'):>7}  {f(r.get('min_sin_angle'), '.4f'):>7}  "
                     f"{str(r.get('mpp_sign_changes', '—')):>5}   {r.get('role', ''):<7} {'; '.join(r.get('reasons') or [])[:160]}")
    return "\n".join(lines)


# --- 選び方 (§4.2) ----------------------------------------------------------------------------------------------------
def select_mknot(eval_fn, M_A: float, M_d: float, *, coarse: float = 0.1, fine: float = 0.01, tie_deg: float = 1e-3,
                 dM_safe: float = 0.05, bisect_tol: float = 1e-3, workers: int = 4, log=None) -> dict:
    """§4.2 の探索で M_knot を選ぶ。eval_fn: M_K → 行 (`evaluate_mknot` の形; `ok`・`reasons`・`theta_max_deg` を見る)。

    粗い格子 M_A + coarse + coarse·k (M_d − coarse まで) → 成立点の θ_max の局所極小を全部 ±coarse を fine で詰める
    (同じ原点の格子) → 境目を二分法で幅 ≤ bisect_tol → 制約と境目からの距離 → 最小 → 同点の規則 (`apply_rule`;
    境目 + ΔM_safe の点は探索の範囲 [M_A + coarse, 格子の上端] の中だけ評価する)。
    戻り: {rule, chosen, chosen_row, table (全評価点), local_minima, windows, boundaries, reason, warnings, ...}。
    成立点が無ければ MknotSelectionError (表つき)。"""
    t0 = time.time()
    ratio = int(round(coarse / fine))
    if ratio < 1 or abs(ratio * fine - coarse) > 1e-12:
        raise ValueError(f"coarse {coarse} は fine {fine} の整数倍にする (同じ原点の格子)")
    M0 = float(M_A) + coarse
    K = int(math.floor((float(M_d) - coarse - M0) / coarse + 1e-9))
    if K < 0:
        raise ValueError(f"探索の範囲が空 (M_A {M_A} + {coarse} > M_d {M_d} − {coarse})")
    J_max = K * ratio

    def grid(j: int) -> float:
        return round(M0 + j * fine, 12)

    warnings = []
    with Evaluator(eval_fn, workers, log) as ev:
        if log:
            log(f"粗い格子: {K + 1} 点 ({grid(0):.6g}〜{grid(J_max):.6g}、{coarse:g} 刻み)")
        crow = ev.many([grid(k * ratio) for k in range(K + 1)], "coarse")
        if not any(_feasible(r) for r in crow):
            raise MknotSelectionError("M_knot を選べない: 粗い格子に成立点が無い (§4.1 の制約を満たす点が 0)", ev.table())
        mins = coarse_local_minima(crow)
        if len(mins) >= 2:
            warnings.append(f"θ_max の局所極小が {len(mins)} 個 (M_K " + ", ".join(f"{crow[i]['M_K']:.6g}" for i in mins)
                            + ") — 全部を詰めた")
        js, windows = set(), []
        for i in mins:
            lo, hi = max(0, i * ratio - ratio), min(J_max, i * ratio + ratio)
            windows.append([grid(lo), grid(hi)])
            js.update(range(lo, hi + 1))
        if log:
            log(f"局所極小 {len(mins)} 個の前後 ±{coarse:g} を {fine:g} で詰める")
        ev.many([grid(j) for j in sorted(js)], "fine")
        boundaries = find_boundaries(ev, bisect_tol)
        rule = None
        for it in range(3):
            rule = apply_rule(ev.table(), boundaries, tie_deg=tie_deg, dM_safe=dM_safe,
                              safe_eval=lambda Ms: ev.many(Ms, "safe"), safe_range=(grid(0), grid(J_max)))
            nb = find_boundaries(ev, bisect_tol)           # 境目 + ΔM_safe の点が新しい切り替わりを作ったら選び直す
            if _same_boundaries(nb, boundaries):
                break
            if it == 2:                                     # 記録 (選択と境目) は選択に使った境目のまま揃える
                warnings.append("境目 + ΔM_safe の点が 3 回の選び直しでも新しい成立の切り替わりを作る — 選択は直前の境目で判定した "
                                f"(新しい境目 {[(b['kind'], round(b['feasible_end'], 6)) for b in nb]})")
                break
            warnings.append("境目 + ΔM_safe の点が新しい成立の切り替わりを作った — 境目を求め直して選び直した")
            boundaries = nb
        table = ev.table()
        n_eval = ev.n_eval
    return {"rule": RULE_NAME, "M_A": float(M_A), "M_d": float(M_d),
            "settings": {"coarse": coarse, "fine": fine, "tie_deg": tie_deg, "dM_safe": dM_safe, "bisect_tol": bisect_tol,
                         "sin_min": SIN_MIN, "margin_min_deg": MARGIN_MIN_DEG, "grid_origin": M0},
            "chosen": float(rule["chosen"]), "chosen_row": rule["chosen_row"], "distance_to_boundary": rule["distance"],
            "reason": rule["reason"], "constrained": rule["constrained"], "theta_min": rule["theta_min"],
            "theta_min_admissible": rule["theta_min_admissible"], "tie_set": rule["tie_set"],
            "excluded_by_distance": rule["excluded_by_distance"], "safe_points": rule["safe_points"],
            "local_minima": [{"M_K": crow[i]["M_K"], "theta_max_deg": crow[i]["theta_max_deg"]} for i in mins],
            "windows": windows, "boundaries": boundaries, "warnings": warnings,
            "table": table, "n_eval": n_eval, "elapsed_s": time.time() - t0}


# --- 解像度のゲート (§4.2) --------------------------------------------------------------------------------------------
def conservative_boundaries(bN: list[dict], b2: list[dict], pair_tol: float) -> list[dict]:
    """N と 2N の境目から両方の解像度の共通の境目を作る (§4.2; 2026-10-14 諮問 mknot-result-v6 の Major)。
    同じ向き (kind) で成立側の端の差が pair_tol 以内の N と 2N の境目を近い順に対にし、成立の範囲を狭める側の成立側の端
    (下の境目は大きい方・上の境目は小さい方) を採る。対にならない境目はそのまま残す (片方の解像度にしか無い境目も
    成立の範囲を狭める側なので、残すのが保守的)。
    戻り: [{kind, feasible_end, lo, hi, width, source, paired, interval_N, interval_2N, shift_2N, reasons_infeasible}]
    (成立側の端の順)。source は採った端の解像度 ('N' / '2N'、同じ値なら 'N=2N')、shift_2N は対のときの
    2N の成立側の端 − N の成立側の端 (下の境目で正・上の境目で負なら 2N が保守側へ動いた)。"""
    pairs = sorted((abs(c["feasible_end"] - b["feasible_end"]), i, j) for i, b in enumerate(bN) for j, c in enumerate(b2)
                   if b["kind"] == c["kind"] and abs(c["feasible_end"] - b["feasible_end"]) <= pair_tol + _EPS)

    def entry(take, src, b, c):
        return {"kind": take["kind"], "feasible_end": float(take["feasible_end"]), "lo": float(take["lo"]),
                "hi": float(take["hi"]), "width": float(take["width"]), "source": src, "paired": b is not None and c is not None,
                "interval_N": (None if b is None else [float(b["lo"]), float(b["hi"])]),
                "interval_2N": (None if c is None else [float(c["lo"]), float(c["hi"])]),
                "shift_2N": (None if b is None or c is None else float(c["feasible_end"] - b["feasible_end"])),
                "reasons_infeasible": list(take.get("reasons_infeasible") or [])}

    usedN, used2, out = set(), set(), []
    for _, i, j in pairs:
        if i in usedN or j in used2:
            continue
        usedN.add(i)
        used2.add(j)
        b, c = bN[i], b2[j]
        sh = c["feasible_end"] - b["feasible_end"]
        take_2N = sh > 0.0 if b["kind"] == "lower" else sh < 0.0
        out.append(entry(c if take_2N else b, "N=2N" if sh == 0.0 else ("2N" if take_2N else "N"), b, c))
    out += [entry(b, "N", b, None) for i, b in enumerate(bN) if i not in usedN]
    out += [entry(c, "2N", None, c) for j, c in enumerate(b2) if j not in used2]
    return sorted(out, key=lambda x: x["feasible_end"])


def reselect_common(evN: Evaluator, ev2: Evaluator, rows_N: list[dict], rows_2N: list[dict], *, tie_deg: float,
                    dM_safe: float, bisect_tol: float, pair_tol: float, safe_range=None, candidate_roles=CANDIDATE_ROLES,
                    max_iter: int = 3, log=None) -> dict:
    """N と 2N の共通の境目 (`conservative_boundaries`) で両方の解像度の選択をそれぞれ選び直す (§4.2; 2N の境目が保守側へ
    動いたときに N の境目のままで判定しない)。evN・ev2: 各解像度の Evaluator (境目を求め直し、境目 + ΔM_safe の点を評価する)。
    rows_N・rows_2N: 各解像度の候補の行 (candidate_roles で絞る)。境目 + ΔM_safe の点が新しい成立の切り替わりを作ったら
    境目を求め直して選び直す (max_iter 回まで)。
    戻り: {boundaries_N, boundaries_2N, common, N: apply_rule の戻り (選べなければ {chosen: None, error}), 2N: 同, warnings}。"""
    warnings = []
    bN, b2 = find_boundaries(evN, bisect_tol), find_boundaries(ev2, bisect_tol)
    res, common = {}, []
    for it in range(max_iter):
        common = conservative_boundaries(bN, b2, pair_tol)
        for tag, ev, rows in (("N", evN, rows_N), ("2N", ev2, rows_2N)):
            try:
                res[tag] = apply_rule(rows, common, tie_deg=tie_deg, dM_safe=dM_safe,
                                      safe_eval=lambda Ms, ev=ev: ev.many(Ms, "safe"),
                                      candidate_roles=candidate_roles, safe_range=safe_range)
            except MknotSelectionError as e:
                res[tag] = {"chosen": None, "error": str(e).splitlines()[0]}
        nbN, nb2 = find_boundaries(evN, bisect_tol), find_boundaries(ev2, bisect_tol)
        if _same_boundaries(nbN, bN) and _same_boundaries(nb2, b2):
            break
        if it == max_iter - 1:                          # 記録 (選択と境目) は選択に使った境目のまま揃える
            warnings.append(f"境目 + ΔM_safe の点が {max_iter} 回の選び直しでも新しい成立の切り替わりを作る — 選択は直前の"
                            "共通の境目で判定した")
            break
        warnings.append("境目 + ΔM_safe の点が新しい成立の切り替わりを作った — 共通の境目を求め直して選び直した")
        bN, b2 = nbN, nb2
    if log:
        log("  共通の境目: " + ", ".join(f"{b['kind']} {b['feasible_end']:.6f} ({b['source']}"
                                         + ("" if b["shift_2N"] is None else f", 2N−N {b['shift_2N']:+.6f}") + ")" for b in common)
            + f" → 選び直し N {res['N'].get('chosen')!r}・2N {res['2N'].get('chosen')!r}")
    return {"boundaries_N": bN, "boundaries_2N": b2, "common": common, "N": res["N"], "2N": res["2N"], "warnings": warnings}


def resolution_gate(problem_dict: dict | None, chosen: float, N: int, *, select_result: dict | None = None,
                    eval_factory=None, base_dir=None, fine: float = 0.01, half_width: float = 0.1,
                    tie_deg: float | None = None, dM_safe: float | None = None, bisect_tol: float | None = None,
                    max_change: float = 0.01, pair_tol: float | None = None, workers: tuple[int, int] = (4, 2),
                    raise_on_fail: bool = True, log=None) -> dict:
    """§4.2 の解像度のゲート。`n_axis_inv` = N と 2N で選んだ点 (chosen) の前後 ±half_width を fine で詰め (chosen を中心の格子)、
    境目を 2N でも求め直す (N の境目のうち chosen から half_width + ΔM_safe 以内のものは区間の両端も 2N で評価する)。
    N と 2N の境目の区間の保守的な側を共通の境目 (`conservative_boundaries`、対にする幅 pair_tol 既定 = half_width) とし、
    その境目で窓の候補から N と 2N の選択をそれぞれ選び直して比べる (`reselect_common`; 境目 + ΔM_safe の点は窓の中だけ評価)。

    判定 (全部で合格):
      change      共通の境目で選び直した N と 2N の選択値の差 ≤ max_change (0.01)
      ok_N・ok_2N 入力の選択点 (chosen) と選び直した 2 点が、その解像度で全部の制約を満たす
                  (選び直した 2 点は両方の解像度で評価する; 高い解像度だけ制約を破る選択点を通さない)
      distance    選び直した 2 点から共通の境目まで ≥ ΔM_safe
      edge_N・edge_2N  窓の中の成立点の θ_max の最小が窓の端に無い
    共通の境目で N を選び直した値 (chosen_N) は入力の chosen と違うことがある (共通の境目が保守側へ動いた、同点の規則の
    境目からの距離が変わった)。その差は moved_N に記録する (合否には入れない)。
    select_result: `select_mknot` の戻り (N の評価の表・境目・設定を再利用)。eval_factory: n_axis → 評価の関数
    (既定は make_eval_fn(problem_dict, n_axis, base_dir); 試験で差し替える)。raise_on_fail なら不合格で MknotUndeterminedError。
    窓の評価は並列 (workers)、選び直しの境目 + ΔM_safe の点と交差の評価は逐次 (少数の点; 2 つの並列を同時に持たない)。"""
    st = (select_result or {}).get("settings") or {}
    tie_deg = st.get("tie_deg", 1e-3) if tie_deg is None else tie_deg
    dM_safe = st.get("dM_safe", 0.05) if dM_safe is None else dM_safe
    bisect_tol = st.get("bisect_tol", 1e-3) if bisect_tol is None else bisect_tol
    pair_tol = half_width if pair_tol is None else pair_tol
    if eval_factory is None:
        def eval_factory(n):
            return make_eval_fn(problem_dict, n_axis_inv=n, base_dir=base_dir)
    t0 = time.time()
    nh = int(round(half_width / fine))
    window = [round(float(chosen) + j * fine, 12) for j in range(-nh, nh + 1)]
    out = {"N": int(N), "N2": int(2 * N), "chosen": float(chosen), "window": [window[0], window[-1]],
           "settings": {"fine": fine, "half_width": half_width, "tie_deg": tie_deg, "dM_safe": dM_safe,
                        "bisect_tol": bisect_tol, "max_change": max_change, "pair_tol": pair_tol}}
    evN = Evaluator(eval_factory(N), workers[0], log)
    if select_result:
        evN.seed(select_result.get("table") or [])
    with evN:
        if log:
            log(f"解像度のゲート: N = {N} の窓 {window[0]:.6g}〜{window[-1]:.6g}")
        rows_N = evN.many(window, "window")
        bN = find_boundaries(evN, bisect_tol)
    near = [b for b in bN if abs(b["feasible_end"] - chosen) <= half_width + dM_safe + _EPS]
    ev2 = Evaluator(eval_factory(2 * N), workers[1], log)
    with ev2:
        if log:
            log(f"解像度のゲート: 2N = {2 * N} の窓と N の境目の区間の両端 ({len(near)} 個)")
        rows_2N = ev2.many(window, "window")
        if near:
            ev2.many([v for b in near for v in (b["lo"], b["hi"])], "bracket")
        find_boundaries(ev2, bisect_tol)            # 2N の二分法は並列の中で (reselect_common は評価済みの点を使う)
    # 共通の境目で両方を選び直す (候補は窓の点だけ; N の窓の点は select_mknot の表の行を再利用するので役割を窓に揃える)
    rc = reselect_common(evN, ev2, [dict(r, role="window") for r in rows_N], [dict(r, role="window") for r in rows_2N],
                         tie_deg=tie_deg, dM_safe=dM_safe, bisect_tol=bisect_tol, pair_tol=pair_tol,
                         safe_range=(window[0], window[-1]), candidate_roles=("window", "safe"), log=log)
    common = rc["common"]
    cN, c2 = rc["N"].get("chosen"), rc["2N"].get("chosen")
    sel_pts = [float(chosen)] + [c for c in (cN, c2) if c is not None]
    at_N = {_key(r["M_K"]): r for r in evN.many(sel_pts, "cross")}          # 逐次 (with の外)
    at_2N = {_key(r["M_K"]): r for r in ev2.many(sel_pts, "cross")}

    def _edge(rows):
        f = [(j, r) for j, r in enumerate(rows) if _feasible(r)]
        if not f:
            return None
        j, _ = min(f, key=lambda t: (t[1]["theta_max_deg"], t[1]["M_K"]))
        return j - nh

    def _bad(at, pts):
        return [m for m in pts if not _feasible(at[_key(m)])]

    eN, e2 = _edge(rows_N), _edge(rows_2N)
    dN = None if cN is None else boundary_distance(cN, common)
    d2 = None if c2 is None else boundary_distance(c2, common)
    badN, bad2 = _bad(at_N, sel_pts), _bad(at_2N, sel_pts)
    checks = {
        "change": cN is not None and c2 is not None and abs(c2 - cN) <= max_change + _EPS,
        "ok_N": cN is not None and c2 is not None and not badN,
        "ok_2N": cN is not None and c2 is not None and not bad2,
        "distance": dN is not None and d2 is not None and min(dN, d2) >= dM_safe - _EPS,
        "edge_N": eN is not None and abs(eN) < nh, "edge_2N": e2 is not None and abs(e2) < nh,
    }
    reasons = []
    if not checks["change"]:
        reasons.append(f"共通の境目で選び直した 2N の値 {c2!r} が N の {cN!r} から {max_change:g} を超えて動く "
                       f"(N: {rc['N'].get('reason') or rc['N'].get('error')}; 2N: {rc['2N'].get('reason') or rc['2N'].get('error')})")
    for tag, n, bad, at in (("N", N, badN, at_N), ("2N", 2 * N, bad2, at_2N)):
        if not checks[f"ok_{tag}"]:
            reasons.append(f"{tag} = {n} で選択点が制約を満たさない: "
                           + ("選び直せない" if cN is None or c2 is None else
                              " / ".join(f"{m:.6g}: {'; '.join(at[_key(m)].get('reasons') or [])}" for m in bad)))
    if not checks["distance"]:
        reasons.append(f"選び直した点から共通の境目まで N {dN!r}・2N {d2!r} < ΔM_safe {dM_safe:g}")
    if not checks["edge_N"]:
        reasons.append(f"N の窓の θ_max の最小が窓の端 (添字 {eN!r} / ±{nh})")
    if not checks["edge_2N"]:
        reasons.append(f"2N の窓の θ_max の最小が窓の端 (添字 {e2!r} / ±{nh})")
    out.update({"chosen_N": cN, "chosen_N_reason": rc["N"].get("reason") or rc["N"].get("error"),
                "chosen_2N": c2, "chosen_2N_reason": rc["2N"].get("reason") or rc["2N"].get("error"),
                "moved_N": (None if cN is None else cN - float(chosen)),
                "common_boundaries": common, "boundaries_N": rc["boundaries_N"], "boundaries_2N": rc["boundaries_2N"],
                "distance_to_conservative_boundary": dN, "distance_2N": d2,
                "distance_input": boundary_distance(chosen, common),
                "rule_N": {k: v for k, v in rc["N"].items() if k != "chosen_row"},
                "rule_2N": {k: v for k, v in rc["2N"].items() if k != "chosen_row"},
                "min_offset_N": eN, "min_offset_2N": e2,
                "row_chosen_N": at_N[_key(chosen)], "row_chosen_2N": at_2N[_key(chosen)],
                "rows_selected": {f"{m:.12g}": {"at_N": at_N[_key(m)], "at_2N": at_2N[_key(m)]} for m in sel_pts},
                "rows_N": evN.table(), "rows_2N": ev2.table(), "checks": checks, "pass": all(checks.values()),
                "warnings": rc["warnings"], "reasons": reasons, "elapsed_s": time.time() - t0})
    if raise_on_fail and not out["pass"]:
        raise MknotUndeterminedError("選択は未確定 (解像度のゲートが不合格): " + " / ".join(reasons), out)
    return out
