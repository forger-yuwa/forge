r"""軸 Mach 則 (knot) の knot のマッハ数 M_knot を、長さを固定したまま選ぶ (最大壁角の最小化)。

規則は 2 つ:
  - `coarse_min_theta_v2` (要求カードの既定、2026-10-15〜; plan plans/active/tooling-nozzle-mknot-fast.md §4.1・§4.2):
    `select_coarse_min_theta`。0.05 の格子 (整数 k、M_K = k/20) を粗い網 (N_s = N_p/4、長さの許容差 0.1 r_t) で 1 回掃き、
    候補の区間ごとの θ_max の局所極小を生産の網 (N_p、カードの L_total_tol) で前後と確かめ (ずらしは最大 3 回)、代表が残った谷を
    一度だけ作る同点の集合 (θ_min + 0.01°) から推定の距離が最大の谷を選ぶ。詳細は下の「規則 coarse_min_theta_v2」の節。
  - `min_theta_max_v1` (旧版、2026-10-14〜15 の既定; plan plans/accepted/tooling-nozzle-mknot-auto.md §4.1・§4.2、§5.1 #2・#2b・#2c・#3b):
    `select_mknot`・`resolution_gate`。カードからは呼ばない (関数と試験は残す)。以下の説明はこの旧版のもの。

設計計算だけで CFD は回さない。

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
  - 境目の探し直しの上限 (2026-10-14 result 段レビュー M1、§5.1 #3d): 境目 + ΔM_safe の点が新しい成立の切り替わりを作ったら、
    評価表 (全部の評価済みの点) から境目を求め直して選び直す。SETTLE_MAX (3) 回選び直しても境目の集合 (向きと成立側の端の列;
    二分法の区間の不成立側の端だけの詰め直しは数えない) が変わるなら、古い境目での選択を返さず `MknotUndeterminedError`
    (表つき; 警告だけで返さない)。
  - 選択値は評価した点の値のまま (格子の原点 M_A + 0.1 からの端数を丸めない; 表示だけ丸める)。
  - 成立点が無ければ例外 (表つき)。黙って 0.45 × M_d には戻さない。
  - 解像度のゲート (`resolution_gate`): n_axis_inv を 2 倍にして選んだ点の ±0.1 を 0.01 で詰め直し、境目も求め直す。
    N と 2N の境目を全部残し (対にしてまとめない; `conservative_boundaries`、2026-10-14 諮問 mknot-gate-output)、全部の境目の
    成立側の端からの距離で N と 2N の選択をそれぞれ選び直して比べる: 選択値の差 ≤ 0.01、両方の解像度で制約と境目からの距離が
    成立、最小が詰めの端に無いこと。満たさなければ「選択は未確定」。出力する値は選び直した N の値 `chosen_N` (§4.2)。
    選び直した値 (同点の規則の後) が詰めの窓の端に来たら、同じ geometry のまま窓を上限つきで広げて (半幅 0.1 → 0.2 → 0.3、
    評価済みの点は使い回す) 選び直す。上限でも端なら「未確定」(2026-10-14 諮問 mknot-card-holes)。θ_max の最小の位置の端の判定
    (同点の規則の前) とは別の判定。
    初回の窓・広げた窓とも、**N と 2N の両方**を ±(半幅 + ΔM_safe) まで評価し (窓の外は候補にしない)、それぞれの解像度で独立に
    境目を求める (2026-10-14 result 段レビュー M2: N だけでは 2N にだけある窓の外の境目を見落とす)。選び直しで境目の集合が
    上限まで収まらなければ未確定 (同 M1)。距離は最後に評価表 (全部の評価済みの点) から求め直した境目で判定する。

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
GATE_HALF_WIDTHS = (0.1, 0.2, 0.3)   # 解像度のゲートの窓の半幅の列 (§4.2 窓の端の選択: 選び直した値が端なら次へ広げる、上限 0.3)
SETTLE_MAX = 3            # 境目の探し直し (境目 + ΔM_safe の点が新しい切り替わりを作ったときの選び直し) の上限 (§5.1 #3d、レビュー M1)
_EPS = 1e-9
# evaluate_mknot が捕捉する逆 MOC の結果の数 (plan tooling-nozzle-mknot-fast §4.2): from_length が採る壁は最後か 1 つ前の評価
# (runner_axismach の from_length は改善が止まった時点で打ち切り、最良は常に最後の 2 回の中) なので最後の 2 つだけを持つ
# (網 levels を全部持つとメモリが評価の回数ぶん増える)。None なら全部 (旧の挙動; 試験の対照だけに使う)
_CAPTURE_KEEP: int | None = 2


class MknotSelectionError(ValueError):
    """選べない (成立点が無い・境目からの距離を満たす点が無い)。`table` に全評価点。"""

    def __init__(self, msg: str, table: list | None = None):
        self.table = list(table or [])
        super().__init__(msg + ("\n" + format_table(self.table) if self.table else ""))


class MknotUndeterminedError(MknotSelectionError):
    """選択は未確定: 解像度のゲートが不合格、または境目の探し直しが上限 (SETTLE_MAX) に達しても境目の集合が変わる
    (`select_mknot`・`reselect_common`; §5.1 #3d)。`result` にゲートの結果 (または探し直しの記録; rows_N・rows_2N が表)。"""

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


def _length_record(raw: dict, p, d: dict | None, caps: list, error: str | None) -> dict:
    """行の長さ合わせの記録 (plan tooling-nozzle-mknot-fast §4.1 の 1: 粗い網の各行に L_c・符号つきの長さの残差・停止の理由)。
    実効の入力 (r_t・L_total_tol・x_end_margin) も残す (どの入力で評価した行かを YAML・card_log と照合するため)。
    成功の行は design_chain の Lc_solve (from_length) から。失敗の行は捕捉した逆 MOC の結果 (最後の 2 つ) のうち |x_F − L_total| が
    最小のもの (from_length の最良は常に最後の 2 回の中) から求め、停止の理由は例外の文から分ける。値は記録だけ (判定には使わない)。"""
    g = (p.geometry if p is not None else raw.get("geometry") or {})
    mode = str(g.get("Lc_mode", "explicit"))
    rt = (raw.get("spec") or {}).get("r_throat")
    out = {"r_t_m": (None if rt is None else float(rt)), "Lc_mode": mode,
           "L_total_tol": (float(g.get("L_total_tol", 0.05)) if mode == "from_length" else None),
           "x_end_margin": float(g.get("x_end_margin", 2.3)), "L_residual_rt": None, "length_stop": None}
    if d is not None:
        sv = d.get("Lc_solve")
        if sv:
            out.update(L_residual_rt=float(sv["xF_residual"]), L_total_tol=float(sv["tol"]),
                       length_stop=f"許容差内 (|x_F − L_total| ≤ {float(sv['tol']):.4g} r_t、逆 MOC {int(sv['n_design_evals'])} 回)")
        else:
            out["length_stop"] = f"長さ合わせなし (Lc_mode {mode})"
        return out
    if mode != "from_length":
        out["length_stop"] = f"設計チェーンが例外 (Lc_mode {mode}、長さ合わせなし)"
        return out
    L_tgt = None
    v = (raw.get("dv") or {}).get("L_total")                 # probdef.dv_value と同じ読み方 ({value: x} か数値)
    try:
        L_tgt = float(v["value"] if isinstance(v, dict) else v)
    except (TypeError, ValueError, KeyError):
        L_tgt = None
    best = None
    for c in caps:
        res, x_E, anc = c[0], c[2], c[3]
        xF = float((res.get("exit") or {}).get("x_F", float("nan")))
        if L_tgt is None or not np.isfinite(xF):
            continue
        f = xF - L_tgt
        if best is None or abs(f) < abs(best[0]):
            best = (f, (None if (x_E is None or anc is None) else float(x_E) - float(anc[0])))
    if best is not None:
        out["L_residual_rt"] = best[0]
        if best[1] is not None:
            out["L_c"] = best[1]
    err = error or ""
    tol = out["L_total_tol"]
    if "from_length 未収束" in err:
        out["length_stop"] = "from_length 未収束 (改善停止か反復上限; 最良の |x_F − L_total| > L_total_tol)"
    elif "実現不能" in err:
        out["length_stop"] = "from_length の L_c が許容窓の端に張り付いた (実現不能)"
    elif "x_F 追跡に失敗" in err:
        out["length_stop"] = "from_length の途中で終端特性線の x_F 追跡に失敗 (場が短い)"
    elif best is not None and tol is not None and abs(best[0]) <= tol:
        out["length_stop"] = "許容差内 (長さ合わせの後の例外で不成立)"
    else:
        out["length_stop"] = "設計チェーンが例外 (長さ合わせの状態は不明)"
    return out


def evaluate_mknot(problem_dict: dict, M_K: float, n_axis_inv: int | None = None, *,
                   base_dir: str | Path | None = None, n_crossing_sample: int = 20,
                   L_total_tol: float | None = None) -> dict:
    """問題 (dict) の `geometry.M_knot` を M_K に (と `geometry.n_axis_inv`・`geometry.L_total_tol` を上書き) して design_chain を回し、
    制約の値と成立の判定を返す。例外は行に記録する (成立しない理由)。

    base_dir: 問題の相対パス (initial_line_run など) の基準。一時の YAML はシステムの一時ディレクトリ ($TMPDIR) に書き、
    読んだ後の `Problem.path` を base_dir の中の仮の名前にする (case の中にファイルを作らない)。
    n_crossing_sample: 特性線の交差の標本の数 (mknot_sweep.py の one() と同じ 20)。
    L_total_tol: from_length の長さの許容差 [r_t] の上書き (規則 coarse_min_theta_v2 の粗い網 L_s = 0.1 r_t; None なら問題の値)。
    逆 MOC の結果の捕捉は最後の _CAPTURE_KEEP (2) 個だけ (plan tooling-nozzle-mknot-fast §4.2; from_length の最良は最後の 2 回の中なので
    数値・判定は全部を捕捉したときと同じ。外を採ったら下の RuntimeError で行の error = 不成立)。
    行には長さ合わせの記録 (`_length_record`: 実効の r_t・L_total_tol・x_end_margin、L_c・符号つきの残差 L_residual_rt・停止の理由
    length_stop) も入れる (記録だけ)。"""
    import forge_design.evaluate.runner_axismach as RA
    from forge_design.geometry.moc_diagnostics import (characteristic_topology_diagnostics, wall_curvature_diagnostics,
                                                        wall_margin_diagnostics)
    from forge_design.probdef import load_problem

    t0 = time.time()
    raw = copy.deepcopy(problem_dict)
    raw.setdefault("geometry", {})["M_knot"] = float(M_K)
    if n_axis_inv:
        raw["geometry"]["n_axis_inv"] = int(n_axis_inv)
    if L_total_tol is not None:
        raw["geometry"]["L_total_tol"] = float(L_total_tol)
    raw.pop("_card", None)
    row = {"M_K": float(M_K), "n_axis": None}
    caps = []
    orig = RA.inverse_design
    p, d = None, None

    def wrapped(*a, **k):
        r = orig(*a, **k)
        caps.append((r, k.get("n_axis"), k.get("x_E"), k.get("axis_anchor")))
        if _CAPTURE_KEEP is not None:
            del caps[:-_CAPTURE_KEEP]           # 最後の _CAPTURE_KEEP 個だけを持つ (§4.2)
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
        hit = [(c[0], c[1]) for c in caps if c[0].get("wall") is w]
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
        d = None
    finally:
        RA.inverse_design = orig
    try:
        lrec = _length_record(raw, p, d, caps, row.get("error"))
    except Exception as e:  # noqa: BLE001 — 記録だけ。判定は変えない
        lrec = {"length_stop": f"長さ合わせの記録に失敗: {type(e).__name__}: {str(e)[:120]}"}
    for k, v in lrec.items():
        row.setdefault(k, v)
    if row.get("n_axis") is None and p is not None:
        row["n_axis_input"] = int(p.geometry.get("n_axis_inv", 500))
    row["reasons"] = feasibility_reasons(row)
    row["ok"] = not row["reasons"]
    row["time_s"] = time.time() - t0
    return row


def make_eval_fn(problem_dict: dict, n_axis_inv: int | None = None, base_dir: str | Path | None = None,
                 L_total_tol: float | None = None):
    """`select_mknot`・`select_coarse_min_theta` に渡す評価 (M_K → 行)。ProcessPoolExecutor に渡せるよう partial で作る。
    L_total_tol: from_length の長さの許容差の上書き (粗い網 L_s; None なら問題の値)。"""
    return functools.partial(evaluate_mknot, problem_dict, n_axis_inv=n_axis_inv, base_dir=base_dir, L_total_tol=L_total_tol)


class _RangeFound(Exception):
    """search_range: knot 則の L_c の許容窓の呼び出しに届いた (引数を持って design_chain を止める)。"""


def search_range(problem_dict: dict, base_dir: str | Path | None = None) -> dict:
    """§4.2 の格子の原点 M_A (軸のアンカーの M) と上端 M_d (軸則・MOC に渡す M_d; 出口較正 Md_moc_offset 込み) を、
    design_chain が knot 則の L_c の許容窓を求める呼び出し `KnotQuinticAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md, M_K)`
    の引数から読む (要求カードで M_knot を選ぶ前に使う; plan §4.3)。その呼び出しで design_chain を止めるので逆 MOC は回さず、
    `geometry.M_knot` の値は使わない (問題から外して渡す)。値は `evaluate_mknot` の行の M_A・M_d_moc と同じ
    (アンカーは初期線から決まり、M_knot に依らない)。axis_law が knot でなければ ValueError、その呼び出しの前の例外はそのまま。
    戻り: {M_A, M_d, x_A, Mp_A, Mpp_A}。"""
    import forge_design.evaluate.runner_axismach as RA
    from forge_design.probdef import load_problem

    raw = copy.deepcopy(problem_dict)
    raw.pop("_card", None)
    g = raw.setdefault("geometry", {})
    if str(g.get("axis_law", "quintic")) != "knot":
        raise ValueError(f"search_range: axis_law が knot でない ({g.get('axis_law', 'quintic')!r}) — M_knot を使わない")
    g.pop("M_knot", None)
    orig = RA.KnotQuinticAxisLaw

    class _Probe(orig):
        @classmethod
        def admissible_Lc_range(cls, x_A, M_A, Mp_A, Mpp_A, M_d, M_K, *a, **k):
            raise _RangeFound(float(x_A), float(M_A), float(Mp_A), float(Mpp_A), float(M_d))

    RA.KnotQuinticAxisLaw = _Probe
    try:
        with tempfile.TemporaryDirectory(prefix="forge_mknot_range_") as td:
            tmp = Path(td) / "problem.yaml"
            tmp.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
            p = load_problem(str(tmp))
            if base_dir is not None:
                p.path = str(Path(base_dir).resolve() / "_mknot_select_eval.yaml")
            RA.design_chain(p)
    except _RangeFound as f:
        x_A, M_A, Mp_A, Mpp_A, M_d = f.args
        return {"M_A": M_A, "M_d": M_d, "x_A": x_A, "Mp_A": Mp_A, "Mpp_A": Mpp_A}
    finally:
        RA.KnotQuinticAxisLaw = orig
    raise ValueError("search_range: design_chain が knot 則の L_c の許容窓を求める前に終わった (M_A・M_d を読めない)")


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
    """境目の集合が同じか。向きと成立側の端 (選び方が使う値: 距離と境目 + ΔM_safe の点) で比べる (§5.1 #3d)。区間の不成立側の
    端だけが詰まった変化 (二分法の幅の中に不成立の点が足された) は選択を変えないので同じとみなす — 探し直しの上限で「境目が
    変わった」と数えて未確定にしない。"""
    return ([(x["kind"], float(x["feasible_end"])) for x in a] == [(x["kind"], float(x["feasible_end"])) for x in b])


def _bounds_text(bs: list[dict]) -> str:
    """境目の列の表示 (例外のメッセージ・ログ用): [(向き, 成立側の端)]。"""
    return str([(b["kind"], round(b["feasible_end"], 6)) + ((b["source"],) if "source" in b else ()) for b in bs])


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
    境目 + ΔM_safe の点は探索の範囲 [M_A + coarse, 格子の上端] の中だけ評価する)。境目 + ΔM_safe の点が新しい成立の切り替わりを
    作ったら評価表から境目を求め直して選び直し、SETTLE_MAX 回でも境目の集合が変わるなら MknotUndeterminedError (表つき; §5.1 #3d)。
    戻り: {rule, chosen, chosen_row, table (全評価点), local_minima, windows, boundaries, reason, warnings, ...}。boundaries は評価表の
    全部の点から求めた境目 (選択に使った境目と同じ)。成立点が無ければ MknotSelectionError (表つき)。"""
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
        for it in range(SETTLE_MAX):
            rule = apply_rule(ev.table(), boundaries, tie_deg=tie_deg, dM_safe=dM_safe,
                              safe_eval=lambda Ms: ev.many(Ms, "safe"), safe_range=(grid(0), grid(J_max)))
            nb = find_boundaries(ev, bisect_tol)           # 評価表 (全部の評価済みの点) から境目を求め直す
            if _same_boundaries(nb, boundaries):
                boundaries = nb                             # 記録は評価表から求めた境目 (成立側の端は選択に使ったものと同じ)
                break
            if it == SETTLE_MAX - 1:
                # 上限でも境目の集合が変わる: 古い境目での選択は返さない (§5.1 #3d、result 段レビュー M1 — 警告だけで返すと、
                # 評価表に既にある境目から ΔM_safe 未満の点を合格として出す反例があった)
                raise MknotUndeterminedError(
                    f"選択は未確定: 境目 + ΔM_safe の点が {SETTLE_MAX} 回の選び直しの後も新しい成立の切り替わりを作る (境目の集合が"
                    f"収まらない) — 古い境目での選択 {rule['chosen']!r} (選択に使った境目から {rule['distance']:.4g}) は返さない。"
                    f"評価表から求め直した境目からは {boundary_distance(rule['chosen'], nb):.4g}; 選択に使った境目 {_bounds_text(boundaries)}、"
                    f"求め直した境目 {_bounds_text(nb)}",
                    {"stage": "select_mknot", "rows_N": ev.table(), "rows_2N": [], "boundaries_used": boundaries,
                     "boundaries_latest": nb, "chosen_stale": rule["chosen"],
                     "distance_stale_latest": boundary_distance(rule["chosen"], nb), "warnings": warnings})
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
def conservative_boundaries(bN: list[dict], b2: list[dict]) -> list[dict]:
    """N と 2N の境目を**全部残して** 1 つの列にする (§4.2; 2026-10-14 諮問 mknot-gate-output の Major)。候補の距離はこの列の
    全部の境目の成立側の端から測る (`boundary_distance`)。同じ向き・近い境目を対にしてまとめない — N と 2N の不成立の区間が
    違う側にあると、まとめた境目からの距離は全部の境目からの距離より大きく出る (反例: N の不成立 (2.00, 2.02)・2N の不成立
    (2.08, 2.10) で、2.05 はまとめると 0.05・全部からは 0.03)。片方の解像度にしか無い境目も残す。
    戻り: [{kind, feasible_end, lo, hi, width, source ('N' / '2N'), reasons_infeasible, nearest_other}] (成立側の端の順、
    同じ端なら N が先)。nearest_other は記録だけ (距離には使わない): もう一方の解像度の同じ向きの境目のうち成立側の端が
    最も近いものの {source, feasible_end, offset (その端 − この端)}。無ければ None。"""
    def entry(b, src, others, osrc):
        same = [c for c in others if c["kind"] == b["kind"]]
        near = min(same, key=lambda c: abs(c["feasible_end"] - b["feasible_end"])) if same else None
        return {"kind": b["kind"], "feasible_end": float(b["feasible_end"]), "lo": float(b["lo"]), "hi": float(b["hi"]),
                "width": float(b["width"]), "source": src, "reasons_infeasible": list(b.get("reasons_infeasible") or []),
                "nearest_other": (None if near is None else
                                  {"source": osrc, "feasible_end": float(near["feasible_end"]),
                                   "offset": float(near["feasible_end"] - b["feasible_end"])})}

    out = [entry(b, "N", b2, "2N") for b in bN] + [entry(c, "2N", bN, "N") for c in b2]
    return sorted(out, key=lambda x: (x["feasible_end"], x["source"] != "N"))


def reselect_common(evN: Evaluator, ev2: Evaluator, rows_N: list[dict], rows_2N: list[dict], *, tie_deg: float,
                    dM_safe: float, bisect_tol: float, safe_range=None, candidate_roles=CANDIDATE_ROLES,
                    max_iter: int = SETTLE_MAX, raise_unsettled: bool = True, log=None) -> dict:
    """N と 2N の**全部の境目** (`conservative_boundaries`; 対にしてまとめない) で両方の解像度の選択をそれぞれ選び直す
    (§4.2; 2N の境目が保守側へ動いたときに N の境目のままで判定しない・まとめた境目で距離を過大に見積もらない)。
    evN・ev2: 各解像度の Evaluator (境目を求め直し、境目 + ΔM_safe の点を評価する)。rows_N・rows_2N: 各解像度の候補の行
    (candidate_roles で絞る)。境目 + ΔM_safe の点は全部の境目について評価する (下の境目は +、上の境目は −)。それが新しい
    成立の切り替わりを作ったら評価表から境目を求め直して選び直す (max_iter 回まで)。
    max_iter 回でも境目の集合が変わるなら選択は未確定 (§5.1 #3d、result 段レビュー M1): raise_unsettled なら
    MknotUndeterminedError (表つき)、そうでなければ settled False を返す (呼び出し側が不合格にする; 解像度のゲートは
    raise_unsettled=False で呼び、判定 settled に入れる)。古い境目での選択を合格として扱わない。
    戻り: {boundaries_N, boundaries_2N, common (選択に使った全部の境目), N: apply_rule の戻り (選べなければ {chosen: None, error}),
    2N: 同, settled, boundaries_N_latest・boundaries_2N_latest・common_latest (評価表から最後に求め直した境目; settled なら
    選択に使った境目と同じ), warnings}。"""
    warnings = []
    bN, b2 = find_boundaries(evN, bisect_tol), find_boundaries(ev2, bisect_tol)
    nbN, nb2 = bN, b2
    res, common, settled = {}, [], True
    for it in range(max_iter):
        common = conservative_boundaries(bN, b2)
        for tag, ev, rows in (("N", evN, rows_N), ("2N", ev2, rows_2N)):
            try:
                res[tag] = apply_rule(rows, common, tie_deg=tie_deg, dM_safe=dM_safe,
                                      safe_eval=lambda Ms, ev=ev: ev.many(Ms, "safe"),
                                      candidate_roles=candidate_roles, safe_range=safe_range)
            except MknotSelectionError as e:
                res[tag] = {"chosen": None, "error": str(e).splitlines()[0]}
        nbN, nb2 = find_boundaries(evN, bisect_tol), find_boundaries(ev2, bisect_tol)
        if _same_boundaries(nbN, bN) and _same_boundaries(nb2, b2):
            bN, b2 = nbN, nb2                           # 記録は評価表から求めた境目 (成立側の端は選択に使ったものと同じ)
            common = conservative_boundaries(bN, b2)
            break
        if it == max_iter - 1:
            # 上限でも境目の集合が変わる: 選択は未確定 (古い境目での選択を合格として返さない; §5.1 #3d)
            settled = False
            latest = conservative_boundaries(nbN, nb2)
            msg = (f"選択は未確定: 境目 + ΔM_safe の点が {max_iter} 回の選び直しの後も新しい成立の切り替わりを作る (全部の境目の集合が"
                   f"収まらない) — 古い境目での選択 N {res['N'].get('chosen')!r}・2N {res['2N'].get('chosen')!r} は使わない "
                   f"(評価表から求め直した境目からの距離 "
                   + "・".join(f"{tag} {('—' if res[tag].get('chosen') is None else format(boundary_distance(res[tag]['chosen'], latest), '.4g'))}"
                               for tag in ("N", "2N"))
                   + f"); 求め直した境目 {_bounds_text(latest)}")
            warnings.append(msg)
            if raise_unsettled:
                raise MknotUndeterminedError(msg, {"stage": "reselect_common", "rows_N": evN.table(), "rows_2N": ev2.table(),
                                                   "common_used": common, "common_latest": latest,
                                                   "N": {k: v for k, v in res["N"].items() if k != "chosen_row"},
                                                   "2N": {k: v for k, v in res["2N"].items() if k != "chosen_row"},
                                                   "warnings": warnings})
            break
        warnings.append("境目 + ΔM_safe の点が新しい成立の切り替わりを作った — 全部の境目を求め直して選び直した")
        bN, b2 = nbN, nb2
    if log:
        log("  全部の境目: " + ", ".join(f"{b['kind']} {b['feasible_end']:.6f} ({b['source']}"
                                         + ("" if b["nearest_other"] is None else
                                            f", {b['nearest_other']['source']} の同じ向き {b['nearest_other']['offset']:+.6f}")
                                         + ")" for b in common)
            + f" → 選び直し N {res['N'].get('chosen')!r}・2N {res['2N'].get('chosen')!r}"
            + ("" if settled else f" — 未確定 (境目の探し直しが上限 {max_iter} 回)"))
    return {"boundaries_N": bN, "boundaries_2N": b2, "common": common, "N": res["N"], "2N": res["2N"], "settled": settled,
            "boundaries_N_latest": nbN, "boundaries_2N_latest": nb2, "common_latest": conservative_boundaries(nbN, nb2),
            "warnings": warnings}


def resolution_gate(problem_dict: dict | None, chosen: float, N: int, *, select_result: dict | None = None,
                    eval_factory=None, base_dir=None, fine: float = 0.01, half_widths=GATE_HALF_WIDTHS,
                    tie_deg: float | None = None, dM_safe: float | None = None, bisect_tol: float | None = None,
                    max_change: float = 0.01, workers: tuple[int, int] = (4, 2), seed_2N: list | None = None,
                    raise_on_fail: bool = True, log=None) -> dict:
    """§4.2 の解像度のゲート。`n_axis_inv` = N と 2N で選んだ点 (chosen) の前後 ±半幅を fine で詰め (chosen を中心の格子)、
    境目を 2N でも求め直す (N の境目のうち chosen から 半幅 + ΔM_safe 以内のものは区間の両端も 2N で評価する)。
    N と 2N の境目を全部残し (`conservative_boundaries`; 対にしてまとめない、2026-10-14 諮問 mknot-gate-output)、全部の境目の
    成立側の端からの距離で窓の候補から N と 2N の選択をそれぞれ選び直して比べる (`reselect_common`; 境目 + ΔM_safe の点は
    全部の境目について窓の中だけ評価)。

    窓の外の文脈 (§5.1 #3d、2026-10-14 result 段レビュー M2): 初回の窓・広げた窓とも、**N と 2N の両方**を ±(半幅 + ΔM_safe) まで
    fine の格子で評価し (役割 context; 候補にしない)、それぞれの解像度で独立に境目を求める (N だけを広げると、2N にだけある窓の
    すぐ外の境目を見落として、そこから ΔM_safe 未満の点を合格にする反例があった)。
    境目の探し直し (同 M1): 選び直しで境目の集合が SETTLE_MAX 回でも収まらなければ判定 settled が不合格 (古い境目での選択を
    合格にしない)。距離の判定は、最後に評価表 (N・2N の全部の評価済みの点; 交差の評価を含む) から求め直した境目で行う
    (選択に使った境目と食い違えば settled も不合格)。

    窓の端の選択 (plan §4.2、2026-10-14 諮問 mknot-card-holes の Major): 選び直した N・2N の値 (同点の規則の後) のどちらかが
    窓の端に来たら、同じ geometry・同じ chosen を中心に半幅を half_widths の次の値に広げ (既定 0.1 → 0.2 → 0.3; 評価済みの点は
    Evaluator が使い回す)、境目・距離・同点の規則を当て直して選び直す。平らな谷では同点の規則が境目から最も遠い点を選ぶので、
    端の選択は窓が選択を縛っている印であって、それだけで「未確定」にはしない (諮問の合成の例: 境目 1.65 / 2.25 で窓 1.95〜2.15
    の選択 1.95 は端、半幅 0.2 で同じ 1.95 が内部に入り距離 0.30)。上限の半幅でも端なら chosen_inside_* が不合格 (未確定)。
    広げた窓でも N と 2N を ±(半幅 + ΔM_safe) まで評価して、窓の外 ΔM_safe 以内の境目を見る (上の「窓の外の文脈」)。
    境目の探し直しが上限で収まらなければ窓は広げない (その窓で未確定)。

    判定 (全部で合格; 最後の窓で):
      settled     選び直しの境目の探し直しが SETTLE_MAX 回以内に収まり、評価表から最後に求め直した境目が選択に使った境目と同じ
      change      全部の境目で選び直した N と 2N の選択値の差 ≤ max_change (0.01)
      ok_N・ok_2N 入力の選択点 (chosen) と選び直した 2 点が、その解像度で全部の制約を満たす
                  (選び直した 2 点は両方の解像度で評価する; 高い解像度だけ制約を破る選択点を通さない)
      distance    選び直した 2 点から全部の境目 (評価表から最後に求め直したもの) まで ≥ ΔM_safe (入力の点の距離は distance_input に
                  記録だけ)
      edge_N・edge_2N  窓の中の成立点の θ_max の最小 (同点の規則の前) が窓の端に無い (窓を広げる契機にはしない; 別の判定)
      chosen_inside_N・chosen_inside_2N  選び直した値 (同点の規則の後) が窓の端に無い (上限の半幅まで広げた後)
    出力する値は選び直した N の値 `chosen_N` (§4.2; 入力の chosen をそのまま出すと境目から ΔM_safe 未満の値が出うる)。
    chosen_N は入力の chosen と違うことがある (2N の境目が保守側へ動いた、同点の規則の境目からの距離が変わった、窓を広げた)。
    その差は moved_N に記録する (合否には入れない)。row_chosen_N・row_chosen_2N は chosen_N を N で・chosen_2N を 2N で評価した行、
    row_input_N・row_input_2N は入力の chosen の行 (rows_selected に全部の組)。
    窓の記録: window・half_width は最後の窓、n_widen は広げた回数、widen は半幅ごとの選び直し (窓・chosen_N・chosen_2N・
    端だった値 edge_values・θ_max の最小の位置 min_offset_*・次の手 action)、edge_values は最後の窓で端だった値 (内部なら空)。
    min_offset_* と edge_* は最後の窓の値。
    select_result: `select_mknot` の戻り (N の評価の表・境目・設定を再利用; 同じ問題・同じ N の評価に限る)。seed_2N: 同じ問題の
    2N の評価済みの行 (同じ問題・同じ 2N に限る)。eval_factory: n_axis → 評価の関数 (既定は make_eval_fn(problem_dict, n_axis,
    base_dir); 試験で差し替える)。half_widths: 窓の半幅の列 (正で増える、fine の整数倍)。raise_on_fail なら不合格で
    MknotUndeterminedError。
    窓の評価は並列 (workers)、選び直しの境目 + ΔM_safe の点と交差の評価は逐次 (少数の点; 2 つの並列を同時に持たない)。"""
    st = (select_result or {}).get("settings") or {}
    tie_deg = st.get("tie_deg", 1e-3) if tie_deg is None else tie_deg
    dM_safe = st.get("dM_safe", 0.05) if dM_safe is None else dM_safe
    bisect_tol = st.get("bisect_tol", 1e-3) if bisect_tol is None else bisect_tol
    hws = [float(h) for h in half_widths]
    nhs = [int(round(h / fine)) for h in hws]
    if (not hws or any(h <= 0.0 for h in hws) or any(b <= a for a, b in zip(hws[:-1], hws[1:]))
            or any(abs(n * fine - h) > 1e-9 for n, h in zip(nhs, hws))):
        raise ValueError(f"half_widths は正で増える列、各値は fine {fine} の整数倍 (受け取った値: {half_widths!r})")
    if eval_factory is None:
        def eval_factory(n):
            return make_eval_fn(problem_dict, n_axis_inv=n, base_dir=base_dir)
    t0 = time.time()
    out = {"N": int(N), "N2": int(2 * N), "chosen": float(chosen),
           "settings": {"fine": fine, "half_widths": hws, "tie_deg": tie_deg, "dM_safe": dM_safe,
                        "bisect_tol": bisect_tol, "max_change": max_change}}
    n_ctx = int(math.ceil(dM_safe / fine - 1e-9)) if dM_safe > 0.0 else 0      # 窓の外の文脈 (N・2N、ΔM_safe ぶん)

    def grid(j: int) -> float:
        return round(float(chosen) + j * fine, 12)

    def _edge(rows, nh):
        f = [(j, r) for j, r in enumerate(rows) if _feasible(r)]
        if not f:
            return None
        j, _ = min(f, key=lambda t: (t[1]["theta_max_deg"], t[1]["M_K"]))
        return j - nh

    def _at_end(c, window):
        return c is not None and (abs(c - window[0]) <= _EPS or abs(c - window[-1]) <= _EPS)

    evN = Evaluator(eval_factory(N), workers[0], log)
    if select_result:
        evN.seed(select_result.get("table") or [])
    ev2 = Evaluator(eval_factory(2 * N), workers[1], log)
    if seed_2N:
        ev2.seed(seed_2N)
    steps = []
    for i, (hw, nh) in enumerate(zip(hws, nhs)):
        window = [grid(j) for j in range(-nh, nh + 1)]
        # 窓の外の文脈 ±(半幅 + ΔM_safe): 初回の窓・広げた窓とも N と 2N の両方で評価する (候補にしない; §5.1 #3d、レビュー M2)
        ctx = [grid(j) for j in range(-nh - n_ctx, nh + n_ctx + 1) if abs(j) > nh]
        with evN:
            if log:
                log(f"解像度のゲート: N = {N} の窓 {window[0]:.6g}〜{window[-1]:.6g}" + ("" if i == 0 else f" (半幅 {hw:g} に広げた)")
                    + (f" と窓の外 ±{hw + dM_safe:g} まで" if ctx else ""))
            rows_N = evN.many(window, "window")
            if ctx:
                evN.many(ctx, "context")
            bN = find_boundaries(evN, bisect_tol)
        near = [b for b in bN if abs(b["feasible_end"] - chosen) <= hw + dM_safe + _EPS]
        with ev2:
            if log:
                log(f"解像度のゲート: 2N = {2 * N} の窓" + (f" と窓の外 ±{hw + dM_safe:g} まで" if ctx else "")
                    + f"、N の境目の区間の両端 ({len(near)} 個)")
            rows_2N = ev2.many(window, "window")
            if ctx:
                ev2.many(ctx, "context")
            if near:
                ev2.many([v for b in near for v in (b["lo"], b["hi"])], "bracket")
            find_boundaries(ev2, bisect_tol)        # 2N の二分法は並列の中で (reselect_common は評価済みの点を使う)
        # 全部の境目で両方を選び直す (候補は窓の点だけ; N の窓の点は select_mknot の表の行を再利用するので役割を窓に揃える)。
        # 境目の探し直しが上限で収まらなければ settled False (下の判定で不合格; 古い境目での選択を合格にしない)
        rc = reselect_common(evN, ev2, [dict(r, role="window") for r in rows_N], [dict(r, role="window") for r in rows_2N],
                             tie_deg=tie_deg, dM_safe=dM_safe, bisect_tol=bisect_tol,
                             safe_range=(window[0], window[-1]), candidate_roles=("window", "safe"), raise_unsettled=False,
                             log=log)
        cN, c2 = rc["N"].get("chosen"), rc["2N"].get("chosen")
        edge_vals = {tag: c for tag, c in (("N", cN), ("2N", c2)) if _at_end(c, window)}
        last = i == len(hws) - 1
        steps.append({"half_width": hw, "window": [window[0], window[-1]], "chosen_N": cN, "chosen_2N": c2,
                      "edge_values": edge_vals, "min_offset_N": _edge(rows_N, nh), "min_offset_2N": _edge(rows_2N, nh),
                      "settled": rc["settled"],
                      "action": (f"境目の探し直しが上限 {SETTLE_MAX} 回で収まらない — 未確定 (窓は広げない)" if not rc["settled"] else
                                 "選び直した値は窓の内部" if not edge_vals else
                                 f"選び直した値が上限の半幅 {hw:g} でも窓の端 — 未確定" if last else
                                 f"選び直した値が窓の端 → 半幅 {hws[i + 1]:g} に広げて選び直す")})
        if log and (edge_vals or not rc["settled"]):
            log(f"  半幅 {hw:g}: " + "・".join(f"{k} {v!r}" for k, v in edge_vals.items())
                + (" が窓の端" if edge_vals else "") + f" — {steps[-1]['action']}")
        if not edge_vals or last or not rc["settled"]:
            break
    sel_pts = [float(chosen)] + [c for c in (cN, c2) if c is not None]
    at_N = {_key(r["M_K"]): r for r in evN.many(sel_pts, "cross")}          # 逐次 (with の外)
    at_2N = {_key(r["M_K"]): r for r in ev2.many(sel_pts, "cross")}
    # 距離は評価表 (N・2N の全部の評価済みの点; 交差の評価を含む) から求め直した境目で判定する (§5.1 #3d、レビュー M1:
    # 選び直しに使った境目の集合だけを信じない)。食い違えば選択に使った境目が古い = settled も不合格
    bN_fin, b2_fin = find_boundaries(evN, bisect_tol), find_boundaries(ev2, bisect_tol)
    common = conservative_boundaries(bN_fin, b2_fin)
    same_fin = _same_boundaries(bN_fin, rc["boundaries_N"]) and _same_boundaries(b2_fin, rc["boundaries_2N"])

    def _bad(at, pts):
        return [m for m in pts if not _feasible(at[_key(m)])]

    eN, e2 = steps[-1]["min_offset_N"], steps[-1]["min_offset_2N"]
    dN = None if cN is None else boundary_distance(cN, common)
    d2 = None if c2 is None else boundary_distance(c2, common)
    badN, bad2 = _bad(at_N, sel_pts), _bad(at_2N, sel_pts)
    checks = {
        "settled": bool(rc["settled"]) and same_fin,
        "change": cN is not None and c2 is not None and abs(c2 - cN) <= max_change + _EPS,
        "ok_N": cN is not None and c2 is not None and not badN,
        "ok_2N": cN is not None and c2 is not None and not bad2,
        "distance": dN is not None and d2 is not None and min(dN, d2) >= dM_safe - _EPS,
        "edge_N": eN is not None and abs(eN) < nh, "edge_2N": e2 is not None and abs(e2) < nh,
        "chosen_inside_N": cN is not None and "N" not in edge_vals,
        "chosen_inside_2N": c2 is not None and "2N" not in edge_vals,
    }
    reasons = []
    if not checks["settled"]:
        reasons.append(f"境目の集合が収まらない — 古い境目での選択 N {cN!r}・2N {c2!r} は使わない ("
                       + (f"選び直しの境目の探し直しが上限 {SETTLE_MAX} 回に達しても境目が変わる" if not rc["settled"] else
                          "選び直しの後の評価 (交差) で評価表の境目が変わった")
                       + f"; 選択に使った境目 {_bounds_text(rc['common'])}、評価表から求め直した境目 {_bounds_text(common)})")
    if not checks["change"]:
        reasons.append(f"全部の境目で選び直した 2N の値 {c2!r} が N の {cN!r} から {max_change:g} を超えて動く "
                       f"(N: {rc['N'].get('reason') or rc['N'].get('error')}; 2N: {rc['2N'].get('reason') or rc['2N'].get('error')})")
    for tag, n, bad, at in (("N", N, badN, at_N), ("2N", 2 * N, bad2, at_2N)):
        if not checks[f"ok_{tag}"]:
            reasons.append(f"{tag} = {n} で選択点が制約を満たさない: "
                           + ("選び直せない" if cN is None or c2 is None else
                              " / ".join(f"{m:.6g}: {'; '.join(at[_key(m)].get('reasons') or [])}" for m in bad)))
    if not checks["distance"]:
        reasons.append(f"選び直した点から全部の境目まで N {dN!r}・2N {d2!r} < ΔM_safe {dM_safe:g}")
    if not checks["edge_N"]:
        reasons.append(f"N の窓の θ_max の最小 (同点の規則の前) が窓の端 (添字 {eN!r} / ±{nh})")
    if not checks["edge_2N"]:
        reasons.append(f"2N の窓の θ_max の最小 (同点の規則の前) が窓の端 (添字 {e2!r} / ±{nh})")
    for tag in ("N", "2N"):
        if tag in edge_vals:
            reasons.append(f"選び直した {tag} の値 {edge_vals[tag]!r} (同点の規則の後) が窓 {window[0]:.6g}〜{window[-1]:.6g} の端 — "
                           f"半幅を上限 {hws[-1]:g} まで広げても端 (半幅の列 {' → '.join(f'{h:g}' for h in hws)})")
    out.update({"window": [window[0], window[-1]], "half_width": hw, "n_widen": len(steps) - 1, "widen": steps,
                "edge_values": edge_vals,
                "chosen_N": cN, "chosen_N_reason": rc["N"].get("reason") or rc["N"].get("error"),
                "chosen_2N": c2, "chosen_2N_reason": rc["2N"].get("reason") or rc["2N"].get("error"),
                "moved_N": (None if cN is None else cN - float(chosen)),
                # 境目は評価表から最後に求め直したもの (距離の判定に使う)。選び直しに使った境目は common_boundaries_used に記録
                "common_boundaries": common, "boundaries_N": bN_fin, "boundaries_2N": b2_fin,
                "common_boundaries_used": rc["common"], "settled": checks["settled"],
                "distance_to_conservative_boundary": dN, "distance_2N": d2,
                "distance_input": boundary_distance(chosen, common),
                "rule_N": {k: v for k, v in rc["N"].items() if k != "chosen_row"},
                "rule_2N": {k: v for k, v in rc["2N"].items() if k != "chosen_row"},
                "min_offset_N": eN, "min_offset_2N": e2,
                # 名前どおりの行 (2026-10-14 諮問 mknot-gate-output: 以前は入力の chosen の行を持っていた)
                "row_chosen_N": (None if cN is None else at_N[_key(cN)]),
                "row_chosen_2N": (None if c2 is None else at_2N[_key(c2)]),
                "row_input_N": at_N[_key(chosen)], "row_input_2N": at_2N[_key(chosen)],
                "rows_selected": {f"{m:.12g}": {"at_N": at_N[_key(m)], "at_2N": at_2N[_key(m)]} for m in sel_pts},
                "rows_N": evN.table(), "rows_2N": ev2.table(), "checks": checks, "pass": all(checks.values()),
                "warnings": rc["warnings"], "reasons": reasons, "elapsed_s": time.time() - t0})
    if raise_on_fail and not out["pass"]:
        raise MknotUndeterminedError("選択は未確定 (解像度のゲートが不合格): " + " / ".join(reasons), out)
    return out


# --- 規則 coarse_min_theta_v2 (plan tooling-nozzle-mknot-fast §4.1、要求カードの既定) ----------------------------------------
#
#   1. 粗い掃引: M_A < 0.05k < M_d かつ両端から 0.05 以上内側の全部の整数 k を、粗い網 (N_s = N_p/4、長さの許容差 L_s = 0.1 r_t) で
#      1 回評価する。制約は feasibility_reasons (θ_max が有限でも ok=False の行はある)。
#   2. 候補の区間: 粗い網で成立し、両隣 (k ± 1) も成立している k を候補にし、k が続く候補を 1 つの区間にまとめる (不成立の点で分かれる)。
#      境目までの推定の距離 = 候補から最も近い粗い網の不成立の格子点までの距離 − 0.05 (記録用と同点の決着用)。
#   3. 谷の代表: 区間ごとに粗い網の θ_max の局所極小を全部拾う (区間の端も内側の隣と比べて極小なら含める。等しい値が続く平らな列は
#      1 つの極小とし、代表は推定の距離が最大の点、同じなら小さい k)。
#   4. 谷ごとの確認とずらし (生産の網 N_p、カードの L_total_tol): 代表 v と v ± 1 を評価する。v が不成立 → 谷を捨てる。v と両隣が
#      成立 → 谷の最終の代表は v。v が成立で片側だけ不成立 → 反対側へ 1 つずらす (移動先が同じ候補の区間にあるときだけ)。移動先で
#      同じ確認をその谷で最大 3 回まで繰り返す (初回の確認 + 移動 3 回)。区間の外・両側が不成立・上限 → 谷を捨てる。
#   5. 谷の比較: 最終の代表が決まった谷の生産の網の θ_max の最小を θ_min とし、θ_max ≤ θ_min + 0.01° の谷を同点の集合として一度だけ
#      作る (逐次比較しない)。その中で推定の距離が最大、同じなら θ_max が小さい、同じなら M_K が小さい谷を選ぶ。最終の代表が 1 つも
#      無ければ例外 (表つき)。
#   6. 粗い網で成立点が無い・粗い網の長さ合わせが全点で未収束・候補が無いときは例外 (表つき)。黙って旧の規則に落ちない。
#   出す値は格子の値 (k / 20)。
RULE_NAME_V2 = "coarse_min_theta_v2"
V2_K_PER_UNIT = 20          # 格子 M_K = k / 20 (= 0.05 k)。整数 k で持ち、浮動小数の比較で隣を決めない
V2_EDGE_STEPS = 1           # 探索の範囲: M_A・M_d から 0.05 (= 1 格子) 以上内側
V2_DIST_OFFSET_STEPS = 1    # 推定の距離 = 最も近い粗い網の不成立の格子点までの距離 − 0.05 (= 1 格子; 旧の ΔM_safe と同じ値)
V2_TIE_DEG = 0.01           # 谷の比較の同点の幅 [°] (θ_min からの差; 設計の許容差)
V2_SHIFT_MAX = 3            # 谷ごとのずらしの上限 (初回の確認 + 移動 3 回)
V2_COARSE_DIV = 4           # 粗い網 N_s = 生産の n_axis_inv / 4 (2400 なら 600)
V2_L_COARSE = 0.1           # 粗い網の長さの許容差 L_s [r_t] (事前に固定; カードの値は 600 の網では全点未収束、plan §3)
_LENGTH_UNCONVERGED = "from_length 未収束"   # runner_axismach.design_chain の from_length の未収束の例外の文 (行の error に入る)


def v2_M(k: int) -> float:
    """格子の値 M_K = k / 20 (0.05 k の最も近い倍精度; 例 41 → 2.05)。"""
    return int(k) / V2_K_PER_UNIT


def v2_k(M: float) -> int:
    """M_K → 格子の整数 k。格子の上に無ければ ValueError (隣を浮動小数の比較で決めない)。"""
    k = int(round(float(M) * V2_K_PER_UNIT))
    if abs(v2_M(k) - float(M)) > 1e-9:
        raise ValueError(f"M_K {M!r} は 0.05 の格子の上に無い (規則 {RULE_NAME_V2})")
    return k


def v2_grid(M_A: float, M_d: float) -> list[int]:
    """§4.1 の 1 の探索の範囲の整数 k: M_A < k/20 < M_d かつ両端から 0.05 以上内側 (k/20 − M_A ≥ 0.05、M_d − k/20 ≥ 0.05)。"""
    edge = V2_EDGE_STEPS / V2_K_PER_UNIT
    lo = math.ceil((float(M_A) + edge) * V2_K_PER_UNIT - 1e-6)
    hi = math.floor((float(M_d) - edge) * V2_K_PER_UNIT + 1e-6)
    ks = [k for k in range(lo, hi + 1) if v2_M(k) - float(M_A) >= edge - _EPS and float(M_d) - v2_M(k) >= edge - _EPS]
    if not ks:
        raise ValueError(f"探索の範囲が空 (M_A {M_A!r} + 0.05 > M_d {M_d!r} − 0.05; 規則 {RULE_NAME_V2})")
    return ks


def _is_length_unconverged(row: dict) -> bool:
    return _LENGTH_UNCONVERGED in str(row.get("error") or "")


def v2_distance_steps(k: int, infeasible_k: list[int]) -> float:
    """推定の距離を格子の数で (最も近い粗い網の不成立の格子点までの |Δk| − 1)。不成立の格子点が無ければ inf。"""
    if not infeasible_k:
        return math.inf
    return float(min(abs(int(k) - j) for j in infeasible_k) - V2_DIST_OFFSET_STEPS)


def v2_distance(k: int, infeasible_k: list[int]) -> float | None:
    """推定の距離 [M_K] (記録用; 不成立の格子点が無ければ None = 無限大)。"""
    s = v2_distance_steps(k, infeasible_k)
    return None if math.isinf(s) else s / V2_K_PER_UNIT


def v2_coarse_analysis(rows: list[dict]) -> dict:
    """§4.1 の 2: 粗い掃引の行 (順不同) → 成立・候補・候補の区間・不成立の格子点。同じ k の行が 2 つあれば ValueError。
    戻り: {by_k, feasible_k, infeasible_k, candidates_k, intervals_k (= [(k_a, k_b)])}。"""
    by_k: dict[int, dict] = {}
    for r in rows:
        k = v2_k(r["M_K"])
        if k in by_k:
            raise ValueError(f"粗い掃引に同じ M_K {v2_M(k)!r} の行が 2 つある")
        by_k[k] = r
    ks = sorted(by_k)
    feas = {k: _feasible(by_k[k]) for k in ks}
    cand = [k for k in ks if feas[k] and feas.get(k - 1) is True and feas.get(k + 1) is True]
    intervals: list[tuple[int, int]] = []
    for k in cand:
        if intervals and k == intervals[-1][1] + 1:
            intervals[-1] = (intervals[-1][0], k)
        else:
            intervals.append((k, k))
    return {"by_k": by_k, "feasible_k": [k for k in ks if feas[k]], "infeasible_k": [k for k in ks if not feas[k]],
            "candidates_k": cand, "intervals_k": intervals}


def v2_valleys(ana: dict) -> list[dict]:
    """§4.1 の 3: 候補の区間ごとに粗い網の θ_max の局所極小を全部拾う。等しい値 (同じ倍精度) が続く平らな列は 1 つの極小で、
    その両外側 (区間の中だけ) の値がどちらも大きければ極小 (区間の端の側は比べない)。代表は列の中で推定の距離が最大の点、同じなら
    小さい k。戻り: [{interval_k, plateau_k, rep_k, theta_coarse}] (区間の順・区間の中の k の順)。"""
    th = {k: float(ana["by_k"][k]["theta_max_deg"]) for k in ana["candidates_k"]}    # 候補は成立 = θ_max が有限
    out = []
    for a, b in ana["intervals_k"]:
        ks = list(range(a, b + 1))
        runs, i = [], 0
        while i < len(ks):
            j = i
            while j + 1 < len(ks) and th[ks[j + 1]] == th[ks[i]]:
                j += 1
            runs.append((i, j))
            i = j + 1
        for i0, i1 in runs:
            t = th[ks[i0]]
            left = th[ks[i0 - 1]] if i0 > 0 else None
            right = th[ks[i1 + 1]] if i1 < len(ks) - 1 else None
            if (left is None or left > t) and (right is None or right > t):
                plateau = ks[i0:i1 + 1]
                rep = min(plateau, key=lambda k: (-v2_distance_steps(k, ana["infeasible_k"]), k))
                out.append({"interval_k": (a, b), "plateau_k": plateau, "rep_k": rep, "theta_coarse": t})
    return out


def v2_confirm(valleys: list[dict], eval_many, *, shift_max: int = V2_SHIFT_MAX, log=None) -> list[dict]:
    """§4.1 の 4: 谷ごとの確認とずらし。eval_many: M_K の列 → 生産の網の行の列 (同じ M_K は 1 回だけ評価する窓口; 評価の順は
    結果を変えない)。各回はまだ確認中の全部の谷の (中心 − 1, 中心, 中心 + 1) をまとめて評価する (並列を活かす; 谷どうしは独立)。
    戻り: 谷の記録の列 (入力の順) — {interval, plateau, start, steps [{center, left_ok, center_ok, right_ok, action}],
    final (格子の値 or None), final_k, discard_reason, n_moves}。"""
    st = [{"v": v, "c": v["rep_k"], "moves": 0, "steps": [], "status": "active", "why": None} for v in valleys]
    while True:
        act = [s for s in st if s["status"] == "active"]
        if not act:
            break
        need = sorted({k for s in act for k in (s["c"] - 1, s["c"], s["c"] + 1)})
        got = {v2_k(r["M_K"]): r for r in eval_many([v2_M(k) for k in need])}
        for s in act:
            c = s["c"]
            okl, okc, okr = (_feasible(got[c - 1]), _feasible(got[c]), _feasible(got[c + 1]))
            step = {"center": v2_M(c), "left": v2_M(c - 1), "right": v2_M(c + 1), "left_ok": okl, "center_ok": okc, "right_ok": okr}
            a, b = s["v"]["interval_k"]
            if not okc:
                s["status"], s["why"] = "discarded", f"中心 {v2_M(c):.2f} が生産の網で不成立 → 谷を捨てる"
            elif okl and okr:
                s["status"] = "final"
                step["action"] = f"中心 {v2_M(c):.2f} と両隣が成立 → 谷の最終の代表"
            elif not okl and not okr:
                s["status"], s["why"] = "discarded", f"中心 {v2_M(c):.2f} の両隣が生産の網で不成立 → 谷を捨てる"
            else:
                dest = c + 1 if not okl else c - 1
                side = "左" if not okl else "右"
                if s["moves"] >= shift_max:
                    s["status"], s["why"] = "discarded", (f"{side}隣 {v2_M(c - 1 if not okl else c + 1):.2f} が不成立だが、ずらしが上限 "
                                                          f"{shift_max} 回 → 谷を捨てる")
                elif not (a <= dest <= b):
                    s["status"], s["why"] = "discarded", (f"{side}隣が不成立、移動先 {v2_M(dest):.2f} が候補の区間 "
                                                          f"[{v2_M(a):.2f}, {v2_M(b):.2f}] の外 → 谷を捨てる")
                else:
                    step["action"] = f"{side}隣が不成立 → 反対側の {v2_M(dest):.2f} へずらす ({s['moves'] + 1} 回目)"
                    s["c"], s["moves"] = dest, s["moves"] + 1
            if s["status"] == "discarded":
                step["action"] = s["why"]
            s["steps"].append(step)
        if log:
            log("  谷の確認: " + "; ".join(f"{v2_M(s['v']['rep_k']):.2f} → {s['status']}"
                                          + (f" ({v2_M(s['c']):.2f})" if s["status"] != "discarded" else "") for s in st))
    out = []
    for s in st:
        v = s["v"]
        fin = s["status"] == "final"
        out.append({"interval": [v2_M(v["interval_k"][0]), v2_M(v["interval_k"][1])], "plateau": [v2_M(k) for k in v["plateau_k"]],
                    "start": v2_M(v["rep_k"]), "theta_coarse": v["theta_coarse"], "steps": s["steps"], "n_moves": s["moves"],
                    "final": (v2_M(s["c"]) if fin else None), "final_k": (s["c"] if fin else None),
                    "discard_reason": (None if fin else s["why"])})
    return out


def v2_pick(finals: list[dict], *, tie_deg: float = V2_TIE_DEG) -> dict:
    """§4.1 の 5: 最終の代表が決まった谷 (各 {M_K, theta_max_deg, distance_steps}) から、θ_min + tie_deg 以内の谷を同点の集合として
    一度だけ作り、推定の距離が最大 → θ_max が小さい → M_K が小さい順で 1 つ選ぶ。入力の順に依らない。戻り: {chosen (入力の dict),
    theta_min, tie_set (M_K の列、昇順)}。finals が空なら ValueError。"""
    if not finals:
        raise ValueError("比べる谷が無い")
    tmin = min(float(f["theta_max_deg"]) for f in finals)
    ties = [f for f in finals if float(f["theta_max_deg"]) <= tmin + tie_deg]
    best = min(ties, key=lambda f: (-float(f["distance_steps"]), float(f["theta_max_deg"]), float(f["M_K"])))
    return {"chosen": best, "theta_min": tmin, "tie_set": sorted(float(f["M_K"]) for f in ties)}


def select_coarse_min_theta(eval_coarse, eval_prod, M_A: float, M_d: float, *, workers: tuple[int, int] = (1, 1),
                            tie_deg: float = V2_TIE_DEG, shift_max: int = V2_SHIFT_MAX, inputs: dict | None = None,
                            log=None) -> dict:
    """規則 coarse_min_theta_v2 (plan tooling-nozzle-mknot-fast §4.1) で M_knot を選ぶ。

    eval_coarse・eval_prod: M_K → 行 (`evaluate_mknot` の形; `ok`・`reasons`・`theta_max_deg` を見る)。粗い網 (N_s、L_s) と生産の網
    (N_p、カードの L_total_tol) の評価を引数で受ける (要求カードは `make_eval_fn(base, N_s, L_total_tol=L_s)`・`make_eval_fn(base, N_p)`;
    試験は合成の関数)。M_A・M_d: 探索の範囲 (`search_range` の戻り)。workers: (粗い網, 生産の網) の並列の数 (2 つの並列を同時に持たない)。
    inputs: 実効の入力の記録 (r_t・N_s・L_s・N_p・L_p・x_end_margin; 記録だけ)。
    例外 (MknotSelectionError、表つき): 粗い網の長さ合わせが全点で未収束・粗い網で成立点が無い・候補が無い・最終の代表が残る谷が無い。
    戻り: {rule, chosen, chosen_row (生産の網の行), theta_max_deg, distance (推定の距離; 無限大なら None), reason, theta_min, tie_set,
    valleys (谷ごとの確認とずらしの記録), final_three (選んだ点と両隣の生産の網の行), coarse_table (粗い掃引の全行), confirm_table
    (生産の網の全行), candidates, intervals, grid, n_eval_coarse, n_eval_prod, settings, inputs, warnings, elapsed_s}。"""
    t0 = time.time()
    ks = v2_grid(M_A, M_d)
    warnings: list[str] = []
    if log:
        log(f"粗い掃引: {len(ks)} 点 ({v2_M(ks[0]):.2f}〜{v2_M(ks[-1]):.2f}、0.05 刻み、並列 {workers[0]})")
    with Evaluator(eval_coarse, workers[0], log) as evc:
        crow = evc.many([v2_M(k) for k in ks], "coarse")
        n_c = evc.n_eval
    ctab = sorted(crow, key=lambda r: r["M_K"])
    if all(_is_length_unconverged(r) for r in ctab):
        raise MknotSelectionError(f"M_knot を選べない ({RULE_NAME_V2}): 粗い網の長さ合わせ (from_length) が全点で未収束", ctab)
    if not any(_feasible(r) for r in ctab):
        raise MknotSelectionError(f"M_knot を選べない ({RULE_NAME_V2}): 粗い網に成立点が無い (制約を満たす点が 0)", ctab)
    ana = v2_coarse_analysis(ctab)
    if not ana["candidates_k"]:
        raise MknotSelectionError(f"M_knot を選べない ({RULE_NAME_V2}): 候補が無い (粗い網で両隣とも成立している成立点が 0)", ctab)
    vs = v2_valleys(ana)
    if len(vs) >= 2:
        warnings.append(f"谷 (候補の区間ごとの θ_max の局所極小) が {len(vs)} 個 (粗い網の代表 "
                        + ", ".join(f"{v2_M(v['rep_k']):.2f}" for v in vs) + ") — 全部を生産の網で確かめた")
    if log:
        log(f"候補の区間 {len(ana['intervals_k'])} 個・谷 {len(vs)} 個 (代表 " + ", ".join(f"{v2_M(v['rep_k']):.2f}" for v in vs)
            + f") を生産の網で確かめる (並列 {workers[1]})")
    with Evaluator(eval_prod, workers[1], log) as evp:
        valleys = v2_confirm(vs, lambda Ms: evp.many(Ms, "confirm"), shift_max=shift_max, log=log)
        prod = dict(evp.rows)
        n_p = evp.n_eval
    ptab = [prod[k] for k in sorted(prod)]
    inf_k = ana["infeasible_k"]
    for v in valleys:
        if v["final"] is not None:
            r = prod[_key(v["final"])]
            v.update(theta_max_deg=float(r["theta_max_deg"]), distance=v2_distance(v["final_k"], inf_k),
                     distance_steps=v2_distance_steps(v["final_k"], inf_k))
        else:
            v.update(theta_max_deg=None, distance=None, distance_steps=None)
    finals = [{"M_K": v["final"], "theta_max_deg": v["theta_max_deg"], "distance_steps": v["distance_steps"], "valley": i}
              for i, v in enumerate(valleys) if v["final"] is not None]
    if not finals:
        raise MknotSelectionError(
            f"M_knot を選べない ({RULE_NAME_V2}): 生産の網で最終の代表が残る谷が無い ("
            + "; ".join(f"谷 {v['start']:.2f}: {v['discard_reason']}" for v in valleys) + ")",
            ctab + [dict(r, role="confirm") for r in ptab])
    pk = v2_pick(finals, tie_deg=tie_deg)
    best = pk["chosen"]
    ck = v2_k(best["M_K"])
    n_disc = sum(v["final"] is None for v in valleys)
    if n_disc:
        warnings.append(f"生産の網の確認で谷を {n_disc} 個捨てた")
    dtxt = "∞ (粗い網に不成立の点なし)" if math.isinf(best["distance_steps"]) else f"{best['distance_steps'] / V2_K_PER_UNIT:.2f}"
    reason = (f"代表が残った谷 {len(finals)} 個 (全 {len(valleys)} 個) の生産の網の θ_max の最小 {pk['theta_min']:.5f}° から "
              f"{tie_deg:g}° 以内の同点の集合 {len(pk['tie_set'])} 個 ({', '.join(f'{m:.2f}' for m in pk['tie_set'])}) の中で "
              + ("唯一の谷" if len(pk["tie_set"]) == 1 else "推定の距離が最大 (同じなら θ_max・M_K の小さい方)")
              + f" → {v2_M(ck):.2f} (θ_max {best['theta_max_deg']:.5f}°、境目までの推定の距離 {dtxt})")
    out = {"rule": RULE_NAME_V2, "M_A": float(M_A), "M_d": float(M_d),
           "settings": {"k_per_unit": V2_K_PER_UNIT, "step": 1.0 / V2_K_PER_UNIT, "edge": V2_EDGE_STEPS / V2_K_PER_UNIT,
                        "dist_offset": V2_DIST_OFFSET_STEPS / V2_K_PER_UNIT, "tie_deg": tie_deg, "shift_max": shift_max,
                        "sin_min": SIN_MIN, "margin_min_deg": MARGIN_MIN_DEG},
           "inputs": dict(inputs or {}),
           "chosen": v2_M(ck), "chosen_row": prod[_key(v2_M(ck))], "theta_max_deg": float(best["theta_max_deg"]),
           "distance": (None if math.isinf(best["distance_steps"]) else best["distance_steps"] / V2_K_PER_UNIT),
           "reason": reason, "theta_min": pk["theta_min"], "tie_set": pk["tie_set"], "chosen_valley": best["valley"],
           "valleys": valleys,
           "final_three": [prod[_key(v2_M(k))] for k in (ck - 1, ck, ck + 1)],
           "coarse_table": ctab, "confirm_table": ptab,
           "candidates": [v2_M(k) for k in ana["candidates_k"]],
           "intervals": [[v2_M(a), v2_M(b)] for a, b in ana["intervals_k"]],
           "grid": [v2_M(ks[0]), v2_M(ks[-1]), len(ks)], "n_eval_coarse": n_c, "n_eval_prod": n_p,
           "warnings": warnings, "elapsed_s": time.time() - t0}
    if log:
        log(f"選択 {out['chosen']:.2f}: {reason}")
    return out
