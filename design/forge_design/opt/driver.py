"""③ベル TOP 幾何 dv のバッチ MOO ドライバ (Phase 2 §5 ステップ 4–5)。(D) 直接モード (`geometry.bell_mode: direct_bspline`、
dv = theta_max_deg・theta_e_deg・L_over_rt) も同じ TP の評価経路で回る (plan bell-two-mode §5.1 #3f; 壁が成立しなければ
INFEASIBLE と理由、前進 MOC ゲートの指標は台帳に記録するだけ)。既存の TP run のゲートのかけ直しは `--gate-only`。

DOE → forge 評価 (VERDICT ゲート) → KRG (f1=−η のみ; f2=L/rt は dv そのもの
なので厳密値) → NSGA-II + EHVI infill → 反復。使い方:

  design/.venv-opt/bin/python -m forge_design.opt.driver \
      case/40.nozzle_design_tool/problem_bell_top.yaml \
      case/40.nozzle_design_tool/run_NNNN_moo_<tag> \
      --seed-res <収束 res_*.h5> [--n-doe 24 --n-iter 13 --batch 2]

設計 (子 plan §4):
- **ゲート**: メッシュ品質 FAIL / 発散 (DIVERGED・非有限メトリクス) / quasisteady
  非 STEADY は評価失敗としてサロゲートに入れない。残差の「NOT CONVERGED
  (stalled)」プラトーはノイズ床挙動として既知 (run_0071) なので、DIVERGED
  でなく ALL STEADY なら採用する (実装定義。子 plan 変更ログ参照)。さらに
  **物理ゲート**: ṁ/ṁ_1D (チョーク流量比 = 実効 Cd) が [0.94, 1.005]、
  η∈(0.5, 1.02) — cross-geometry warm start は発散せず**不始動流という偽
  アトラクタ** (machmax~1.4, CF~7) に落ちることがあり (run_0073 smoke で実測)、
  これは quasisteady だけでは弾けないため。
- **warm seed**: 正規化 dv 距離が最近傍の PASS 評価の最終 res から interp
  warm start。初期 DOE は --seed-res (基準 run_0072 等) から。
- **段階起動を常用**: cross-geometry warm start を陰解法 cfl4 に直投入すると
  発散 (run_0019) または不始動化 (run_0073 smoke) するため、**全評価を
  soft 段 (convMethod 0 + cfl 0.5, 3000 step) → 本設定 (12000 step)** の
  2 段で回す (run_0019/0024/0034 の実証レシピの常用化)。
- **決定性**: LHS/pymoo/KRG とも seed 固定。台帳は ledger.jsonl (逐次追記)。

TP 経路 (2026-10-11, plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2b; 問題が `gas.model: semiperfect` + `gas.exhaust`):
- 物理ゲートの分母は同じ凍結 TP の ṁ_1D・F_1D (`runner.ideal_1d_tp`; 定比熱の式は使わない)。
- 評価は `runner.run_staged_tp` (soft → mid → 本段を同じ run・同じメッシュで、段の間は restart_field.py、StageManifest)。
  **全点を冷間の段階起動** (準 1D の等エントロピー IC → soft → mid → 本段) で回し、cross-mesh の warm start は使わない
  (`--seed-res` も無視)。`interp_field.py` は生座標の最近傍なので、形状が変わると新旧の壁の間の帯に旧い壁の値 (u≈0・回復温度) が入り、
  本段で出口の壁の角から逆流が育って発散した (case/68 run_0010 doe_001 → 同じ形の冷間 run_0011 は全ゲート PASS、plan §4.1)。
  段の長さは同じなので費用は変わらず、評価が評価順に依存しない。
  停止規則: 同じ fail_class (SUSPECT または RETRYABLE) が 3 評価続く、または DOE で SUSPECT + RETRYABLE が DOE 点数の 20 % を超えたら止める。
- 受理は `tp_gates` の (i)〜(v) を**すべて**満たすとき (2026-10-11 組み替え、plan §4.1・§5.1 #2c'、diagnostician の記録
  notes/reviews/2026-10-11-bell-2b-gates-diagnostician.md): (i) 本段の区間の `check_convergence.py --segment` の VERDICT は記録し、
  不合格は DIVERGED と有意な末尾上昇 (rising) と判定不能だけ (NOT CONVERGED (stalled) は GATES.json に併記)、(ii) 推力系 4 列
  (両推力会計・η の出口面と壁) が drift・osc 1.5e-4、ṁ が 1e-3 で `check_quasisteady.py --series-csv` の STEADY、単調な列は
  漸近値の最終比も閾値以内、末尾窓に 8 点以上 (cv_residual は列にしない)、(iii) Phase 2 の物理ゲート (ṁ/ṁ_1D・η)、
  (iv) 最終スナップショットの |CV 残差| ≤ 1e-4、(v) メッシュ品質 PASS (AR > 1000 の格子は AR ≤ 5000 の例外の条件も) と
  局所 y1+ (`check_wall_resolution.py`) PASS。η は出口面会計を主、壁会計を照合 (台帳には両方と cv)。
  Euler (`evaluate.model: euler`、slip 壁; plan bell-two-mode §5.1 #6 (W1)) の run では局所 y1+ を回さず、GATES.json に
  「対象外」の理由を書いて ok 扱い ((i)〜(iv)・メッシュ品質は同じ)。
- 評価ごとの問題 YAML は campaign_dir に書くので、問題 YAML 基準の相対パス (`gas.exhaust`・`gas.species_db`・
  `geometry.rao_wall_file`) は書き出す前に絶対パスへ直す (以前は切れていた; 引き継ぎ (b))。
CPG 経路 (既定) の挙動は変えない。
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

from ..evaluate import config_builder as CB
from ..evaluate import runner
from ..geometry.wall import NozzleWall
from ..probdef import load_problem
from .doe import lhs
from .ehvi import hypervolume2d, nondominated_mask
from .moo import propose_infill
from .surrogate import KrigingSet

DV_ORDER = ("theta_a_deg", "theta_e_deg", "L_over_rt")
# (D) 直接モード (geometry.bell_mode: direct_bspline; plan bell-two-mode §4.3・§5.1 #3f) の設計変数 (L・θ_e・θ_max、x_m は内部決定)
DV_ORDER_DIRECT = ("theta_max_deg", "theta_e_deg", "L_over_rt")


def dv_order_of(raw: dict) -> tuple:
    """問題 YAML の dict → 設計変数の並び (bell_mode 無し = TOP の DV_ORDER、direct_bspline = DV_ORDER_DIRECT)。"""
    mode = (raw.get("geometry") or {}).get("bell_mode")
    if mode is None:
        return DV_ORDER
    if mode == "direct_bspline":
        return DV_ORDER_DIRECT
    if mode == "moc_rao":
        # (M) の摂動の探索 (plan bell-two-mode §5.1 #8c): pert_a1..K と固定の L
        return runner.rao_perturb_dv(raw) + ("L_over_rt",)
    raise ValueError(f"geometry.bell_mode '{mode}' の設計変数の並びは未定義")

# 問題 YAML の所在を基準に解決する相対パスのキー (probdef: gas.species_db・gas.exhaust・geometry.rao_wall_file [(M) の壁ファイル、
# plan bell-two-mode §5.1 #6])
PROBLEM_PATH_KEYS = (("gas", "exhaust"), ("gas", "species_db"), ("geometry", "rao_wall_file"), ("geometry", "rao_offset_file"))


def absolutize_problem_paths(raw: dict, base_dir) -> dict:
    """問題 YAML の dict の相対パス (`PROBLEM_PATH_KEYS`) を `base_dir` 基準の絶対パスに直した写しを返す。
    評価ごとの YAML を別ディレクトリに書くと相対パスが切れる (plan bell-two-mode §5.1 #2b 引き継ぎ (b))。
    文字列でない値は ValueError (黙って書き出さない)。"""
    out = json.loads(json.dumps(raw))      # deep copy
    for path in PROBLEM_PATH_KEYS:
        cur = out
        for k in path[:-1]:
            cur = cur.get(k) if isinstance(cur, dict) else None
        if not isinstance(cur, dict) or cur.get(path[-1]) in (None, ""):
            continue
        v = cur[path[-1]]
        if not isinstance(v, str):
            raise ValueError(f"{'.'.join(path)} はパスの文字列であること (値 {v!r})")
        pth = Path(v)
        if not pth.is_absolute():
            cur[path[-1]] = str((Path(base_dir) / pth).resolve())
    return out


# ---------------------------------------------------------------- TP の受理ゲート (plan bell-two-mode §4.1 ゲート (i)〜(v))
# 2026-10-11 組み替え (diagnostician [記録 notes/reviews/2026-10-11-bell-2b-gates-diagnostician.md]、§5.1 #2c の継続 run が A):
# 残差は記録 (FAIL は DIVERGED と有意な末尾上昇だけ)、推力系 STEADY 1.5e-4・ṁ 1e-3、最終 |cv| ≤ 1e-4、η は出口面会計を主。
TP_PHYS_BOUNDS = {"mdot_ratio_1d": (0.94, 1.005), "eta": (0.5, 1.02)}   # (iii) Phase 2 の物理ゲート (CPG 経路と同じ値)
TP_CV_TOL = 1.0e-4                                                     # (iv) 最終スナップショットの |CV 残差| ≤ 1e-4 (F_1D 正規化)
FAIL_STOP_N = 3             # 同じ fail_class (SUSPECT・RETRYABLE) が続いたらキャンペーンを止める評価数 (`_stop_on_repeated_failures`)
DOE_FAIL_STOP_FRAC = 0.20   # DOE の SUSPECT + RETRYABLE がこの割合を超えたら止める (`_stop_on_doe_failure_rate`)
# (ii) 推力系 4 列 (runner.THRUST_QS_COLS) の閾値: drift・osc とも 1.5e-4 (相対、末尾窓 4 割)。問題 YAML の evaluate.quasisteady で上書き
TP_QS_DEFAULT = {"tail": 0.4, "drift": 1.5e-4, "osc": 1.5e-4, "min_snaps": 4}
TP_QS_MDOT = {"drift": 1.0e-3, "osc": 1.0e-3}                          # (ii) ṁ (runner.THRUST_QS_MDOT_COLS) の閾値
TP_QS_MIN_TAIL = 8          # (ii) 末尾窓の点数の下限 (本段の出力間隔は本段長の 1/20 以下 = runner.TP_MAIN_OUT_DIV)
# 単調な列は増分の幾何減衰から外挿した漸近値と最終値の差 (最終比) も列の drift の閾値以内 (classify の実効の単調許容は drift/2)
_ASYM_RE = r"最終比\s*([+-]?[0-9.]+(?:[eE][+-]?\d+)?)\s*%"


def _tool(name: str) -> str:
    return str(runner.FORGE_TOOLS / name)


def gate_convergence(run_dir) -> dict:
    """(i) 本段の区間 (stage_manifest.json の最後の区間) の `check_convergence.py --segment` の VERDICT を**記録**する。
    不合格は DIVERGED (NaN/Inf) と有意な末尾上昇 (列の trend が rising — check_convergence の `_tail_rise` の有意性込み) と判定不能
    (manifest 無し・必須列の欠け・末尾窓の代表値 0・VERDICT 行が無い) だけ。`NOT CONVERGED` (stalled / still converging) は合格とし
    `class` と `verdict` に併記する (plan §4.1 (i)、⑤ SERN §4.39 と同型: 基準ケースのプラトーは ω 残差の高周波ジッタと遅い単調
    トランジェントで、推力には 1e-4 以下しか効かない)。出力を CONVERGENCE_SEGMENT_VERDICT.txt に残す。"""
    import re
    r = subprocess.run([sys.executable, _tool("check_convergence.py"), str(run_dir), "--segment"],
                       capture_output=True, text=True, env=runner._ENV)
    txt = r.stdout + r.stderr
    Path(run_dir, "CONVERGENCE_SEGMENT_VERDICT.txt").write_text(txt)
    heads = [l for l in txt.splitlines() if l.startswith("=== ") and "->" in l]
    verdict = heads[-1].split("->", 1)[1].strip(" =") if heads else "判定不能 (VERDICT 行が無い)"
    seg = [l.strip() for l in txt.splitlines() if "[segment] 判定区間" in l]
    trends, rising, nan_cols = {}, [], []
    for line in txt.splitlines():
        m = re.match(r"\s+(rms_\w+)\s*:\s*(.*)$", line)
        if not m:
            continue
        col, msg = m.group(1), m.group(2)
        t = re.search(r"drop=\s*[-\d.]+dec\s+(rising|falling|flat)", msg)
        if t:
            trends[col] = t.group(1)
            if t.group(1) == "rising":
                rising.append(col)
        elif "NaN/Inf" in msg:
            nan_cols.append(col)
    undecidable = (not heads) or ("判定不能" in txt) or (not seg)
    diverged = "DIVERGED" in verdict or bool(nan_cols)
    if diverged:
        cls = "DIVERGED"
    elif undecidable:
        cls = "UNDECIDABLE"
    elif rising:
        cls = "RISING"
    elif verdict.startswith("PASS") and r.returncode == 0:
        cls = "PASS"
    elif "stalled" in verdict:
        cls = "STALLED"
    elif "still converging" in verdict:
        cls = "CONVERGING"
    else:
        cls = "NOT CONVERGED"
    ok = cls in ("PASS", "STALLED", "CONVERGING", "NOT CONVERGED")
    return {"ok": ok, "verdict": verdict, "class": cls, "segment": seg[-1] if seg else None, "trends": trends,
            "rising": rising, "file": "CONVERGENCE_SEGMENT_VERDICT.txt",
            "note": ("" if cls == "PASS" else f"残差 {verdict} ({'併記' if ok else '不合格'})")}


def _qs_run(series_csv: Path, cols, q: dict) -> tuple:
    """`check_quasisteady.py --series-csv` を 1 回回して (出力, 列ごとの VERDICT, 列ごとの detail, OVERALL 行, rc) を返す。"""
    r = subprocess.run([sys.executable, _tool("check_quasisteady.py"), "--series-csv", str(series_csv),
                        "--series-cols", ",".join(cols), "--tail", repr(float(q["tail"])),
                        "--drift", repr(float(q["drift"])), "--osc", repr(float(q["osc"])),
                        "--min-snaps", str(int(q["min_snaps"]))],
                       capture_output=True, text=True, env=runner._ENV)
    txt = r.stdout + r.stderr
    per, det = {}, {}
    for line in txt.splitlines():
        s = line.strip()
        for c in cols:
            if s.startswith(c + " ") or s.startswith(c + ":"):
                per[c] = s.split()[-1]
                det[c] = s.split(":", 1)[1].rsplit(None, 1)[0].strip() if ":" in s else ""
    overall = [l for l in txt.splitlines() if l.startswith("OVERALL:")]
    return txt, per, det, (overall[-1] if overall else None), r.returncode


def gate_quasisteady(series_csv, cols=None, qs: dict | None = None, mdot_cols=None) -> dict:
    """(ii) 推力系 4 列 (両推力会計・η の出口面と壁; `runner.THRUST_QS_COLS`) を drift・osc 1.5e-4、ṁ (`runner.THRUST_QS_MDOT_COLS`)
    を 1e-3 で `check_quasisteady.py --series-csv` が STEADY。加えて (a) 末尾窓 (tail 4 割) に `TP_QS_MIN_TAIL` 点以上、(b) 単調な列は
    漸近値の最終比 (check_quasisteady の外挿) が列の drift の閾値以内。cv_residual は 0 を横切るので列にしない (iv で見る)。
    出力 (2 回分) を series_csv の隣の QUASISTEADY_THRUST_VERDICT.txt に残す。`qs` は推力系の閾値の上書き (tail・drift・osc・
    min_snaps)、`qs["mdot_drift"]`・`qs["mdot_osc"]` で ṁ の閾値を上書き。"""
    import math
    import re
    cols = tuple(cols or runner.THRUST_QS_COLS)
    mcols = tuple(runner.THRUST_QS_MDOT_COLS if mdot_cols is None else mdot_cols)
    qs = dict(qs or {})
    q = {**TP_QS_DEFAULT, **{k: v for k, v in qs.items() if k in TP_QS_DEFAULT}}
    qm = {**q, "drift": float(qs.get("mdot_drift", TP_QS_MDOT["drift"])), "osc": float(qs.get("mdot_osc", TP_QS_MDOT["osc"]))}
    series_csv = Path(series_csv)
    per, det, asym, problems, txts = {}, {}, {}, [], []
    for group, gcols, gq in (("thrust", cols, q), ("mdot", mcols, qm)):
        if not gcols:
            continue
        txt, p_, d_, ov, rc = _qs_run(series_csv, gcols, gq)
        txts.append(f"##### {group}: cols {','.join(gcols)} drift {gq['drift']:g} osc {gq['osc']:g} tail {gq['tail']:g}\n{txt}")
        per.update(p_)
        det.update(d_)
        if ov is None or len(p_) != len(gcols):
            problems.append(f"{group}: 判定不能 (OVERALL 行または列の VERDICT が無い)")
        for c in gcols:
            v = p_.get(c)
            if v is not None and v != "STEADY":
                problems.append(f"{c}: {v}")
            m = re.search(_ASYM_RE, d_.get(c, ""))
            if m:
                rel = float(m.group(1)) / 100.0
                asym[c] = {"final_ratio": rel, "tol": float(gq["drift"]), "ok": abs(rel) <= float(gq["drift"])}
                if not asym[c]["ok"]:
                    problems.append(f"{c}: 単調で漸近値の最終比 {rel:+.2e} が閾値 {gq['drift']:g} を超える")
    n_rows = 0
    try:
        import csv as _csv
        with open(series_csv) as fh:
            n_rows = sum(1 for _ in _csv.DictReader(fh))
    except OSError:
        problems.append("時系列 CSV が読めない")
    n_tail = max(3, int(math.ceil(float(q["tail"]) * n_rows))) if n_rows else 0
    if n_tail < TP_QS_MIN_TAIL:
        problems.append(f"末尾窓の点数 {n_tail} < {TP_QS_MIN_TAIL} (全 {n_rows} 点; 本段の出力間隔を本段長の 1/20 以下に)")
    (series_csv.parent / "QUASISTEADY_THRUST_VERDICT.txt").write_text("\n".join(txts))
    ok = not problems and all(per.get(c) == "STEADY" for c in cols + mcols)
    return {"ok": ok, "per_column": per, "detail": det, "asymptote": asym, "n_rows": n_rows, "n_tail": n_tail,
            "problems": problems, "overall": "ALL STEADY" if ok else "NOT STEADY: " + "; ".join(problems),
            "thresholds": {"thrust": q, "mdot": {"drift": qm["drift"], "osc": qm["osc"]}, "min_tail": TP_QS_MIN_TAIL},
            "file": "QUASISTEADY_THRUST_VERDICT.txt"}


def gate_physics(metrics: dict, bounds: dict | None = None) -> dict:
    """(iii) Phase 2 の物理ゲート: ṁ/ṁ_1D (実効 Cd) と η (出口面会計) が範囲内 (開区間)。欠損・非有限は不合格。"""
    b = {**TP_PHYS_BOUNDS, **(bounds or {})}
    eta, mrel = metrics.get("eta_cf"), metrics.get("mdot_ratio_1d")
    fin = all(v is not None and np.isfinite(v) for v in (eta, mrel))
    ok = fin and b["eta"][0] < eta < b["eta"][1] and b["mdot_ratio_1d"][0] < mrel < b["mdot_ratio_1d"][1]
    return {"ok": bool(ok), "eta": eta, "mdot_ratio_1d": mrel, "bounds": b, "missing": not fin}


def gate_cv(metrics: dict, tol: float = TP_CV_TOL) -> dict:
    """(iv) 収束 (ユーザ決定 2026-10-11): **最終スナップショット**の |CV 残差| = |F_wall − F_exit|/F_1D ≤ tol (既定 1e-4; 旧 1e-3 の
    精度ゲートを吸収)。(ii) の推力系 STEADY と組で使う (単一スナップショットの |cv| は両会計の交差点で偽合格する — run_0003 の
    −5.2e-7 は交差点)。開境界の粘性応力・壁の摩擦を評価できなかった (cv_closable False) ものは不合格。"""
    cv = metrics.get("cv_residual")
    ok = bool(metrics.get("cv_closable")) and cv is not None and np.isfinite(cv) and abs(cv) <= tol
    return {"ok": ok, "cv_residual": cv, "tol": tol, "closable": bool(metrics.get("cv_closable"))}


def gate_mesh_quality(run_dir) -> dict:
    """(v-a) メッシュ品質 (prepare が書く MESH_QUALITY.txt の VERDICT が PASS)。VERDICT の AR の上限が 1000 を超える (AGENTS.md の
    AR ≤ 5000 の例外を使う) 格子は、例外の条件 (AR > 1000 の全セルでスキュー ≤ 0.1、入口・出口面は平面; prepare が書く
    MESH_AR_RELAX.json の ok) も要る。記録が無ければ判定不能 = 不合格 (plan §4.1 (v))。"""
    import re
    f = Path(run_dir) / "MESH_QUALITY.txt"
    txt = f.read_text() if f.exists() else ""
    v = [l.strip() for l in txt.splitlines() if l.strip().startswith("VERDICT")]
    ok = bool(v) and v[-1].startswith("VERDICT: PASS")
    out = {"ok": ok, "verdict": v[-1] if v else "判定不能 (MESH_QUALITY.txt に VERDICT が無い)"}
    m = re.search(r"AR<=\s*([0-9.]+)", v[-1]) if v else None
    ar_lim = float(m.group(1)) if m else None
    out["ar_limit"] = ar_lim
    if ok and (ar_lim is None or ar_lim > runner.AR_RELAX_AR):
        rf = Path(run_dir) / "MESH_AR_RELAX.json"
        if not rf.exists():
            out.update(ok=False, ar_relax="判定不能 (AR の上限 > 1000 なのに MESH_AR_RELAX.json が無い)")
        else:
            c = json.loads(rf.read_text())
            out["ar_relax"] = {k: c.get(k) for k in ("ok", "n_high_ar", "max_skew_high_ar", "n_violations", "inlet_planar",
                                                     "outlet_planar", "skew_max", "ar_threshold")}
            out["ok"] = bool(c.get("ok"))
    return out


def gate_wall_resolution(run_dir, groups: str = "wall", target: float = 1.0) -> dict:
    """(v-b) 局所 y1+ (`check_wall_resolution.py --groups wall`; 第一内部ノード距離と接線 traction、ソルバの ypls は使わない)
    が PASS (y1+ > target の点の割合が既定の許容 2 % 以下)。出力を WALL_RESOLUTION_VERDICT.txt に残す。"""
    r = subprocess.run([sys.executable, _tool("check_wall_resolution.py"), str(run_dir), "--groups", groups,
                        "--target", repr(float(target))], capture_output=True, text=True, env=runner._ENV)
    txt = r.stdout + r.stderr
    Path(run_dir, "WALL_RESOLUTION_VERDICT.txt").write_text(txt)
    v = [l.strip() for l in txt.splitlines() if l.strip().startswith("VERDICT")]
    ok = r.returncode == 0 and bool(v) and v[-1].startswith("VERDICT: PASS")
    over = [l.strip() for l in txt.splitlines() if "y1+ >" in l]
    return {"ok": ok, "verdict": v[-1] if v else "判定不能 (VERDICT 行が無い)", "detail": over[-1] if over else None,
            "file": "WALL_RESOLUTION_VERDICT.txt"}


def combine_tp_gates(g: dict) -> tuple:
    """ゲート結果 {名前: {"ok": ...}} → (PASS/FAIL, fail_class, note)。分類は CPG 経路と同じ語彙:
    RETRYABLE = 未収束・未定常 (段や CFL を変えた再投入で回復し得る)、SUSPECT = 計測・格子の問題 (CV 残差・メトリクス欠損・
    Cd > 1・メッシュ品質・y1+; 要人間)、INFEASIBLE = 収束・定常なのに物理ゲート外 (不始動流などの偽解)。"""
    failed = [k for k in ("convergence", "quasisteady", "physics", "cv", "mesh_quality", "wall_resolution")
              if not g.get(k, {}).get("ok", False)]
    if not failed:
        return "PASS", None, ""
    note = "gate: " + ", ".join(failed)
    if "convergence" in failed or "quasisteady" in failed:
        return "FAIL", "RETRYABLE", note
    ph = g.get("physics", {})
    if ph.get("missing") or "cv" in failed or "mesh_quality" in failed or "wall_resolution" in failed:
        return "FAIL", "SUSPECT", note
    if ph.get("mdot_ratio_1d") is not None and ph["mdot_ratio_1d"] >= ph["bounds"]["mdot_ratio_1d"][1]:
        return "FAIL", "SUSPECT", note + " (Cd > 1 は物理的にあり得ない: 計測を疑う)"
    return "FAIL", "INFEASIBLE", note


# Euler (slip 壁) の run で壁解像ゲート (v-b) を対象外にする理由 (plan bell-two-mode §5.1 #6 (W1)、GATES.json に残す)
EULER_WALL_RESOLUTION_SKIP = ("対象外: Euler (evaluate.model euler、slip 壁) は粘性・境界層が無く、y1+ (第一内部ノード距離 × 壁の接線応力) "
                              "が定義されない (slip の壁ダンプには twall も無い)")


def gate_wall_resolution_skipped(reason: str = EULER_WALL_RESOLUTION_SKIP) -> dict:
    """(v-b) を回さない run の記録 (ok 扱い・`skipped` True・理由)。check_wall_resolution.py は起動しない。"""
    return {"ok": True, "skipped": True, "verdict": reason, "detail": None, "file": None}


def tp_gates(run_dir, metrics: dict, qs: dict | None = None, cv_tol: float = TP_CV_TOL,
             phys_bounds: dict | None = None, flow_model: str = "ns") -> dict:
    """TP の評価 1 点の受理判定 (i)〜(v)。各ゲートの結果と総合 (status・fail_class・note) を返し、GATES.json に残す。
    `flow_model` (`runner.flow_model` の値): `"euler"` (slip 壁) では (v-b) 壁解像を回さず、理由つきで ok 扱いにする
    (`gate_wall_resolution_skipped`)。(i)〜(iv)・(v-a) は NS と同じ。metrics の `wall_model` (推力の会計の内壁の扱い) と
    `flow_model` が食い違えば ValueError (Euler の run を NS として、またはその逆で判定しない)。推力の会計が無い metrics
    (最終 step の会計が出せなかった run) は照合せず、(iii)(iv) が欠損で落とす。"""
    if flow_model not in runner.FLOW_MODELS:
        raise ValueError(f"tp_gates: flow_model '{flow_model}' は未知 ({list(runner.FLOW_MODELS)})")
    wm = metrics.get("wall_model", "viscous" if "eta_exit" in metrics else None)
    if wm is not None and (flow_model == "euler") != (wm == "slip"):
        raise ValueError(f"tp_gates: flow_model {flow_model} と推力の会計の内壁の扱い wall_model {wm} が食い違う")
    run_dir = Path(run_dir)
    g = {"convergence": gate_convergence(run_dir),
         "quasisteady": gate_quasisteady(run_dir / "thrust_history.csv", qs=qs),
         "physics": gate_physics(metrics, phys_bounds),
         "cv": gate_cv(metrics, cv_tol),
         "mesh_quality": gate_mesh_quality(run_dir),
         "wall_resolution": (gate_wall_resolution_skipped() if flow_model == "euler" else gate_wall_resolution(run_dir))}
    status, fc, note = combine_tp_gates(g)
    # 残差の VERDICT は合格でも併記する (NOT CONVERGED (stalled) を GATES.json の上位に残す; plan §4.1 (i))
    rnote = g["convergence"].get("note") or ""
    out = {"status": status, "fail_class": fc, "note": "; ".join(x for x in (note, rnote) if x),
           "residual_verdict": g["convergence"]["verdict"], "residual_class": g["convergence"]["class"],
           "eta_primary": "exit (出口面会計; 壁会計は照合)", "gates": g,
           "thresholds": {"cv_tol": cv_tol, "quasisteady": g["quasisteady"]["thresholds"]}}
    if flow_model != "ns":
        out["flow_model"] = flow_model                       # NS の GATES.json の中身は従来どおり
    (run_dir / "GATES.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    return out


class _F1KrgF2Exact:
    """f1=−η は KRG、f2=L/rt は dv そのもの (厳密・σ≈0) の predict アダプタ。"""

    def __init__(self, krg_f1: KrigingSet, l_index: int) -> None:
        self.krg = krg_f1
        self.l_index = l_index

    def predict(self, X):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        mu1, sg1 = self.krg.predict(X)
        L = X[:, self.l_index]
        mu = np.column_stack([mu1[:, 0], L])
        sg = np.column_stack([sg1[:, 0], np.full(L.shape, 1e-9)])
        return mu, sg


class Campaign:
    def __init__(self, base_yaml, campaign_dir, seed_res, ref=(-0.90, 13.0),
                 seed: int = 0) -> None:
        self.base_yaml = Path(base_yaml)
        self.dir = Path(campaign_dir)
        self.tp = runner._is_tp(load_problem(self.base_yaml))   # TP 経路 (gas.model: semiperfect + gas.exhaust)
        if seed_res is None and not self.tp:
            raise ValueError("CPG 経路は --seed-res (収束 res_*.h5) が必須 (node の冷間 IC は発散)")
        self.dir.mkdir(parents=True, exist_ok=True)
        self.seed_res = Path(seed_res).resolve() if seed_res is not None else None
        self.ref = tuple(ref)
        self.seed = int(seed)
        self.base_raw = yaml.safe_load(self.base_yaml.read_text())
        self.dv_order = dv_order_of(self.base_raw)
        # 固定の設計変数 (`fixed: true`、min・max なし可) は [value, value] (単目的の長さ固定 — plan bell-two-mode §5.1 #8b)。
        # MOO の経路 (`run`) は固定の変数を持たない問題だけ (幅 0 の変数を正規化できない)
        dvs = self.base_raw["dv"]
        self.fixed_mask = np.array([bool(dvs[k].get("fixed", False)) for k in self.dv_order])
        self.bounds = np.array([[float(dvs[k]["value"]), float(dvs[k]["value"])] if dvs[k].get("fixed", False) else
                                [float(dvs[k]["min"]), float(dvs[k]["max"])] for k in self.dv_order])
        self.ledger = self.dir / "ledger.jsonl"
        self.rows: list = []
        if self.ledger.exists():  # 再開: 既存台帳を読み戻す (評価済みはスキップ)
            self.rows = [json.loads(l) for l in self.ledger.read_text().splitlines() if l]
        if self.tp:
            # TP: 同じ凍結 TP の 1D 理想値 (ṁ_1D・F_1D(ε; p_a)) を物理ゲートと η の分母に (定比熱の式は使わない)
            self.ideal_1d = runner.ideal_1d_tp(load_problem(self.base_yaml))
            self.mdot_1d = float(self.ideal_1d["mdot"])
            self.F_1d = float(self.ideal_1d["F"])
            return
        # 1D チョーク流量 (Cd=1): 物理ゲートの分母
        g, cp = float(self.base_raw["gas"]["gamma"]), float(self.base_raw["gas"]["cp"])
        R = cp * (g - 1.0) / g
        s = self.base_raw["spec"]
        At = np.pi * float(s["r_throat"]) ** 2
        self.mdot_1d = (float(s["Pt"]) * At * np.sqrt(g / (R * float(s["Tt"])))
                        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0))))

    # -- 幾何整合の事前フィルタ (CFD なし) ------------------------------------
    def feasible(self, x) -> bool:
        try:
            self._wall(x)
            return True
        except ValueError:
            return False

    def _wall(self, x) -> NozzleWall:
        p = load_problem(self._write_problem(x, self.dir / "_feas_probe.yaml"))
        wall, _ = runner.build_wall(p)
        return wall

    def _write_problem(self, x, path: Path) -> Path:
        # 相対パス (gas.exhaust・gas.species_db) は元の問題 YAML の所在で解決した絶対パスにする (引き継ぎ (b))
        raw = absolutize_problem_paths(self.base_raw, self.base_yaml.resolve().parent)   # deep copy
        for k, v in zip(self.dv_order, np.asarray(x, dtype=float)):
            raw["dv"][k]["value"] = float(v)
        path.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
        return path

    # -- warm seed: 正規化 dv 距離の最近傍 PASS -------------------------------
    def _warm_for(self, x) -> Path:
        lo, span = self.bounds[:, 0], self.bounds[:, 1] - self.bounds[:, 0]
        best, best_d = self.seed_res, None
        xn = (np.asarray(x, dtype=float) - lo) / span
        for r in self.rows:
            if r["status"] != "PASS":
                continue
            dn = float(np.linalg.norm((np.array(r["x"]) - lo) / span - xn))
            if best_d is None or dn < best_d:
                best, best_d = Path(r["run_dir"]) / r["res_file"], dn
        return best

    # -- 1 点評価 (ゲート込み) --------------------------------------------------
    def evaluate(self, x, tag: str) -> dict:
        if self.tp:
            row = self._evaluate_tp(x, tag)
            self._stop_on_repeated_failures()
            return row
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time()
        run_dir = self.dir / tag
        prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "run_dir": str(run_dir), "status": "FAIL",
               "fail_class": None, "warm_from": "", "note": ""}
        try:
            warm = self._warm_for(x)
            row["warm_from"] = str(warm)
            # soft 段 (常用): convMethod 0 + cfl 0.5 で 3000 step 流れを再確立
            stage = self.dir / (tag + "_stgA")
            n_stage = 3000
            runner.prepare(prob, stage, nsteps=n_stage, warm_from=warm)
            cfg = (stage / "solverConfig.yaml").read_text()
            # cfl 0.5・1 次。収束場を必ずダンプさせる (outStepInterval > n_stage だと res_0 しか
            # 出ず、本段が seed 直投入相当になって発散する — smoke で実測したバグ)。値は YAML 上の位置で書き換える
            # (2026-10-15、plan tooling-generator-defaults-and-launch-gate §5.1 #11; 旧: 文字列置換・正規表現、書く文字列は同じ)
            cfg = CB.set_scalars(cfg, {"time.deltaT.cfl": "0.5", "time.deltaT.cfl_pseudo": "0.5", "space.convMethod": "0",
                                       "time.outStepInterval": str(n_stage)})
            # 明示すべきキーを足して書き、段の宣言 RECIPE.json を書く (機種 bell_cpg の soft; config_builder)
            CB.launch_stage(stage, cfg, (stage / "bcondConfig.yaml").read_text(), "soft", "bell_cpg", problem=prob,
                            generator="opt.driver.Campaign.evaluate (CPG の soft)")
            runner.run_forge(stage)
            sres = sorted(stage.glob("res_[0-9]*.h5"),
                          key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
            if (not sres
                    or int("".join(c for c in sres[-1].stem if c.isdigit())) < n_stage
                    or self._diverged(stage, {"eta_cf": 0.9})):
                raise RuntimeError("soft 段が発散または最終場未出力")
            # 本段 (12000 step, 全域 2 次 + cfl4)。soft 段の変換済みの格子を写して使い、soft の最終場を段終了ゲート +
            # restart_field.py で初期場にする (2026-10-15 まで同じ幾何から格子を作り直して interp_field.py; plan
            # tooling-generator-defaults-and-launch-gate §5.1 #19)。ゲート落ちは soft 段の失敗として扱う (RETRYABLE)
            try:
                runner.continue_cpg(prob, stage, run_dir)
            except RuntimeError as e:
                if "段終了ゲート" in str(e):
                    raise RuntimeError(f"soft 段の最終場が段終了ゲートで不合格: {e}") from None
                raise
            CB.launch_stage(run_dir, (run_dir / "solverConfig.yaml").read_text(), (run_dir / "bcondConfig.yaml").read_text(),
                            "main", "bell_cpg", problem=prob, generator="opt.driver.Campaign.evaluate (CPG の本段)")
            runner.run_forge(run_dir)
            out = runner.collect(prob, run_dir)
            verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text()
            qs = self._quasisteady(run_dir)
            eta = out.get("eta_cf")
            mrel = (out.get("mdot_kg_s") / self.mdot_1d
                    if out.get("mdot_kg_s") is not None else None)
            row.update(eta_cf=eta, mdot=out.get("mdot_kg_s"),
                       mdot_ratio_1d=mrel,
                       thrust_N=out.get("thrust_N"), res_file=out.get("res_file", ""),
                       conv="DIVERGED" if "DIVERGED" in verdict else
                            ("PASS" if "PASS" in verdict else "STALLED"),
                       quasisteady=qs)
            n_main = int(self.base_raw["evaluate"].get("nStepOuter", 12000))
            res_step = int("".join(c for c in Path(out.get("res_file", "res_-1"))
                                   .stem if c.isdigit()) or -1)
            row["res_step"] = res_step
            ok = (row["conv"] != "DIVERGED" and qs == "ALL STEADY"
                  and res_step >= n_main                    # 最終場が最後まで到達
                  and eta is not None and np.isfinite(eta)
                  and 0.5 < eta < 1.02                      # 物理ゲート (η)
                  and mrel is not None and 0.94 < mrel < 1.005)  # 物理ゲート (Cd)
            row["status"] = "PASS" if ok else "FAIL"
            if not ok:
                # ゲート落ちの分類 (例外なし経路)。RETRYABLE=条件を変えた再投入で
                # 回復し得る / INFEASIBLE=偽アトラクタ (不始動流) 系で再投入無意味
                # / SUSPECT=メトリクス欠損 (原因不明・要人間)。
                if row["conv"] == "DIVERGED":
                    row["fail_class"] = "RETRYABLE"      # 本段発散 → recipe 弱めて再投入
                    row["note"] = "gate: main-stage DIVERGED"
                elif qs != "ALL STEADY" or res_step < n_main:
                    row["fail_class"] = "RETRYABLE"      # 未定常/未到達 → 延長で回復見込み
                    row["note"] = f"gate: qs={qs} res_step={res_step}/{n_main}"
                elif eta is None or mrel is None or not np.isfinite(eta):
                    row["fail_class"] = "SUSPECT"        # 完走したのにメトリクス欠損
                    row["note"] = f"gate: metrics missing (eta={eta} mrel={mrel})"
                elif mrel > 1.005:
                    # Cd 上限超えは物理的にあり得ない (チョーク流量は超えられない)
                    # → 設計点でなく計測 (積分面/メッシュ) が怪しい。要人間確認
                    row["fail_class"] = "SUSPECT"
                    row["note"] = f"gate: Cd>1 unphysical (mdot_ratio={mrel:.4f} eta={eta:.4f})"
                else:
                    row["fail_class"] = "INFEASIBLE"     # 不始動流等の偽解 (Cd 下限/η 逸脱)
                    row["note"] = f"gate: physics eta={eta:.4f} mdot_ratio={mrel:.4f}"
        except ValueError as e:
            # 幾何構築失敗 (fold, 無効パラメータ等) — 設計空間外 (再投入不可)
            row["fail_class"] = "INFEASIBLE"
            row["note"] = f"geometry: {str(e)[:250]}"
        except RuntimeError as e:
            # 段階起動失敗 (soft 段 NaN/未出力) — 再投入で回復可能性あり
            err_msg = str(e)
            if "soft 段" in err_msg:
                row["fail_class"] = "RETRYABLE"
                row["note"] = f"soft_stage: {err_msg[:200]}"
            else:
                row["fail_class"] = "SUSPECT"
                row["note"] = f"runtime: {err_msg[:250]}"
        except Exception as e:
            # I/O・OOM・その他予期不測 — 一律 SUSPECT (人間が見る)
            row["fail_class"] = "SUSPECT"
            row["note"] = f"{type(e).__name__}: {str(e)[:250]}"
        row["wall_s"] = round(time.time() - t0, 1)
        self.rows.append(row)
        with open(self.ledger, "a") as f:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
        fail_marker = f" [{row['fail_class']}]" if row['fail_class'] else ""
        print(f"[{tag}] {row['status']}{fail_marker} x={np.round(x,3)} eta={row.get('eta_cf')}"
              f" ({row['wall_s']}s) {row['note']}", flush=True)
        return row

    def _evaluate_tp(self, x, tag: str) -> dict:
        """TP の 1 点評価: 段階起動 (`runner.run_staged_tp`; warm は cross-mesh、無ければ冷間) → 推力の CV 会計
        (`runner.collect`) → 受理ゲート (i)〜(v) (`tp_gates`)。台帳の行は CPG 経路と同じ列 + gates。"""
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time()
        run_dir = self.dir / tag
        prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "run_dir": str(run_dir), "status": "FAIL",
               "fail_class": None, "warm_from": "", "note": "", "path": "tp"}
        try:
            # TP は全点冷間 (warm start は使わない — 冒頭の docstring・plan §4.1)。warm_from は空のまま台帳に残す
            runner.run_staged_tp(prob, run_dir, warm_from=None)
            out = runner.collect(prob, run_dir)
            qs = (self.base_raw.get("evaluate") or {}).get("quasisteady")
            g = tp_gates(run_dir, out, qs=qs, flow_model=runner.flow_model(load_problem(prob)))
            row.update(eta_cf=out.get("eta_cf"), eta_wall=out.get("eta_wall"), mdot=out.get("mdot_kg_s"),
                       mdot_ratio_1d=out.get("mdot_ratio_1d"), thrust_N=out.get("thrust_N"),
                       thrust_wall_N=out.get("thrust_wall_N"), cv_residual=out.get("cv_residual"),
                       res_file=out.get("res_file", ""), conv=g["gates"]["convergence"]["class"],
                       residual_verdict=g["residual_verdict"],
                       quasisteady=g["gates"]["quasisteady"]["overall"],
                       gates={k: bool(v["ok"]) for k, v in g["gates"].items()})
            for k in ("moc_gate", "wall_design"):           # (D): 前進 MOC ゲートの指標 (記録だけ) と壁の設計 (x_m 等)
                if k in out:
                    row[k] = out[k]
            row["status"], row["fail_class"], row["note"] = g["status"], g["fail_class"], g["note"]
        except ValueError as e:
            row["fail_class"] = "INFEASIBLE"                # 幾何構築失敗 (設計空間外・(D) の成立区間が空)
            row["note"] = f"geometry: {str(e)[:250]}"
            if isinstance(e, runner.WallInfeasible):
                row["wall_design"] = e.record
        except RuntimeError as e:
            row["fail_class"] = "RETRYABLE" if "段" in str(e) else "SUSPECT"   # 段の失敗 (発散・段終了ゲート) / メッシュ品質等
            nan = sorted(run_dir.glob("res_nan_*.h5"))                       # 最初の NaN の step (ソルバの detectNaN の出力名)
            row["note"] = f"runtime: {str(e)[:250]}" + (f" [NaN step {nan[0].stem.split('_')[-1]}]" if nan else "")
        except subprocess.CalledProcessError as e:  # 下請けツール (interp_field 等) の失敗: 理由は stderr にある
            row["fail_class"] = "SUSPECT"
            tool = Path(str(e.cmd[1] if isinstance(e.cmd, (list, tuple)) and len(e.cmd) > 1 else e.cmd)).name
            err = (e.stderr or e.stdout or "").strip().replace("\n", " ")
            row["note"] = f"CalledProcessError rc={e.returncode} {tool}: {err[-400:]}"
        except Exception as e:  # noqa: BLE001 — I/O・OOM・その他は一律 SUSPECT (人間が見る)
            row["fail_class"] = "SUSPECT"
            row["note"] = f"{type(e).__name__}: {str(e)[:250]}"
        row["wall_s"] = round(time.time() - t0, 1)
        self.rows.append(row)
        with open(self.ledger, "a") as f:
            f.write(json.dumps(row, ensure_ascii=False, default=float) + "\n")
        fail_marker = f" [{row['fail_class']}]" if row['fail_class'] else ""
        print(f"[{tag}] {row['status']}{fail_marker} x={np.round(x, 3)} eta={row.get('eta_cf')}"
              f" ({row['wall_s']}s) {row['note']}", flush=True)
        return row

    def _stop_on_repeated_failures(self, n: int = FAIL_STOP_N) -> None:
        """直近 n 評価が続けて同じ fail_class (SUSPECT = 計測・下請けツール・I/O、RETRYABLE = 段の発散・未定常) なら止める。
        同じ原因で DOE 点を空費しないため (2026-10-11 case/68 run_0010: warm start の interp_field が化学種署名を解決できず
        17 点が 10 秒ずつ落ちた)。失敗点はサロゲートに入らないので、infill が同じ失敗領域を繰り返し提案するのもここで捕まえる。"""
        tail = self.rows[-n:]
        for fc in ("SUSPECT", "RETRYABLE"):
            if len(tail) == n and all(r.get("fail_class") == fc for r in tail):
                raise SystemExit(f"[driver] 直近 {n} 評価が続けて {fc} — 停止 (台帳 {self.ledger} の note を見る。"
                                 "原因を直すか上位に諮ってから、失敗行の扱いを決めて再開する)")

    def _stop_on_doe_failure_rate(self, n_doe: int, frac: float = DOE_FAIL_STOP_FRAC) -> None:
        """DOE の評価で SUSPECT + RETRYABLE が DOE 点数の frac を超えたら止める (INFEASIBLE と壁解像などの FAIL は数えない)。"""
        bad = [r for r in self.rows if r["tag"].startswith("doe_") and r.get("fail_class") in ("SUSPECT", "RETRYABLE")]
        if len(bad) > frac * n_doe:
            raise SystemExit(f"[driver] DOE の SUSPECT + RETRYABLE が {len(bad)} 点 (> {frac:.0%} × {n_doe}) — 停止して上位に諮る")

    def _diverged(self, run_dir: Path, out: dict) -> bool:
        v = (run_dir / "CONVERGENCE_VERDICT.txt")
        txt = v.read_text() if v.exists() else ""
        eta = out.get("eta_cf")
        return "DIVERGED" in txt or eta is None or not np.isfinite(eta)

    def _quasisteady(self, run_dir: Path) -> str:
        r = subprocess.run([sys.executable,
                            str(runner.FORGE_TOOLS / "check_quasisteady.py"), str(run_dir)],
                           capture_output=True, text=True, env=runner._ENV)
        txt = r.stdout + r.stderr
        (run_dir / "QUASISTEADY_VERDICT.txt").write_text(txt)
        overall = [l for l in txt.splitlines() if "OVERALL:" in l]
        if not overall:
            return "UNKNOWN"
        line = overall[-1]
        if "NOT" in line:  # "NOT ALL STEADY" は "ALL STEADY" を部分文字列に含む
            return "NOT ALL STEADY"
        return "ALL STEADY" if "ALL STEADY" in line else line.split("OVERALL:")[-1].strip()

    # -- 目的ベクトル (最小化): f1=−η, f2=L/rt ---------------------------------
    def _XF(self) -> tuple:
        X = np.array([r["x"] for r in self.rows if r["status"] == "PASS"])
        F = np.array([[-r["eta_cf"], r["x"][self.dv_order.index("L_over_rt")]]
                      for r in self.rows if r["status"] == "PASS"])
        return X, F

    # -- メインループ -----------------------------------------------------------
    def run(self, n_doe: int, n_iter: int, batch: int) -> None:
        if self.fixed_mask.any():
            raise ValueError("MOO の経路は固定の設計変数を持たない問題だけ (長さ固定の単目的は --soo)")
        done = {r["tag"] for r in self.rows}
        cand = lhs(self.bounds, 3 * n_doe, seed=self.seed)
        feas = [x for x in cand if self.feasible(x)]
        doe = feas[:n_doe]
        print(f"DOE {len(doe)} 点 (LHS {3*n_doe} 点のうち成立 {len(feas)} 点の先頭)", flush=True)
        for i, x in enumerate(doe):
            tag = f"doe_{i:03d}"
            if tag not in done:
                self.evaluate(x, tag)
                if self.tp:
                    self._stop_on_doe_failure_rate(len(doe))
        for it in range(n_iter):
            X, F = self._XF()
            if len(X) < 4:
                raise RuntimeError("PASS 評価が少なすぎる (<4) — ゲート/レシピを確認")
            hv = hypervolume2d(F, self.ref)
            print(f"--- iter {it}: PASS {len(X)} 点, HV={hv:.5f} ---", flush=True)
            krg = KrigingSet(self.bounds).fit(X, F[:, :1])
            surro = _F1KrgF2Exact(krg, self.dv_order.index("L_over_rt"))
            xs = propose_infill(surro, X, F, self.ref, self.bounds, n_infill=batch,
                                seed=self.seed + 1000 + it, feasible=self.feasible)
            for j, x in enumerate(xs):
                tag = f"inf_{it:02d}_{j}"
                if tag not in done:
                    self.evaluate(x, tag)
        X, F = self._XF()
        mask = nondominated_mask(F)
        pareto = [{"x": X[i].tolist(), "eta_cf": float(-F[i, 0]),
                   "L_over_rt": float(F[i, 1])} for i in np.where(mask)[0]]
        pareto.sort(key=lambda r: r["L_over_rt"])
        summary = {"n_eval": len(self.rows),
                   "n_pass": int(len(X)), "hv": hypervolume2d(F, self.ref),
                   "ref": self.ref, "pareto": pareto}
        (self.dir / "pareto.json").write_text(
            json.dumps(summary, indent=1, ensure_ascii=False))
        print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)


    # -- 単目的 (長さ固定) の再最適化: plan bell-two-mode §5.1 #8b (codex 2026-10-13 の改訂の事前登録) ------
    def _soo_full(self, xa) -> np.ndarray:
        x = self.bounds[:, 0].copy()
        x[~self.fixed_mask] = np.asarray(xa, dtype=float)
        return x

    def _soo_best_history(self) -> list:
        """評価順の各行に、その時点までの PASS の最良 η (b) と PASS の通し番号 K を付けた表。"""
        out, best, K = [], None, 0
        for r in self.rows:
            if r["tag"].startswith("rep_"):        # 再現性の回し直しは探索の履歴に入れない
                continue
            if r["status"] == "PASS" and r.get("eta_cf") is not None:
                K += 1
                best = r["eta_cf"] if best is None else max(best, r["eta_cf"])
            out.append({"n": len(out) + 1, "tag": r["tag"], "status": r["status"], "eta_cf": r.get("eta_cf"), "K": K, "best": best,
                        "wall_s": r.get("wall_s")})
        return out

    def _soo_propose(self, seed: int, n_lhs: int = 2048, n_local: int = 512, min_dist: float = 5e-3):
        """Kriging (−η) の期待改善 (EI) 最大の、成立・未評価の点 (能動変数だけ)。評価済み (失敗を含む) から正規化距離 min_dist 以内は除く。"""
        from scipy.stats import norm
        act = ~self.fixed_mask
        ab = self.bounds[act]
        lo, span = ab[:, 0], ab[:, 1] - ab[:, 0]
        P = [(np.asarray(r["x"], dtype=float)[act], -float(r["eta_cf"])) for r in self.rows
             if r["status"] == "PASS" and r.get("eta_cf") is not None]
        X = np.array([p[0] for p in P]); Y = np.array([p[1] for p in P])
        krg = KrigingSet(ab).fit(X, Y[:, None])
        rng = np.random.default_rng(seed)
        xb = (X[int(np.argmin(Y))] - lo) / span
        local = np.clip(xb + 0.05 * rng.standard_normal((n_local, len(lo))), 0.0, 1.0) * span + lo
        cand = np.vstack([lhs(ab, n_lhs, seed=seed), local])
        mu, sg = krg.predict(cand)
        mu, sg = mu[:, 0], np.maximum(sg[:, 0], 1e-15)
        z = (Y.min() - mu) / sg
        ei = (Y.min() - mu) * norm.cdf(z) + sg * norm.pdf(z)
        Xall = np.array([np.asarray(r["x"], dtype=float)[act] for r in self.rows])
        Xn_all = (Xall - lo) / span
        for i in np.argsort(-ei):
            xn = (cand[i] - lo) / span
            if np.min(np.linalg.norm(Xn_all - xn, axis=1)) < min_dist:
                continue
            if not self.feasible(self._soo_full(cand[i])):
                continue
            return cand[i], float(ei[i]), float(-mu[i]), float(sg[i])
        raise RuntimeError("EI の候補に成立・未評価の点が無い")

    def run_soo(self, n_doe: int = 12, n_max: int = 40, n_confirm: int = 8, window: int = 8, tol: float = 1e-4,
                repeat_best: bool = True, with_start: bool = False) -> dict:
        """§5.1 #8b: LHS n_doe → EI 1 点ずつ。**停止候補** = EI による PASS が window 点以上かつ b_K − b_{K−window} < tol、
        **予算打切り** = 評価 n_max 回 (失敗込み)。その後 **確認** n_confirm 評価 (最大 n_max + n_confirm)、G = 確認後の最良 − 停止時の最良。
        最後に最良形を同じ設定で回し直す (再現性)。要約は `soo_summary.json` (評価ごとに書き直す)。"""
        act = ~self.fixed_mask
        if act.sum() < 1:
            raise ValueError("単目的の経路: 能動の設計変数が無い")
        done = {r["tag"] for r in self.rows}
        cand = lhs(self.bounds[act], 3 * n_doe, seed=self.seed)
        feas = [x for x in cand if self.feasible(self._soo_full(x))]
        doe = feas[:n_doe]
        print(f"[soo] DOE {len(doe)} 点 (LHS {3*n_doe} 点のうち成立 {len(feas)} 点の先頭)"
              + (" + 出発点 (問題 YAML の dv.value)" if with_start else ""), flush=True)
        if with_start and "doe_start" not in done:
            # 出発点 (plan bell-two-mode §5.1 #8c: a = 0 = 補正済み Rao) を DOE の先頭に必ず 1 点
            x0 = np.array([float(self.base_raw["dv"][k]["value"]) for k in self.dv_order])
            if not self.feasible(x0):
                raise ValueError("出発点 (dv.value) が幾何の成立を満たさない")
            self.evaluate(x0, "doe_start")
        for i, xa in enumerate(doe):
            tag = f"doe_{i:03d}"
            if tag not in done:
                self.evaluate(self._soo_full(xa), tag)
                if self.tp:
                    self._stop_on_doe_failure_rate(len(doe))
        n_doe_eval = sum(1 for r in self.rows if r["tag"].startswith("doe_"))
        summ = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #8b", "seed": self.seed, "n_doe": n_doe, "n_max": n_max,
                "n_confirm": n_confirm, "window": window, "tol": tol, "fixed": {k: float(self.bounds[i, 0]) for i, k in
                                                                                 enumerate(self.dv_order) if self.fixed_mask[i]}}
        it = 0
        announced = False
        while True:
            hist = self._soo_best_history()
            n = len(hist)
            stop = soo_stop(hist, n_doe_eval, window, tol, n_max)
            if stop is not None and not announced:
                announced = True
                print(f"[soo] 停止 ({stop['reason']}) n {stop['n']}・最良 η {stop['best']} — 確認の {n_confirm} 評価へ", flush=True)
            if stop is not None and n >= stop["n"] + n_confirm:
                break
            if n >= n_max + n_confirm:
                break
            tag = f"ei_{it:03d}" if stop is None else f"conf_{n - stop['n']:03d}"
            it += 1
            if tag in done:
                continue
            xa, ei, mu, sg = self._soo_propose(seed=self.seed + 1000 + n)
            print(f"[soo] {tag}: EI {ei:.3e} 予測 η {mu:.6f} ± {sg:.1e} x={np.round(xa, 4)}", flush=True)
            self.evaluate(self._soo_full(xa), tag)
            summ["history"] = self._soo_best_history()
            (self.dir / "soo_summary.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
        hist = self._soo_best_history()
        stop = soo_stop(hist, n_doe_eval, window, tol, n_max)
        best_row = max((r for r in self.rows if r["status"] == "PASS" and r.get("eta_cf") is not None), key=lambda r: r["eta_cf"])
        summ.update({"history": hist, "stop": stop, "best": {"tag": best_row["tag"], "x": best_row["x"], "eta_cf": best_row["eta_cf"]},
                     "G": (hist[-1]["best"] - stop["best"]) if stop else None,
                     "confirm": ("A (確認で 1e-4 超の改善 — 停止規則は改善を取り逃した)" if stop and hist[-1]["best"] - stop["best"] > tol else
                                 "B (確認の区間では反例なし)") if stop else None})
        if repeat_best and "rep_best" not in done:
            r2 = self.evaluate(np.asarray(best_row["x"], dtype=float), "rep_best")
            summ["repeat"] = {"tag": "rep_best", "of": best_row["tag"], "status": r2["status"], "eta_cf": r2.get("eta_cf"),
                              "diff": (r2["eta_cf"] - best_row["eta_cf"]) if r2.get("eta_cf") is not None else None}
        elif "rep_best" in done:
            r2 = next(r for r in self.rows if r["tag"] == "rep_best")
            summ["repeat"] = {"tag": "rep_best", "status": r2["status"], "eta_cf": r2.get("eta_cf")}
        (self.dir / "soo_summary.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
        print(json.dumps({k: summ[k] for k in ("stop", "best", "G", "confirm")}, indent=1, ensure_ascii=False, default=float), flush=True)
        return summ


def soo_stop(hist: list, n_doe_eval: int, window: int, tol: float, n_max: int) -> dict | None:
    """単目的の停止 (plan bell-two-mode §5.1 #8b): 履歴の先頭から順に、**最初に**条件を満たした評価の位置を返す (再開しても同じ位置)。
    停止候補 = DOE の後の PASS が window 点以上かつ b_K − b_{K−window} < tol (b_K = PASS の K 点目までの最良)、
    予算打切り = 評価 n_max 回 (失敗込み)。どちらも無ければ None。"""
    for m in range(n_doe_eval + 1, len(hist) + 1):
        h = hist[:m]
        pass_ei = sum(1 for x in h[n_doe_eval:] if x["status"] == "PASS")
        K, bK = h[-1]["K"], h[-1]["best"]
        b_prev = next((x["best"] for x in reversed(h) if x["K"] == K - window), None) if K > window else None
        if pass_ei >= window and bK is not None and b_prev is not None and bK - b_prev < tol:
            return {"reason": "stop_candidate", "n": m, "K": K, "best": bK, "best_K_minus_window": b_prev}
        if m >= n_max:
            return {"reason": "budget_cut", "n": m, "K": K, "best": bK}
    return None


def gate_run(problem, run_dir) -> dict:
    """TP の run 1 本に受理ゲート (i)〜(v) をかけ直す (`runner.collect` → `tp_gates`; GATES.json・各 VERDICT を run に書く)。
    段階起動・継続を CLI で回した run の判定と、ゲートの定数を変えたあとの既存 run の再判定に使う。"""
    p = load_problem(problem)
    if not runner._is_tp(p):
        raise ValueError("gate_run は TP 経路 (gas.model: semiperfect + gas.exhaust) の run 専用")
    out = runner.collect(problem, run_dir)
    return tp_gates(run_dir, out, qs=(p.raw.get("evaluate") or {}).get("quasisteady"), flow_model=runner.flow_model(p))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="③ベルの MOO キャンペーン (TOP dv / (D) 直接モード)")
    ap.add_argument("problem")
    ap.add_argument("campaign_dir", help="キャンペーンのディレクトリ (--gate-only では判定する run)")
    ap.add_argument("--gate-only", action="store_true",
                    help="TP: campaign_dir を既存の run とみなし、受理ゲート (i)〜(v) だけをかけ直す (GATES.json を書く)")
    ap.add_argument("--seed-res", default=None,
                    help="初期 warm start の収束 res_*.h5 (CPG は必須。TP は使わない — 全点を冷間の段階起動)")
    ap.add_argument("--n-doe", type=int, default=24)
    ap.add_argument("--n-iter", type=int, default=13)
    ap.add_argument("--batch", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref", type=float, nargs=2, default=(-0.90, 13.0))
    ap.add_argument("--soo", action="store_true",
                    help="長さ固定の単目的 (η 最大) の再最適化 (plan bell-two-mode §5.1 #8b): --n-doe・--n-max・--n-confirm・--window・--tol")
    ap.add_argument("--n-max", type=int, default=40)
    ap.add_argument("--n-confirm", type=int, default=8)
    ap.add_argument("--window", type=int, default=8)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--with-start", action="store_true", help="--soo: 問題 YAML の dv.value (出発点) を DOE の先頭に必ず評価")
    a = ap.parse_args(argv)
    if a.gate_only:
        g = gate_run(a.problem, a.campaign_dir)
        print(json.dumps({"status": g["status"], "fail_class": g["fail_class"], "note": g["note"],
                          "gates": {k: {kk: v.get(kk) for kk in ("ok", "verdict", "class", "overall", "cv_residual")
                                        if kk in v} for k, v in g["gates"].items()}},
                         indent=1, ensure_ascii=False, default=float))
        return 0 if g["status"] == "PASS" else 1
    c = Campaign(a.problem, a.campaign_dir, a.seed_res, ref=tuple(a.ref), seed=a.seed)
    if a.soo:
        c.run_soo(a.n_doe, a.n_max, a.n_confirm, a.window, a.tol, with_start=a.with_start)
    else:
        c.run(a.n_doe, a.n_iter, a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
