"""plan discretization-moc-axis-limit-and-corrector §6 V5 (事前登録 2026-10-07) の判定。
腕 B (現行の単調壁) run_0143〜0145 対 腕 M (MOC の軸処理を変えた壁) run_0150〜0152 の許容判定、出口較正、IC 依存 (run_0153 対 腕 M)。
入力: eval_wallfit_euler.py --series (tag v5) が書いた各 run の wallfit_series_v5.csv・QUASISTEADY_wallfit_v5.txt と
  _band_ab/wallfit_series_v5.json (7 本を同じ呼び出し・同じ X_E・X_F で評価した記録)、各 run の prepare_info.json・IC_MAP.json・
  IC_INSPECTION.json (腕 M)・CONVERGENCE_VERDICT_segment.txt・nozzle.h5 (壁の証拠)。
量と Δq (§6 V5): M 波・P 波・オーバーシュート・出口規格化オーバーシュート・|P 傾き| (いずれも r/r_w = 0.1) と |出口コア M − 6|。
  Δq = 0.001 / 0.010 / 0.003 / 0.003 / 0.03 %pt、1.8e-4。出口コア M は腕 B 対 腕 M の判定では記録用 (出口較正と IC 依存の出口 M 条件に使う)。
窓: 本段 step 6000〜18000 の 13 枚 (1000 ごと) に固定する。条件を満たさない run があっても、窓を動かしたり延長したりしない (追加の検証は別登録)。
§6 V5 の追加登録 (2026-10-07、commit 9a09c0a5、V5 の結果を誰も読む前。諮問 notes/reviews/2026-10-07-moc-v5-eval-interpretation-diagnose.md):
  準定常: 全 7 run の判定対象 6 量と出口コア M に check_quasisteady (classify_series、drift 0.05、osc 0.10) を窓の末尾 5 枚と
    全 13 枚の両方で当て、両方 STEADY を原則とする (加えて eval_wallfit_euler が書いた記録の準定常ファイルも STEADY であること)。
  例外は |出口コア M − 6| (exit_M_dev) だけ: 符号付き M_exit − 6 が窓の末尾 10 枚で §6 E4 の「絶対許容内」
    (throat_mono_judge.abs_tolerance_verdict、許容 1.8e-5) を満たし、かつ全 13 枚の最大絶対値と幅がともに 1.8e-5 以下なら比較に使う。
    元の VERDICT は残し、表示は「絶対許容内」(STEADY と書かない)。それ以外の非 STEADY・欠損・非有限は保留。
統計 (§6 E′ と同じ定義): run i の窓平均 m_i、時間変動の標準誤差 s_i = sd_i/√13 (自己相関は補正しない; 限界として記録)。
  腕の平均 M = mean(m_i)、腕の標準誤差 SE = max(sd(m_i)/√k, √(Σ s_i²)/k) (k = 腕の run 数; k = 1 のときは s_1)。
  D = M_M − M_B、SE_D = √(SE_B² + SE_M²)。
許容判定 (大きいほど悪い量の片側): 全量で D + 2·SE_D ≤ Δq → 「許容幅内」/ D − 2·SE_D ≥ Δq の量があれば「悪化」/ それ以外は保留。
  別欄で差の検出の有無 (|D| > 2·SE_D) を記録する (判定には使わない)。
出口較正: 出口コア M の D について 3 区分。|D| + 2·SE_D ≤ 1e-4 → 「据え置き」(Md_moc_offset を変えない。旧較正から変える必要が
  ないという意味で、出口 M = 6 の達成の証明ではない) / |D| − 2·SE_D > 1e-4 → 「較正をやり直す」/ それ以外 → 「保留 (判別不能)」
  (較正は変えない。V5 の較正判断も未完了)。
IC 依存 (§6 V5「2 本の評価量の差が許容判定の条件を満たすことを前提にする」): D_ic = q_ISEN − M_M、SE_ic = √(SE_M² + s_ISEN²)。
  IC の違いには良し悪しの向きが無いので両側で判定する (片側だと ISEN が「良い」向きの大差を通してしまう):
  全量で |D_ic| + 2·SE_ic ≤ Δq → 「IC 依存は許容幅内」/ |D_ic| − 2·SE_ic ≥ Δq の量があれば「IC 依存あり」/ それ以外は保留。
  6 量に加え、出口コア M そのもの (符号つき) の差に |D_ic,M| + 2·SE_ic ≤ 1e-4 を要求する (出口誤差 |M − 6| の差が 0 でも、
  5.9998 と 6.0002 は出口 M が 4e-4 違う)。前提 (下の前提検査) が 1 つでも不成立なら、IC の差は参考値として出し、
  「IC 依存あり/なし」と確定しない。前提なので、「IC 依存は許容幅内」でなければ総合は保留。
前提検査 (どれか不成立なら総合は「保留 (前提不成立)」。証拠が読めないこと・欠損も不成立。既定値で通さない):
  全 run: 時系列の記録 (wallfit_series_v5.json) にあり status ok、CSV・準定常ファイルの sha256 が記録と同じ (古い成果物でない)、
    記録の X_E・X_F の基準が腕 B の r1 で、各 run 自身の X_E・X_F との差 ≤ 1e-6 r_t (座標の整合の許容差であって、評価誤差の
    保証ではない)、記録に全 run 共通の評価座標 (X_E・X_F・WIN_T・WIN_O) がある (標本と選択範囲が全 run で同じ)、
    窓の 13 枚がそろい有限、判定区間 (CONVERGENCE_VERDICT_segment.txt) の残差判定が pass か既知の plateau、
    上の準定常 (例外は exit_M_dev だけ)、乾式確認 (DRY) の準備から作った run でない、wall_fit.mono_r2 = [0.0, 1.5]、
    壁節点が自腕の当てはめ後 spline に一致 (M4)、保存 spline が [0, 1.5] で単調 (形状ゲート S1 と同じ許容差)。
  記録だけ (判定に使わない): 各 run 自身の座標・±1e-6 r_t で動かしたときの x の標本数と、窓の最後の枚 (step 18000) の評価量の差
    (評価座標の感度。eval_wallfit_euler の quantities をそのまま使う)。
  腕 B: prepare_info に moc キーが無ければ legacy・fixed2 とみなす (2026-10-07 の MOC の実装より前に準備した run。例外として記録)。
    moc キーがあれば legacy・fixed2 であること。IC_MAP.json が VERDICT OK・mode index・IC は run_0114 の res_6000、3 本の写像後 sha256 が同じ。
  腕 M: moc が analytic・converge でゲート合格 (applicable かつ pass)。IC_MAP.json が VERDICT OK・mode index・IC は run_0114 の res_6000、
    IC_INSPECTION.json が VERDICT OK・全条件合格で、C1〜C5 の 5 キー (moc_v5_euler.inspect_ic の条件名) がすべて存在し個別に合格
    (欠けたら保留。幾何の対応の合格を IC 非依存の証明に代用しない)、上限が IC_MAP.json と prepare_info の上限と同じ、写像の移動の最大 ≤ 上限、
    3 本の写像後 sha256 が同じ。腕 B と見分けられる壁節点があり、そのすべてが自腕の当てはめに一致 (x の float32 の丸めは許す)。
  ISEN: moc のゲート合格、腕 M と同じ壁 (保存 spline が完全一致) と同じ格子 (IC 記録の格子ハッシュが腕 M と同じ)、
    IC が等エントロピー (IC_MAP.json の mode isentropic・VERDICT OK)。
  同じ実効設定: 7 本すべてで、本段の solverConfig.yaml・bcondConfig.yaml (YAML として)、forge のバイナリ (RUN_PROVENANCE.txt の
    forge_sha256)、起動時の実効値の行、段の並び (soft → main) と hard キー (stage_manifest.json) が腕 B の r1 と同じ。
SE は自己相関を補正しない実務の指標であり、保証された信頼限界・厳密な非劣化・IC 独立性を主張しない。
出力 JSON に評価器自身の sha256 と、plan の追加登録の commit (9a09c0a5) を記録する。
usage: python3 moc_v5_euler_eval.py [case_dir]   → _band_ab/moc_v5_euler_eval.json と標準出力の表
"""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PLAN = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5"
PLAN_REG_COMMIT = "9a09c0a5"                      # §6 V5 の追加登録 (2026-10-07、結果を読む前) の commit
ARMS = {"B": [f"run_0{n}_euler_wallfit_monoG1_r{k}" for n, k in ((143, 1), (144, 2), (145, 3))],
        "M": [f"run_0{n}_euler_wallfit_mocG1_r{k}" for n, k in ((150, 1), (151, 2), (152, 3))]}
ISEN = "run_0153_euler_icdep_mocG1_isen"
TAG = "v5"
GEOM_REF = ARMS["B"][0]
GEOM_TOL = 1e-6                                   # [r_t] 座標の整合の許容差 (各 run 自身の X_E・X_F と評価の基準との差の上限。評価誤差の保証ではない)
GEOM_TOL_LABEL = "座標の整合の許容差"
SERIES_H = 0.05                                   # [r_t] eval_wallfit_euler の時系列の評価刻み (quantities(F, 0.05))
WIN = (6000, 18000)
WIN_STEPS = list(range(WIN[0], WIN[1] + 1, 1000))   # 13 枚
DQ = {"M_wave_eta0.1": 0.001, "P_wave_eta0.1": 0.010, "overshoot_eta0.1": 0.003, "overshoot_exitnorm_eta0.1": 0.003,
      "P_slope_abs_eta0.1": 0.03, "exit_M_dev": 0.00018}
RECORD_ONLY = ("exit_core_M",)
QS_COLS = list(DQ) + list(RECORD_ONLY)            # 準定常を見る量 (判定対象 6 量と出口コア M)
QS_DRIFT, QS_OSC, QS_MIN_SNAPS = 0.05, 0.10, 4     # check_quasisteady の drift・osc (追加登録) と最小枚数 (CLI の既定)
QS_TAIL = 5                                       # 窓の末尾 5 枚 (もう 1 つは全 13 枚)
M_TARGET = 6.0
EXIT_EXC_QTY = "exit_M_dev"                       # 絶対許容内の例外を認める唯一の量
EXIT_EXC_FULL_TOL = 1.8e-5                        # 例外の追加条件: 全 13 枚の最大絶対値と幅 (末尾 10 枚の E4 は throat_mono_judge の Δq/10)
EXIT_CAL_TOL = 1e-4
IC_EXIT_M_TOL = 1e-4                              # IC 依存: 出口コア M そのもの (符号つき) の差の許容
# 腕 M の IC の検査 (moc_v5_euler.inspect_ic が IC_INSPECTION.json の conditions に書く条件名)。5 キーを個別に要求する
IC_INSPECTION_KEYS = ("C1_ic_index_map_checks", "C2_columns_x_unchanged", "C3_radial_scaling_of_wall_shift",
                      "C4_wall_shift_explained_by_design", "C5_mesh_quality")
EXIT_CAL_NOTE = "「据え置き」は旧較正から変える必要がないという意味で、出口 M = 6 の達成の証明ではない"
IC_RUN_NAME, IC_RES = "run_0114_euler_pin_G1_recal_ext6k", "res_6000.h5"
MONO_R2 = [0.0, 1.5]
MOC_B = {"axis_limit": "legacy", "corrector": "fixed2"}
MOC_M = {"axis_limit": "analytic", "corrector": "converge"}
WIN_LABEL, WORSE_LABEL, HOLD_LABEL = "許容幅内", "悪化", "保留"


def _sha(p: Path):
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def _json(p: Path) -> dict:
    return json.loads(p.read_text())


def load_window(case: Path, run: str) -> dict:
    """wallfit_series_v5.csv の窓 13 枚 → {量: ndarray}。欠けた枚・非有限・列の欠損は ValueError (0 で埋めない)。"""
    f = case / run / f"wallfit_series_{TAG}.csv"
    if not f.is_file():
        raise FileNotFoundError(f"{f.name} が無い (eval_wallfit_euler.py --series を先に回すこと)")
    rows = list(csv.DictReader(open(f)))
    win = [r for r in rows if WIN[0] <= int(float(r["step"])) <= WIN[1]]
    steps = [int(float(r["step"])) for r in win]
    if steps != WIN_STEPS:
        raise ValueError(f"窓 {WIN} の 1000 step ごとの 13 枚がそろっていない (実際 {steps})")
    out = {}
    for k in list(DQ) + list(RECORD_ONLY):
        if k not in win[0]:
            raise ValueError(f"列 {k} が無い")
        v = np.array([float(r[k]) if r[k] not in ("", None) else np.nan for r in win])
        if not np.all(np.isfinite(v)):
            raise ValueError(f"{k} に非有限値・空欄")
        out[k] = v
    return out


def quasisteady_window(vals: dict) -> dict:
    """窓 13 枚の check_quasisteady (classify_series、drift 0.05・osc 0.10) を、末尾 5 枚と全 13 枚の両方で (§6 V5 追加登録)。
    vals = load_window の戻り値。戻り値 {量: {"tail5", "full13" (VERDICT), "tail5_detail", "full13_detail"}}。"""
    sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda/tools"))
    from check_quasisteady import classify_series
    n = len(WIN_STEPS)
    out = {}
    for c in QS_COLS:
        v = [float(x) for x in vals[c]]
        if len(v) != n:
            raise ValueError(f"{c}: 窓の枚数 {len(v)} が {n} でない")
        # check_quasisteady は末尾 k = max(3, ceil(tail·n)) 枚 → (5 − 0.5)/13 でちょうど 5 枚、1.0 で 13 枚 (eval_wallfit_euler と同じ取り方)
        v5, d5, _ = classify_series(WIN_STEPS, v, (QS_TAIL - 0.5) / n, QS_DRIFT, QS_OSC, QS_MIN_SNAPS)
        v13, d13, _ = classify_series(WIN_STEPS, v, 1.0, QS_DRIFT, QS_OSC, QS_MIN_SNAPS)
        out[c] = {"tail5": v5, "full13": v13, "tail5_detail": d5, "full13_detail": d13}
    return out


def exit_exception(signed) -> dict:
    """exit_M_dev だけの例外 (§6 V5 追加登録): 符号付き M_exit − 6 (窓の 13 枚) が、末尾 10 枚で §6 E4 の「絶対許容内」
    (throat_mono_judge.abs_tolerance_verdict、許容 Δq/10 = 1.8e-5) を満たし、かつ全 13 枚の最大絶対値と幅がともに 1.8e-5 以下。
    戻り値の applies が True のときだけ比較に使う (表示は「絶対許容内」。元の VERDICT は別に残す)。"""
    from throat_mono_judge import ABS_TOL_N, ABS_TOL_WITHIN, DELTA_Q, abs_tolerance_verdict
    v = np.asarray(signed, dtype=float)
    rec = {"tail_n": ABS_TOL_N, "tail_tol": DELTA_Q["exit_M_dev"] / 10.0, "full_tol": EXIT_EXC_FULL_TOL, "applies": False}
    if v.shape != (len(WIN_STEPS),) or not np.all(np.isfinite(v)):
        rec["reason"] = f"窓の 13 枚がそろわない・非有限 (形 {v.shape})"
        return rec
    e4 = abs_tolerance_verdict(WIN_STEPS, [float(x) for x in v], DELTA_Q["exit_M_dev"] / 10.0, n=ABS_TOL_N)
    maxabs, width = float(np.abs(v).max()), float(np.ptp(v))
    why = []
    if e4.get("status") != ABS_TOL_WITHIN:
        why.append(f"末尾 {ABS_TOL_N} 枚の E4 が {e4.get('status')} ({e4.get('reason')})")
    if not maxabs <= EXIT_EXC_FULL_TOL:
        why.append(f"全 13 枚の最大絶対値 {maxabs:.3g} > {EXIT_EXC_FULL_TOL:g}")
    if not width <= EXIT_EXC_FULL_TOL:
        why.append(f"全 13 枚の幅 {width:.3g} > {EXIT_EXC_FULL_TOL:g}")
    rec.update(e4_tail=e4, full_max_abs=maxabs, full_width=width, applies=not why, reason=("; ".join(why) or None),
               label=(ABS_TOL_WITHIN if not why else None))
    return rec


def judge_window_quasisteady(qsw: dict, qfile: dict, exc: dict) -> tuple:
    """1 run の準定常の判定 (§6 V5 追加登録)。qsw = quasisteady_window、qfile = 記録の準定常ファイルの VERDICT ({量: VERDICT})、
    exc = exit_exception。窓の末尾 5 枚・全 13 枚・記録のファイルがすべて STEADY なら STEADY。exit_M_dev だけは例外が成り立てば
    「絶対許容内」(STEADY と表示しない)。それ以外は保留。戻り値 (不成立の理由 [文], {量: STEADY / 絶対許容内 / 保留})。"""
    from throat_mono_judge import ABS_TOL_WITHIN
    status, bad = {}, []
    for c in QS_COLS:
        q = qsw.get(c) or {}
        src = {"窓の末尾 5 枚": q.get("tail5", "UNKNOWN"), "窓の全 13 枚": q.get("full13", "UNKNOWN"),
               "記録の準定常ファイル": qfile.get(c, "UNKNOWN")}
        if all(v == "STEADY" for v in src.values()):
            status[c] = "STEADY"
        elif c == EXIT_EXC_QTY and exc.get("applies") is True:
            status[c] = ABS_TOL_WITHIN
        else:
            status[c] = HOLD_LABEL
            bad.append(", ".join(f"{c} {v} ({k})" for k, v in src.items() if v != "STEADY")
                       + (f" [絶対許容内の例外に当たらない: {exc.get('reason')}]" if c == EXIT_EXC_QTY else ""))
    why = (["準定常でない量 (STEADY でなく、exit_M_dev の絶対許容内の例外にも当たらない): " + "; ".join(bad)] if bad else [])
    return why, status


def window_conditions(case: Path, run: str, vals: dict) -> tuple:
    """1 run の窓の前提 (準定常と exit_M_dev の例外)。vals = load_window の戻り値 (13 枚そろい有限)。戻り値 (不成立の理由, 記録)。"""
    qp = case / run / f"QUASISTEADY_wallfit_{TAG}.txt"
    qfile = parse_quasisteady(qp.read_text() if qp.is_file() else "", QS_COLS)
    qsw = quasisteady_window(vals)
    signed = np.asarray(vals["exit_core_M"], dtype=float) - M_TARGET
    exc = exit_exception(signed)
    why, status = judge_window_quasisteady(qsw, qfile, exc)
    rec = {"status": status, "window": qsw, "file_tail5": qfile, "exit_M_dev_exception": exc,
           # 記録: CSV の exit_M_dev と |exit_core_M − 6| の差 (どちらも %.10g で書かれている)
           "exit_M_dev_vs_signed_maxabs": float(np.abs(np.asarray(vals["exit_M_dev"], dtype=float) - np.abs(signed)).max())}
    return why, rec


def parse_quasisteady(text: str, cols) -> dict:
    """check_quasisteady --series-csv の出力 → {列: VERDICT} (行が無ければ "UNKNOWN")。eval_wallfit_euler._verdict_lines と同じ読み方。"""
    out = {}
    for c in cols:
        m = re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", text, re.M)
        out[c] = m.group(1) if m else "UNKNOWN"
    return out


def arm_stats(per_run: list) -> dict:
    """§6 E′ の腕の統計。per_run = [窓の値の配列 (13 枚)]。"""
    m = np.array([v.mean() for v in per_run])
    s = np.array([v.std(ddof=1) / np.sqrt(len(v)) for v in per_run])
    k = len(per_run)
    se_between = float(m.std(ddof=1) / np.sqrt(k)) if k >= 2 else None
    se_time = float(np.sqrt((s ** 2).sum())) / k
    se = max(se_between, se_time) if se_between is not None else se_time
    return {"mean": float(m.mean()), "run_means": [float(x) for x in m], "temporal_se": [float(x) for x in s],
            "se_between_runs": se_between, "se_temporal": se_time, "se": float(se), "k": k}


def judge_one_sided(D: float, SE: float, dq: float) -> str:
    if D + 2.0 * SE <= dq:
        return WIN_LABEL
    if D - 2.0 * SE >= dq:
        return WORSE_LABEL
    return HOLD_LABEL


def judge_two_sided(D: float, SE: float, dq: float) -> str:
    if abs(D) + 2.0 * SE <= dq:
        return "IC 依存は許容幅内"
    if abs(D) - 2.0 * SE >= dq:
        return "IC 依存あり"
    return HOLD_LABEL


def judge_exit_calibration(D: float, SE: float, tol: float = EXIT_CAL_TOL) -> str:
    """出口較正の 3 区分 (据え置き / やり直す / 保留)。「据え置き」は出口 M = 6 の達成の証明ではない (EXIT_CAL_NOTE)。"""
    if abs(D) + 2.0 * SE <= tol:
        return f"据え置き (Md_moc_offset を変えない。{EXIT_CAL_NOTE})"
    if abs(D) - 2.0 * SE > tol:
        return "較正をやり直す"
    return "保留 (判別不能; 較正は変えない。V5 の較正判断は未完了)"


# --- 壁の照合 (nozzle.h5 を読む。試験では差し替える) -------------------------------------------------------------
def _ulp(a):
    return np.spacing(np.abs(np.asarray(a, dtype=np.float64)).astype(np.float32)).astype(np.float64)


def _x_last_design(rd: Path) -> float:
    return float(np.loadtxt(rd / "wall_design.csv", delimiter=",", skiprows=1)[-1, 0])


def wall_vs_other(rd: Path, od: Path) -> dict:
    """自腕と他腕の当てはめ後 spline を自腕の壁節点の x で評価し、見分けられる節点 (2 本の差 > 2·float32 の許容) が自腕に一致するか。
    throat_mono_ab.wall_evidence の他腕照合と同じ方法で、壁節点の x の差は float32 の丸め + 設計の出口 x の差まで許す
    (腕 M と腕 B は X_F が 1e-8 r_t 違い、x の配置が 1 ulp 動く列がある)。"""
    import throat_mono_ab as TM
    info, io = _json(rd / "prepare_info.json"), _json(od / "prepare_info.json")
    S = float(info["scale_m"])
    spl, spo = TM.wall_spline(info), TM.wall_spline(io)
    if spl is None or spo is None:
        return {"status": "missing", "reason": "wall_fit.spline が無い"}
    xw, rw, _ = TM.wall_nodes(rd, info)
    xo, _, _ = TM.wall_nodes(od, io)
    if xw.shape != xo.shape:
        return {"status": "shape_differs"}
    tolx = _ulp(np.maximum(np.abs(xw), np.abs(xo))) + abs(_x_last_design(rd) - _x_last_design(od))
    if not np.all(np.abs(xw - xo) <= tolx):
        return {"status": "x_differs", "max_dx_m": float(np.abs(xw - xo).max())}
    m = (xw / S >= max(spl.t[0], spo.t[0])) & (xw / S <= min(spl.t[-1], spo.t[-1]))
    xs = xw[m] / S
    r_own, r_oth = S * spl(xs), S * spo(xs)
    tol = _ulp(rw[m]) + np.abs(spl(xs, 1)) * _ulp(xw[m])
    disc = np.abs(r_own - r_oth) > 2.0 * tol
    own_better = np.abs(rw[m] - r_own) < np.abs(rw[m] - r_oth)
    return {"status": "ok", "other": od.name, "n_discriminable": int(disc.sum()),
            "n_discriminable_matching_own": int((disc & own_better).sum()),
            "analytic_dr_max_m": float(np.abs(r_own - r_oth).max()), "x_analytic_dr_max_rt": float(xs[np.argmax(np.abs(r_own - r_oth))]),
            "x_range_discriminable_rt": ([float(xs[disc].min()), float(xs[disc].max())] if disc.any() else None)}


def default_wall_check(case: Path, run: str, other: str | None) -> dict:
    """{"own": throat_mono_ab.wall_evidence (M4), "vs_other": wall_vs_other (other があるとき)}。"""
    import throat_mono_ab as TM
    out = {"own": TM.wall_evidence(case / run)}
    if other is not None:
        out["vs_other"] = wall_vs_other(case / run, case / other)
    return out


def spline_shape_failures(info: dict) -> list:
    """保存 spline が [0, 1.5] で単調 (r‴ 最大 ≤ 1e-6、r″ の最大増加 ≤ 1e-7; throat_mono_practical_eval と同じ)。"""
    from scipy.interpolate import BSpline
    sys.path.insert(0, str(HERE.parents[1] / "design"))
    from forge_design.geometry.wall_axismach import r3_piecewise_exact
    sp = (info.get("wall_fit") or {}).get("spline")
    if not sp:
        return ["prepare_info の wall_fit.spline が無い (形状を検査できない)"]
    spl = BSpline(np.asarray(sp["t"], dtype=float), np.asarray(sp["c"], dtype=float), int(sp["k"]))
    sh = r3_piecewise_exact(spl, None, 0.0, 1.5)
    if not (sh["r3_max"] <= 1e-6 and sh["r2_max_increase"] <= 1e-7):
        return [f"保存 spline が [0, 1.5] で単調でない (r‴ 最大 {sh['r3_max']:.3g}、r″ の最大増加 {sh['r2_max_increase']:.3g})"]
    return []


# --- 評価座標 (§6 V5 追加登録: 座標の整合の許容差、標本数・選択範囲、感度の記録) -----------------------------------------
def geometry_failures(geom) -> list:
    """時系列の記録の共通の評価座標 (X_E・X_F・WIN_T・WIN_O) がそろい有限で、区間が正しい向きか (全 run を同じ標本・選択範囲で評価した証拠)。"""
    if not isinstance(geom, dict):
        return [f"時系列の記録に共通の評価座標 (geometry) が無い ({geom!r})"]
    bad = []
    for k in ("X_E", "X_F"):
        v = geom.get(k)
        if not (isinstance(v, (int, float)) and not isinstance(v, bool) and np.isfinite(v)):
            bad.append(f"評価座標 {k} = {v!r} が有限の数でない")
    for k in ("WIN_T", "WIN_O"):
        v = geom.get(k)
        if not (isinstance(v, (list, tuple)) and len(v) == 2 and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                                                                     and np.isfinite(x) for x in v) and v[0] < v[1]):
            bad.append(f"評価の選択範囲 {k} = {v!r} が有限の [下限, 上限] でない")
    if not bad and not geom["X_E"] < geom["X_F"]:
        bad.append(f"評価座標 X_E {geom['X_E']!r} < X_F {geom['X_F']!r} でない")
    return bad


def _wallfit_defs() -> dict:
    """eval_wallfit_euler.py の grid・pspline_fixed・quantities と定数 P_REF・MD を、ソースから定義だけ取り出す (あのモジュールは
    import すると旧 A/B の評価が走るので import しない)。評価座標 (X_E・X_F・WIN_T・WIN_O) は _set_geometry で入れる。"""
    import ast
    from scipy.interpolate import BSpline
    p = HERE / "eval_wallfit_euler.py"
    src = p.read_text()
    want = {"grid", "pspline_fixed", "quantities"}
    body = []
    for n in ast.parse(src).body:
        if isinstance(n, ast.FunctionDef) and n.name in want:
            body.append(n)
        elif isinstance(n, ast.Assign) and {t.id for t in ast.walk(n.targets[0]) if isinstance(t, ast.Name)} & {"P_REF", "MD"}:
            body.append(n)
    got = {n.name for n in body if isinstance(n, ast.FunctionDef)}
    if got != want:
        raise ValueError(f"eval_wallfit_euler.py に {sorted(want - got)} が無い")
    ns = {"np": np, "BSpline": BSpline}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(p), "exec"), ns)   # noqa: S102 — リポジトリ内の評価器の定義だけ
    ns["_source_sha256"] = hashlib.sha256(src.encode()).hexdigest()
    return ns


def _set_geometry(ns: dict, geom: dict, dE: float = 0.0, dF: float = 0.0) -> dict:
    """評価座標を (X_E + dE, X_F + dF) にする。窓は記録の WIN_T・WIN_O を、下端は dE、上端は dF だけ動かす
    (run_series の取り方 WIN_T = (X_E + 2, X_F − 1)・WIN_O = (X_E − 15, X_F) と同じ依存。0 のずれでは記録の値そのもの)。"""
    xe, xf = float(geom["X_E"]) + dE, float(geom["X_F"]) + dF
    wt = (float(geom["WIN_T"][0]) + dE, float(geom["WIN_T"][1]) + dF)
    wo = (float(geom["WIN_O"][0]) + dE, float(geom["WIN_O"][1]) + dF)
    ns.update(X_E=xe, X_F=xf, WIN_T=wt, WIN_O=wo)
    return {"X_E": xe, "X_F": xf, "WIN_T": list(wt), "WIN_O": list(wo)}


def sample_sets(ns: dict, geom: dict, dE: float = 0.0, dF: float = 0.0) -> dict:
    """記録用: 評価座標 (X_E + dE, X_F + dF) での x の標本数と選択範囲 (eval_wallfit_euler.grid(0.05) の標本、窓 WIN_T・WIN_O、
    当てはめ区間 [X_E, X_F])。"""
    g = _set_geometry(ns, geom, dE, dF)
    xq = ns["grid"](SERIES_H)

    def n_in(lo, hi):
        return int(((xq >= lo) & (xq <= hi)).sum())
    return {"n_grid": int(len(xq)), "n_WIN_T": n_in(*g["WIN_T"]), "n_WIN_O": n_in(*g["WIN_O"]), "n_fit": n_in(g["X_E"], g["X_F"])}


def default_coord_sensitivity(case: Path, run: str, geom: dict, shifts: dict) -> dict:
    """記録用 (判定に使わない): 窓の最後の枚 (step 18000) で、評価座標を shifts {名前: (dX_E, dX_F)} だけ動かしたときの
    評価量の差 (eval_wallfit_euler.quantities をそのまま使う)。res_18000.h5 を読む (試験では差し替える)。"""
    sys.path.insert(0, str(HERE.parents[1] / "design"))
    from forge_design.report.nozzle_report import eta_line, load_field
    ns = _wallfit_defs()
    ns["eta_line"] = eta_line
    F = load_field(case / run, f"res_{WIN[1]}.h5")

    def q_at(dE, dF):
        _set_geometry(ns, geom, dE, dF)
        q, _ = ns["quantities"](F, SERIES_H)
        return {c: float(q[c]) for c in QS_COLS}
    base = q_at(0.0, 0.0)
    diff = {}
    for name, (dE, dF) in shifts.items():
        q = q_at(dE, dF)
        diff[name] = {c: q[c] - base[c] for c in QS_COLS}
    return {"step": WIN[1], "base": base, "diff": diff, "eval_wallfit_euler_sha256": ns["_source_sha256"]}


def coordinate_record(case: Path, srec: dict, runs: list, S: dict, coord_sensitivity) -> dict:
    """記録用 (判定に使わない): 評価座標の感度。各 run 自身の座標 (dX_E, dX_F) と ±1e-6 r_t (座標の整合の許容差) で動かしたときの
    x の標本数と、窓の最後の枚の評価量の差 (量ごとの最大絶対値と Δq に対する比)。"""
    geom = srec.get("geometry")
    out = {"tolerance_label": GEOM_TOL_LABEL, "tolerance_r_t": GEOM_TOL, "used_in_judgement": False,
           "note": f"{GEOM_TOL:g} r_t は{GEOM_TOL_LABEL}であって評価誤差の保証ではない", "geometry": geom, "samples": {}, "sensitivity": {}}
    if geometry_failures(geom):
        out["error"] = "共通の評価座標が無い・不正 (前提検査で保留)"
        return out
    tol_shifts = {"+tol_X_E": (GEOM_TOL, 0.0), "-tol_X_E": (-GEOM_TOL, 0.0), "+tol_X_F": (0.0, GEOM_TOL), "-tol_X_F": (0.0, -GEOM_TOL)}
    try:
        ns = _wallfit_defs()
        out["eval_wallfit_euler_sha256"] = ns["_source_sha256"]
        out["series_evaluator_sha256"] = srec.get("evaluator_sha256")
        common = sample_sets(ns, geom)
        out["samples"]["common"] = common
        out["samples"]["tol_shifts"] = {k: sample_sets(ns, geom, *v) for k, v in tol_shifts.items()}
        out["samples"]["tol_shifts_same_as_common"] = all(v == common for v in out["samples"]["tol_shifts"].values())
    except Exception as e:  # noqa: BLE001 — 記録用 (判定に使わない)
        out["samples"]["error"] = f"{type(e).__name__}: {e}"
        ns, common = None, None
    for run in runs:
        r = ((srec.get("runs") or {}).get(run)) or {}
        dE, dF = r.get("dX_E"), r.get("dX_F")
        own_ok = all(isinstance(v, (int, float)) and not isinstance(v, bool) and np.isfinite(v) for v in (dE, dF))
        shifts = dict(tol_shifts)
        if own_ok:
            shifts = {"own": (float(dE), float(dF)), **tol_shifts}
            if ns is not None:
                try:
                    own = sample_sets(ns, geom, float(dE), float(dF))
                    out["samples"].setdefault("own", {})[run] = {"counts": own, "same_as_common": own == common}
                except Exception as e:  # noqa: BLE001
                    out["samples"].setdefault("own", {})[run] = {"error": f"{type(e).__name__}: {e}"}
        if run not in S:
            out["sensitivity"][run] = {"error": "窓の評価量が無い (欠損)"}
            continue
        try:
            sv = coord_sensitivity(case, run, geom, shifts)
            base = sv.get("base") or {}
            sv["base_minus_csv_at_step"] = {c: (float(base[c]) - float(S[run][c][-1]) if c in base else None) for c in QS_COLS}
            mx = {c: max((abs(float(d[c])) for d in (sv.get("diff") or {}).values() if c in d), default=None) for c in QS_COLS}
            sv["max_abs_diff"] = mx
            sv["max_abs_diff_over_dq"] = {c: (mx[c] / DQ[c] if mx[c] is not None else None) for c in DQ}
            out["sensitivity"][run] = sv
        except Exception as e:  # noqa: BLE001 — 記録用 (判定に使わない)
            out["sensitivity"][run] = {"error": f"{type(e).__name__}: {e}"}
    return out


# --- 前提検査 -----------------------------------------------------------------------------------------------------
def _ic_src_ok(rec: dict) -> bool:
    p = Path(str(rec.get("src_res")))
    return p.name == IC_RES and p.parent.name == IC_RUN_NAME


def run_preconditions(case: Path, run: str, role: str, srec: dict, wall_check) -> tuple:
    """1 run の前提。戻り値 (不成立の理由 [文], 例外の記録 [文], 記録 dict, prepare_info, IC_MAP.json)。role: B / M / ISEN。"""
    from throat_mono_judge import SEGMENT_VERDICT_FILE, mono_r2_matches, parse_segment_verdict
    rd = case / run
    why, exc, meta = [], [], {}
    # 時系列の記録
    r = (srec.get("runs") or {}).get(run)
    if r is None:
        why.append(f"時系列の記録 (wallfit_series_{TAG}.json) に無い")
    else:
        if r.get("status") != "ok":
            why.append(f"時系列の評価が {r.get('status')!r} ({r.get('reason')})")
        for key, name in (("csv_sha256", f"wallfit_series_{TAG}.csv"), ("quasisteady_sha256", f"QUASISTEADY_wallfit_{TAG}.txt")):
            got = _sha(rd / name)
            if not r.get(key) or got != r.get(key):
                why.append(f"{name} が時系列の記録と違う (sha256 記録 {str(r.get(key))[:12]} / 実際 {str(got)[:12]}; 古い成果物か別の呼び出し)")
        for k in ("dX_E", "dX_F"):
            v = r.get(k)
            if not (isinstance(v, (int, float)) and not isinstance(v, bool) and np.isfinite(v) and abs(v) <= GEOM_TOL):
                why.append(f"自身の {k[1:]} と評価の基準の差 {v!r} が{GEOM_TOL_LABEL} {GEOM_TOL:g} r_t を超える・欠損")
        meta["series"] = {k: r.get(k) for k in ("dX_E", "dX_F", "n_snaps", "last_step", "status")}
    # 判定区間の残差
    p = rd / SEGMENT_VERDICT_FILE
    seg = parse_segment_verdict(p.read_text() if p.is_file() else None)
    meta["segment"] = seg
    if seg["status"] not in ("pass", "plateau"):
        why.append(f"判定区間の残差判定が {seg['status']} ({seg.get('reason') or seg.get('line')})")
    # 準定常 (窓の末尾 5 枚・全 13 枚・記録の準定常ファイル) は窓の値が要るので evaluate の window_conditions で見る
    # prepare_info・壁
    info = _json(rd / "prepare_info.json")
    if info.get("DRY") or info.get("DRY_NO_IC"):
        why.append("乾式確認 (DRY) の準備から作った run")
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf or not mono_r2_matches(wf.get("mono_r2"), MONO_R2):
        why.append(f"wall_fit.mono_r2 {wf.get('mono_r2', '(記録なし)')!r} が {MONO_R2} と一致しない")
    why += spline_shape_failures(info)
    moc = info.get("moc")
    meta["moc"] = ({k: (moc or {}).get(k) for k in ("axis_limit", "corrector")} | {"gate": (moc or {}).get("gate")}) if moc else None
    if role == "B":
        if moc is None:
            exc.append(f"{run}: prepare_info に moc キーが無い — MOC の実装 (2026-10-07) より前に準備した run なので legacy・fixed2 とみなした")
        elif {k: moc.get(k) for k in MOC_B} != MOC_B:
            why.append(f"腕 B の moc が legacy・fixed2 でない ({moc.get('axis_limit')}, {moc.get('corrector')})")
    else:
        gate = (moc or {}).get("gate") or {}
        if moc is None or {k: moc.get(k) for k in MOC_M} != MOC_M:
            why.append(f"moc が analytic・converge でない ({None if moc is None else (moc.get('axis_limit'), moc.get('corrector'))})")
        if gate.get("applicable") is not True or gate.get("pass") is not True:
            why.append(f"MOC のゲートが合格でない (applicable {gate.get('applicable')!r}, pass {gate.get('pass')!r})")
    # IC
    rec = _json(rd / "IC_MAP.json")
    meta["ic"] = {k: rec.get(k) for k in ("VERDICT", "mode", "src_res", "dst_sha256_after", "dst_mesh_digest")}
    if role in ("B", "M"):
        if rec.get("VERDICT") != "OK" or rec.get("mode") != "index":
            why.append(f"IC_MAP.json の VERDICT {rec.get('VERDICT')!r}・mode {rec.get('mode')!r} (OK・index であること)")
        if not _ic_src_ok(rec):
            why.append(f"IC が {IC_RUN_NAME}/{IC_RES} でない ({rec.get('src_res')})")
    if role == "M":
        lim = (((rec.get("checks") or {}).get("displacement") or {}).get("detail") or {}).get("limit_m")
        dmax = (((rec.get("checks") or {}).get("displacement") or {}).get("detail") or {}).get("max_m")
        ins = _json(rd / "IC_INSPECTION.json")
        lim_info = ((info.get("ic") or {}).get("inspection") or {}).get("limit_m")
        conds = ins.get("conditions") if isinstance(ins.get("conditions"), dict) else {}
        meta["ic"].update(limit_m=lim, max_m=dmax, inspection_VERDICT=ins.get("VERDICT"), inspection_limit_m=ins.get("limit_m"),
                          inspection_conditions={k: (conds[k].get("ok") if isinstance(conds.get(k), dict) else None)
                                                 for k in IC_INSPECTION_KEYS})
        if (ins.get("VERDICT") != "OK" or not conds
                or not all(isinstance(v, dict) and v.get("ok") is True for v in conds.values())):
            why.append(f"IC の検査 (IC_INSPECTION.json) が OK・全条件合格でない ({ins.get('VERDICT')!r})")
        # C1〜C5 の 5 キーを個別に要求する (存在する条件が全部合格でも、欠けた条件は合格にしない)
        miss = [k for k in IC_INSPECTION_KEYS if k not in conds]
        nok = [k for k in IC_INSPECTION_KEYS if k in conds and not (isinstance(conds[k], dict) and conds[k].get("ok") is True)]
        if miss:
            why.append(f"IC の検査 (IC_INSPECTION.json) に条件 {', '.join(miss)} が無い (C1〜C5 の 5 キーを個別に要求)")
        if nok:
            why.append(f"IC の検査の条件 {', '.join(nok)} が合格 (ok: true) でない")
        if not (isinstance(lim, (int, float)) and lim == ins.get("limit_m") and lim == lim_info):
            why.append(f"IC 写像の上限 {lim!r} が IC の検査の上限 {ins.get('limit_m')!r}・prepare_info の上限 {lim_info!r} と一致しない")
        if not (isinstance(dmax, (int, float)) and isinstance(lim, (int, float)) and dmax <= lim):
            why.append(f"写像の移動の最大 {dmax!r} が上限 {lim!r} 以下でない")
    if role == "ISEN":
        if rec.get("VERDICT") != "OK" or rec.get("mode") != "isentropic":
            why.append(f"IC_MAP.json の VERDICT {rec.get('VERDICT')!r}・mode {rec.get('mode')!r} (OK・isentropic であること)")
    # 壁の証拠
    other = None if role == "ISEN" else (ARMS["B"][0] if role == "M" else None)
    try:
        wc = wall_check(case, run, other)
    except Exception as e:  # noqa: BLE001 — 証拠が取れないことは不成立
        wc = {"own": {"status": "error", "reason": f"{type(e).__name__}: {e}"}}
    meta["wall"] = wc
    if (wc.get("own") or {}).get("status") != "consistent":
        why.append(f"壁節点が自腕の当てはめ後 spline に一致しない ({(wc.get('own') or {}).get('status')}: {(wc.get('own') or {}).get('reason', '')})")
    if other is not None:
        vo = wc.get("vs_other") or {}
        nd, nown = vo.get("n_discriminable"), vo.get("n_discriminable_matching_own")
        if vo.get("status") != "ok":
            why.append(f"腕 B との壁の照合が {vo.get('status')!r}")
        elif not (isinstance(nd, int) and nd > 0):
            why.append("腕 B と見分けられる壁節点が無い (同じ壁の可能性)")
        elif nown != nd:
            why.append(f"見分けられる壁節点 {nd} のうち自腕の当てはめに一致するのは {nown!r} (壁の取り違えの可能性)")
    return why, exc, meta, info, rec


def preconditions(case: Path, srec: dict, wall_check) -> tuple:
    """全 run の前提と腕をまたぐ照合。戻り値 (不成立 {run: [理由]}, 例外 [文], 記録 {run: meta})。"""
    bad, exc, metas, infos, recs = {}, [], {}, {}, {}
    common = []
    if srec.get("tag") != TAG:
        common.append(f"時系列の記録の tag {srec.get('tag')!r} が {TAG!r} でない")
    if srec.get("geom_ref") != GEOM_REF:
        common.append(f"時系列の記録の X_E・X_F の基準 {srec.get('geom_ref')!r} が {GEOM_REF} (腕 B の r1) でない")
    # 標本・選択範囲が全 run で同じ: 記録に全 run 共通の評価座標があり、それが腕 B の r1 自身の座標 (差 0) であること
    common += geometry_failures(srec.get("geometry"))
    gr = ((srec.get("runs") or {}).get(GEOM_REF)) or {}
    if gr and not (gr.get("dX_E") == 0.0 and gr.get("dX_F") == 0.0):
        common.append(f"評価座標が基準 {GEOM_REF} 自身の座標でない (差 {gr.get('dX_E')!r}, {gr.get('dX_F')!r})")
    for role, runs in (("B", ARMS["B"]), ("M", ARMS["M"]), ("ISEN", [ISEN])):
        for run in runs:
            if not (case / run).is_dir():
                bad[run] = ["run dir が無い (欠損)"]
                continue
            try:
                why, ex, meta, info, rec = run_preconditions(case, run, role, srec, wall_check)
                infos[run], recs[run] = info, rec
            except Exception as e:  # noqa: BLE001 — 証拠が読めないことは不成立
                why, ex, meta = [f"証拠を読めない: {type(e).__name__}: {e}"], [], {}
            metas[run] = meta
            exc += ex
            if why:
                bad[run] = why
    # 腕の中の 3 本が同じ準備から作られたこと (写像後 sha256)
    for arm in ("B", "M"):
        shas = [recs.get(r, {}).get("dst_sha256_after") for r in ARMS[arm]]
        if any(s is None for s in shas) or len(set(shas)) != 1:
            common.append(f"腕 {arm} の 3 本の IC 写像後 sha256 がそろわない ({[str(s)[:12] for s in shas]})")
    # ISEN が腕 M と同じ壁・同じ格子
    sp_m = [((infos.get(r) or {}).get("wall_fit") or {}).get("spline") for r in ARMS["M"]]
    sp_i = ((infos.get(ISEN) or {}).get("wall_fit") or {}).get("spline")
    if sp_i is None or any(s is None or s != sp_i for s in sp_m):
        common.append("ISEN の保存 spline が腕 M の 3 本と完全一致しない (同じ壁でない・欠損)")
    md_m = {recs.get(r, {}).get("dst_mesh_digest") for r in ARMS["M"]}
    md_i = recs.get(ISEN, {}).get("dst_mesh_digest")
    if md_i is None or md_m != {md_i}:
        common.append("ISEN の格子ハッシュが腕 M と同じでない (同じ格子でない・欠損)")
    # 同じ実効設定 (§6 V5「同じ評価器・同じ判定区間・同じ実効設定」): 腕 B の r1 を基準に、本段の solverConfig・bcondConfig
    # (YAML として)、forge のバイナリ (sha256)、起動時の実効値の行、段の並びと hard キー
    eff = {}
    for run in ARMS["B"] + ARMS["M"] + [ISEN]:
        if (case / run).is_dir():
            try:
                eff[run] = effective_settings(case / run)
            except Exception as e:  # noqa: BLE001 — 読めないことは不成立
                eff[run] = {"error": f"{type(e).__name__}: {e}"}
    ref = eff.get(GEOM_REF)
    for run, e in eff.items():
        if "error" in e:
            common.append(f"{run}: 実効設定を読めない ({e['error']})")
            continue
        if ref is None or "error" in ref:
            continue
        diff = [k for k in e if e[k] != ref.get(k)]
        if diff:
            common.append(f"{run}: 実効設定が {GEOM_REF} と違う ({', '.join(diff)})")
    if ref is None:
        common.append(f"実効設定の基準 {GEOM_REF} が無い")
    for run in metas:
        if run in eff and "error" not in eff[run]:
            metas[run]["effective"] = {"forge_sha256": eff[run]["forge_sha256"], "stages": [s[0] for s in eff[run]["stages"]]}
    if common:
        bad["共通"] = common
    return bad, exc, metas


def effective_settings(rd: Path) -> dict:
    """本段の実効設定: solverConfig.yaml・bcondConfig.yaml (YAML として読んだもの。本段の後は本段の設定が残る)、
    RUN_PROVENANCE.txt の forge_sha256 と起動時の実効値の行、stage_manifest.json の段の並びと hard キー。欠損は例外 (既定値で埋めない)。"""
    import yaml

    def _y(p):
        t = p.read_text()
        try:
            return yaml.safe_load(t)
        except yaml.YAMLError:
            return " ".join(t.split())                 # YAML として読めなければ空白を正規化した全文で比べる
    prov = (rd / "RUN_PROVENANCE.txt").read_text().splitlines()
    sha = [l.split(":", 1)[1].strip() for l in prov if l.startswith("forge_sha256")]
    if len(sha) != 1 or not re.fullmatch(r"[0-9a-f]{64}", sha[0]):
        raise ValueError(f"RUN_PROVENANCE.txt の forge_sha256 を読めない ({sha})")
    man = _json(rd / "stage_manifest.json")
    stages = [(s["tag"], s["key"]) for s in man["stages"]]
    if [s[0] for s in stages] != ["soft", "main"]:
        raise ValueError(f"段の並びが soft → main でない ({[s[0] for s in stages]})")
    return {"solverConfig": _y(rd / "solverConfig.yaml"), "bcondConfig": _y(rd / "bcondConfig.yaml"), "forge_sha256": sha[0],
            "effective_lines": sorted(l.strip() for l in prov if l.startswith("'") and " effective" in l),
            "stages": stages}


IC_WIN_LABEL, IC_DEP_LABEL = "IC 依存は許容幅内", "IC 依存あり"
IC_REF_LABEL = "参考値 (前提未達; IC 依存あり/なしを確定しない)"


def _ic_overall(ic_verdicts: dict) -> str:
    if all(v == IC_WIN_LABEL for v in ic_verdicts.values()):
        return IC_WIN_LABEL
    if any(v == IC_DEP_LABEL for v in ic_verdicts.values()):
        return f"{IC_DEP_LABEL}: " + ", ".join(k for k, v in ic_verdicts.items() if v == IC_DEP_LABEL)
    return "保留 (判別不能: " + ", ".join(k for k, v in ic_verdicts.items() if v != IC_WIN_LABEL) + ")"


def evaluate(case: Path, wall_check=None, coord_sensitivity=None) -> dict:
    """§6 V5 の判定。wall_check・coord_sensitivity は nozzle.h5・res_*.h5 を読む部分 (試験では差し替える)。"""
    case = Path(case).resolve()
    sys.path.insert(0, str(HERE))
    wall_check = wall_check or default_wall_check
    coord_sensitivity = coord_sensitivity or default_coord_sensitivity
    sp = case / f"_band_ab/wallfit_series_{TAG}.json"
    try:
        srec = _json(sp)
    except (OSError, ValueError) as e:
        srec = {}
        srec_err = f"時系列の記録 {sp.name} を読めない ({type(e).__name__})"
    else:
        srec_err = None
    pre_bad, pre_exc, metas = preconditions(case, srec, wall_check)
    if srec_err:
        pre_bad.setdefault("共通", []).insert(0, srec_err)
    runs_all = ARMS["B"] + ARMS["M"] + [ISEN]
    # 評価量 (読めない run は欠損。窓は 6000〜18000 に固定し、動かさない)
    S, missing = {}, {}
    for run in runs_all:
        try:
            S[run] = load_window(case, run)
        except Exception as e:  # noqa: BLE001 — 欠損は欠損として記録 (0 で埋めない)
            missing[run] = f"{type(e).__name__}: {e}"
    # 窓の準定常と exit_M_dev の例外 (§6 V5 追加登録)。全 7 run に同じ条件 (腕 B・ISEN も免除しない)
    qs_status = {}
    for run in runs_all:
        if run not in S:
            continue
        try:
            why, wrec = window_conditions(case, run, S[run])
        except Exception as e:  # noqa: BLE001 — 判定できないことは不成立
            why, wrec = [f"準定常を判定できない: {type(e).__name__}: {e}"], {"error": f"{type(e).__name__}: {e}"}
        metas.setdefault(run, {})["window_quasisteady"] = wrec
        qs_status[run] = wrec.get("status")
        if why:
            pre_bad.setdefault(run, []).extend(why)
    complete = not missing
    pre_ok = not pre_bad
    rows, verdicts, ic_rows, ic_verdicts, exit_cal = [], {}, [], {}, None
    if complete:
        for k in QS_COLS:
            st = {arm: arm_stats([S[r][k] for r in runs]) for arm, runs in ARMS.items()}
            si = arm_stats([S[ISEN][k]])
            D = st["M"]["mean"] - st["B"]["mean"]
            SE = float(np.hypot(st["B"]["se"], st["M"]["se"]))
            Dic = si["mean"] - st["M"]["mean"]
            SEic = float(np.hypot(st["M"]["se"], si["se"]))
            row = {"qty": k, "B": st["B"], "M": st["M"], "D_M_minus_B": float(D), "SE_D": SE, "detected": bool(abs(D) > 2 * SE)}
            irow = {"qty": k, "ISEN": si, "D_ISEN_minus_M": float(Dic), "SE_ic": SEic}
            if k in DQ:
                dq = DQ[k]
                row.update(dq=dq, verdict=judge_one_sided(D, SE, dq), D_plus_2SE_over_dq=float((D + 2 * SE) / dq))
                irow.update(dq=dq, verdict=judge_two_sided(Dic, SEic, dq), absD_plus_2SE_over_dq=float((abs(Dic) + 2 * SEic) / dq))
                verdicts[k] = row["verdict"]
            else:
                # 腕 B 対 腕 M では記録のみ (出口較正に使う)。IC 依存では出口コア M そのもの (符号つき) の差を 1e-4 で判定する
                row["verdict"] = "記録のみ"
                irow.update(dq=IC_EXIT_M_TOL, signed=True, verdict=judge_two_sided(Dic, SEic, IC_EXIT_M_TOL),
                            absD_plus_2SE_over_dq=float((abs(Dic) + 2 * SEic) / IC_EXIT_M_TOL))
            ic_verdicts[k] = irow["verdict"]
            rows.append(row)
            ic_rows.append(irow)
        ex = next(r for r in rows if r["qty"] == "exit_core_M")
        exit_cal = {"qty": "exit_core_M", "D_M_minus_B": ex["D_M_minus_B"], "SE_D": ex["SE_D"], "tol": EXIT_CAL_TOL,
                    "absD_plus_2SE": abs(ex["D_M_minus_B"]) + 2 * ex["SE_D"], "absD_minus_2SE": abs(ex["D_M_minus_B"]) - 2 * ex["SE_D"],
                    "verdict": judge_exit_calibration(ex["D_M_minus_B"], ex["SE_D"]), "note": EXIT_CAL_NOTE}
    # IC 依存の総合 (前提が 1 つでも不成立なら、差は参考値で「IC 依存あり/なし」と確定しない)
    ic_overall_if_held = _ic_overall(ic_verdicts) if complete else None
    if not complete:
        ic_overall = "保留 (欠損)"
    elif not pre_ok:
        ic_overall = "保留 (前提未達: IC の差は参考値で、IC 依存あり/なしを確定しない)"
        for irow in ic_rows:
            irow["verdict_if_preconditions_held"] = irow["verdict"]
            irow["verdict"] = IC_REF_LABEL
    else:
        ic_overall = ic_overall_if_held
    # 総合 (どれか 1 run でも条件を満たさなければ保留。窓の移動・延長はしない)
    if pre_bad:
        overall = "保留 (前提不成立): " + "; ".join(f"{r}: {', '.join(w)}" for r, w in pre_bad.items())
    elif not complete:
        overall = "保留 (欠損): " + "; ".join(f"{r}: {w}" for r, w in missing.items())
    elif ic_overall != IC_WIN_LABEL:
        overall = f"保留 (IC 依存の前提が成り立たない: {ic_overall})"
    elif all(v == WIN_LABEL for v in verdicts.values()):
        det = [r["qty"] for r in rows if r.get("detected") and r["qty"] in DQ]
        overall = "許容幅内 (Euler で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし")
    elif any(v == WORSE_LABEL for v in verdicts.values()):
        overall = "悪化 (D − 2·SE ≥ Δq の量: " + ", ".join(k for k, v in verdicts.items() if v == WORSE_LABEL) + ")"
    else:
        overall = "保留 (ユーザ判断: " + ", ".join(k for k, v in verdicts.items() if v != WIN_LABEL) + ")"
    if exit_cal is not None:
        if pre_bad or ic_overall != IC_WIN_LABEL:
            exit_cal["verdict_if_preconditions_held"] = exit_cal["verdict"]
            exit_cal["verdict"] = "保留 (前提不成立・IC 依存の前提; 較正は変えない。V5 の較正判断は未完了)"
        exit_cal["decision_complete"] = not exit_cal["verdict"].startswith("保留")
    # 記録用 (判定に使わない): 評価座標の感度 (座標の整合の許容差 1e-6 r_t と各 run 自身の座標)
    coord = coordinate_record(case, srec, runs_all, S, coord_sensitivity)
    # 記録用: 腕 B の v5 の時系列が E′ の時系列 (wallfit_series_e3.csv) と同じか (同じ評価器で再計算したことの確認)
    recompute = {}
    for run in ARMS["B"]:
        f3 = case / run / "wallfit_series_e3.csv"
        if run in S and f3.is_file():
            try:
                rows3 = [r for r in csv.DictReader(open(f3)) if WIN[0] <= int(float(r["step"])) <= WIN[1]]
                recompute[run] = max(float(np.abs(np.array([float(r[k]) for r in rows3]) - S[run][k]).max())
                                     for k in list(DQ) + list(RECORD_ONLY))
            except Exception as e:  # noqa: BLE001
                recompute[run] = f"{type(e).__name__}: {e}"
        else:
            recompute[run] = None
    out = {"plan": PLAN, "plan_registration_commit": PLAN_REG_COMMIT,
           "evaluator": Path(__file__).name, "evaluator_sha256": _sha(Path(__file__).resolve()),
           "window_steps": WIN, "n_per_run": len(WIN_STEPS), "window_fixed": True, "arms": ARMS, "isen": ISEN, "tag": TAG,
           "series_record": str(sp.relative_to(case)) if sp.is_file() else None,
           "series_geometry": srec.get("geometry"), "series_geom_ref": srec.get("geom_ref"),
           "series_evaluator_sha256": srec.get("evaluator_sha256"),
           "quasisteady": {"drift": QS_DRIFT, "osc": QS_OSC, "tail": QS_TAIL, "full": len(WIN_STEPS), "status": qs_status,
                           "exception_qty": EXIT_EXC_QTY, "exception_full_tol": EXIT_EXC_FULL_TOL},
           "rows": rows, "overall": overall,
           "ic_dependence": {"rows": ic_rows, "overall": ic_overall, "overall_if_preconditions_held": ic_overall_if_held,
                             "two_sided": True, "exit_core_M_tol": IC_EXIT_M_TOL},
           "exit_calibration": exit_cal, "preconditions_failed": pre_bad, "precondition_exceptions": pre_exc,
           "missing": missing, "run_meta": metas, "evaluation_coordinates": coord,
           "armB_v5_vs_e3_series_maxabs_record_only": recompute,
           "limits": ["自己相関は補正していない (SE は実務の指標で、保証された信頼限界・厳密な非劣化・IC 独立性を主張しない)",
                      "窓 6000〜18000 は E′ と同じ (E′ では結果を見た後に決めた。V5 では結果を見る前に固定し、動かさない・延長しない)",
                      "IC 依存は 1 本 (ISEN) と腕 M の比較で、ISEN の SE は時間変動だけ", "軸 (η0) の量は判定対象外",
                      f"「絶対許容内」は {EXIT_EXC_QTY} だけの例外で、STEADY ではない (元の VERDICT は run_meta に残す)",
                      f"{GEOM_TOL:g} r_t は{GEOM_TOL_LABEL}であって評価誤差の保証ではない (感度は evaluation_coordinates に記録)",
                      "IC の検査 C1〜C5 の合格は幾何の対応であって IC 非依存の証明ではない (IC 依存は ISEN で見る)",
                      "腕 B は MOC の実装より前に準備したので prepare_info に moc が無く、legacy・fixed2 とみなした"]}
    (case / "_band_ab").mkdir(exist_ok=True)
    (case / "_band_ab/moc_v5_euler_eval.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    case = Path(argv[0]).resolve() if argv else HERE
    out = evaluate(case)
    print(f"評価器 {out['evaluator']} sha256 {out['evaluator_sha256']}  (plan の追加登録 {out['plan_registration_commit']})")
    print(f"{'量':28s} {'Δq':>8s} {'B 平均':>11s} {'M 平均':>11s} {'D=M−B':>10s} {'2·SE_D':>9s} {'(D+2SE)/Δq':>10s} {'検出':>4s}  判定"
          f"   | {'ISEN−M':>10s} {'2·SE_ic':>9s}  IC 依存")
    ic = {r["qty"]: r for r in out["ic_dependence"]["rows"]}
    for r in out["rows"]:
        i = ic[r["qty"]]
        print(f"{r['qty']:28s} {r.get('dq', float('nan')):8.2e} {r['B']['mean']:11.6g} {r['M']['mean']:11.6g} {r['D_M_minus_B']:+10.2e} "
              f"{2 * r['SE_D']:9.2e} {r.get('D_plus_2SE_over_dq', float('nan')):10.2f} {'あり' if r['detected'] else 'なし':>4s}  {r['verdict']}"
              f"   | {i['D_ISEN_minus_M']:+10.2e} {2 * i['SE_ic']:9.2e}  {i['verdict']} (Δq {i.get('dq', float('nan')):.2g})")
    for run, st in out["quasisteady"]["status"].items():
        ns = {c: v for c, v in (st or {}).items() if v != "STEADY"}
        print(f"準定常 {run}: " + ("全量 STEADY (末尾 5 枚・全 13 枚・記録のファイル)" if st and not ns else
                                  ", ".join(f"{c} {v}" for c, v in ns.items()) if ns else "判定できない"))
    for e in out["precondition_exceptions"]:
        print("前提の例外:", e)
    for r, w in out["preconditions_failed"].items():
        print(f"前提不成立 {r}: " + " / ".join(w))
    if out["exit_calibration"]:
        ec = out["exit_calibration"]
        print(f"出口較正: D {ec['D_M_minus_B']:+.3e}  |D| + 2SE {ec['absD_plus_2SE']:.3e} (≤ {ec['tol']:g})  → {ec['verdict']}")
        print(f"  注: {ec['note']}")
    print(f"評価座標: {GEOM_TOL_LABEL} {GEOM_TOL:g} r_t (評価誤差の保証ではない。感度は JSON の evaluation_coordinates に記録)")
    print("IC 依存:", out["ic_dependence"]["overall"])
    print("総合:", out["overall"])
    print("->", case / "_band_ab/moc_v5_euler_eval.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
