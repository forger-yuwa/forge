"""plan discretization-moc-axis-limit-and-corrector §6 V5 (事前登録 2026-10-07) の判定。
腕 B (現行の単調壁) run_0143〜0145 対 腕 M (MOC の軸処理を変えた壁) run_0150〜0152 の許容判定、出口較正、IC 依存 (run_0153 対 腕 M)。
入力: eval_wallfit_euler.py --series (tag v5) が書いた各 run の wallfit_series_v5.csv・QUASISTEADY_wallfit_v5.txt と
  _band_ab/wallfit_series_v5.json (7 本を同じ呼び出し・同じ X_E・X_F で評価した記録)、各 run の prepare_info.json・IC_MAP.json・
  IC_INSPECTION.json (腕 M)・CONVERGENCE_VERDICT_segment.txt・nozzle.h5 (壁の証拠)。
量と Δq (§6 V5): M 波・P 波・オーバーシュート・出口規格化オーバーシュート・|P 傾き| (いずれも r/r_w = 0.1) と |出口コア M − 6|。
  Δq = 0.001 / 0.010 / 0.003 / 0.003 / 0.03 %pt、1.8e-4。出口コア M は記録用 (出口較正に使う)。
窓: 本段 step 6000〜18000 の 13 枚 (1000 ごと)。
統計 (§6 E′ と同じ定義): run i の窓平均 m_i、時間変動の標準誤差 s_i = sd_i/√13 (自己相関は補正しない; 限界として記録)。
  腕の平均 M = mean(m_i)、腕の標準誤差 SE = max(sd(m_i)/√k, √(Σ s_i²)/k) (k = 腕の run 数; k = 1 のときは s_1)。
  D = M_M − M_B、SE_D = √(SE_B² + SE_M²)。
許容判定 (大きいほど悪い量の片側): 全量で D + 2·SE_D ≤ Δq → 「許容幅内」/ D − 2·SE_D ≥ Δq の量があれば「悪化」/ それ以外は保留。
  別欄で差の検出の有無 (|D| > 2·SE_D) を記録する (判定には使わない)。
出口較正: 出口コア M の D について |D| + 2·SE_D ≤ 1e-4 → 「据え置き」(Md_moc_offset を変えない) / |D| − 2·SE_D > 1e-4 → 「較正をやり直す」/
  それ以外 → 「保留 (判別不能)」(較正は変えない)。
IC 依存 (§6 V5「2 本の評価量の差が許容判定の条件を満たすことを前提にする」): D_ic = q_ISEN − M_M、SE_ic = √(SE_M² + s_ISEN²)。
  IC の違いには良し悪しの向きが無いので両側で判定する (片側だと ISEN が「良い」向きの大差を通してしまう):
  全量で |D_ic| + 2·SE_ic ≤ Δq → 「IC 依存は許容幅内」/ |D_ic| − 2·SE_ic ≥ Δq の量があれば「IC 依存あり」/ それ以外は保留。
  前提なので、「IC 依存は許容幅内」でなければ総合は保留。
前提検査 (どれか不成立なら総合は「保留 (前提不成立)」。証拠が読めないこと・欠損も不成立。既定値で通さない):
  全 run: 時系列の記録 (wallfit_series_v5.json) にあり status ok、CSV・準定常ファイルの sha256 が記録と同じ (古い成果物でない)、
    記録の X_E・X_F の基準が腕 B の r1 で、各 run 自身の X_E・X_F との差 ≤ 1e-6 r_t (評価の刻み 0.05 r_t の 2e-5)、
    窓の 13 枚がそろい有限、判定区間 (CONVERGENCE_VERDICT_segment.txt) の残差判定が pass か既知の plateau、
    判定する 6 量と出口コア M の check_quasisteady (末尾 5 枚) が STEADY (§6 E4 と同じ扱い。近零量の「絶対許容内」は V5 に
    登録していないので適用せず、記録だけする)、乾式確認 (DRY) の準備から作った run でない、wall_fit.mono_r2 = [0.0, 1.5]、
    壁節点が自腕の当てはめ後 spline に一致 (M4)、保存 spline が [0, 1.5] で単調 (形状ゲート S1 と同じ許容差)。
  腕 B: prepare_info に moc キーが無ければ legacy・fixed2 とみなす (2026-10-07 の MOC の実装より前に準備した run。例外として記録)。
    moc キーがあれば legacy・fixed2 であること。IC_MAP.json が VERDICT OK・mode index・IC は run_0114 の res_6000、3 本の写像後 sha256 が同じ。
  腕 M: moc が analytic・converge でゲート合格 (applicable かつ pass)。IC_MAP.json が VERDICT OK・mode index・IC は run_0114 の res_6000、
    IC_INSPECTION.json が VERDICT OK・全条件合格で、上限が IC_MAP.json と prepare_info の上限と同じ、写像の移動の最大 ≤ 上限、
    3 本の写像後 sha256 が同じ。腕 B と見分けられる壁節点があり、そのすべてが自腕の当てはめに一致 (x の float32 の丸めは許す)。
  ISEN: moc のゲート合格、腕 M と同じ壁 (保存 spline が完全一致) と同じ格子 (IC 記録の格子ハッシュが腕 M と同じ)、
    IC が等エントロピー (IC_MAP.json の mode isentropic・VERDICT OK)。
  同じ実効設定: 7 本すべてで、本段の solverConfig.yaml・bcondConfig.yaml (YAML として)、forge のバイナリ (RUN_PROVENANCE.txt の
    forge_sha256)、起動時の実効値の行、段の並び (soft → main) と hard キー (stage_manifest.json) が腕 B の r1 と同じ。
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
ARMS = {"B": [f"run_0{n}_euler_wallfit_monoG1_r{k}" for n, k in ((143, 1), (144, 2), (145, 3))],
        "M": [f"run_0{n}_euler_wallfit_mocG1_r{k}" for n, k in ((150, 1), (151, 2), (152, 3))]}
ISEN = "run_0153_euler_icdep_mocG1_isen"
TAG = "v5"
GEOM_REF = ARMS["B"][0]
GEOM_TOL = 1e-6                                   # [r_t] 各 run 自身の X_E・X_F と評価の基準との差の上限
WIN = (6000, 18000)
WIN_STEPS = list(range(WIN[0], WIN[1] + 1, 1000))   # 13 枚
DQ = {"M_wave_eta0.1": 0.001, "P_wave_eta0.1": 0.010, "overshoot_eta0.1": 0.003, "overshoot_exitnorm_eta0.1": 0.003,
      "P_slope_abs_eta0.1": 0.03, "exit_M_dev": 0.00018}
RECORD_ONLY = ("exit_core_M",)
EXIT_CAL_TOL = 1e-4
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


def load_signed_exit_tail(case: Path, run: str, n: int):
    """記録用: 末尾 n 枚の (step, 出口コア M − 6) (近零量の「絶対許容内」の参考値)。"""
    rows = list(csv.DictReader(open(case / run / f"wallfit_series_{TAG}.csv")))
    return [int(float(r["step"])) for r in rows][-n:], [float(r["exit_core_M"]) - 6.0 for r in rows][-n:]


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
    if abs(D) + 2.0 * SE <= tol:
        return "据え置き (Md_moc_offset を変えない)"
    if abs(D) - 2.0 * SE > tol:
        return "較正をやり直す"
    return "保留 (判別不能; 較正は変えない)"


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
            if not (isinstance(v, (int, float)) and np.isfinite(v) and abs(v) <= GEOM_TOL):
                why.append(f"自身の {k[1:]} と評価の基準の差 {v!r} が {GEOM_TOL:g} r_t を超える・欠損")
        meta["series"] = {k: r.get(k) for k in ("dX_E", "dX_F", "n_snaps", "last_step", "status")}
    # 判定区間の残差
    p = rd / SEGMENT_VERDICT_FILE
    seg = parse_segment_verdict(p.read_text() if p.is_file() else None)
    meta["segment"] = seg
    if seg["status"] not in ("pass", "plateau"):
        why.append(f"判定区間の残差判定が {seg['status']} ({seg.get('reason') or seg.get('line')})")
    # 準定常 (末尾 5 枚)
    qp = rd / f"QUASISTEADY_wallfit_{TAG}.txt"
    qs = parse_quasisteady(qp.read_text() if qp.is_file() else "", list(DQ) + list(RECORD_ONLY))
    meta["quasisteady"] = qs
    ns = [c for c, v in qs.items() if v != "STEADY"]
    if ns:
        why.append("末尾 5 枚が STEADY でない量: " + ", ".join(f"{c} {qs[c]}" for c in ns))
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
        meta["ic"].update(limit_m=lim, max_m=dmax, inspection_VERDICT=ins.get("VERDICT"), inspection_limit_m=ins.get("limit_m"))
        if ins.get("VERDICT") != "OK" or not ins.get("conditions") or not all(v.get("ok") is True for v in ins["conditions"].values()):
            why.append(f"IC の検査 (IC_INSPECTION.json) が OK・全条件合格でない ({ins.get('VERDICT')!r})")
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


def evaluate(case: Path, wall_check=None) -> dict:
    case = Path(case).resolve()
    sys.path.insert(0, str(HERE))
    wall_check = wall_check or default_wall_check
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
    # 評価量 (読めない run は欠損)
    S, missing = {}, {}
    for run in ARMS["B"] + ARMS["M"] + [ISEN]:
        try:
            S[run] = load_window(case, run)
        except Exception as e:  # noqa: BLE001 — 欠損は欠損として記録 (0 で埋めない)
            missing[run] = f"{type(e).__name__}: {e}"
    complete = not missing
    rows, verdicts, ic_rows, ic_verdicts, exit_cal = [], {}, [], {}, None
    if complete:
        for k in list(DQ) + list(RECORD_ONLY):
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
                verdicts[k], ic_verdicts[k] = row["verdict"], irow["verdict"]
            else:
                row["verdict"] = irow["verdict"] = "記録のみ"
            rows.append(row)
            ic_rows.append(irow)
        ex = next(r for r in rows if r["qty"] == "exit_core_M")
        exit_cal = {"qty": "exit_core_M", "D_M_minus_B": ex["D_M_minus_B"], "SE_D": ex["SE_D"], "tol": EXIT_CAL_TOL,
                    "absD_plus_2SE": abs(ex["D_M_minus_B"]) + 2 * ex["SE_D"],
                    "verdict": judge_exit_calibration(ex["D_M_minus_B"], ex["SE_D"])}
    # IC 依存の総合
    if not complete:
        ic_overall = "保留 (欠損)"
    elif all(v == "IC 依存は許容幅内" for v in ic_verdicts.values()):
        ic_overall = "IC 依存は許容幅内"
    elif any(v == "IC 依存あり" for v in ic_verdicts.values()):
        ic_overall = "IC 依存あり: " + ", ".join(k for k, v in ic_verdicts.items() if v == "IC 依存あり")
    else:
        ic_overall = "保留 (判別不能: " + ", ".join(k for k, v in ic_verdicts.items() if v != "IC 依存は許容幅内") + ")"
    # 総合
    if pre_bad:
        overall = "保留 (前提不成立): " + "; ".join(f"{r}: {', '.join(w)}" for r, w in pre_bad.items())
    elif not complete:
        overall = "保留 (欠損): " + "; ".join(f"{r}: {w}" for r, w in missing.items())
    elif ic_overall != "IC 依存は許容幅内":
        overall = f"保留 (IC 依存の前提が成り立たない: {ic_overall})"
    elif all(v == WIN_LABEL for v in verdicts.values()):
        det = [r["qty"] for r in rows if r.get("detected") and r["qty"] in DQ]
        overall = "許容幅内 (Euler で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし")
    elif any(v == WORSE_LABEL for v in verdicts.values()):
        overall = "悪化 (D − 2·SE ≥ Δq の量: " + ", ".join(k for k, v in verdicts.items() if v == WORSE_LABEL) + ")"
    else:
        overall = "保留 (ユーザ判断: " + ", ".join(k for k, v in verdicts.items() if v != WIN_LABEL) + ")"
    if exit_cal is not None and (pre_bad or ic_overall != "IC 依存は許容幅内"):
        exit_cal["verdict_if_preconditions_held"] = exit_cal["verdict"]
        exit_cal["verdict"] = "保留 (前提不成立・IC 依存の前提; 較正は変えない)"
    # 記録用: 近零量の「絶対許容内」(§6 E4 の方法; V5 では判定に使わない)
    abs_tol = {}
    try:
        from throat_mono_judge import ABS_TOL_N, DELTA_Q, abs_tolerance_verdict
        for run in S:
            st_, sv_ = load_signed_exit_tail(case, run, ABS_TOL_N)
            abs_tol[run] = abs_tolerance_verdict(st_, sv_, DELTA_Q["exit_M_dev"] / 10.0).get("status")
    except Exception as e:  # noqa: BLE001 — 記録用 (判定に使わない)
        abs_tol["error"] = f"{type(e).__name__}: {e}"
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
    out = {"plan": PLAN, "window_steps": WIN, "n_per_run": len(WIN_STEPS), "arms": ARMS, "isen": ISEN, "tag": TAG,
           "series_record": str(sp.relative_to(case)) if sp.is_file() else None,
           "series_geometry": srec.get("geometry"), "series_geom_ref": srec.get("geom_ref"),
           "series_evaluator_sha256": srec.get("evaluator_sha256"),
           "rows": rows, "overall": overall, "ic_dependence": {"rows": ic_rows, "overall": ic_overall, "two_sided": True},
           "exit_calibration": exit_cal, "preconditions_failed": pre_bad, "precondition_exceptions": pre_exc,
           "missing": missing, "run_meta": metas, "abs_tol_exit_M_dev_record_only": abs_tol,
           "armB_v5_vs_e3_series_maxabs_record_only": recompute,
           "limits": ["自己相関は補正していない (2SE は保証された信頼限界でない)", "窓 6000〜18000 は E′ と同じ (E′ では結果を見た後に決めた)",
                      "IC 依存は 1 本 (ISEN) と腕 M の比較で、ISEN の SE は時間変動だけ", "軸 (η0) の量は判定対象外",
                      "腕 B は MOC の実装より前に準備したので prepare_info に moc が無く、legacy・fixed2 とみなした"]}
    (case / "_band_ab").mkdir(exist_ok=True)
    (case / "_band_ab/moc_v5_euler_eval.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    case = Path(argv[0]).resolve() if argv else HERE
    out = evaluate(case)
    print(f"{'量':28s} {'Δq':>8s} {'B 平均':>11s} {'M 平均':>11s} {'D=M−B':>10s} {'2·SE_D':>9s} {'(D+2SE)/Δq':>10s} {'検出':>4s}  判定"
          f"   | {'ISEN−M':>10s} {'2·SE_ic':>9s}  IC 依存")
    ic = {r["qty"]: r for r in out["ic_dependence"]["rows"]}
    for r in out["rows"]:
        i = ic[r["qty"]]
        print(f"{r['qty']:28s} {r.get('dq', float('nan')):8.2e} {r['B']['mean']:11.6g} {r['M']['mean']:11.6g} {r['D_M_minus_B']:+10.2e} "
              f"{2 * r['SE_D']:9.2e} {r.get('D_plus_2SE_over_dq', float('nan')):10.2f} {'あり' if r['detected'] else 'なし':>4s}  {r['verdict']}"
              f"   | {i['D_ISEN_minus_M']:+10.2e} {2 * i['SE_ic']:9.2e}  {i['verdict']}")
    for e in out["precondition_exceptions"]:
        print("前提の例外:", e)
    for r, w in out["preconditions_failed"].items():
        print(f"前提不成立 {r}: " + " / ".join(w))
    if out["exit_calibration"]:
        ec = out["exit_calibration"]
        print(f"出口較正: D {ec['D_M_minus_B']:+.3e}  |D| + 2SE {ec['absD_plus_2SE']:.3e} (≤ {ec['tol']:g})  → {ec['verdict']}")
    print("IC 依存:", out["ic_dependence"]["overall"])
    print("総合:", out["overall"])
    print("->", case / "_band_ab/moc_v5_euler_eval.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
