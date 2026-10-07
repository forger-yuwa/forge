"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「NS の 3 条件」(N0・N1・N2) の評価器。
合否は plan tooling-nozzle-throat-monotone-r2 §6 N・K の転記 (U4 の登録):
  dry: 出口コア M 5.997〜6.003 (6.000 ± 0.05 %、2026-10-07 ユーザ決定。旧 ± 0.02 %)・波 η0.1 ≤ 0.01 %・オーバーシュート η0.1 ≤ +0.035 %・δ_E/δ_C = 1 ± 0.5 %・出口半径 0.775 m ± 0.1 mm、
       4 量 (出口コア M・オーバーシュート・δ_E/δ_C・波) が check_quasisteady --series-csv で STEADY、壁解像 PASS (y1+ > 1 の面積 ≤ 5 %)、
       check_convergence (--segment) は全残差で RISING なし (plateau の NOT CONVERGED は既知として記録)。
  凝縮: 凝縮 4 量 (cond_series.py: 軸の凝縮開始・S_max・出口コア g・出口コア M) が末尾 5 枚 STEADY、残差 RISING なし。
判定窓: dry は本段 60000〜80000 の 5 枚 (5000 ごと。延長したら連結した系列の末尾 5 枚 = 80000〜100000)、凝縮は 14000〜18000 の 5 枚。
記録 (判定なし): 質量流量・音速線の位置・スロート付近の壁圧・排除厚の補正量 (ns_n012.py record-series)、N1 − N0 (U4) と
  N2 − N1 (V5′) の差の表、生産 (run_0147 + 0149 / run_0148) との差。差は窓の平均の差で、窓の幅を並べる (差の信頼区間ではない)。

登録に無く、ここで決めた細部 (報告に列挙):
  - 各量の合否は判定窓 5 枚の平均で判定する (run_0147 + 0149 の報告と同じ)。5 枚の最小・最大と、5 枚すべてが範囲内かも記録する。
  - 出口半径は prepare_info.json の sizing.exit_radius_m (物理壁の x_e での値) を使い、report.json の値と 1e-9 m 以内で一致することを
    前提にする。
  - 準定常の VERDICT は check_quasisteady を本段の全系列 (延長は連結) に、末尾 5 枚 (--tail 4.5/n) で当てる
    (throat_mono_ns_verdicts.py と同じ呼び方。drift 0.05・osc 0.10 は既定)。
  - 残差: --segment の判定が pass・plateau・still converging なら「RISING なし」、rising・diverged は未達、それ以外は判定不能。延長したら
    本段と延長の両方に課す。
  - 前提 (満たさなければその条件は判定不能): 準備の記録 (較正値・問題の sha256・k_f)・IC の検査 OK・3 条件の固定の検査 OK・
    forge の sha256 が run_0147 と同じ・段の記録 (soft → mid → 本段)・RUN_RC 0・早期停止なし・NAN_SCAN CLEAN。
  - r_t (plan §6 U4「r_t の解き直し」): 渡した run (dry 3 本・延長・凝縮) の実効の r_t (prepare_info.json の scale_m) がすべて同じで、
    --r-throat を渡したときはその値とビット一致。どれかが違えば全条件を判定不能にする (r_t が違う比較は単独の変更にならない)。

usage (case dir):
  python3 ns_n012_eval.py eval --md-offset V [--r-throat R] --set N0=<dry>[+<ext>] --set N1=... --set N2=... [--cond N0=<cond> ...]
                                [--out JSON] [--dry]
  python3 ns_n012_eval.py gate-dry --md-offset V [--r-throat R] --cond-name N0 --run <dry> [--ext <ext>]   (終了コード 0 = dry の全ゲート合格)
出力: _band_ab/ns_n012_eval.json (既定) と、nozzle_report --verdicts 用の _band_ab/verdicts_<run_NNNN>_n012.json。
"""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parents[1] / "solver_density_cuda/tools"
sys.path.insert(0, str(HERE))

PLAN = "plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §6 U4 (NS の 3 条件); 合否は plans/accepted/tooling-nozzle-throat-monotone-r2.md §6 N・K"
CONDS = ("N0", "N1", "N2")
RECORD = "NS_N012.json"
OUT_JSON = "_band_ab/ns_n012_eval.json"
SEG_FILE = "CONVERGENCE_VERDICT_segment.txt"
MAIN_STEPS, OUT_INT, EXT_STEPS = 80000, 5000, 20000
COND_STEPS, COND_OUT = 18000, 1000
WIN = 5
# dry のゲート (monotone plan §6 N の転記)。値は判定窓 5 枚の平均
# 出口コア M: 6.000 ± 0.05 % (2026-10-07 ユーザ決定「プラマイ 0.05% にしようか」。NS の 3 条件の結果を見た後の決定で、
# 旧の ± 0.02 % [5.9988〜6.0012] は 3e-4 ≈ スロート半径 10 µm 相当で物理的に細かすぎるため。plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4)
EXIT_M_LO, EXIT_M_HI = 5.997, 6.003
WAVE_MAX = 0.01            # [%]
OVERSHOOT_MAX = 0.035      # [%]
# オーバーシュートの準定常: 平均が零に近く check_quasisteady の相対の drift が大きく出るので、窓 5 枚の絶対の幅が
# ゲートの 1/10 (0.0035 %) 以下なら準定常とみなす (2026-10-07 ユーザ決定「絶対の幅で可とする」。延長の後も 3 条件とも
# DRIFTING、窓の動きは約 0.0005 %。plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4)。判定の VERDICT はそのまま残す
OVERSHOOT_QS_ABS = OVERSHOOT_MAX / 10.0
DE_DC_TOL = 0.005          # |δ_E/δ_C − 1|
EXIT_R_M, EXIT_R_TOL = 0.775, 1.0e-4
WALL_OVER_PCT = 5.0
EXIT_R_AGREE = 1e-9        # prepare_info の出口半径と report.json の値の一致 [m]
QS_DRY = {"exitM_A": "exit_core_M", "overshoot01": "overshoot_eta0.1", "wave01": "wave_eta0.1", "dE_over_dC": "dE_over_dC"}
QS_COND = {"onset_x_axis": "cond_onset_x_axis", "S_max": "cond_S_max", "exit_g_core": "cond_exit_g_core", "exit_core_M": "exit_core_M"}
QS_VERDICTS = ("STEADY", "OSCILLATING", "TRANSIENT-UNSETTLED", "DRIFTING", "NONFINITE")
SEG_OK = ("pass", "plateau", "converging")      # 「RISING なし」(plateau・still converging は記録)
SEG_FAIL = ("rising", "diverged")
PASS, FAIL, UNDET = "合格", "未達", "判定不能"
REF_DRY = ("run_0147_ns_mono_final", "run_0149_ns_mono_final_ext")    # 生産 (比較元、記録のみ)
REF_COND = "run_0148_ns_mono_final_cond"
REF_PROV_RUN = "run_0147_ns_mono_final"
REF_RECORD_DIR = "_band_ab/ns_n012_ref"
RECORD_CSV = "ns_n012_record_series.csv"
STAGE_TAGS_DRY = ["S1_soft", "S2_mid", "main"]


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(p) -> str | None:
    p = Path(p)
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def jload(p):
    return json.loads(Path(p).read_text())


class SeriesError(ValueError):
    """時系列を判定に使えない (欠落・重複・逆順・窓の不一致・列の欠損・非数)。"""


def read_rows(p) -> list:
    p = Path(p)
    if not p.is_file():
        raise SeriesError(f"{p} が無い")
    with open(p) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SeriesError(f"{p} が空")
    return rows


def int_steps(rows, who) -> list:
    out = []
    for r in rows:
        try:
            x = float(r.get("step"))
        except (TypeError, ValueError):
            raise SeriesError(f"{who}: step {r.get('step')!r} が数でない") from None
        if not math.isfinite(x) or x != int(x):
            raise SeriesError(f"{who}: step {r.get('step')!r} が有限の整数でない")
        out.append(int(x))
    return out


def expect_steps(rows, interval, end, who) -> None:
    s = int_steps(rows, who)
    want = list(range(interval, end + 1, interval))
    if s != want:
        lack, extra = sorted(set(want) - set(s)), sorted(set(s) - set(want))
        raise SeriesError(f"{who}: step が {interval} ごとの {interval}〜{end} でない (欠落 {lack[:5]}、余分 {extra[:5]}、"
                          f"{'順序違い・重複' if not lack and not extra else ''})")


def join_rows(base_rows, ext_rows, parent_end, who="") -> list:
    """本段 (step 5000〜parent_end) と延長 (local 5000〜EXT_STEPS) を「通算 = local + parent_end」で連結する。"""
    expect_steps(base_rows, OUT_INT, parent_end, f"{who}本段")
    expect_steps(ext_rows, OUT_INT, EXT_STEPS, f"{who}延長")
    if set(base_rows[0]) != set(ext_rows[0]):
        raise SeriesError(f"{who}本段と延長の列が違う")
    out = [dict(r) for r in base_rows]
    for r in ext_rows:
        q = dict(r)
        q["step"] = str(int(float(r["step"])) + parent_end)
        out.append(q)
    return out


def column(rows, col, who) -> list:
    if col not in rows[0]:
        raise SeriesError(f"{who}: 列 {col} が無い")
    out = []
    for r in rows:
        try:
            v = float(r[col])
        except (TypeError, ValueError):
            v = float("nan")
        out.append(v)
    return out


def window_stats(rows, col, want_steps, who) -> dict:
    """末尾 WIN 枚の平均・最小・最大・幅。窓の step が want_steps と違えば SeriesError。非有限があれば平均は nan。"""
    tail = rows[-WIN:]
    st = int_steps(tail, who)
    if st != list(want_steps):
        raise SeriesError(f"{who}: 判定窓の step が {st} ({list(want_steps)} であること)")
    v = column(tail, col, who)
    fin = all(math.isfinite(x) for x in v)
    return {"mean": (sum(v) / len(v)) if fin else float("nan"), "min": min(v) if fin else float("nan"),
            "max": max(v) if fin else float("nan"), "range": (max(v) - min(v)) if fin else float("nan"),
            "values": v, "steps": st, "finite": fin}


def parse_quasisteady(text, cols) -> dict:
    """check_quasisteady --series-csv の出力 → {列: {verdict, detail}}。列の行が無い・読めない列は verdict None。"""
    out = {c: {"verdict": None, "detail": None} for c in cols}
    if "ERROR" in text:
        for c in cols:
            out[c]["detail"] = text.strip().splitlines()[-1] if text.strip() else "ERROR"
        return out
    for line in text.splitlines():
        m = re.match(r"^\s{2}(\S+)\s*: (.*?)\s+(" + "|".join(QS_VERDICTS) + r")\s*$", line)
        if m and m.group(1) in out:
            out[m.group(1)] = {"verdict": m.group(3), "detail": m.group(2).strip()}
    return out


def quasisteady(csv_path: Path, cols, n_rows: int) -> dict:
    """check_quasisteady.py --series-csv を末尾 5 枚 (--tail 4.5/n) で回す (throat_mono_ns_verdicts.py と同じ呼び方)。"""
    tail = 4.5 / n_rows
    q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(csv_path),
                        "--series-cols", ",".join(cols), "--tail", f"{tail:.6f}"], capture_output=True, text=True)
    res = parse_quasisteady(q.stdout + q.stderr, cols)
    for c in cols:
        res[c]["tail_arg"] = f"{tail:.6f}"
    return res


def write_rows(p: Path, rows) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def segment(run: Path) -> dict:
    from throat_mono_judge import parse_segment_verdict
    p = run / SEG_FILE
    st = parse_segment_verdict(p.read_text() if p.is_file() else None)
    st["run"] = run.name
    st["gate"] = PASS if st["status"] in SEG_OK else (FAIL if st["status"] in SEG_FAIL else UNDET)
    return st


def item(key, value, threshold, ok, detail=None) -> dict:
    """ゲート 1 項目。ok は True / False / None (判定不能)。"""
    return {"key": key, "value": value, "threshold": threshold, "status": PASS if ok is True else (FAIL if ok is False else UNDET),
            "detail": detail}


def overall(items) -> str:
    st = [i["status"] for i in items]
    if UNDET in st:
        return UNDET
    return PASS if all(s == PASS for s in st) else FAIL


# --- 前提 -----------------------------------------------------------------------------------------------------------
def provenance_sha(run: Path) -> str | None:
    p = run / "RUN_PROVENANCE.txt"
    if not p.is_file():
        return None
    m = re.search(r"^forge_sha256:\s*(\S+)", p.read_text(), re.M)
    return m.group(1) if m else None


def scale_of(run: Path):
    """実効の r_t (prepare_info.json の scale_m)。読めなければ None。"""
    try:
        return jload(run / "prepare_info.json").get("scale_m")
    except (OSError, ValueError):
        return None


def r_throat_failures(runs, r_throat: float | None) -> list:
    """渡した run の実効の r_t がすべて同じで、r_throat (指定時) と一致するか。不成立の理由の list。"""
    got = {Path(r).name: scale_of(Path(r)) for r in runs}
    why = []
    bad = {k: v for k, v in got.items() if not isinstance(v, float)}
    if bad:
        why.append(f"実効の r_t (scale_m) を読めない run: {bad}")
    vals = {v for v in got.values() if isinstance(v, float)}
    if len(vals) > 1:
        why.append(f"実効の r_t が run で違う: {got}")
    if r_throat is not None:
        off = {k: v for k, v in got.items() if v != r_throat}
        if off:
            why.append(f"実効の r_t が指定値 {r_throat!r} と違う run: {off}")
    return why


def run_preconditions(case: Path, run: Path, role: str, cond: str, md_offset: float, dry: bool, ref_sha: str | None) -> list:
    """1 本の run の前提。不成立の理由の list (空 = 成立)。"""
    why = []
    try:
        r = jload(run / RECORD)
    except (OSError, ValueError) as e:
        return [f"{run.name}: 準備の記録 {RECORD} を読めない ({e})"]
    if r.get("role") != role or r.get("condition") != cond:
        why.append(f"{run.name}: 記録の役割・条件 {r.get('role')}/{r.get('condition')} ({role}/{cond} であること)")
    if r.get("md_offset") != md_offset:
        why.append(f"{run.name}: 記録の較正値 {r.get('md_offset')!r} ({md_offset!r} であること)")
    if bool(r.get("dry")) != bool(dry):
        why.append(f"{run.name}: 乾式の印 {r.get('dry')!r} が要求 ({dry}) と違う")
    prob = case / str(r.get("problem"))
    if sha256_file(prob) != r.get("problem_sha256"):
        why.append(f"{run.name}: 問題 {r.get('problem')} の sha256 が準備のときと違う")
    if role == "dry":
        ic = jload(run / "IC_CHECK.json") if (run / "IC_CHECK.json").is_file() else {}
        if ic.get("VERDICT") != "OK":
            why.append(f"{run.name}: IC の検査が OK でない ({ic.get('VERDICT')})")
    if not dry:
        want_tags = STAGE_TAGS_DRY if role == "dry" else ["main"]      # 延長・凝縮は段なし (本段だけ)
        try:
            sm = jload(run / "stage_manifest.json")
            tags = [s.get("tag") for s in sm.get("stages", [])]
            if tags != want_tags:
                why.append(f"{run.name}: 段の記録 {tags} ({want_tags} であること)")
        except (OSError, ValueError) as e:
            why.append(f"{run.name}: stage_manifest.json を読めない ({e})")
    if role == "ext":
        if r.get("parent_end") != MAIN_STEPS or r.get("ext_steps") != EXT_STEPS:
            why.append(f"{run.name}: 延長の記録 (親の終了 {r.get('parent_end')}、{r.get('ext_steps')} step) が登録と違う")
    if not dry:
        rc = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
        if rc != "0":
            why.append(f"{run.name}: RUN_RC が 0 でない ({rc})")
        if (run / "EARLY_STOP.txt").exists():
            why.append(f"{run.name}: 早期停止 (EARLY_STOP.txt)")
        ns = jload(run / "NAN_SCAN.json") if (run / "NAN_SCAN.json").is_file() else {}
        if ns.get("VERDICT") != "CLEAN":
            why.append(f"{run.name}: NAN_SCAN が CLEAN でない ({ns.get('VERDICT')})")
        got = provenance_sha(run)
        if ref_sha is None or got != ref_sha:
            why.append(f"{run.name}: forge の sha256 {got} が {REF_PROV_RUN} ({ref_sha}) と違う・読めない")
    return why


# --- dry ------------------------------------------------------------------------------------------------------------
def dry_series(base: Path, ext: Path | None, name: str, who: str) -> tuple:
    """(行, 窓の step, 系列の説明)。延長があれば連結する (延長の親が base であること)。"""
    b = read_rows(base / name)
    if ext is None:
        expect_steps(b, OUT_INT, MAIN_STEPS, f"{who}{base.name}")
        return b, list(range(MAIN_STEPS - (WIN - 1) * OUT_INT, MAIN_STEPS + 1, OUT_INT)), f"{base.name} 5000〜{MAIN_STEPS}"
    er = jload(ext / RECORD)
    if er.get("parent") != base.name or er.get("parent_end") != MAIN_STEPS:
        raise SeriesError(f"{who}延長 {ext.name} の親が {er.get('parent')} (終了 {er.get('parent_end')}) — {base.name} の {MAIN_STEPS} であること")
    rows = join_rows(b, read_rows(ext / name), MAIN_STEPS, who)
    end = MAIN_STEPS + EXT_STEPS
    return rows, list(range(end - (WIN - 1) * OUT_INT, end + 1, OUT_INT)), f"{base.name} 5000〜{MAIN_STEPS} + {ext.name} (+{MAIN_STEPS})"


def gate_dry(case: Path, cond: str, base: Path, ext: Path | None, out_dir: Path) -> dict:
    """dry 1 条件のゲート (monotone plan §6 N の転記)。"""
    res = {"condition": cond, "runs": [base.name] + ([ext.name] if ext else []), "items": [], "quasisteady": {}}
    final = ext or base
    try:
        rows, wsteps, desc = dry_series(base, ext, "quantities_series.csv", "")
        res["series"] = desc
        res["window_steps"] = wsteps
        csvp = out_dir / f"ns_n012_series_{base.name[:8]}{'_joined' if ext else ''}.csv"
        write_rows(csvp, rows)
        res["series_csv"] = str(csvp.relative_to(case)) if csvp.is_relative_to(case) else str(csvp)
        qs = quasisteady(csvp, list(QS_DRY), len(rows))
        win = {c: window_stats(rows, c, wsteps, base.name) for c in QS_DRY}
        res["window"] = {QS_DRY[c]: {k: win[c][k] for k in ("mean", "min", "max", "range", "steps", "finite")} for c in QS_DRY}
        for c, key in QS_DRY.items():
            res["quasisteady"][key] = qs[c]
    except SeriesError as e:
        res["items"].append(item("時系列 (quantities_series.csv)", None, "5000 ごと・窓 5 枚", None, str(e)))
        win, qs = None, None
    if win is not None:
        m = win["exitM_A"]
        res["items"].append(item("出口コア M", m["mean"], f"{EXIT_M_LO}〜{EXIT_M_HI}",
                                 (EXIT_M_LO <= m["mean"] <= EXIT_M_HI) if m["finite"] else None,
                                 f"窓 5 枚の最小 {m['min']:.7f}・最大 {m['max']:.7f}、5 枚とも範囲内 "
                                 f"{all(EXIT_M_LO <= v <= EXIT_M_HI for v in m['values'])}"))
        w = win["wave01"]
        res["items"].append(item("Mach 波 η0.1 [%]", w["mean"], f"≤ {WAVE_MAX}", (w["mean"] <= WAVE_MAX) if w["finite"] else None,
                                 f"5 枚の最大 {w['max']:.6g}"))
        o = win["overshoot01"]
        res["items"].append(item("オーバーシュート η0.1 [%]", o["mean"], f"≤ +{OVERSHOOT_MAX}",
                                 (o["mean"] <= OVERSHOOT_MAX) if o["finite"] else None, f"5 枚の最大 {o['max']:.6g}"))
        d = win["dE_over_dC"]
        res["items"].append(item("δ_E/δ_C (x_F)", d["mean"], f"1 ± {DE_DC_TOL}",
                                 (abs(d["mean"] - 1.0) <= DE_DC_TOL) if d["finite"] else None, None))
        for c, key in QS_DRY.items():
            v = qs[c]["verdict"]
            if c == "overshoot01" and v is not None and v != "STEADY" and o["finite"]:
                ok = o["range"] <= OVERSHOOT_QS_ABS
                res["items"].append(item(f"準定常 {key}", v, f"STEADY または窓の幅 ≤ {OVERSHOOT_QS_ABS:g} %", ok,
                                         f"{qs[c]['detail']}; 窓 5 枚の絶対の幅 {o['range']:.3g} % "
                                         f"({'≤' if ok else '>'} {OVERSHOOT_QS_ABS:g} %、2026-10-07 ユーザ決定の絶対の幅の条件)"))
            else:
                res["items"].append(item(f"準定常 {key}", v, "STEADY", None if v is None else (v == "STEADY"), qs[c]["detail"]))
    # 出口半径 (幾何。prepare_info の sizing と report.json の一致を前提に)
    try:
        info = jload(base / "prepare_info.json")
        r_pi = (info.get("sizing") or {}).get("exit_radius_m")
        rep = jload(base / "report/report.json")["metrics"]
        r_rep = rep["wall_shape"]["exit_radius_m"]
        agree = r_pi is not None and abs(r_pi - r_rep) <= EXIT_R_AGREE
        res["items"].append(item("出口半径 [m]", r_pi, f"{EXIT_R_M} ± {EXIT_R_TOL}",
                                 (abs(r_pi - EXIT_R_M) <= EXIT_R_TOL) if agree else None,
                                 f"report.json {r_rep!r}" + ("" if agree else f" — prepare_info と {EXIT_R_AGREE:g} m 以内で一致しない・無い")))
    except (OSError, ValueError, KeyError, TypeError) as e:
        res["items"].append(item("出口半径 [m]", None, f"{EXIT_R_M} ± {EXIT_R_TOL}", None, f"読めない: {type(e).__name__}: {e}"))
    # 壁解像 (最終の run の最終場)
    try:
        wr = jload(final / "report/report.json")["metrics"]["wall_resolution"]
        ok = None
        if wr.get("over_area_allow_pct") == WALL_OVER_PCT and wr.get("verdict") in ("PASS", "FAIL"):
            ok = (wr["verdict"] == "PASS" and wr.get("over_area_pct") is not None and wr["over_area_pct"] <= WALL_OVER_PCT)
        res["items"].append(item("壁解像 (y1+ > 1 の面積 %)", wr.get("over_area_pct"), f"PASS・≤ {WALL_OVER_PCT}", ok,
                                 f"{final.name}: {wr.get('verdict_line')} (許容 {wr.get('over_area_allow_pct')} %、y1+ 最大 {wr.get('y1p_max')})"))
    except (OSError, ValueError, KeyError, TypeError) as e:
        res["items"].append(item("壁解像", None, f"PASS・≤ {WALL_OVER_PCT}", None, f"読めない: {type(e).__name__}: {e}"))
    # 残差 (本段と延長の両方)
    segs = [segment(base)] + ([segment(ext)] if ext else [])
    res["convergence"] = segs
    st = [s["gate"] for s in segs]
    res["items"].append(item("残差 RISING なし (check_convergence --segment)", "; ".join(f"{s['run']}: {s['status']}" for s in segs),
                             "RISING・DIVERGED なし", None if UNDET in st else all(x == PASS for x in st),
                             "; ".join(f"{s['run']}: {s['line']}" for s in segs)))
    res["overall"] = overall(res["items"])
    res["unmet"] = [i["key"] for i in res["items"] if i["status"] == FAIL]
    res["undetermined"] = [i["key"] for i in res["items"] if i["status"] == UNDET]
    if res["overall"] == FAIL and ext is None:
        res["next"] = "未達の量がある — 登録では延長 1 回 (20000 step) の対象 (延長するかは主セッションが判断)"
    elif res["overall"] == FAIL:
        res["next"] = "延長後も未達 — 未達のままユーザに判断を仰ぐ (登録)"
    return res


# --- 凝縮 ------------------------------------------------------------------------------------------------------------
def gate_cond(case: Path, cond: str, run: Path, out_dir: Path) -> dict:
    res = {"condition": cond, "run": run.name, "items": [], "quasisteady": {}}
    wsteps = list(range(COND_STEPS - (WIN - 1) * COND_OUT, COND_STEPS + 1, COND_OUT))
    try:
        rows = read_rows(run / "cond_series.csv")
        expect_steps(rows, COND_OUT, COND_STEPS, run.name)
        qs = quasisteady(run / "cond_series.csv", list(QS_COND), len(rows))
        win = {c: window_stats(rows, c, wsteps, run.name) for c in QS_COND}
        res["window_steps"] = wsteps
        res["window"] = {QS_COND[c]: {k: win[c][k] for k in ("mean", "min", "max", "range", "finite")} for c in QS_COND}
        for c, key in QS_COND.items():
            v = qs[c]["verdict"]
            res["quasisteady"][key] = qs[c]
            res["items"].append(item(f"準定常 {key}", v, "STEADY", None if v is None else (v == "STEADY"), qs[c]["detail"]))
    except SeriesError as e:
        res["items"].append(item("時系列 (cond_series.csv)", None, "1000 ごと・窓 5 枚", None, str(e)))
    s = segment(run)
    res["convergence"] = s
    res["items"].append(item("残差 RISING なし (check_convergence --segment)", s["status"], "RISING・DIVERGED なし",
                             None if s["gate"] == UNDET else s["gate"] == PASS, s["line"]))
    res["overall"] = overall(res["items"])
    res["unmet"] = [i["key"] for i in res["items"] if i["status"] == FAIL]
    if res["overall"] == FAIL:
        res["next"] = "未達 — 登録では延長 1 回、なお未達ならユーザ判断 (延長するかは主セッションが判断)"
    return res


# --- 記録 (判定なし) と差の表 -------------------------------------------------------------------------------------------
def record_window(case: Path, base: Path, ext: Path | None, csv_name: str = RECORD_CSV, base_csv: Path | None = None,
                  ext_csv: Path | None = None, parent_end: int = MAIN_STEPS) -> dict:
    """記録の時系列 (ns_n012.py record-series) の窓の平均と幅。列ごと。"""
    b = read_rows(base_csv or (base / csv_name))
    if ext is not None or ext_csv is not None:
        e = read_rows(ext_csv or (ext / csv_name))
        expect_steps(e, OUT_INT, EXT_STEPS, "延長の記録")
        bs = int_steps(b, "本段の記録")
        if bs[-1] != parent_end or bs != list(range(OUT_INT, parent_end + 1, OUT_INT)):
            raise SeriesError(f"本段の記録の step が 5000〜{parent_end} でない")
        rows = [dict(r) for r in b] + [dict(r, step=str(int(float(r["step"])) + parent_end)) for r in e]
        end = parent_end + EXT_STEPS
    else:
        expect_steps(b, OUT_INT, MAIN_STEPS, "記録")
        rows, end = b, MAIN_STEPS
    want = list(range(end - (WIN - 1) * OUT_INT, end + 1, OUT_INT))
    out = {}
    for c in rows[0]:
        if c == "step":
            continue
        w = window_stats(rows, c, want, "記録")
        out[c] = {"mean": w["mean"], "range": w["range"]}
    return {"window_steps": want, "values": out}


def diff_table(a: dict, b: dict) -> dict:
    """b − a (窓の平均の差) と、両方の窓の幅の大きい方 (差の信頼区間ではない)。a, b = {量: {mean, range}}。"""
    out = {}
    for k in a:
        if k in b:
            ma, mb = a[k]["mean"], b[k]["mean"]
            out[k] = {"diff": mb - ma, "a": ma, "b": mb, "range_max": max(a[k]["range"], b[k]["range"])}
    return out


def window_of(gate: dict) -> dict:
    return {k: {"mean": v["mean"], "range": v["range"]} for k, v in (gate.get("window") or {}).items()}


def reference_dry(case: Path) -> dict:
    """生産 (run_0147 + 延長 run_0149) の窓 60000〜80000 (記録のみ)。"""
    b, e = case / REF_DRY[0], case / REF_DRY[1]
    out = {"runs": list(REF_DRY)}
    try:
        base_rows = read_rows(b / "quantities_series.csv")
        ext_p = e / "quantities_series_extonly.csv"
        rows = join_rows(base_rows, read_rows(ext_p), 60000, "生産") if ext_p.is_file() else read_rows(e / "quantities_series.csv")
        want = list(range(60000, 80001, OUT_INT))
        out["window"] = {QS_DRY[c]: {k: window_stats(rows, c, want, "生産")[k] for k in ("mean", "range")} for c in QS_DRY}
    except SeriesError as ex:
        out["window_error"] = str(ex)
    try:
        rb, re_ = case / REF_RECORD_DIR / f"{REF_DRY[0]}_record.csv", case / REF_RECORD_DIR / f"{REF_DRY[1]}_record.csv"
        out["record"] = record_window(case, b, e, base_csv=rb, ext_csv=re_, parent_end=60000)
    except SeriesError as ex:
        out["record_error"] = str(ex)
    return out


def reference_cond(case: Path) -> dict:
    run = case / REF_COND
    out = {"run": REF_COND}
    try:
        rows = read_rows(run / "cond_series.csv")
        want = list(range(COND_STEPS - (WIN - 1) * COND_OUT, COND_STEPS + 1, COND_OUT))
        out["window"] = {QS_COND[c]: {k: window_stats(rows, c, want, REF_COND)[k] for k in ("mean", "range")} for c in QS_COND}
    except SeriesError as ex:
        out["window_error"] = str(ex)
    return out


def verdicts_json(case: Path, res: dict, kind: str) -> Path | None:
    """nozzle_report --verdicts 用 ({量: {verdict, window, note}})。"""
    qs = res.get("quasisteady") or {}
    if not qs or any(v.get("verdict") is None for v in qs.values()):
        return None
    win = res.get("window") or {}
    run = res["runs"][-1] if kind == "dry" else res["run"]
    out = {}
    for k, v in qs.items():
        w = win.get(k, {})
        out[k] = {"verdict": v["verdict"], "window": f"{res.get('series', run)} 窓 {res['window_steps'][0]}〜{res['window_steps'][-1]}",
                  "note": f"窓 5 枚の平均 {w.get('mean', float('nan')):.7g}・幅 {w.get('range', float('nan')):.3g} (ns_n012_eval.py)"}
    m = re.match(r"(run_\d{4})", Path(run).name)
    tag = m.group(1) if m else Path(run).name
    p = case / "_band_ab" / f"verdicts_{tag}_n012.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    return p


def parse_set(spec: str) -> tuple:
    """'N0=run_a' / 'N0=run_a+run_b' → (N0, run_a, run_b|None)。"""
    m = re.fullmatch(r"(N[012])=([^+=\s]+)(?:\+([^+=\s]+))?", spec)
    if not m:
        raise ValueError(f"{spec!r} は N0=<run>[+<延長の run>] の形でない")
    return m.group(1), m.group(2), m.group(3)


def evaluate(case: Path, md_offset: float, sets: dict, conds: dict, dry: bool, out_path: Path, r_throat: float | None = None) -> dict:
    """sets = {N0: (base, ext|None), ...}、conds = {N0: run, ...} (省略可)。r_throat = --r-throat (None = 同じことだけを要求)。"""
    ref_sha = provenance_sha(case / REF_PROV_RUN)
    all_runs = [case / x for c in sets for x in sets[c] if x] + [case / conds[c] for c in conds]
    rt_why = r_throat_failures(all_runs, r_throat)
    out = {"plan": PLAN, "tool": "ns_n012_eval.py", "tool_sha256": sha256_file(Path(__file__)), "date": now(), "dry": bool(dry),
           "md_offset": md_offset, "forge_sha256_ref": {"run": REF_PROV_RUN, "sha256": ref_sha},
           "r_throat": {"requested": r_throat, "effective": {r.name: scale_of(r) for r in all_runs}, "failures": rt_why},
           "inputs": {c: {"dry": [sets[c][0]] + ([sets[c][1]] if sets[c][1] else []), "cond": conds.get(c)} for c in sets},
           "preconditions": {}, "dry_gates": {}, "cond_gates": {}, "record": {}, "diffs": {}, "cond_diffs": {}}
    band = case / "_band_ab"
    band.mkdir(exist_ok=True)
    # 3 条件の固定 (準備時の検査)
    sc = case / ("_band_ab/ns_n012_prep_check_dry.json" if dry else "_band_ab/ns_n012_prep_check.json")
    set_ok = sc.is_file() and jload(sc).get("VERDICT") == "OK" and \
        [jload(sc)["runs"][c]["run"] for c in CONDS if c in sets] == [sets[c][0] for c in CONDS if c in sets]
    for c, (b, e) in sets.items():
        why = [] if set_ok else [f"3 条件の固定の検査 ({sc.name}) が OK でない・run の組が違う"]
        why += rt_why
        why += run_preconditions(case, case / b, "dry", c, md_offset, dry, ref_sha)
        if e:
            why += run_preconditions(case, case / e, "ext", c, md_offset, dry, ref_sha)
        if c in conds:
            why += run_preconditions(case, case / conds[c], "cond", c, md_offset, dry, ref_sha)
            try:
                kr = jload(case / conds[c] / RECORD)
                if kr.get("parent") not in (b, e):
                    why.append(f"{conds[c]}: 凝縮の親 {kr.get('parent')} が {b}{' / ' + e if e else ''} でない")
            except (OSError, ValueError):
                pass
        out["preconditions"][c] = {"ok": not why, "failures": why}
        g = gate_dry(case, c, case / b, (case / e) if e else None, band)
        if why:
            g["overall_raw"] = g["overall"]
            g["overall"] = UNDET
        out["dry_gates"][c] = g
        vp = verdicts_json(case, g, "dry")
        g["verdicts_json"] = str(vp.relative_to(case)) if vp else None
        if c in conds:
            k = gate_cond(case, c, case / conds[c], band)
            if why:
                k["overall_raw"] = k["overall"]
                k["overall"] = UNDET
            out["cond_gates"][c] = k
            vp = verdicts_json(case, k, "cond")
            k["verdicts_json"] = str(vp.relative_to(case)) if vp else None
        try:
            out["record"][c] = record_window(case, case / b, (case / e) if e else None)
        except SeriesError as ex:
            out["record"][c] = {"error": str(ex)}
    # 差の表 (判定なし)
    ref = reference_dry(case)
    out["reference_dry"] = ref
    pairs = [("U4 (N1 − N0)", "N0", "N1"), ("V5′ (N2 − N1)", "N1", "N2")]
    for label, a, b in pairs:
        if a in out["dry_gates"] and b in out["dry_gates"]:
            out["diffs"][label] = {"gates": diff_table(window_of(out["dry_gates"][a]), window_of(out["dry_gates"][b]))}
            ra, rb = out["record"].get(a, {}), out["record"].get(b, {})
            if "values" in ra and "values" in rb:
                out["diffs"][label]["record"] = diff_table(ra["values"], rb["values"])
        if a in out["cond_gates"] and b in out["cond_gates"]:
            out["cond_diffs"][label] = diff_table(window_of(out["cond_gates"][a]), window_of(out["cond_gates"][b]))
    for c in out["dry_gates"]:
        if "window" in ref:
            out["diffs"][f"{c} − 生産 (run_0147+0149)"] = {"gates": diff_table(ref["window"], window_of(out["dry_gates"][c]))}
            if "record" in ref and "values" in out["record"].get(c, {}):
                out["diffs"][f"{c} − 生産 (run_0147+0149)"]["record"] = diff_table(ref["record"]["values"], out["record"][c]["values"])
    refc = reference_cond(case)
    out["reference_cond"] = refc
    for c in out["cond_gates"]:
        if "window" in refc:
            out["cond_diffs"][f"{c} − 生産 (run_0148)"] = diff_table(refc["window"], window_of(out["cond_gates"][c]))
    out["summary"] = {c: {"dry": out["dry_gates"][c]["overall"], "dry_unmet": out["dry_gates"][c]["unmet"],
                          "cond": out["cond_gates"].get(c, {}).get("overall")} for c in out["dry_gates"]}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    return out


def print_summary(out: dict) -> None:
    for c, g in out["dry_gates"].items():
        pre = out["preconditions"][c]
        print(f"== {c} dry {'+'.join(g['runs'])}: {g['overall']}" + ("" if pre["ok"] else f"  (前提不成立: {pre['failures']})"))
        for i in g["items"]:
            v = i["value"]
            vs = f"{v:.7g}" if isinstance(v, float) else str(v)
            print(f"   {i['status']:4s} {i['key']:34s} {vs:>16s}  [{i['threshold']}]")
        if g.get("next"):
            print("   → " + g["next"])
        k = out["cond_gates"].get(c)
        if k:
            print(f"== {c} 凝縮 {k['run']}: {k['overall']}")
            for i in k["items"]:
                print(f"   {i['status']:4s} {i['key']:34s} {str(i['value']):>16s}  [{i['threshold']}]")
    for label, d in out["diffs"].items():
        print(f"-- 差 {label} (判定なし; 窓の平均の差、幅は両窓の大きい方)")
        for q, v in d.get("gates", {}).items():
            print(f"   {q:22s} {v['diff']:+.3e}  (幅 {v['range_max']:.2e})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("eval")
    e.add_argument("--md-offset", required=True)
    e.add_argument("--set", action="append", default=[], help="N0=<dry>[+<延長>] (3 条件ぶん)")
    e.add_argument("--cond", action="append", default=[], help="N0=<凝縮の run>")
    e.add_argument("--out", default=OUT_JSON)
    e.add_argument("--dry", action="store_true")
    g = sub.add_parser("gate-dry")
    g.add_argument("--md-offset", required=True)
    g.add_argument("--cond-name", required=True, choices=CONDS)
    g.add_argument("--run", required=True)
    g.add_argument("--ext")
    g.add_argument("--r-throat", default=None)
    e.add_argument("--r-throat", default=None, help="r_t [m] (make_ns_n012_problems.py --r-throat と同じ値)")
    a = ap.parse_args(argv)
    md = float(a.md_offset)
    rt = float(a.r_throat) if a.r_throat is not None else None
    if a.cmd == "gate-dry":
        case = HERE
        res = gate_dry(case, a.cond_name, case / a.run, (case / a.ext) if a.ext else None, case / "_band_ab")
        why = run_preconditions(case, case / a.run, "dry", a.cond_name, md, False, provenance_sha(case / REF_PROV_RUN))
        if a.ext:
            why += run_preconditions(case, case / a.ext, "ext", a.cond_name, md, False, provenance_sha(case / REF_PROV_RUN))
        why += r_throat_failures([case / a.run] + ([case / a.ext] if a.ext else []), rt)
        st = UNDET if why else res["overall"]
        print(f"[gate-dry] {a.cond_name} {a.run}{'+' + a.ext if a.ext else ''}: {st}"
              + (f" 前提不成立 {why}" if why else "") + (f" 未達 {res['unmet']}" if res["unmet"] else "")
              + (f" 判定不能 {res['undetermined']}" if res["undetermined"] else ""))
        return 0 if st == PASS else 1
    sets = {}
    for s in a.set:
        c, b, x = parse_set(s)
        sets[c] = (b, x)
    conds = {}
    for s in a.cond:
        c, r, x = parse_set(s)
        if x:
            raise SystemExit(f"--cond {s}: 凝縮は 1 本 (+ は付けない)")
        conds[c] = r
    if sorted(sets) != list(CONDS):
        raise SystemExit(f"--set は N0・N1・N2 の 3 つ (受け取った {sorted(sets)})")
    out = evaluate(HERE, md, sets, conds, a.dry, HERE / a.out, rt)
    print_summary(out)
    print(f"→ {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
