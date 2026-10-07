"""ns_n012_eval.py・make_ns_n012_problems.py・ns_n012.py (純粋な部分)・run_ns_n012.sh (引数の検査) の試験
(plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「NS の 3 条件」)。forge・AWS・変換器は使わない。合成した系列と、一時ディレクトリの
模擬 case (dry 3 本・延長・凝縮 3 本・生産の比較元の CSV と JSON) だけを使う。
  評価器: 全ゲート合格 / 出口コア M・波・オーバーシュート・δ_E/δ_C・出口半径・壁解像・残差 RISING の各未達 / 境界値 /
    準定常 DRIFTING / 非有限 / 欠落・窓のずれ・列の欠損 → 判定不能 / prepare_info と report の出口半径の食い違い → 判定不能 /
    壁解像の許容の違い → 判定不能 / 延長の連結 (窓 80000〜100000)・親の違い → 判定不能 / 前提 (forge の sha256・IC の検査・較正値・
    RUN_RC・段・早期停止・NAN_SCAN・3 条件の固定の検査) の不成立 → 判定不能 / 凝縮の合格・未達 / 差の表 (N1 − N0・N2 − N1・生産との差)、
    check_quasisteady の出力の読み取り、nozzle_report --verdicts で読める JSON。
  問題の生成: 6 本の差が登録どおり、N1・N2 に pw_ramp が無い、較正値のトークンの往復、不正な較正値、改竄の検出。
  ns_n012: 壁の層の対応の記録 (壁が動いた列だけで層を取り違える)、interp_field と同じ転送の式、roe の無い res の拒否、nan-scan、
    同一格子の検査。
  起動スクリプト: 不正な run 名・既存の dir・不正な COND・引数不足で止まる (forge を起動しない)。
usage: python3 test_ns_n012_eval.py
"""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import ns_n012_eval as E  # noqa: E402
import make_ns_n012_problems as MK  # noqa: E402
import ns_n012 as NS  # noqa: E402

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    if not cond:
        fails += 1


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


MD = 3.91e-4
FSHA = "ab" * 32


def write_csv(p, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


SEG_PLATEAU = """  [segment] 判定区間 = main  (80000 行) -> x/residual_history_segment.csv

=== x  [last step 79999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===
  rms_ro      : init=3.90e-07 fin=3.34e-07 drop= 0.1dec flat     <-- STALLED (plateau)

OVERALL: CHECK FAILURES ABOVE
"""
SEG_RISING = SEG_PLATEAU.replace("<-- STALLED (plateau)", "<-- RISING (divergent)")


def dry_values(step, **over):
    """dry の合格する系列 (定常 + わずかなジッタ)。over で量ごとに step の関数を差し替える。"""
    j = 1e-7 * ((step // 5000) % 3 - 1)
    v = {"exitM_A": 5.99990 + j, "exitM_B": 5.99990 + j, "dM_BA": 0.0, "overshoot01": 0.0080 + 10 * j, "wave01": 0.0065 + 10 * j,
         "delta_E": 0.7327, "n_core_A": 12, "n_core_B": 12, "dE_over_dC": 0.99985 + j}
    for k, f in over.items():
        v[k] = f(step)
    return v


def record_values(step, base=0.0):
    return {"mdot_throat_kgs": 100.0 + base, "mdot_median_kgs": 100.0 + base, "mdot_spread_rel": 1e-4, "mdot_exit_kgs": 100.0 + base,
            "x_sonic_eta0": 0.1 + base, "x_sonic_eta0.5": 0.05, "x_sonic_eta0.9": -0.02, "pw_over_pt_x0": 0.55, "deltaE_rt_x0": 0.0014}


def make_problem(case, cond):
    p = case / MK.out_name(cond, "dry")
    if not p.exists():
        p.write_text(f"name: {cond}\n")
    pk = case / MK.out_name(cond, "cond")
    if not pk.exists():
        pk.write_text(f"name: {cond}_cond\n")
    return p, pk


def make_dry(case, name, cond, md=MD, steps_end=80000, over=None, exit_r=0.7750001, rep_r=None, wall=("PASS", 3.5, 5.0),
             seg=SEG_PLATEAU, ic="OK", tags=("S1_soft", "S2_mid", "main"), rc="0", nan="CLEAN", fsha=FSHA, role="dry", parent=None,
             rec_base=0.0, interval=5000, scale=0.0766539):
    d = case / name
    d.mkdir(parents=True)
    p, _ = make_problem(case, cond)
    rec = {"role": role, "condition": cond, "md_offset": md, "dry": False, "problem": p.name, "problem_sha256": sha(p), "k_f": 1.05}
    if role == "ext":
        rec.update(parent=parent, parent_end=80000, ext_steps=20000)
    (d / E.RECORD).write_text(json.dumps(rec))
    (d / "prepare_info.json").write_text(json.dumps({"scale_m": scale, "sizing": {"exit_radius_m": exit_r}}))
    rows = [{"step": s, **dry_values(s, **(over or {}))} for s in range(interval, steps_end + 1, interval)]
    write_csv(d / "quantities_series.csv", rows)
    write_csv(d / E.RECORD_CSV, [{"step": s, **record_values(s, rec_base)} for s in range(interval, steps_end + 1, interval)])
    (d / "report").mkdir()
    (d / "report/report.json").write_text(json.dumps({"metrics": {
        "wall_shape": {"exit_radius_m": exit_r if rep_r is None else rep_r},
        "wall_resolution": {"verdict": wall[0], "over_area_pct": wall[1], "over_area_allow_pct": wall[2], "verdict_line": f"VERDICT: {wall[0]}",
                            "y1p_max": 1.5}}}))
    if seg is not None:
        (d / E.SEG_FILE).write_text(seg)
    if ic is not None and role == "dry":
        (d / "IC_CHECK.json").write_text(json.dumps({"VERDICT": ic}))
    (d / "stage_manifest.json").write_text(json.dumps({"stages": [{"tag": t} for t in tags]}))
    (d / "RUN_RC").write_text(rc + "\n")
    (d / "NAN_SCAN.json").write_text(json.dumps({"VERDICT": nan}))
    (d / "RUN_PROVENANCE.txt").write_text(f"forge_sha256: {fsha}\n")
    return d


def make_cond(case, name, cond, parent, over=None, seg=SEG_PLATEAU, md=MD, scale=0.0766539):
    d = case / name
    d.mkdir(parents=True)
    (d / "prepare_info.json").write_text(json.dumps({"scale_m": scale}))
    _, pk = make_problem(case, cond)
    (d / E.RECORD).write_text(json.dumps({"role": "cond", "condition": cond, "md_offset": md, "dry": False, "problem": pk.name,
                                          "problem_sha256": sha(pk), "parent": parent}))
    rows = []
    for s in range(1000, 18001, 1000):
        j = 1e-6 * ((s // 1000) % 3 - 1)
        v = {"step": s, "onset_x_axis": 57.96, "S_max": 16.86 + j, "exit_g_core": 3.128e-4 + 1e-9 * j, "exit_core_M": 5.98639 + j}
        for k, f in (over or {}).items():
            v[k] = f(s)
        rows.append(v)
    write_csv(d / "cond_series.csv", rows)
    (d / E.SEG_FILE).write_text(seg.replace("80000 行", "18000 行"))
    (d / "stage_manifest.json").write_text(json.dumps({"stages": [{"tag": "main"}]}))
    (d / "RUN_RC").write_text("0\n")
    (d / "NAN_SCAN.json").write_text(json.dumps({"VERDICT": "CLEAN"}))
    (d / "RUN_PROVENANCE.txt").write_text(f"forge_sha256: {FSHA}\n")
    return d


def make_case(tmp, dry_kw=None, names=("run_0165_n0", "run_0166_n1", "run_0167_n2"), set_ok=True):
    case = Path(tmp) / "case"
    case.mkdir()
    (case / "_band_ab").mkdir()
    ref = case / E.REF_PROV_RUN
    ref.mkdir()
    (ref / "RUN_PROVENANCE.txt").write_text(f"forge_sha256: {FSHA}\n")
    dry_kw = dry_kw or {}
    for i, (c, n) in enumerate(zip(E.CONDS, names)):
        make_dry(case, n, c, rec_base=0.01 * i, **dry_kw.get(c, {}))
    (case / "_band_ab/ns_n012_prep_check.json").write_text(json.dumps(
        {"VERDICT": "OK" if set_ok else "NG", "runs": {c: {"run": n} for c, n in zip(E.CONDS, names)}}))
    return case


def run_eval(case, sets=None, conds=None, r_throat=None):
    sets = sets or {"N0": ("run_0165_n0", None), "N1": ("run_0166_n1", None), "N2": ("run_0167_n2", None)}
    return E.evaluate(case, MD, sets, conds or {}, False, case / "_band_ab/ns_n012_eval.json", r_throat)


def status_of(out, cond, key):
    for i in out["dry_gates"][cond]["items"]:
        if i["key"] == key:
            return i["status"]
    return None


# ===== 1. 評価器 =====================================================================================================
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    out = run_eval(case)
    check("全ゲート合格 → 3 条件とも合格", all(out["dry_gates"][c]["overall"] == E.PASS for c in E.CONDS))
    check("窓は 60000〜80000 の 5 枚", out["dry_gates"]["N0"]["window_steps"] == [60000, 65000, 70000, 75000, 80000])
    check("準定常 4 量が STEADY", all(v["verdict"] == "STEADY" for v in out["dry_gates"]["N0"]["quasisteady"].values()))
    d = out["diffs"]["U4 (N1 − N0)"]["record"]["mdot_median_kgs"]["diff"]
    check("差の表 (記録): N1 − N0 の質量流量の差 = 0.01", abs(d - 0.01) < 1e-12)
    check("差の表 (ゲートの量) に 4 量", sorted(out["diffs"]["V5′ (N2 − N1)"]["gates"]) == sorted(E.QS_DRY.values()))
    vj = out["dry_gates"]["N1"]["verdicts_json"]
    sys.path.insert(0, str(C.parents[1] / "design"))
    try:
        from forge_design.report.nozzle_report import load_verdicts
        lv = load_verdicts(case / vj)
        check("verdicts JSON は nozzle_report --verdicts で読める", set(lv["items"]) == set(E.QS_DRY.values()))
    except Exception as e:  # noqa: BLE001
        check(f"verdicts JSON は nozzle_report --verdicts で読める ({type(e).__name__}: {e})", False)
    check("生産の比較元が無ければ差の表に生産は出ない (エラーを記録)", "window_error" in out["reference_dry"])

cases = [
    ("出口コア M 5.9969 → 未達", {"N1": {"over": {"exitM_A": lambda s: 5.9969}}}, "N1", "出口コア M", E.FAIL),
    ("出口コア M 6.0031 → 未達", {"N1": {"over": {"exitM_A": lambda s: 6.0031}}}, "N1", "出口コア M", E.FAIL),
    ("出口コア M 5.997 (下限ちょうど) → 合格", {"N1": {"over": {"exitM_A": lambda s: 5.997}}}, "N1", "出口コア M", E.PASS),
    ("出口コア M 5.9985 (旧ゲートでは未達、± 0.05 % では合格) → 合格", {"N1": {"over": {"exitM_A": lambda s: 5.9985}}}, "N1", "出口コア M", E.PASS),
    ("波 0.011 → 未達", {"N0": {"over": {"wave01": lambda s: 0.011}}}, "N0", "Mach 波 η0.1 [%]", E.FAIL),
    ("オーバーシュート 0.036 → 未達", {"N2": {"over": {"overshoot01": lambda s: 0.036}}}, "N2", "オーバーシュート η0.1 [%]", E.FAIL),
    ("δ_E/δ_C 1.006 → 未達", {"N0": {"over": {"dE_over_dC": lambda s: 1.006}}}, "N0", "δ_E/δ_C (x_F)", E.FAIL),
    ("δ_E/δ_C 0.996 → 合格", {"N0": {"over": {"dE_over_dC": lambda s: 0.996}}}, "N0", "δ_E/δ_C (x_F)", E.PASS),
    ("出口半径 0.7752 → 未達", {"N1": {"exit_r": 0.7752}}, "N1", "出口半径 [m]", E.FAIL),
    ("出口半径: prepare_info と report が食い違う → 判定不能", {"N1": {"rep_r": 0.77501}}, "N1", "出口半径 [m]", E.UNDET),
    ("壁解像 FAIL → 未達", {"N2": {"wall": ("FAIL", 7.0, 5.0)}}, "N2", "壁解像 (y1+ > 1 の面積 %)", E.FAIL),
    ("壁解像の許容が 5 % でない → 判定不能", {"N2": {"wall": ("PASS", 3.0, 10.0)}}, "N2", "壁解像 (y1+ > 1 の面積 %)", E.UNDET),
    ("残差 RISING → 未達", {"N0": {"seg": SEG_RISING}}, "N0", "残差 RISING なし (check_convergence --segment)", E.FAIL),
    ("残差の判定ファイルが無い → 判定不能", {"N0": {"seg": None}}, "N0", "残差 RISING なし (check_convergence --segment)", E.UNDET),
    ("オーバーシュートが窓で単調に増える (窓の幅 0.006 % > 0.0035 %) → DRIFTING → 未達", {"N0": {"over": {"overshoot01": lambda s: 0.008 + 3e-7 * s}}}, "N0",
     "準定常 overshoot_eta0.1", E.FAIL),
    ("零に近いオーバーシュートが窓で動く (DRIFTING、窓の幅 0.0005 % ≤ 0.0035 %) → 絶対の幅の条件で合格 (2026-10-07 ユーザ決定)",
     {"N0": {"over": {"overshoot01": lambda s: 0.0005 + 2.5e-8 * s}}}, "N0", "準定常 overshoot_eta0.1", E.PASS),
    ("窓に非有限 → 波の値は判定不能", {"N0": {"over": {"wave01": lambda s: float("nan") if s == 70000 else 0.0065}}}, "N0",
     "Mach 波 η0.1 [%]", E.UNDET),
    ("窓に非有限 → 準定常 NONFINITE → 未達", {"N0": {"over": {"wave01": lambda s: float("nan") if s == 70000 else 0.0065}}}, "N0",
     "準定常 wave_eta0.1", E.FAIL),
    ("本段が 75000 で終わる (窓のずれ) → 判定不能", {"N2": {"steps_end": 75000}}, "N2", "時系列 (quantities_series.csv)", E.UNDET),
]
for label, kw, cond, key, want in cases:
    with tempfile.TemporaryDirectory() as tmp:
        case = make_case(tmp, kw)
        out = run_eval(case)
        got = status_of(out, cond, key)
        check(f"{label} (得た {got})", got == want)
        if want == E.FAIL and key != "準定常 wave_eta0.1":
            check(f"  {label}: 総合は未達・延長の候補", out["dry_gates"][cond]["overall"] == E.FAIL and "延長" in out["dry_gates"][cond].get("next", ""))

# 欠落・重複の step
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    p = case / "run_0165_n0/quantities_series.csv"
    rows = list(csv.DictReader(open(p)))
    write_csv(p, rows[:5] + rows[6:])
    out = run_eval(case)
    check("step の欠落 → 判定不能", out["dry_gates"]["N0"]["overall"] == E.UNDET)
    rows2 = rows[:-1] + [rows[-2]]
    write_csv(p, rows2)
    out = run_eval(case)
    check("step の重複・逆順 → 判定不能", out["dry_gates"]["N0"]["overall"] == E.UNDET)

# 前提の不成立 → 判定不能
pre = [("forge の sha256 が run_0147 と違う", {"N1": {"fsha": "cd" * 32}}, "N1"),
       ("IC の検査が NG", {"N2": {"ic": "NG"}}, "N2"),
       ("較正値が違う", {"N0": {"md": 3.92e-4}}, "N0"),
       ("RUN_RC 1", {"N0": {"rc": "1"}}, "N0"),
       ("段の記録が本段だけ (段階起動なし)", {"N1": {"tags": ("main",)}}, "N1"),
       ("NAN_SCAN が NAN", {"N2": {"nan": "NAN"}}, "N2")]
for label, kw, cond in pre:
    with tempfile.TemporaryDirectory() as tmp:
        case = make_case(tmp, kw)
        out = run_eval(case)
        check(f"前提: {label} → 判定不能", out["dry_gates"][cond]["overall"] == E.UNDET and not out["preconditions"][cond]["ok"])
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    (case / "run_0166_n1/EARLY_STOP.txt").write_text("x")
    out = run_eval(case)
    check("前提: 早期停止 → 判定不能", out["dry_gates"]["N1"]["overall"] == E.UNDET)
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp, set_ok=False)
    out = run_eval(case)
    check("前提: 3 条件の固定の検査が NG → 3 条件とも判定不能", all(out["dry_gates"][c]["overall"] == E.UNDET for c in E.CONDS))
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    (case / MK.out_name("N1", "dry")).write_text("name: changed\n")
    out = run_eval(case)
    check("前提: 問題 YAML が準備の後に変わった → 判定不能", out["dry_gates"]["N1"]["overall"] == E.UNDET)

# 延長の連結
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp, {"N0": {"over": {"overshoot01": lambda s: 0.008 + (2e-7 * s if s <= 80000 else 0.016)}}})
    make_dry(case, "run_0171_n0_ext", "N0", steps_end=20000, role="ext", parent="run_0165_n0",
             over={"overshoot01": lambda s: 0.024 + 1e-9 * (s % 7)}, ic=None, tags=("main",))
    out0 = run_eval(case)
    check("延長前: オーバーシュートが DRIFTING で未達", status_of(out0, "N0", "準定常 overshoot_eta0.1") == E.FAIL)
    out = run_eval(case, sets={"N0": ("run_0165_n0", "run_0171_n0_ext"), "N1": ("run_0166_n1", None), "N2": ("run_0167_n2", None)})
    g = out["dry_gates"]["N0"]
    check("延長: 窓は 80000〜100000", g["window_steps"] == [80000, 85000, 90000, 95000, 100000])
    check(f"延長: 連結した系列でオーバーシュートが STEADY ({status_of(out, 'N0', '準定常 overshoot_eta0.1')})",
          status_of(out, "N0", "準定常 overshoot_eta0.1") == E.PASS)
    check("延長: 残差は本段と延長の両方", len(g["convergence"]) == 2)
    check("延長: 記録の連結", out["record"]["N0"]["window_steps"] == [80000, 85000, 90000, 95000, 100000])
    shutil.rmtree(case / "run_0171_n0_ext")
    make_dry(case, "run_0171_n0_ext", "N0", steps_end=20000, role="ext", parent="run_0999_other", ic=None, tags=("main",))
    out = run_eval(case, sets={"N0": ("run_0165_n0", "run_0171_n0_ext"), "N1": ("run_0166_n1", None), "N2": ("run_0167_n2", None)})
    check("延長: 親が違う → 判定不能", out["dry_gates"]["N0"]["overall"] == E.UNDET)
    try:
        E.join_rows([{"step": s, "a": 1} for s in range(5000, 80001, 5000)], [{"step": s, "a": 1} for s in range(0, 20001, 5000)], 80000)
        check("延長: local step 0 を含む延長は拒否", False)
    except E.SeriesError:
        check("延長: local step 0 を含む延長は拒否", True)

# 凝縮
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    make_cond(case, "run_0168_k0", "N0", "run_0165_n0")
    make_cond(case, "run_0169_k1", "N1", "run_0166_n1", over={"S_max": lambda s: 16.0 + 1e-3 * s})
    make_cond(case, "run_0170_k2", "N2", "run_0167_n2", seg=SEG_RISING)
    out = run_eval(case, conds={"N0": "run_0168_k0", "N1": "run_0169_k1", "N2": "run_0170_k2"})
    check("凝縮: 4 量 STEADY・RISING なし → 合格", out["cond_gates"]["N0"]["overall"] == E.PASS)
    check("凝縮: S_max が DRIFTING → 未達", out["cond_gates"]["N1"]["overall"] == E.FAIL
          and "準定常 cond_S_max" in out["cond_gates"]["N1"]["unmet"])
    check("凝縮: 残差 RISING → 未達", out["cond_gates"]["N2"]["overall"] == E.FAIL)
    check("凝縮: 窓は 14000〜18000", out["cond_gates"]["N0"]["window_steps"] == [14000, 15000, 16000, 17000, 18000])
    check("凝縮: 差の表 (U4)", "U4 (N1 − N0)" in out["cond_diffs"])
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    make_cond(case, "run_0168_k0", "N0", "run_0166_n1")
    out = run_eval(case, conds={"N0": "run_0168_k0"})
    check("凝縮: 親が同じ条件の dry でない → 判定不能", out["cond_gates"]["N0"]["overall"] == E.UNDET)

# 生産の比較元 (run_0147 + 延長 run_0149) との差
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    a, b = case / E.REF_DRY[0], case / E.REF_DRY[1]
    write_csv(a / "quantities_series.csv", [{"step": s, **dry_values(s)} for s in range(5000, 60001, 5000)])
    b.mkdir()
    write_csv(b / "quantities_series_extonly.csv", [{"step": s, **dry_values(s, exitM_A=lambda s: 5.99880)} for s in range(5000, 20001, 5000)])
    (case / E.REF_RECORD_DIR).mkdir(parents=True)
    write_csv(case / E.REF_RECORD_DIR / f"{E.REF_DRY[0]}_record.csv", [{"step": s, **record_values(s, -1.0)} for s in range(5000, 60001, 5000)])
    write_csv(case / E.REF_RECORD_DIR / f"{E.REF_DRY[1]}_record.csv", [{"step": s, **record_values(s, -1.0)} for s in range(5000, 20001, 5000)])
    rc_ = case / E.REF_COND
    rc_.mkdir()
    shutil.copy(make_cond(case, "run_0168_k0", "N0", "run_0165_n0") / "cond_series.csv", rc_ / "cond_series.csv")
    out = run_eval(case, conds={"N0": "run_0168_k0"})
    dd = out["diffs"]["N0 − 生産 (run_0147+0149)"]
    w = [5.99880, 5.99880, 5.99880, 5.99880, 5.99880]
    check("生産との差: 出口コア M (窓は 60000 = run_0147 の最終 + run_0149 の 4 枚)",
          abs(dd["gates"]["exit_core_M"]["a"] - (5.99990 + 0 * 1e-7 + 4 * 5.99880) / 5) < 2e-7)
    check("生産との差: 記録 (質量流量 +1.0)", abs(dd["record"]["mdot_median_kgs"]["diff"] - 1.0) < 1e-9)
    check("凝縮の生産との差 (run_0148)", "N0 − 生産 (run_0148)" in out["cond_diffs"])

# check_quasisteady の出力の読み取り
txt = """
=== x.csv  [16 rows, steps 5000..80000]  -> DRIFTING ===
  exitM_A         : tail mean=6  drift=0.0%/tail  fluct=0.0%                      STEADY
  overshoot01     : tail mean=0.0161  drift=12.3%/tail  fluct=20.1%  (extremum at tail-end) DRIFTING
  wave01          : tail mean=0.006  drift=3.0%/tail  fluct=1.0%  (still trending at tail-end) TRANSIENT-UNSETTLED
"""
r = E.parse_quasisteady(txt, ["exitM_A", "overshoot01", "wave01", "dE_over_dC"])
check("check_quasisteady の読み取り (STEADY・DRIFTING・TRANSIENT-UNSETTLED・行なし)",
      [r[c]["verdict"] for c in ("exitM_A", "overshoot01", "wave01", "dE_over_dC")] == ["STEADY", "DRIFTING", "TRANSIENT-UNSETTLED", None])
r = E.parse_quasisteady("=== x.csv  -> ERROR: missing column(s) ['wave01'] ===", ["wave01"])
check("check_quasisteady の ERROR → verdict なし", r["wave01"]["verdict"] is None)
try:
    E.parse_set("N3=run_x")
    check("parse_set: N3 は拒否", False)
except ValueError:
    check("parse_set: N3 は拒否", True)
check("parse_set: 延長つき", E.parse_set("N1=run_0166_a+run_0171_b") == ("N1", "run_0166_a", "run_0171_b"))

# ===== 2. 問題の生成 =================================================================================================
texts, rec = MK.build(MD, C)
check("生成: 6 本", sorted(texts) == sorted(MK.out_name(c, k) for c in MK.CONDS for k in ("dry", "cond")))
docs = {n: MK.yaml_load(t) for n, t in texts.items()}
base = MK.yaml_load((C / MK.BASE["dry"]).read_text())
for c in MK.CONDS:
    g = docs[MK.out_name(c, "dry")]["geometry"]
    check(f"生成 {c}: 較正値・pw_upstream・MOC", g["Md_moc_offset"] == MD and g["pw_upstream"] == MK.SPEC[c]["pw_upstream"]
          and g["moc_axis_limit"] == MK.SPEC[c]["moc_axis_limit"] and g["moc_corrector"] == MK.SPEC[c]["moc_corrector"])
    check(f"生成 {c}: pw_ramp は N0 だけ", ("pw_ramp" in g) == (c == "N0"))
    dd = MK.diff_paths(base, docs[MK.out_name(c, "dry")])
    check(f"生成 {c}: 元との差は許したパスだけ", set(dd) <= MK.ALLOWED)
    check(f"生成 {c}: 格子・r_t・k_f は元のまま", docs[MK.out_name(c, "dry")]["mesh"] == base["mesh"]
          and docs[MK.out_name(c, "dry")]["spec"] == base["spec"]
          and docs[MK.out_name(c, "dry")]["deltastar_initializer"] == base["deltastar_initializer"])
for v in (3.91e-4, 1e-5, -2.5e-6, 0.0, 1.234567890123e-4):
    t = MK.offset_token(v)
    check(f"較正値のトークン {t!r} は float として往復", MK.yaml_load(f"v: {t}\n")["v"] == v and isinstance(MK.yaml_load(f"v: {t}\n")["v"], float))
for bad in ("abc", "nan", "inf", "0.5"):
    try:
        MK.parse_offset(bad)
        check(f"不正な較正値 {bad!r} は拒否", False)
    except ValueError:
        check(f"不正な較正値 {bad!r} は拒否", True)
t = texts[MK.out_name("N1", "dry")].replace("  ni: 2000", "  ni: 2001")
try:
    MK.check_generated(t, base, "N1", "dry", MD)
    check("生成物の改竄 (格子) を検出", False)
except ValueError:
    check("生成物の改竄 (格子) を検出", True)
t = texts[MK.out_name("N2", "dry")].replace("moc_corrector: converge", "moc_corrector: fixed2")
try:
    MK.check_generated(t, base, "N2", "dry", MD)
    check("生成物の改竄 (MOC) を検出", False)
except ValueError:
    check("生成物の改竄 (MOC) を検出", True)
with tempfile.TemporaryDirectory() as tmp:
    rc = MK.main(["--md-offset", "3.91e-4", "--out-dir", tmp])
    rc2 = MK.main(["--md-offset", "0.000391", "--out-dir", tmp, "--check"])
    rc3 = MK.main(["--md-offset", "0.000392", "--out-dir", tmp])
    rc4 = MK.main(["--md-offset", "0.000392", "--out-dir", tmp, "--check"])
    check("生成 → 同じ値で --check は OK、違う値の上書きは拒否、違う値の --check は NG", (rc, rc2, rc3, rc4) == (0, 0, 2, 2))

# ===== 3. ns_n012 の純粋な部分 ==========================================================================================
ni, nj = 6, 8
x = np.repeat(np.linspace(0.0, 5.0, ni), nj)
eta = np.tile(1.0 - (np.linspace(1.0, 0.0, nj) ** 3), ni)      # 壁に寄せた配点 (j = nj−1 が壁)
rw_s = np.full(ni, 1.0)
rw_d = rw_s.copy()
rw_d[3:] -= 0.05                                                # 列 3〜5 で壁が動く (第 1 セル厚より大きい)
cs = np.c_[x, eta * np.repeat(rw_s, nj), np.zeros(ni * nj)]
cd = np.c_[x, eta * np.repeat(rw_d, nj), np.zeros(ni * nj)]
from scipy.spatial import cKDTree  # noqa: E402
_, idx = cKDTree(cs[:, :2]).query(cd[:, :2])
rec = NS.wall_layer_record(cd, cs, idx, ni, nj, ni, nj, 1.0)
check("壁の層: 壁の動かない列では取り違えなし", rec["near_wall"]["n_in_columns_wall_moved_le_half_h1"] == 0)
check("壁の層: 壁が動いた列で取り違えを記録", rec["near_wall"]["n_in_columns_wall_moved_gt_half_h1"] > 0
      and rec["layers"]["0"]["n_donor_other_layer"] == 3)
check("壁の層: 壁の移動の最大 0.05 (= 5e4 µm)", abs(rec["wall_shift"]["max_abs_um"] - 5e4) < 1e-6)

with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp) / "res.h5"
    rng = np.random.default_rng(1)
    with h5py.File(p, "w") as f:
        for k in ("ro", "P", "Ux", "Uy", "Uz", "k", "omega", "roe", "Y0", "Y1"):
            f.create_dataset(f"VALUE/{k}", data=rng.random(7).astype(np.float32) + 0.5)
    fl = NS.interp_like_fields(p)
    with h5py.File(p) as f:
        ro, ux = np.array(f["VALUE/ro"]), np.array(f["VALUE/Ux"])
    check("interp_field と同じ式: roUx = ro·Ux (float32)、roY0・roY1・roK・roOmega", np.array_equal(fl["roUx"], ro * ux)
          and fl["roUx"].dtype == np.float32 and {"roY0", "roY1", "roK", "roOmega", "roe"} <= set(fl))
    with h5py.File(p, "r+") as f:
        del f["VALUE/roe"]
    try:
        NS.interp_like_fields(p)
        check("roe の無い res は拒否 (CPG 式の組み直しを使わない)", False)
    except ValueError:
        check("roe の無い res は拒否 (CPG 式の組み直しを使わない)", True)

with tempfile.TemporaryDirectory() as tmp:
    run = Path(tmp) / "run_0999_x"
    run.mkdir()
    write_csv(run / "residual_history_S1_soft.csv", [{"step": s, "rms_ro": 1e-3, "rms_dq_ro": "nan"} for s in range(5)])
    write_csv(run / "residual_history_main.csv", [{"step": s, "rms_ro": (1e-3 if s != 3 else "nan"), "rms_dq_ro": 0} for s in range(5)])
    with h5py.File(run / "res_80000.h5", "w") as f:
        f.create_dataset("VALUE/ro", data=np.ones(4)); f.create_dataset("VALUE/T", data=np.ones(4)); f.create_dataset("VALUE/P", data=np.ones(4))
    o = NS.nan_scan(run)
    check("nan-scan: rms_dq_* は見ない・本段の NaN の最初の行", o["VERDICT"] == "NAN" and o["first_nonfinite"]["file"] == "residual_history_main.csv"
          and o["first_nonfinite"]["step"] == "3")
    (run / "residual_history_main.csv").unlink()
    o = NS.nan_scan(run)
    check("nan-scan: NaN なし → CLEAN", o["VERDICT"] == "CLEAN")
    with h5py.File(run / "res_80000.h5", "r+") as f:
        f["VALUE/T"][1] = -1.0
    o = NS.nan_scan(run)
    check("nan-scan: 最終 res の T ≤ 0 → NAN", o["VERDICT"] == "NAN" and "T" in o["final_res"]["nonpositive"])


def mini_mesh(p, coord, conne):
    with h5py.File(p, "w") as f:
        f.create_dataset("MESH/COORD", data=coord)
        f.create_dataset("MESH/CONNE", data=conne)
        f.create_dataset("BCONDS/1/vizBfaceNodes", data=np.array([0, 1]))


with tempfile.TemporaryDirectory() as tmp:
    a, b = Path(tmp) / "a.h5", Path(tmp) / "b.h5"
    mini_mesh(a, np.arange(9.0), np.array([1, 2, 3]))
    mini_mesh(b, np.arange(9.0), np.array([1, 2, 3]))
    check("同一格子の検査: 同じなら空", NS.same_mesh(a, b) == [])
    mini_mesh(b, np.arange(9.0) + 1e-9, np.array([1, 2, 3]))
    check("同一格子の検査: 座標が違えば検出", any("COORD" in w for w in NS.same_mesh(a, b)))

# ===== 3b. r_t の解き直し (--r-throat; plan §6 U4「r_t の解き直し」) =====================================================
RT = 0.07671234567891234
texts_rt, rec_rt = MK.build(MD, C, RT)
base_c = MK.yaml_load((C / MK.BASE["cond"]).read_text())
docs_rt = {n: MK.yaml_load(t) for n, t in texts_rt.items()}
check("r_t: 6 本すべての spec.r_throat が指定値の float", all(type(d["spec"]["r_throat"]) is float and d["spec"]["r_throat"] == RT
                                                           for d in docs_rt.values()))
for c in MK.CONDS:
    dd = set(MK.diff_paths(base, docs_rt[MK.out_name(c, "dry")])) | set(MK.diff_paths(base_c, docs_rt[MK.out_name(c, "cond")]))
    check(f"r_t {c}: 元との差に spec.r_throat が入り、許したパスだけ", MK.RT_PATH in dd and dd <= MK.ALLOWED | {MK.RT_PATH})
check("r_t: 記録に r_throat・トークン・元の値", rec_rt["r_throat"] == RT and rec_rt["r_throat_token"] == repr(RT)
      and rec_rt["base_r_throat"] == base["spec"]["r_throat"])
texts0, rec0 = MK.build(MD, C)
check("r_t 省略: 6 本の r_throat は元のまま・記録は None", all(MK.yaml_load(t)["spec"]["r_throat"] == base["spec"]["r_throat"]
                                                         for t in texts0.values()) and rec0["r_throat"] is None)
for v in (RT, 0.0766539, 0.05, 1.234e-2):
    t = MK.float_token(v, "r_throat")
    check(f"r_t のトークン {t!r} は float として往復", MK.yaml_load(f"v: {t}\n")["v"] == v and type(MK.yaml_load(f"v: {t}\n")["v"]) is float)
for bad in ("abc", "76.7", "-0.0766", "nan", "0", "inf"):
    try:
        MK.parse_r_throat(bad)
        check(f"不正な r_t {bad!r} は拒否", False)
    except ValueError:
        check(f"不正な r_t {bad!r} は拒否", True)
t = texts_rt[MK.out_name("N1", "cond")].replace(f"r_throat: {RT!r}", "r_throat: 0.0767")
try:
    MK.check_generated(t, base_c, "N1", "cond", MD, RT)
    check("r_t: 1 本だけ違う r_t を検出", False)
except ValueError:
    check("r_t: 1 本だけ違う r_t を検出", True)
t = texts0[MK.out_name("N0", "dry")].replace("r_throat: 0.0766539", "r_throat: 0.0767")
try:
    MK.check_generated(t, base, "N0", "dry", MD)
    check("r_t 省略時: r_throat の書き換えを検出 (差を許さない)", False)
except ValueError:
    check("r_t 省略時: r_throat の書き換えを検出 (差を許さない)", True)
with tempfile.TemporaryDirectory() as tmp:
    rcs = (MK.main(["--md-offset", repr(MD), "--r-throat", repr(RT), "--out-dir", tmp]),
           MK.main(["--md-offset", repr(MD), "--r-throat", repr(RT), "--out-dir", tmp, "--check"]),
           MK.main(["--md-offset", repr(MD), "--r-throat", "0.0767", "--out-dir", tmp, "--check"]),
           MK.main(["--md-offset", repr(MD), "--out-dir", tmp, "--check"]),
           MK.main(["--md-offset", repr(MD), "--out-dir", tmp]))
    check(f"r_t: 生成 → 同じ r_t の --check OK、違う r_t・省略の --check は NG、省略での上書きは拒否 ({rcs})", rcs == (0, 0, 2, 2, 2))
    rj = json.loads((Path(tmp) / MK.RECORD).read_text())
    check("r_t: 記録の実効値 (runner の読み取り) が 6 本とも指定値", all(v["r_throat"] == RT for v in rj["effective"].values()))
with tempfile.TemporaryDirectory() as tmp:
    rcs = (MK.main(["--md-offset", repr(MD), "--out-dir", tmp]), MK.main(["--md-offset", repr(MD), "--r-throat", repr(RT), "--out-dir", tmp, "--check"]))
    check("r_t: 元の r_t で作った 6 本に --r-throat の --check は NG", rcs == (0, 2))

# verify-set の r_t (3 条件で同じ・指定値と一致・問題の記録と一致)。乾式の印つきの模擬 run で (forge・変換器は使わない)


def fake_prep(root, name, cond, scale):
    d = Path(root) / name
    d.mkdir(parents=True)
    sp = MK.SPEC[cond]
    (d / NS.RECORD).write_text(json.dumps({"condition": cond, "role": "dry", "dry": True, "ic": {"src_run": "x"}, "stages": "full",
                                            "main_steps": 80000, "out_interval": 5000, "cfl_main": 1.0, "implicit_relax": 0.7,
                                            "k_f": 1.05, "md_offset": MD}))
    (d / "prepare_info.json").write_text(json.dumps({"DRY": True, "pw_upstream": {"value": sp["pw_upstream"]},
                                                     "moc": {"axis_limit": sp["moc_axis_limit"], "corrector": sp["moc_corrector"]},
                                                     "Md_moc_offset": MD, "scale_m": scale, "initializer": {"cf_scale": 1.05},
                                                     "mesh": {"ni": 2000, "nj": 97}}))
    (d / NS.IC_CHECK).write_text(json.dumps({"VERDICT": "OK"}))
    for f in ("solverConfig.yaml", "bcondConfig.yaml", "species_meta.yaml", "probe.yaml"):
        (d / f).write_text("a: 1\n")
    return d


def vset(scales, rt_arg, rec_rt_val):
    with tempfile.TemporaryDirectory() as tmp:
        runs = [fake_prep(tmp, f"run_990{i}_x", c, sc) for i, (c, sc) in enumerate(zip(MK.CONDS, scales))]
        recp = Path(tmp) / "rec.json"
        recp.write_text(json.dumps({"r_throat": rec_rt_val, "base_r_throat": 0.0766539}))
        return NS.verify_set(runs, repr(MD), True, Path(tmp) / "out.json", rt_arg, recp)


check("verify-set: 3 本とも指定の r_t・記録も同じ → OK", vset([RT] * 3, repr(RT), RT)["VERDICT"] == "OK")
check("verify-set: r_t 省略・3 本とも元の r_t → OK", vset([0.0766539] * 3, None, None)["VERDICT"] == "OK")
check("verify-set: 1 本だけ r_t が違う → NG", vset([RT, RT, 0.0767], repr(RT), RT)["VERDICT"] == "NG")
check("verify-set: 3 本同じだが指定値と違う → NG", vset([0.0767] * 3, repr(RT), RT)["VERDICT"] == "NG")
check("verify-set: 問題の記録の r_t が指定値と違う → NG", vset([RT] * 3, repr(RT), None)["VERDICT"] == "NG")
check("verify-set: r_t 省略なのに 3 本が元の r_t でない → NG", vset([RT] * 3, None, None)["VERDICT"] == "NG")

# 評価器の前提の r_t
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp, {c: {"scale": RT} for c in E.CONDS})
    make_cond(case, "run_0168_k0", "N0", "run_0165_n0", scale=RT)
    out = run_eval(case, conds={"N0": "run_0168_k0"}, r_throat=RT)
    check("評価器: 6 本の r_t が指定値と同じ → 合格", all(out["dry_gates"][c]["overall"] == E.PASS for c in E.CONDS)
          and out["cond_gates"]["N0"]["overall"] == E.PASS)
    out = run_eval(case, conds={"N0": "run_0168_k0"}, r_throat=0.0766539)
    check("評価器: r_t が指定値と違う → 全条件判定不能", all(out["dry_gates"][c]["overall"] == E.UNDET for c in E.CONDS))
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp, {"N2": {"scale": RT}})
    out = run_eval(case)
    check("評価器: r_t 省略でも 3 本の r_t が違えば全条件判定不能", all(out["dry_gates"][c]["overall"] == E.UNDET for c in E.CONDS)
          and out["r_throat"]["failures"])
with tempfile.TemporaryDirectory() as tmp:
    case = make_case(tmp)
    make_cond(case, "run_0168_k0", "N0", "run_0165_n0", scale=RT)
    out = run_eval(case, conds={"N0": "run_0168_k0"})
    check("評価器: 凝縮だけ r_t が違う → 判定不能", out["cond_gates"]["N0"]["overall"] == E.UNDET)

# ===== 4. 起動スクリプトの引数の検査 (forge を起動しない) ==================================================================
sh = C / "run_ns_n012.sh"
dry_dir_before = (C / "_dry_ns_n012").exists()
env = dict(os.environ, DRY="1", FORGE_BIN="/nonexistent")


def sh_run(*args, **kw):
    e = dict(env, **kw)
    return subprocess.run(["bash", str(sh), *args], capture_output=True, text=True, env=e, cwd=str(C))


r = sh_run("main", "3.91e-4", "bad_name", "run_9001_b", "run_9002_c", "run_9003_d", "run_9004_e", "run_9005_f")
check("起動: 不正な run 名で止まる", r.returncode == 2 and "run_NNNN" in r.stdout)
r = sh_run("main", "3.91e-4", "run_9001_a", "run_9001_a", "run_9002_c", "run_9003_d", "run_9004_e", "run_9005_f")
check("起動: 重複した run 名で止まる", r.returncode == 2 and "重複" in r.stdout)
r = sh_run("main", "3.91e-4", "run_9001_a")
check("起動: 引数不足で止まる", r.returncode == 2)
r = sh_run("main", "3.91e-4", "run_9001_a", "run_9002_b", "run_9003_c", "run_9004_d", "run_9005_e", "run_9006_f", COND="maybe")
check("起動: 不正な COND で止まる", r.returncode == 2 and "COND" in r.stdout)
r = sh_run("ext", "3.91e-4", "N5", "run_9001_a", "run_9002_b")
check("起動: 不正な条件で止まる", r.returncode == 2 and "N0 / N1 / N2" in r.stdout)
r = sh_run("bogus")
check("起動: 不明なモードは使い方を出して止まる", r.returncode == 2 and "usage" in r.stdout)
check("起動: どの検査でも _dry_ns_n012 を作っていない", (C / "_dry_ns_n012").exists() == dry_dir_before)

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
