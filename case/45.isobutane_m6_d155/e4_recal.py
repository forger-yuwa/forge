"""plan verification-case45-euler-total-enthalpy §5.1 #5b・§6 E4 (2026-10-07 登録 99431498) の run 準備・実行: 出口較正のやり直し。
  段 1 (d0) = problem_d155_euler_e4_recal_d0.yaml (Md_moc_offset δ₀ = +3.770e-4) → run_0163_euler_e4_recal_d0
  段 2 (d1) = problem_d155_euler_e4_recal_d1.yaml (δ₁ = 段 1 の評価器が出した値; make-d1 が d0 から Md_moc_offset と name だけを変えて作る)
             → run_0164_euler_e4_recal_d1。段 2 は主セッションが段 1 の判定を見てから起動する (段 1 が「更新」のときだけ)。
固定 (登録): 単調壁 [0, 1.5]・legacy + fixed2 (MOC のキー無し)・凍結の初期線 run_0062 (res_6000)・r_t・BC・熱物性・バイナリ。格子は E3 の
  mesh_euler (2000 × 97・全域 0.005)。IC は新しい格子の上の等エントロピー IC (RA.prepare; G1 の場を移さない)。
準備 (prep): 問題の検査 (check-problem) → RA.prepare → 検査 (step 数の同期・MOC legacy/fixed2・mono_r2・格子の採用元と実効値・
  Md_moc_offset・凍結の初期線・r_t・壁の証拠 M4・メッシュ品質 PASS) → 起動前の IC の検査 (全節点 |T0 − 1600| ≤ 1 K; E2 と同じ
  euler_t0_e2.ic_check) → 実効のメッシュの記録 (E4_MESH.json・E4_MESH_sections.csv; E2 と同じ euler_t0_e2.mesh_record) →
  段 2 は段 1 の prep (run_0163 の E4_PREP.json) と bcondConfig・solverConfig・probe・species_meta が同一か → E4_PREP.json。
  不成立はすべて例外で止める (投入スクリプトは forge を起動しない)。
実行 (run): runner_axismach.run_staged(stages="soft", mid_stage=False, cfl_main=2) — soft (1 次・cfl 0.5・3000 step) → 本段 (2 次・
  cfl 2・implicitRelax 0.7・54000 step・1000 ごと出力)。本段の step 数・出力間隔は評価器の定数 (e4_recal_eval.MAIN_NSTEPS・
  OUT_INTERVAL) を読み、prep の後と起動の直前に solverConfig.yaml の実効値と照合する (食い違えば止める)。soft 段の出力は
  runner が段の終わりに消す前に `_soft_stage/` に退避する (E2 の install_soft_saver と同じ方式)。design/ は変えない (import だけ)。
usage: [CASE_RUNS=<run_0062・run_0114・run_0163 のある case dir>] python3 e4_recal.py check-problem d0|d1 [--delta X]
       python3 e4_recal.py make-d1 --delta X        (段 1 の評価 _band_ab/e4_recal_eval.json が「更新」で δ₁ の repr が X のときだけ)
       python3 e4_recal.py prep d0|d1 [--out-root DIR] [--dry]
       python3 e4_recal.py verify-prep <prep_dir> <run_dir>
       python3 e4_recal.py run <run_dir>
--dry: ローカルの乾式確認 (forge を起動しない)。FORGE_BIN を存在しない道に向け、FORGE_ALLOW_UNVERIFIED_SPECIES=1 にする (化学種の属性は
  付かない)。起動前の IC の検査の熱物性は CASE_RUNS の run_0114 の解決済み記録で代用する (E2 の乾式と同じ)。段 2 の段 1 との照合は
  run_0163 が無ければ「未確認 (DRY)」。その prep は prepare_info に DRY の印が付き、run・verify-prep が拒否する。
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import re
import shutil
import sys
from pathlib import Path

import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import throat_mono_ab as TM  # noqa: E402  (RUNS・問題の読み込み・壁の証拠を共有する)
import e4_recal_eval as EV  # noqa: E402  (段・run 名・step 数・登録の値を評価器と共有する)
import euler_t0_e2 as E2RUN  # noqa: E402  (起動前の IC の検査・実効のメッシュの記録を E2 と同じ関数で)

ROOT = C.parents[1]
PLAN = EV.PLAN
PLAN_REG_COMMIT = EV.PLAN_REG_COMMIT
BASE_PROBLEM = "problem_d155_euler_pin_G1_recal_mono.yaml"      # 単調壁・legacy + fixed2 の問題 (E3 で mesh_euler を足した)
PREPS = {"d0": "_prep_e4_d0", "d1": "_prep_e4_d1"}
MESH_JSON, MESH_CSV = "E4_MESH.json", "E4_MESH_sections.csv"
SOFT_SAVE_EXTRA = ("forge_run.log", "RUN_PROVENANCE.txt", "run_case_stdout.log", "forge_launches.jsonl")
DESIGN_FILES = ("design/forge_design/evaluate/runner_axismach.py", "design/forge_design/evaluate/runner.py",
                "design/forge_design/meshing/mesh2d.py", "design/forge_design/evaluate/ic.py", "design/forge_design/probdef.py")
MOC_KEYS = ("moc_axis_limit", "moc_corrector")


def _sha(p) -> str | None:
    return EV._sha(p)


# --- 問題の検査 ----------------------------------------------------------------------------------------------------------------
def check_problem(stage: str, delta: float | None = None, case: Path = C) -> list:
    """登録どおりか。d0 = 元の問題 (BASE_PROBLEM) と name を除いて同一で、Md_moc_offset = δ₀・MOC のキー無し・mono_r2 [0, 1.5]・
    凍結の初期線 run_0062/res_6000・mesh_euler が 2000 × 97・全域 0.005 (それ以外のキー無し)。d1 = d0 と name・Md_moc_offset を除いて同一で、
    Md_moc_offset = delta (float の完全一致)。不成立の理由のリスト。"""
    load = lambda n: yaml.safe_load((case / n).read_text())  # noqa: E731
    bad = []
    if not (case / EV.PROBLEMS[stage]).is_file():
        return [f"{EV.PROBLEMS[stage]} が無い"]
    doc = load(EV.PROBLEMS[stage])
    base = load(BASE_PROBLEM) if stage == "d0" else load(EV.PROBLEMS["d0"])
    ign = {("name",)} if stage == "d0" else {("name",), ("geometry", "Md_moc_offset")}
    bad += [f"{EV.PROBLEMS[stage]} と {BASE_PROBLEM if stage == 'd0' else EV.PROBLEMS['d0']} の違い: {d}"
            for d in E2RUN._diff_paths(base, doc, ignore=ign)]
    g = doc.get("geometry") or {}
    want_d = EV.DELTA0 if stage == "d0" else delta
    v = g.get("Md_moc_offset")
    if want_d is None or not isinstance(v, float) or v != float(want_d):
        bad.append(f"geometry.Md_moc_offset = {v!r} (登録 {want_d!r}、float であること)")
    if any(k in g for k in MOC_KEYS):
        bad.append(f"MOC のキー {[k for k in MOC_KEYS if k in g]} がある (legacy + fixed2 = キー無し)")
    if g.get("wall_fit_mono_r2") != EV.MONO_R2:
        bad.append(f"geometry.wall_fit_mono_r2 = {g.get('wall_fit_mono_r2')!r} ({EV.MONO_R2} でない)")
    if (g.get("initial_line"), g.get("initial_line_run"), g.get("initial_line_res")) != ("cfd", *EV.INITIAL_LINE):
        bad.append(f"凍結の初期線 {(g.get('initial_line'), g.get('initial_line_run'), g.get('initial_line_res'))} が登録と違う")
    if g.get("wall_repr") != "joint":
        bad.append(f"geometry.wall_repr = {g.get('wall_repr')!r} (joint でない)")
    me = doc.get("mesh_euler")
    want_me = {"ni": 2000, "nj": 97, "throat_refine": 4.0, "throat_width": 3.0, "wall_first_frac": 0.005}
    if me != want_me or any(type(me[k]) is not type(want_me[k]) for k in want_me):
        bad.append(f"mesh_euler = {me!r} (登録 {want_me})")
    if float((doc.get("spec") or {}).get("r_throat", -1)) != EV.SCALE_M or float((doc.get("spec") or {}).get("Tt", -1)) != EV.TT_REG:
        bad.append("spec.r_throat・spec.Tt が登録と違う")
    if stage == "d1" and doc.get("name") == base.get("name"):
        bad.append("d1 の name が d0 と同じ")
    return bad


def yaml_float(s: str) -> str:
    """float の repr を PyYAML (YAML 1.1) が float と読む書き方にする (仮数に小数点が無い `1e-05` は文字列と読まれるので `1.0e-05` に)。
    読み直して元の float と完全一致することを確かめる。"""
    v = float(s)
    t = repr(v)
    if ("e" in t or "E" in t) and "." not in t.split("e")[0].split("E")[0]:
        m, e = re.split(r"[eE]", t)
        t = f"{m}.0e{e}"
    got = yaml.safe_load(f"x: {t}")["x"]
    if not isinstance(got, float) or got != v:
        raise ValueError(f"{s!r} を YAML の float として書けない ({t!r} → {got!r})")
    return t


def make_d1(delta_s: str, eval_json: Path = C / EV.OUT_JSON, case: Path = C) -> Path:
    """段 1 の評価が「更新」で、その δ₁ の repr が delta_s と同じときだけ、d0 から d1 の問題を作る (name と Md_moc_offset の行だけを書き換える)。"""
    ev = json.loads(Path(eval_json).read_text())
    j1 = ((ev.get("stages") or {}).get("d0") or {}).get("judgment") or {}
    if j1.get("verdict") != EV.LBL_UPDATE:
        raise SystemExit(f"段 1 の判定が「更新」でない ({j1.get('verdict')!r}) — 段 2 は登録外。作らない")
    if j1.get("delta1_repr") != delta_s or float(delta_s) != float(j1.get("delta1")):
        raise SystemExit(f"--delta {delta_s} が段 1 の δ₁ ({j1.get('delta1_repr')}) と一致しない — 作らない")
    if ev.get("evaluator_sha256") != _sha(C / "e4_recal_eval.py"):
        raise SystemExit("段 1 の評価器の sha256 が今の e4_recal_eval.py と違う — 作らない")
    src = (case / EV.PROBLEMS["d0"]).read_text()
    lines = src.splitlines(keepends=True)
    i_name = [i for i, l in enumerate(lines) if re.match(r"^name:\s", l)]
    i_off = [i for i, l in enumerate(lines) if re.match(r"^  Md_moc_offset:\s", l)]
    if len(i_name) != 1 or len(i_off) != 1:
        raise SystemExit("d0 の name・Md_moc_offset の行が 1 本ずつでない")
    lines[i_name[0]] = "name: isobutane_m6_d155_euler_e4_recal_d1\n"
    lines[i_off[0]] = (f"  Md_moc_offset: {yaml_float(delta_s)}   # δ₁ = δ₀ − (平均 M_common − 6) (段 1 run_0163 の評価; "
                       f"評価器 {ev.get('evaluator_sha256', '')[:16]})\n")
    hdr = (f"# plan verification-case45-euler-total-enthalpy §6 E4 段 2 (δ₁): e4_recal.py make-d1 が {EV.PROBLEMS['d0']} から name と\n"
           f"#   geometry.Md_moc_offset だけを書き換えて作った (段 1 の判定 {j1.get('verdict')[:2]}、δ₁ = {delta_s})。\n")
    out = case / EV.PROBLEMS["d1"]
    text = hdr + "".join(lines)
    if out.exists() and out.read_text() != text:
        raise SystemExit(f"{out.name} が既にあり内容が違う — 消してから作る")
    out.write_text(text)
    bad = check_problem("d1", float(delta_s), case=case)
    if bad:
        out.unlink()
        raise SystemExit("作った d1 の検査が不成立 (消した):\n  " + "\n  ".join(bad))
    return out


# --- 準備 --------------------------------------------------------------------------------------------------------------------------
def _input_hashes(d: Path) -> dict:
    return {p.name: _sha(p) for p in sorted(d.iterdir()) if p.is_file() and p.name != EV.PREP_JSON}


def _binaries() -> dict:
    return {k: {"path": os.environ.get(k), "sha256": (_sha(os.environ[k]) if os.environ.get(k) else None)}
            for k in ("FORGE_BIN", "REAL_CONVERTER", "FORGE_CONVERTER")}


def _git(*args) -> str:
    return E2RUN._git(*args)


def prep_checks(d: Path, info: dict, stage: str, expected_delta: float) -> tuple:
    """prepare の後の検査 (不成立の理由のリスト, 壁の証拠, Mesh2DParams)。"""
    bad = [f"step 数の同期: {x}" for x in EV.step_problems(EV.config_steps((d / "solverConfig.yaml").read_text()))]
    bad += EV.prepare_info_problems(info, expected_delta)
    gate = (info.get("moc") or {}).get("gate") or {}
    if gate.get("applicable") not in (False, None):
        bad.append(f"MOC のゲートが legacy + fixed2 の「適用外」でない ({gate})")
    ev = TM.wall_evidence(d)
    if ev.get("status") != "consistent":
        bad.append(f"壁の証拠 (M4) が {ev.get('status')}")
    q = (d / "MESH_QUALITY.txt").read_text() if (d / "MESH_QUALITY.txt").is_file() else ""
    if not any(l.strip().startswith("VERDICT: PASS") for l in q.splitlines()):
        bad.append("メッシュ品質が PASS でない")
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.evaluate import runner_axismach as RA
    mp = RA.mesh_params_euler(RA.load_problem(C / EV.PROBLEMS[stage]), float(info["scale_m"]))   # prepare と同じ呼び方 (mesh_euler)
    if dataclasses.asdict(mp) != (info.get("mesh") or {}).get("params"):
        bad.append("prepare_info の mesh.params が問題の mesh_euler の解決値と違う")
    return bad, ev, mp


def compare_with_d0(d: Path, dry: bool) -> dict:
    """段 2 の prep の bcondConfig・solverConfig・probe・species_meta が段 1 の prep (run_0163 の E4_PREP.json の記録) と同一か。"""
    p0 = TM.RUNS / EV.RUNS["d0"] / EV.PREP_JSON
    rec = {"ref": str(p0), "failures": []}
    if not p0.is_file():
        rec["status"] = "未確認 (DRY: run_0163 が無い)" if dry else "run_0163 の E4_PREP.json が無い"
        if not dry:
            rec["failures"].append(f"{p0} が無い (段 1 の prep と照合できない)")
        return rec
    f0 = json.loads(p0.read_text()).get("files_sha256") or {}
    rec["files"] = {n: {"d0": f0.get(n), "d1": _sha(d / n)} for n in EV.SAME_AS_D0}
    diff = [n for n, v in rec["files"].items() if (v["d0"] or v["d1"]) and v["d0"] != v["d1"]]
    if diff:
        rec["failures"].append(f"段 1 の prep と違う: {diff}")
    rec["status"] = "OK" if not rec["failures"] else "FAIL"
    return rec


def prep(stage: str, out_root: Path, dry: bool = False) -> dict:
    if stage not in EV.STAGES:
        raise SystemExit(f"段は {EV.STAGES} のどれか ({stage!r})")
    if dry:
        # 乾式確認: forge を起動しない (--resolve-species も含めて)。runner の _ENV は import 時に環境を写すので import より前に設定する
        os.environ["FORGE_BIN"] = str(Path("/nonexistent/forge-dry-run-guard"))
        os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
    expected = EV.DELTA0
    if stage == "d1":
        doc = yaml.safe_load((C / EV.PROBLEMS["d1"]).read_text()) if (C / EV.PROBLEMS["d1"]).is_file() else {}
        expected = (doc.get("geometry") or {}).get("Md_moc_offset")
    bad = check_problem(stage, expected)
    if bad:
        raise SystemExit("問題の検査が不成立 — 止める:\n  " + "\n  ".join(bad))
    out_root = Path(out_root).resolve()
    d = out_root / PREPS[stage]
    if d.exists():
        raise SystemExit(f"{d} が既にある (投入スクリプトが消してから作る)")
    RA = TM._load_problem_with_runs()
    prob = C / EV.PROBLEMS[stage]
    info = RA.prepare(prob, d, nsteps=EV.MAIN_NSTEPS, ic_from=None, cfl_main=EV.CFL_MAIN, implicit_relax=EV.RELAX)
    fails, ev, mp = prep_checks(d, info, stage, expected)
    print(stage, prob.name, "| mesh", json.dumps({k: info["mesh"].get(k) for k in ("ni", "nj", "wall_first_frac", "source")}),
          "| M4", ev.get("status"), "|", (d / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
    same = compare_with_d0(d, dry) if stage == "d1" else None
    if same:
        fails += same["failures"]
    thermo = TM.RUNS / E2RUN.DRY_THERMO_RUN if dry else d
    try:
        ic = E2RUN.ic_check(d, thermo, dry)
    except Exception as e:  # noqa: BLE001 — 検査できないことも不成立として記録する
        ic = {"VERDICT": "FAIL", "problems": [f"IC の検査を完了できない: {type(e).__name__}: {e}"]}
    fails += [f"起動前の IC の検査: {p}" for p in ic["problems"]]
    rec, rows = E2RUN.mesh_record(d, prob, mp, float(info["scale_m"]))
    rec.update(plan=PLAN, plan_reg_commit=PLAN_REG_COMMIT, stage=stage,
               note="prepare_info.json の mesh 欄 (全 Mesh2DParams・採用元・座標と接続のハッシュ) に加えた、断面ごとの実測の記録 (E2 と同じ形式)",
               prepare_info_mesh_hashes=(info.get("mesh") or {}).get("hashes"))
    E2RUN.write_sections(d / MESH_CSV, rows)
    rec["hashes"]["sections_csv"] = _sha(d / MESH_CSV)
    rec.update(binaries=_binaries(), base_problem_sha256=_sha(C / BASE_PROBLEM), design_sha256={f: _sha(ROOT / f) for f in DESIGN_FILES},
               git={"head": _git("rev-parse", "HEAD"), "status_design": _git("status", "--porcelain", "--", "design")})
    (d / MESH_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    summary = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "tool": "e4_recal.prep", "stage": stage, "dry": dry,
               "problem": prob.name, "problem_sha256": _sha(prob), "Md_moc_offset": expected,
               "steps": {"evaluator": {"MAIN_NSTEPS": EV.MAIN_NSTEPS, "OUT_INTERVAL": EV.OUT_INTERVAL,
                                       "windows": {k: list(v) for k, v in EV.WINDOWS.items()}},
                         "config": EV.config_steps((d / "solverConfig.yaml").read_text())},
               "wall_evidence": ev, "same_as_d0": same, "ic_check": ic, "mesh": rec, "failures": fails,
               "VERDICT": "OK" if not fails else "REFUSED"}
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / f"_band_ab/e4_recal_prep_{stage}{'_dry' if dry else ''}.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    if fails:
        raise RuntimeError("E4 の準備の検査が不成立 — 止める (forge を起動しない):\n  " + "\n  ".join(fails))
    info = json.loads((d / "prepare_info.json").read_text())
    info.update(stages="soft", plan=PLAN,
                e4={"stage": stage, "run": EV.RUNS[stage], "problem": prob.name, "Md_moc_offset": expected, "plan_reg_commit": PLAN_REG_COMMIT},
                ic={"mode": "isentropic", "tool": "paste_isentropic_ic (runner_axismach.prepare)", "VERDICT": ("OK" if not dry else "DRY"),
                    "check": {k: ic.get(k) for k in ("max_abs_dev_K", "T0_min", "T0_max", "thermo_source", "VERDICT")}})
    if dry:
        info["DRY"] = True
    (d / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    prec = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "stage": stage, "run": EV.RUNS[stage], "problem": prob.name,
            "Md_moc_offset": expected, "dry": dry, "ic_check": ic, "same_as_d0": same,
            "steps": {"config": summary["steps"]["config"], "MAIN_NSTEPS": EV.MAIN_NSTEPS, "OUT_INTERVAL": EV.OUT_INTERVAL},
            "files_sha256": _input_hashes(d)}
    (d / EV.PREP_JSON).write_text(json.dumps(prec, indent=1, ensure_ascii=False, default=float))
    tb = rec["summary_bands_rt"]
    print(f"{stage} 実効のメッシュ: 壁側 (/r_w) スロート {tb['-4<=x<=1'].get('wall_gap_rw')} 下流 x≥17 {tb['x>=17'].get('wall_gap_rw')} | "
          f"軸側 (/r_w) スロート {tb['-4<=x<=1'].get('axis_gap_rw')} 下流 {tb['x>=17'].get('axis_gap_rw')} | "
          f"IC max|T0−1600| {ic['max_abs_dev_K']:.4g} K ({ic['VERDICT']}) | Md_moc_offset {expected!r}")
    return summary


def verify_prep(prep_dir: Path, run_dir: Path) -> list:
    """run_dir が prep_dir の複製のまま (E4_PREP.json の入力の sha256 が一致) か。戻り = 不成立の理由。"""
    bad = []
    try:
        rec = json.loads((prep_dir / EV.PREP_JSON).read_text())
        info = json.loads((run_dir / "prepare_info.json").read_text())
    except (OSError, ValueError) as e:
        return [f"{EV.PREP_JSON} / prepare_info.json を読めない ({e})"]
    if rec.get("dry") or info.get("DRY"):
        bad.append("乾式確認 (DRY) の prep")
    if (rec.get("ic_check") or {}).get("VERDICT") != "OK":
        bad.append(f"起動前の IC の検査が OK でない ({(rec.get('ic_check') or {}).get('VERDICT')!r})")
    if run_dir.name != rec.get("run"):
        bad.append(f"run の名前 {run_dir.name} が prep の記録 {rec.get('run')} と違う")
    for n, sha in (rec.get("files_sha256") or {}).items():
        if _sha(run_dir / n) != sha:
            bad.append(f"{run_dir.name}/{n} が prep の記録と違う")
    if not rec.get("files_sha256"):
        bad.append("prep の入力の sha256 の記録が無い")
    return bad


# --- 実行 (soft 段の出力の退避つき) ---------------------------------------------------------------------------------------------
def install_soft_saver(RA, rd: Path) -> list:
    """RA._restart_same_mesh を包む: 呼ばれる前に run 直下の res_* と段のログを rd/_soft_stage/ に写す (1 回だけ。soft 段の引き継ぎ)。
    euler_t0_e2.install_soft_saver と同じ方式 (記録の plan・tool の名前だけが違う)。戻り値は呼ばれた記録のリスト。"""
    orig = RA._restart_same_mesh
    calls = []

    def wrapped(res_h5, mesh_h5):
        res_h5 = Path(res_h5)
        if calls:
            raise RuntimeError("段の引き継ぎが 2 回目 — 登録の段は soft → 本段だけ (mid 段は無い)")
        if res_h5.parent.resolve() != rd.resolve() or res_h5.name != f"res_{EV.SOFT_STEPS}.h5":
            raise RuntimeError(f"soft 段の最終場が {res_h5} (期待 {rd}/res_{EV.SOFT_STEPS}.h5)")
        dst = rd / EV.SOFT_DIR
        dst.mkdir(exist_ok=False)
        rec = {"plan": PLAN, "tool": "e4_recal.install_soft_saver", "stage": "soft", "restart_src": res_h5.name, "files": {}}
        srcs = sorted(p for p in rd.glob("res_*") if p.is_file()) + [rd / n for n in SOFT_SAVE_EXTRA if (rd / n).is_file()]
        for p in srcs:
            shutil.copy2(p, dst / p.name)
            a, b = _sha(p), _sha(dst / p.name)
            if a != b:
                raise RuntimeError(f"soft 段の出力の写し {p.name} が元と一致しない")
            rec["files"][p.name] = a
        rec["nozzle_h5_sha256_before_restart"] = _sha(mesh_h5)
        calls.append(str(res_h5))
        try:
            orig(res_h5, mesh_h5)
        finally:
            rec["nozzle_h5_sha256_after_restart"] = _sha(mesh_h5)
            (dst / "SOFT_STAGE_SAVED.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    RA._restart_same_mesh = wrapped
    return calls


def check_run_steps(rd: Path) -> list:
    """起動の直前: run の solverConfig.yaml の nStepOuter・outStepInterval が評価器の定数と同じで、窓がその定数から作られた形か。"""
    try:
        return EV.step_problems(EV.config_steps((rd / "solverConfig.yaml").read_text()))
    except Exception as e:  # noqa: BLE001 — 読めないことも食い違いとして止める
        return [f"solverConfig.yaml を読めない: {type(e).__name__}: {e}"]


def run(rd: Path) -> int:
    from forge_design.evaluate import runner_axismach as RA
    info = json.loads((rd / "prepare_info.json").read_text())
    if info.get("DRY"):
        raise SystemExit(f"{rd} は乾式確認 (--dry) の prep から作られている — 回さない")
    e4 = info.get("e4") or {}
    stage = e4.get("stage")
    if stage not in EV.STAGES or rd.name != EV.RUNS[stage]:
        raise SystemExit(f"{rd}: E4 の段 {stage!r} とその run 名でない — 回さない")
    if (info.get("ic") or {}).get("VERDICT") != "OK":
        raise SystemExit(f"{rd}: IC の VERDICT が OK でない ({(info.get('ic') or {}).get('VERDICT')!r}) — 回さない")
    bad = EV.prepare_info_problems(info, e4.get("Md_moc_offset"))
    if bad or (stage == "d0" and float(e4.get("Md_moc_offset")) != EV.DELTA0):
        raise SystemExit(f"{rd}: prepare_info の固定の条件が不成立 — 回さない: {bad}")
    if (rd / EV.SOFT_DIR).exists() or any(rd.glob("res_*")):
        raise SystemExit(f"{rd}: 既に出力がある (_soft_stage または res_*) — 回さない")
    sp = check_run_steps(rd)
    if sp:
        raise SystemExit(f"{rd}: 本段の step 数が評価器と同期していない — 回さない:\n  " + "\n  ".join(sp))
    calls = install_soft_saver(RA, rd)
    rc = RA.run_staged(rd, cfl_main=EV.CFL_MAIN, mid_stage=False, stages="soft")
    if len(calls) != 1:
        raise RuntimeError(f"soft 段の出力の退避が {len(calls)} 回 (1 回のはず)")
    print(f"forge exit={rc}")
    return rc


def main(argv) -> int:
    if argv and argv[0] == "check-problem":
        ap = argparse.ArgumentParser(prog="e4_recal.py check-problem")
        ap.add_argument("stage", choices=EV.STAGES)
        ap.add_argument("--delta", default=None, help="d1 の Md_moc_offset (段 1 の δ₁ の repr)")
        a = ap.parse_args(argv[1:])
        if a.stage == "d1" and a.delta is None:
            print("CHECK-PROBLEM: FAIL\n  d1 は --delta が必要")
            return 1
        bad = check_problem(a.stage, None if a.delta is None else float(a.delta))
        print("CHECK-PROBLEM: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if argv and argv[0] == "make-d1":
        ap = argparse.ArgumentParser(prog="e4_recal.py make-d1")
        ap.add_argument("--delta", required=True)
        ap.add_argument("--eval-json", default=str(C / EV.OUT_JSON))
        a = ap.parse_args(argv[1:])
        print("MAKE-D1:", make_d1(a.delta, Path(a.eval_json)).name)
        return 0
    if argv and argv[0] == "prep":
        ap = argparse.ArgumentParser(prog="e4_recal.py prep")
        ap.add_argument("stage", choices=EV.STAGES)
        ap.add_argument("--out-root", default=str(C), help="_prep_e4_* を作る場所 (既定は case dir)")
        ap.add_argument("--dry", action="store_true", help="乾式確認 (forge を起動しない)")
        a = ap.parse_args(argv[1:])
        prep(a.stage, Path(a.out_root), dry=a.dry)
        return 0
    if len(argv) == 3 and argv[0] == "verify-prep":
        bad = verify_prep(Path(argv[1]).resolve(), Path(argv[2]).resolve())
        print("VERIFY-PREP: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if len(argv) == 2 and argv[0] == "run":
        return run(Path(argv[1]).resolve())
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
