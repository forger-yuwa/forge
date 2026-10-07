"""plan discretization-moc-axis-limit-and-corrector §6 V5d (2026-10-07 登録 99431498) の問題の生成・run の準備・実行。
新しい Euler の格子 (E3 の mesh_euler、2000 × 97・全域 0.005) で、腕 B = 単調壁・legacy + fixed2 と 腕 M = 単調壁・analytic + converge を
比べる。両腕とも同じ出口較正の値 (plan verification-case45-euler-total-enthalpy の E4/E4V の結果) を使い、等エントロピー IC から各 3 本。
腕・run 名・窓・判定は評価器 moc_v5d_eval.py と共有する (起動は run_moc_v5d.sh)。design/ は変えない (import だけ)。

問題 (make-problems): 腕 B = problem_d155_euler_e4_recal_d0.yaml (E4 の問題: 単調壁・legacy + fixed2・mesh_euler・凍結の初期線 run_0062) の
  写しで name と geometry.Md_moc_offset (= --md-offset) だけを変えたもの → problem_d155_euler_v5d_B.yaml。
  腕 M = 腕 B に geometry.moc_axis_limit: analytic・geometry.moc_corrector: converge の 2 行を足し name を変えたもの
  → problem_d155_euler_v5d_M.yaml。検査 (check-problems): B と E4 の d0 の差が name・Md_moc_offset だけ、M と B の差が name・MOC の 2 キー
  だけ、M と V5 の腕 M の問題 (problem_d155_euler_pin_G1_recal_mono_moc.yaml) の差が name・Md_moc_offset だけ、Md_moc_offset が float で
  --md-offset と完全一致。共用 (--share-e4v) なら B と E4V の問題 (problem_d155_euler_e4_recal_d1.yaml) の差が name だけ。
  本番では E4 の評価 (_band_ab/e4_recal_eval.json) が較正値を採用し、その値が --md-offset と同じことを要求する (乾式は記録だけ)。
  既存の問題と中身が違う場合、V5d の run dir が 1 つも無ければ書き直し、あれば止める。
準備 (prep B|M): RA.prepare (等エントロピー IC; ic_from なし) → 検査 (step 数の同期・MOC・単調壁・格子の採用元と実効値・Md_moc_offset・
  凍結の初期線・r_t・壁の証拠 M4・メッシュ品質 PASS) → 起動前の IC の検査 (全節点 |T0 − 1600| ≤ 1 K; E2・E4 と同じ euler_t0_e2.ic_check)
  → 実効のメッシュの記録 (V5D_MESH.json・V5D_MESH_sections.csv) → 共用なら腕 B の prep と run_0164 の照合 (設定ファイル・保存 spline・
  格子のハッシュ・Md_moc_offset・MOC・x_E) → V5D_PREP.json。不成立はすべて例外で止める (投入スクリプトは forge を起動しない)。
照合 (cross-check): 腕 B と腕 M の prep の設定ファイル (bcondConfig・solverConfig・probe・species_meta) と格子の実効値が同じで、
  壁は腕 M の壁節点が腕 B と見分けられ、そのすべてが自腕の spline に一致すること。
実行 (run): runner_axismach.run_staged(stages="soft", mid_stage=False, cfl_main=2) — soft (1 次・cfl 0.5・3000 step) → 本段 (2 次・cfl 2・
  implicitRelax 0.7・54000 step・1000 ごと出力)。soft 段の出力は runner が消す前に <run>/_soft_stage/ に退避する (E4 と同じ方式)。
usage: [CASE_RUNS=<run_0062・run_0114 (・run_0164) のある case dir>] python3 moc_v5d.py make-problems --md-offset X [--share-e4v] [--dry]
       python3 moc_v5d.py check-problems --md-offset X [--share-e4v] [--dry]
       python3 moc_v5d.py check-share --md-offset X [--dry]          (共用の前提: run_0164 の完走・較正値・E4 の参照 run)
       python3 moc_v5d.py prep B|M --md-offset X [--share-e4v] [--out-root DIR] [--dry]
       python3 moc_v5d.py cross-check <prep_B> <prep_M> --md-offset X [--dry]
       python3 moc_v5d.py verify-prep <prep_dir> <run_dir> [<run_dir> ...]
       python3 moc_v5d.py run <run_dir>
       python3 moc_v5d.py plan-runs [--share-e4v]                  (投入スクリプト用: 新しい run・prep・起動の順)
       python3 moc_v5d.py launch-record --md-offset X [--share-e4v] [--dry] [--extra JSON]
--dry: ローカルの乾式確認 (forge を起動しない)。FORGE_BIN を存在しない道に向け、FORGE_ALLOW_UNVERIFIED_SPECIES=1 にする。
  起動前の IC の検査の熱物性は CASE_RUNS の run_0114 の解決済み記録で代用する (E2・E4 の乾式と同じ)。その prep は prepare_info に DRY の印が
  付き、run・verify-prep が拒否する。E4 の採用・run_0164 の照合は「未確認 (DRY)」として記録だけする。
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import throat_mono_ab as TM  # noqa: E402  (RUNS・問題の読み込み・壁の証拠)
import e4_recal as E4R  # noqa: E402  (YAML の float の書き方)
import e4_recal_eval as EV4  # noqa: E402
import euler_t0_e2 as E2RUN  # noqa: E402  (起動前の IC の検査・実効のメッシュの記録・YAML の木の差)
import euler_t0_e2_eval as EV2  # noqa: E402
import moc_v5d_eval as EV  # noqa: E402  (腕・run 名・窓・固定の条件の検査を評価器と共有する)
import moc_v5_euler_eval as V5  # noqa: E402  (壁の他腕との照合)

ROOT = C.parents[1]
PLAN = EV.PLAN
BASE_B = EV4.PROBLEMS["d0"]                                   # E4 の段 1 の問題 (単調壁・legacy + fixed2・mesh_euler)
E4V_PROBLEM = EV4.PROBLEMS["d1"]                              # E4V の問題 (共用の照合)
V5_M_PROBLEM = "problem_d155_euler_pin_G1_recal_mono_moc.yaml"  # V5 の腕 M の問題 (MOC の 2 キーの書き方の照合)
MOC_LINES = ("  moc_axis_limit: analytic       # 軸上の端点の sinθ/r を解析極限 θ_r に (plan discretization-moc-axis-limit-and-corrector §4.1)\n",
             "  moc_corrector: converge        # 予測修正を収束まで (更新量 ≤ 1e-12、上限 50 回; 同 plan §4.2)\n")
PROBLEMS_RECORD = "_band_ab/moc_v5d_problems.json"
MESH_JSON, MESH_CSV = "V5D_MESH.json", "V5D_MESH_sections.csv"
SOFT_SAVE_EXTRA = ("forge_run.log", "RUN_PROVENANCE.txt", "run_case_stdout.log", "forge_launches.jsonl")
DESIGN_FILES = ("design/forge_design/evaluate/runner_axismach.py", "design/forge_design/evaluate/runner.py",
                "design/forge_design/meshing/mesh2d.py", "design/forge_design/evaluate/ic.py", "design/forge_design/probdef.py",
                "design/forge_design/geometry/moc_kernel.py", "design/forge_design/geometry/moc_inverse.py")
IC_PLACEHOLDER = {"mode": "isentropic", "VERDICT": "OK"}       # prepare の直後の検査では IC は別に検査する (下の ic_check)


def _sha(p) -> str | None:
    return EV._sha(p)


def _dry_env() -> None:
    """乾式確認: forge を起動しない (--resolve-species も含めて)。runner の _ENV は import 時に環境を写すので import より前に設定する。"""
    os.environ["FORGE_BIN"] = str(Path("/nonexistent/forge-dry-run-guard"))
    os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"


def _load_yaml(p: Path) -> dict:
    sys.path.insert(0, str(EV.TOOLS))
    import yaml_strict
    return yaml_strict.load(Path(p).read_text())


def any_v5d_run_exists(case: Path = C) -> list:
    return [r for arm in ("B", "M") for r in EV.RUNS_NEW[arm] if (case / r).exists()]


# --- 問題 ----------------------------------------------------------------------------------------------------------------------
def problem_texts(md_repr: str, case: Path = C) -> dict:
    """{腕: YAML の本文}。E4 の d0 の name と Md_moc_offset の行だけを書き換え、腕 M は Md_moc_offset の直後に MOC の 2 行を足す。"""
    src = (case / BASE_B).read_text()
    lines = src.splitlines(keepends=True)
    i_name = [i for i, l in enumerate(lines) if re.match(r"^name:\s", l)]
    i_off = [i for i, l in enumerate(lines) if re.match(r"^  Md_moc_offset:\s", l)]
    if len(i_name) != 1 or len(i_off) != 1:
        raise SystemExit(f"{BASE_B} の name・Md_moc_offset の行が 1 本ずつでない")
    if any(re.match(r"^  moc_(axis_limit|corrector):", l) for l in lines):
        raise SystemExit(f"{BASE_B} に MOC のキーがある (legacy + fixed2 = キー無しのはず)")
    tok = E4R.yaml_float(md_repr)
    out = {}
    for arm in ("B", "M"):
        ls = list(lines)
        ls[i_name[0]] = f"name: {EV.PROBLEM_NAMES[arm]}\n"
        ls[i_off[0]] = (f"  Md_moc_offset: {tok}   # V5d の両腕に共通の出口較正 (E4 の採用値; moc_v5d.py make-problems --md-offset {md_repr})\n")
        if arm == "M":
            ls[i_off[0] + 1:i_off[0] + 1] = list(MOC_LINES)
        role = ("腕 B (単調壁・legacy + fixed2 [MOC のキー無し])" if arm == "B"
                else "腕 M (単調壁・analytic + converge: 腕 B に geometry.moc_axis_limit・moc_corrector の 2 行だけを足した)")
        hdr = (f"# plan discretization-moc-axis-limit-and-corrector §6 V5d (登録 {EV.PLAN_REG_COMMIT}): {role}。\n"
               f"#   生成: moc_v5d.py make-problems --md-offset {md_repr} (手で編集しない)。元: {BASE_B} (sha256 {_sha(case / BASE_B)[:16]}…) の\n"
               f"#   name と geometry.Md_moc_offset だけを書き換えた。格子は mesh_euler (2000 × 97・全域 0.005; E3)。両腕で Md_moc_offset は同じ値。\n"
               f"#   以下は元の問題のコメント (E4 の段 1 の説明を含む) をそのまま残している。\n")
        out[arm] = hdr + "".join(ls)
    return out


def check_problems(md: float, share: bool, dry: bool = False, case: Path = C) -> tuple:
    """登録どおりか。戻り (不成立の理由, 記録)。"""
    bad, rec = [], {"files": {}}
    try:
        doc = {arm: _load_yaml(case / EV.PROBLEMS[arm]) for arm in ("B", "M")}
        d0 = _load_yaml(case / BASE_B)
        v5m = _load_yaml(case / V5_M_PROBLEM)
    except Exception as e:  # noqa: BLE001 — 読めない・重複キーは不成立
        return [f"問題を読めない: {type(e).__name__}: {e}"], rec
    rec["files"] = {n: _sha(case / n) for n in (EV.PROBLEMS["B"], EV.PROBLEMS["M"], BASE_B, V5_M_PROBLEM)}
    D = E2RUN._diff_paths
    bad += [f"腕 B と {BASE_B} の違い (name・Md_moc_offset 以外): {d}" for d in D(d0, doc["B"], ignore={("name",), ("geometry", "Md_moc_offset")})]
    bad += [f"腕 M と腕 B の違い (name・MOC の 2 キー以外): {d}"
            for d in D(doc["B"], doc["M"], ignore={("name",), ("geometry", "moc_axis_limit"), ("geometry", "moc_corrector")})]
    bad += [f"腕 M と V5 の腕 M の問題 {V5_M_PROBLEM} の違い (name・Md_moc_offset 以外): {d}"
            for d in D(v5m, doc["M"], ignore={("name",), ("geometry", "Md_moc_offset")})]
    gB, gM = doc["B"].get("geometry") or {}, doc["M"].get("geometry") or {}
    if any(k in gB for k in ("moc_axis_limit", "moc_corrector")):
        bad.append("腕 B に MOC のキーがある (legacy + fixed2 = キー無し)")
    if (gM.get("moc_axis_limit"), gM.get("moc_corrector")) != ("analytic", "converge"):
        bad.append(f"腕 M の MOC のキーが analytic・converge でない ({gM.get('moc_axis_limit')!r}, {gM.get('moc_corrector')!r})")
    for arm, g in (("B", gB), ("M", gM)):
        v = g.get("Md_moc_offset")
        if not isinstance(v, float) or v != md:
            bad.append(f"腕 {arm} の geometry.Md_moc_offset = {v!r} が --md-offset {md!r} と完全一致する float でない")
        if doc[arm].get("name") != EV.PROBLEM_NAMES[arm]:
            bad.append(f"腕 {arm} の name {doc[arm].get('name')!r} が {EV.PROBLEM_NAMES[arm]} でない")
    if share:
        p1 = case / E4V_PROBLEM
        if p1.is_file():
            bad += [f"共用: 腕 B と E4V の問題 {E4V_PROBLEM} の違い (name 以外): {d}"
                    for d in D(_load_yaml(p1), doc["B"], ignore={("name",)})]
            rec["e4v_problem"] = {"file": E4V_PROBLEM, "sha256": _sha(p1)}
        elif dry:
            rec["e4v_problem"] = {"status": f"未確認 (DRY: {E4V_PROBLEM} が無い)"}
        else:
            bad.append(f"共用: E4V の問題 {E4V_PROBLEM} が無い (照合できない)")
    return bad, rec


def make_problems(md_text: str, share: bool, dry: bool = False, case: Path = C) -> dict:
    md, md_repr = EV.parse_md(md_text)
    e4 = EV.e4_adoption(case, md, share)
    if e4["problems"] and not dry:
        raise SystemExit("E4 の結果と合わない — 問題を作らない:\n  " + "\n  ".join(e4["problems"]))
    texts = problem_texts(md_repr, case)
    existing = any_v5d_run_exists(case)
    written = []
    for arm, text in texts.items():
        p = case / EV.PROBLEMS[arm]
        if p.exists() and p.read_text() == text:
            continue
        if p.exists() and existing:
            raise SystemExit(f"{p.name} が既にあり中身が違う。V5d の run ({existing}) があるので書き直さない — 止める")
        p.write_text(text)
        written.append(p.name)
    bad, prec = check_problems(md, share, dry, case)
    rec = {"plan": PLAN, "plan_reg_commit": EV.PLAN_REG_COMMIT, "tool": "moc_v5d.make_problems", "date": datetime.now(timezone.utc).isoformat(),
           "md_offset": md, "md_offset_repr": md_repr, "md_offset_yaml_token": E4R.yaml_float(md_repr), "share_e4v": share, "dry": dry,
           "e4_adoption": ({**e4, "status": "記録だけ (DRY)"} if dry and e4["problems"] else e4), "written": written,
           "check": {"problems": bad, **prec}}
    (case / "_band_ab").mkdir(exist_ok=True)
    (case / (PROBLEMS_RECORD.replace(".json", "_dry.json") if dry else PROBLEMS_RECORD)).write_text(
        json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    if bad:
        raise SystemExit("作った問題の検査が不成立:\n  " + "\n  ".join(bad))
    return rec


def check_share(md: float, dry: bool = False) -> list:
    """共用の前提: run_0164 が完走 (RUN_RC 0・早期停止なし)・本段の窓の res と残差判定があり、Md_moc_offset が --md-offset と同じ。"""
    rd = TM.RUNS / EV.E4V_RUN
    if not rd.is_dir():
        return [] if dry else [f"{rd} が無い"]
    bad = []
    rc = (rd / "RUN_RC").read_text().strip() if (rd / "RUN_RC").is_file() else None
    if rc != "0":
        bad.append(f"{EV.E4V_RUN} の RUN_RC が 0 でない ({rc!r})")
    if (rd / "EARLY_STOP.txt").exists():
        bad.append(f"{EV.E4V_RUN} に EARLY_STOP.txt がある")
    miss = [s for s in EV.WIN13 if not (rd / f"res_{s}.h5").is_file()]
    if miss:
        bad.append(f"{EV.E4V_RUN} に判定窓の res が無い: {miss}")
    if not (rd / EV4.SEGMENT_VERDICT_FILE).is_file():
        bad.append(f"{EV.E4V_RUN} に {EV4.SEGMENT_VERDICT_FILE} が無い")
    try:
        v = json.loads((rd / "prepare_info.json").read_text()).get("Md_moc_offset")
    except (OSError, ValueError):
        v = None
    if not (isinstance(v, float) and v == md):
        bad.append(f"{EV.E4V_RUN} の Md_moc_offset {v!r} が --md-offset {md!r} と違う")
    return bad


# --- 準備 --------------------------------------------------------------------------------------------------------------------------
def _input_hashes(d: Path) -> dict:
    return {p.name: _sha(p) for p in sorted(d.iterdir()) if p.is_file() and p.name != EV.PREP_JSON}


def _binaries() -> dict:
    return {k: {"path": os.environ.get(k), "sha256": (_sha(os.environ[k]) if os.environ.get(k) and Path(os.environ[k]).is_file() else None)}
            for k in ("FORGE_BIN", "REAL_CONVERTER", "FORGE_CONVERTER")}


def _initial_line_key(info: dict):
    """凍結の初期線の照合用の値: run の名前・snapshot・内容のハッシュ (sha256_16)。run の絶対パスは起動の仕方 (CASE_RUNS) で変わるので比べない。"""
    il = info.get("initial_line") or {}
    if not isinstance(il, dict) or not il.get("run"):
        return None
    return (Path(str(il["run"])).name, il.get("res"), il.get("sha256_16"))


def compare_with_e4v(d: Path, info: dict, dry: bool) -> dict:
    """共用: 腕 B の prep が run_0164 (E4V) と同じ条件か (設定ファイル・保存 spline・格子の実効値とハッシュ・Md_moc_offset・MOC・x_E)。"""
    rd = TM.RUNS / EV.E4V_RUN
    rec = {"ref": str(rd), "failures": []}
    if not rd.is_dir():
        rec["status"] = "未確認 (DRY: run_0164 が無い)" if dry else "run_0164 が無い"
        if not dry:
            rec["failures"].append(f"{rd} が無い (照合できない)")
        return rec
    try:
        pi = json.loads((rd / "prepare_info.json").read_text())
        pr = json.loads((rd / EV4.PREP_JSON).read_text())
    except (OSError, ValueError) as e:
        rec["failures"].append(f"run_0164 の prepare_info.json・{EV4.PREP_JSON} を読めない ({type(e).__name__})")
        rec["status"] = "FAIL"
        return rec
    f0 = pr.get("files_sha256") or {}
    rec["files"] = {n: {"e4v": f0.get(n), "v5d_B": _sha(d / n)} for n in EV.SAME_FILES}
    diff = [n for n, v in rec["files"].items() if (v["e4v"] or v["v5d_B"]) and v["e4v"] != v["v5d_B"]]
    if diff:
        rec["failures"].append(f"run_0164 の prep と違う設定ファイル: {diff}")
    for key, get in (("wall_fit.spline", lambda i: (i.get("wall_fit") or {}).get("spline")),
                     ("mesh.params", lambda i: (i.get("mesh") or {}).get("params")),
                     ("mesh.hashes", lambda i: (i.get("mesh") or {}).get("hashes")),
                     ("Md_moc_offset", lambda i: i.get("Md_moc_offset")), ("x_E", lambda i: i.get("x_E")),
                     ("moc", lambda i: {k: (i.get("moc") or {}).get(k) for k in ("axis_limit", "corrector")}),
                     ("initial_line", _initial_line_key), ("scale_m", lambda i: i.get("scale_m"))):
        a, b = get(pi), get(info)
        if a is None or a != b:
            rec["failures"].append(f"run_0164 と {key} が違う・欠ける")
    rec["status"] = "OK" if not rec["failures"] else "FAIL"
    return rec


def prep(arm: str, md_text: str, share: bool, out_root: Path, dry: bool = False) -> dict:
    if arm not in ("B", "M"):
        raise SystemExit(f"腕は B か M ({arm!r})")
    if dry:
        _dry_env()
    md, md_repr = EV.parse_md(md_text)
    bad, _ = check_problems(md, share, dry)
    if bad:
        raise SystemExit("問題の検査が不成立 — 止める:\n  " + "\n  ".join(bad))
    out_root = Path(out_root).resolve()
    d = out_root / EV.PREPS[arm]
    if d.exists():
        raise SystemExit(f"{d} が既にある (投入スクリプトが消してから作る)")
    RA = TM._load_problem_with_runs()
    prob = C / EV.PROBLEMS[arm]
    info = RA.prepare(prob, d, nsteps=EV.MAIN_NSTEPS, ic_from=None, cfl_main=EV.CFL_MAIN, implicit_relax=EV.RELAX)
    fails = [f"step 数の同期: {x}" for x in EV.step_problems(EV2.config_steps((d / "solverConfig.yaml").read_text()))]
    fails += EV.prepare_info_problems({**info, "ic": IC_PLACEHOLDER}, arm, md)
    ev = TM.wall_evidence(d)
    if ev.get("status") != "consistent":
        fails.append(f"壁の証拠 (M4) が {ev.get('status')}")
    q = (d / "MESH_QUALITY.txt").read_text() if (d / "MESH_QUALITY.txt").is_file() else ""
    if not any(l.strip().startswith("VERDICT: PASS") for l in q.splitlines()):
        fails.append("メッシュ品質が PASS でない")
    mp = RA.mesh_params_euler(RA.load_problem(prob), float(info["scale_m"]))
    if dataclasses.asdict(mp) != (info.get("mesh") or {}).get("params"):
        fails.append("prepare_info の mesh.params が問題の mesh_euler の解決値と違う")
    print(arm, prob.name, "| mesh", json.dumps({k: info["mesh"].get(k) for k in ("ni", "nj", "wall_first_frac", "source")}),
          "| moc", {k: (info.get("moc") or {}).get(k) for k in ("axis_limit", "corrector")}, "gate", ((info.get("moc") or {}).get("gate") or {}).get("pass"),
          "| M4", ev.get("status"), "|", q.strip().splitlines()[-1] if q.strip() else "(品質の記録なし)")
    same = compare_with_e4v(d, info, dry) if (share and arm == "B") else None
    if same:
        fails += same["failures"]
    thermo = TM.RUNS / E2RUN.DRY_THERMO_RUN if dry else d
    try:
        ic = E2RUN.ic_check(d, thermo, dry)
    except Exception as e:  # noqa: BLE001 — 検査できないことも不成立として記録する
        ic = {"VERDICT": "FAIL", "problems": [f"IC の検査を完了できない: {type(e).__name__}: {e}"]}
    fails += [f"起動前の IC の検査: {p}" for p in ic["problems"]]
    rec, rows = E2RUN.mesh_record(d, prob, mp, float(info["scale_m"]))
    rec.update(plan=PLAN, plan_reg_commit=EV.PLAN_REG_COMMIT, arm=arm,
               note="prepare_info.json の mesh 欄に加えた、断面ごとの実測の記録 (E2・E4 と同じ形式)",
               prepare_info_mesh_hashes=(info.get("mesh") or {}).get("hashes"))
    E2RUN.write_sections(d / MESH_CSV, rows)
    rec["hashes"]["sections_csv"] = _sha(d / MESH_CSV)
    rec.update(binaries=_binaries(), design_sha256={f: _sha(ROOT / f) for f in DESIGN_FILES},
               git={"head": E2RUN._git("rev-parse", "HEAD"), "status_design": E2RUN._git("status", "--porcelain", "--", "design")})
    (d / MESH_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    runs = EV.new_runs(share)[arm]
    summary = {"plan": PLAN, "plan_reg_commit": EV.PLAN_REG_COMMIT, "tool": "moc_v5d.prep", "arm": arm, "runs": runs, "dry": dry,
               "problem": prob.name, "problem_sha256": _sha(prob), "md_offset": md, "md_offset_repr": md_repr, "share_e4v": share,
               "moc": info.get("moc"), "x_E": info.get("x_E"), "wall_evidence": ev, "same_as_e4v": same, "ic_check": ic, "mesh": rec,
               "steps": {"MAIN_NSTEPS": EV.MAIN_NSTEPS, "OUT_INTERVAL": EV.OUT_INTERVAL,
                         "config": EV2.config_steps((d / "solverConfig.yaml").read_text())},
               "failures": fails, "VERDICT": "OK" if not fails else "REFUSED"}
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / f"_band_ab/moc_v5d_prep_{arm}{'_dry' if dry else ''}.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    if fails:
        raise RuntimeError("V5d の準備の検査が不成立 — 止める (forge を起動しない):\n  " + "\n  ".join(fails))
    info = json.loads((d / "prepare_info.json").read_text())
    info.update(stages="soft", plan=PLAN,
                v5d={"arm": arm, "runs": runs, "problem": prob.name, "md_offset": md, "md_offset_repr": md_repr, "share_e4v": share,
                     "plan_reg_commit": EV.PLAN_REG_COMMIT},
                ic={"mode": "isentropic", "tool": "paste_isentropic_ic (runner_axismach.prepare)", "VERDICT": ("OK" if not dry else "DRY"),
                    "check": {k: ic.get(k) for k in ("max_abs_dev_K", "T0_min", "T0_max", "thermo_source", "VERDICT")}})
    if dry:
        info["DRY"] = True
    (d / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    prec = {"plan": PLAN, "plan_reg_commit": EV.PLAN_REG_COMMIT, "arm": arm, "runs": runs, "problem": prob.name,
            "md_offset": md, "md_offset_repr": md_repr, "share_e4v": share, "dry": dry, "ic_check": ic, "same_as_e4v": same,
            "steps": summary["steps"], "files_sha256": _input_hashes(d)}
    (d / EV.PREP_JSON).write_text(json.dumps(prec, indent=1, ensure_ascii=False, default=float))
    tb = rec["summary_bands_rt"]
    print(f"{arm} 実効のメッシュ: 壁側 (/r_w) スロート {tb['-4<=x<=1'].get('wall_gap_rw')} 下流 x≥17 {tb['x>=17'].get('wall_gap_rw')} | "
          f"IC max|T0−1600| {ic['max_abs_dev_K']:.4g} K ({ic['VERDICT']}) | Md_moc_offset {md_repr} | x_E {info.get('x_E')!r}"
          + (f" | run_0164 との照合 {same['status']}" if same else ""))
    return summary


def cross_check(dB: Path, dM: Path, md_text: str, dry: bool = False) -> dict:
    """腕 B と腕 M の prep: 設定ファイル・格子の実効値・Md_moc_offset が同じで、壁は腕 M が腕 B と見分けられ自腕に一致する。"""
    md, md_repr = EV.parse_md(md_text)
    fails, rec = [], {"prep_B": str(dB), "prep_M": str(dM)}
    iB, iM = (json.loads((x / "prepare_info.json").read_text()) for x in (dB, dM))
    rec["files"] = {n: {"B": _sha(dB / n), "M": _sha(dM / n)} for n in EV.SAME_FILES}
    diff = [n for n, v in rec["files"].items() if (v["B"] or v["M"]) and v["B"] != v["M"]]
    if diff:
        fails.append(f"腕 B と腕 M の prep で違う設定ファイル: {diff}")
    for key, get in (("mesh.params", lambda i: (i.get("mesh") or {}).get("params")), ("Md_moc_offset", lambda i: i.get("Md_moc_offset")),
                     ("initial_line", _initial_line_key), ("scale_m", lambda i: i.get("scale_m")),
                     ("wall_fit.mono_r2", lambda i: (i.get("wall_fit") or {}).get("mono_r2"))):
        if get(iB) is None or get(iB) != get(iM):
            fails.append(f"腕 B と腕 M の {key} が違う・欠ける")
    if iB.get("Md_moc_offset") != md:
        fails.append(f"Md_moc_offset {iB.get('Md_moc_offset')!r} が --md-offset {md!r} と違う")
    for arm, info in (("B", iB), ("M", iM)):
        fails += [f"腕 {arm}: {p}" for p in EV.prepare_info_problems({**info, "DRY": None, "ic": IC_PLACEHOLDER}, arm, md)]
    try:
        vo = V5.wall_vs_other(dM, dB)
    except Exception as e:  # noqa: BLE001
        vo = {"status": "error", "reason": f"{type(e).__name__}: {e}"}
    rec["wall_M_vs_B"] = vo
    nd, nown = vo.get("n_discriminable"), vo.get("n_discriminable_matching_own")
    if vo.get("status") != "ok" or not (isinstance(nd, int) and nd > 0) or nown != nd:
        fails.append(f"腕 M の壁節点が腕 B と見分けられ、すべて自腕に一致する、が成り立たない ({vo.get('status')}, {nd}, {nown})")
    rec["x_E"] = {"B": iB.get("x_E"), "M": iM.get("x_E")}
    rec["X_F_design_m"] = {arm: V5._x_last_design(x) for arm, x in (("B", dB), ("M", dM))}
    rec.update(plan=PLAN, plan_reg_commit=EV.PLAN_REG_COMMIT, md_offset_repr=md_repr, dry=dry, failures=fails,
               VERDICT="OK" if not fails else "REFUSED")
    (C / f"_band_ab/moc_v5d_prep_cross{'_dry' if dry else ''}.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"cross-check: {rec['VERDICT']} | 壁 M − B 最大 {vo.get('analytic_dr_max_m')!r} m (x {vo.get('x_analytic_dr_max_rt')!r} r_t) "
          f"見分けられる節点 {nd} (自腕 {nown}) | x_E {rec['x_E']}")
    if fails:
        raise SystemExit("腕 B と腕 M の照合が不成立 — 止める:\n  " + "\n  ".join(fails))
    return rec


def verify_prep(prep_dir: Path, run_dirs) -> list:
    """run_dirs が prep_dir の複製のまま (V5D_PREP.json の入力の sha256 が一致) で、run 名が prep の腕の run か。"""
    bad = []
    try:
        rec = json.loads((prep_dir / EV.PREP_JSON).read_text())
    except (OSError, ValueError) as e:
        return [f"{EV.PREP_JSON} を読めない ({e})"]
    if rec.get("dry"):
        bad.append("乾式確認 (DRY) の prep")
    if (rec.get("ic_check") or {}).get("VERDICT") != "OK":
        bad.append(f"起動前の IC の検査が OK でない ({(rec.get('ic_check') or {}).get('VERDICT')!r})")
    if not rec.get("files_sha256"):
        bad.append("prep の入力の sha256 の記録が無い")
    for rd in run_dirs:
        rd = Path(rd)
        try:
            info = json.loads((rd / "prepare_info.json").read_text())
        except (OSError, ValueError) as e:
            bad.append(f"{rd.name}: prepare_info.json を読めない ({e})")
            continue
        if info.get("DRY"):
            bad.append(f"{rd.name}: 乾式確認の prepare_info")
        if rd.name not in (rec.get("runs") or []):
            bad.append(f"run の名前 {rd.name} が prep の腕の run {rec.get('runs')} に無い")
        for n, sha in (rec.get("files_sha256") or {}).items():
            if _sha(rd / n) != sha:
                bad.append(f"{rd.name}/{n} が prep の記録と違う")
    return bad


# --- 実行 (soft 段の出力の退避つき) ---------------------------------------------------------------------------------------------
def install_soft_saver(RA, rd: Path) -> list:
    """RA._restart_same_mesh を包む: 呼ばれる前に run 直下の res_* と段のログを rd/_soft_stage/ に写す (1 回だけ。soft 段の引き継ぎ)。
    e4_recal.install_soft_saver と同じ方式 (記録の plan・tool の名前だけが違う)。"""
    orig = RA._restart_same_mesh
    calls = []

    def wrapped(res_h5, mesh_h5):
        res_h5 = Path(res_h5)
        if calls:
            raise RuntimeError("段の引き継ぎが 2 回目 — 登録の段は soft → 本段だけ (mid 段は無い)")
        if res_h5.parent.resolve() != rd.resolve() or res_h5.name != f"res_{EV.SOFT_STEPS}.h5":
            raise RuntimeError(f"soft 段の最終場が {res_h5} (期待 {rd}/res_{EV.SOFT_STEPS}.h5)")
        dst = rd / EV4.SOFT_DIR
        dst.mkdir(exist_ok=False)
        rec = {"plan": PLAN, "tool": "moc_v5d.install_soft_saver", "stage": "soft", "restart_src": res_h5.name, "files": {}}
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
    try:
        return EV.step_problems(EV2.config_steps((rd / "solverConfig.yaml").read_text()))
    except Exception as e:  # noqa: BLE001 — 読めないことも食い違いとして止める
        return [f"solverConfig.yaml を読めない: {type(e).__name__}: {e}"]


def run(rd: Path) -> int:
    info = json.loads((rd / "prepare_info.json").read_text())
    if info.get("DRY"):
        raise SystemExit(f"{rd} は乾式確認 (--dry) の prep から作られている — 回さない")
    v = info.get("v5d") or {}
    arm, share = v.get("arm"), v.get("share_e4v")
    if arm not in ("B", "M") or not isinstance(share, bool) or rd.name not in EV.new_runs(share)[arm]:
        raise SystemExit(f"{rd}: V5d の腕 {arm!r}・共用 {share!r} の run 名でない — 回さない")
    if (info.get("ic") or {}).get("VERDICT") != "OK":
        raise SystemExit(f"{rd}: IC の VERDICT が OK でない ({(info.get('ic') or {}).get('VERDICT')!r}) — 回さない")
    md, md_repr = EV.parse_md(v.get("md_offset_repr"))
    bad = EV.prepare_info_problems(info, arm, md)
    if bad or md_repr != v.get("md_offset_repr"):
        raise SystemExit(f"{rd}: prepare_info の固定の条件が不成立 — 回さない: {bad}")
    if (rd / EV4.SOFT_DIR).exists() or any(rd.glob("res_*")):
        raise SystemExit(f"{rd}: 既に出力がある (_soft_stage または res_*) — 回さない")
    sp = check_run_steps(rd)
    if sp:
        raise SystemExit(f"{rd}: 本段の step 数が評価器と同期していない — 回さない:\n  " + "\n  ".join(sp))
    from forge_design.evaluate import runner_axismach as RA
    calls = install_soft_saver(RA, rd)
    rc = RA.run_staged(rd, cfl_main=EV.CFL_MAIN, mid_stage=False, stages="soft")
    if len(calls) != 1:
        raise RuntimeError(f"soft 段の出力の退避が {len(calls)} 回 (1 回のはず)")
    print(f"forge exit={rc}")
    return rc


# --- 投入スクリプト用 -----------------------------------------------------------------------------------------------------------
def launch_order(share: bool) -> list:
    """起動の順: 腕 B と腕 M を交互に (並列の 1 束に両腕が入るように。時間帯・GPU の混み具合が片方の腕に偏らない)。"""
    nr = EV.new_runs(share)
    out = []
    for i in range(max(len(nr["B"]), len(nr["M"]))):
        out += [x[i] for x in (nr["B"], nr["M"]) if i < len(x)]
    return out


def plan_runs(share: bool) -> str:
    a = EV.arms(share)
    nr = EV.new_runs(share)
    lines = [f"ORDER {' '.join(launch_order(share))}", f"ALL {','.join(a['B'] + a['M'])}", f"GEOM_REF {EV.GEOM_REF}",
             f"SHARED {EV.E4V_RUN if share else '-'}"]
    lines += [f"RUN {r} {arm} {EV.PREPS[arm]}" for arm in ("B", "M") for r in nr[arm]]
    return "\n".join(lines)


def launch_record(md_text: str, share: bool, dry: bool, extra: dict | None = None) -> Path:
    md, md_repr = EV.parse_md(md_text)
    rec = {"plan": PLAN, "plan_reg_commit": EV.PLAN_REG_COMMIT, "date": datetime.now(timezone.utc).isoformat(), "dry": dry,
           "md_offset": md, "md_offset_repr": md_repr, "share_e4v": share, "arms": EV.arms(share), "new_runs": EV.new_runs(share),
           "order": launch_order(share), "geom_ref": EV.GEOM_REF,
           "steps": {"soft": EV.SOFT_STEPS, "main": EV.MAIN_NSTEPS, "out": EV.OUT_INTERVAL, "win13": [EV.WIN13[0], EV.WIN13[-1]],
                     "tail5": [EV.TAIL5[0], EV.TAIL5[-1]]},
           "scripts_sha256": {n: _sha(C / n) for n in ("moc_v5d.py", "moc_v5d_eval.py", "run_moc_v5d.sh", "eval_wallfit_euler.py",
                                                      "e4_recal_eval.py", "moc_v5_euler_eval.py", "euler_t0_e2_eval.py", "euler_t0_e2.py")},
           "problems_sha256": {n: _sha(C / n) for n in EV.PROBLEMS.values()},
           "git": {"head": E2RUN._git("rev-parse", "HEAD"),
                   "status": E2RUN._git("status", "--porcelain", "--", "design", str(C / "moc_v5d.py"), str(C / "moc_v5d_eval.py"),
                                        str(C / "run_moc_v5d.sh"))},
           **(extra or {})}
    p = C / (EV.LAUNCH_JSON.replace(".json", "_dry.json") if dry else EV.LAUNCH_JSON)
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    return p


def main(argv) -> int:
    cmd = argv[0] if argv else ""
    ap = argparse.ArgumentParser(prog=f"moc_v5d.py {cmd}")
    if cmd in ("make-problems", "check-problems", "check-share", "prep", "cross-check", "launch-record"):
        ap.add_argument("--md-offset", required=True)
        ap.add_argument("--dry", action="store_true")
    if cmd in ("make-problems", "check-problems", "prep", "plan-runs", "launch-record"):
        ap.add_argument("--share-e4v", action="store_true")
    if cmd == "prep":
        ap.add_argument("arm", choices=("B", "M"))
        ap.add_argument("--out-root", default=str(C))
    if cmd == "cross-check":
        ap.add_argument("prep_B")
        ap.add_argument("prep_M")
    if cmd == "launch-record":
        ap.add_argument("--extra", default=None, help="起動の記録に足す JSON (投入スクリプトのバイナリの sha256 など)")
    if cmd == "make-problems":
        a = ap.parse_args(argv[1:])
        rec = make_problems(a.md_offset, a.share_e4v, a.dry)
        print(f"MAKE-PROBLEMS: OK ({rec['md_offset_repr']}、書いた {rec['written'] or 'なし (同じ中身)'}; E4 の採用 "
              f"{'OK' if rec['e4_adoption'].get('ok') else rec['e4_adoption'].get('status', '不成立')})")
        return 0
    if cmd == "check-problems":
        a = ap.parse_args(argv[1:])
        md, _ = EV.parse_md(a.md_offset)
        bad, _ = check_problems(md, a.share_e4v, a.dry)
        print("CHECK-PROBLEMS: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if cmd == "check-share":
        a = ap.parse_args(argv[1:])
        md, _ = EV.parse_md(a.md_offset)
        bad = check_share(md, a.dry)
        print("CHECK-SHARE: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad))
              + (" (DRY: run_0164 が無ければ照合しない)" if a.dry and not bad else ""))
        return 0 if not bad else 1
    if cmd == "prep":
        a = ap.parse_args(argv[1:])
        prep(a.arm, a.md_offset, a.share_e4v, Path(a.out_root), dry=a.dry)
        return 0
    if cmd == "cross-check":
        a = ap.parse_args(argv[1:])
        cross_check(Path(a.prep_B).resolve(), Path(a.prep_M).resolve(), a.md_offset, a.dry)
        return 0
    if cmd == "verify-prep" and len(argv) >= 3:
        bad = verify_prep(Path(argv[1]).resolve(), [Path(x).resolve() for x in argv[2:]])
        print("VERIFY-PREP: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if cmd == "run" and len(argv) == 2:
        return run(Path(argv[1]).resolve())
    if cmd == "plan-runs":
        a = ap.parse_args(argv[1:])
        print(plan_runs(a.share_e4v))
        return 0
    if cmd == "launch-record":
        a = ap.parse_args(argv[1:])
        print("LAUNCH-RECORD:", launch_record(a.md_offset, a.share_e4v, a.dry, json.loads(a.extra) if a.extra else None))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
