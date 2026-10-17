#!/usr/bin/env python3
"""承認の記録 (`forge_design.intake.approval`) と判定の正式な入口 (`wt_qualify --approval`) の試験
(plan tooling-design-orchestrator §4.7.1 (1)〜(4)、§5.1 #2a)。正例と負例を並べる (全部止める実装を合格にしない)。

固定の入力: design/tests/fixtures/approval_trial_a/ (手順書の試走 A のカード・問題 YAML・NS p0 の停止規則の基準)。
試験はカードを版 2 に書き換えて一時ディレクトリで使う。

  正例
  (P1) 版 2 のカード (出所・予算つき) の preview に問題が無く、write で承認の記録ができる
  (P2) 承認の原本のカードと問題 YAML の check が OK (original)
  (P3) 較正の係数だけを変えた試行の入力 (card_p1 + approval の参照、problem_p1) の check が OK (trial、再承認なし)
  (P4) run の停止規則の基準の check_run が OK、承認の値の公差・窓 (approved_tolerances) が基準の線と同じ
  負例
  (N1) 出所の無いユーザ欄 (N2) 出所の未知のパス (N3) recipe_default なのに既定と違う値 (N4) 予算の欠落 (N5) 出所 proposal
  (N6) requirements_without_gate の削除 (N7) CLI の公差 (値が同じでも)・--optional-parts・反復の上限の指定、--iteration の欠落
  (N8) 承認の後の line_n の変更 (N9) 停止規則の閾値の変更 (N10) 物性 (輸送) の変更 (N11) 問題 YAML の qualification
  (N12) 承認の原本に approval の参照 (N13) 試行の入力に参照が無い・sha256 が違う (N14) write の拒否 (既存・問題あり・改ざん)
  (N15) 承認の版より古い生成器の停止規則の基準 (実物の run_0003; residual_finite.required が無い) は check_run が MISMATCH
  判別 A/B (codex 2026-10-17)
  (i)  承認値 0.001 に CLI 0.1 → 拒否 (承認なしの今の経路は 0.1 を採用する、を併記)
  run の物理の照合 (§5.1 #3): (R0) 試走 A の本段の設定が承認と一致 (正例)、(R1〜R14) 入口の全温・全圧・組成・分布、壁温・断熱壁、
  出口の静圧、軸対称、H2O の輸送、粘性のモデル、熱力学の基準温度、壁関数・圧縮性補正・エネルギーの k を変えると MISMATCH、
  (R+) 数値の設定 (CFL・リミッタ) だけの変更は止めない
  (ii) line_n 2001 を承認し、run の基準 (今のコードで作り直した) の line_mach.n だけを A 2001 / B 2 → A は OK、B は MISMATCH
       (A も止まれば試験の不成立)

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_approval_tests.py
所要: 数秒。
"""
import copy
import json
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.feedback import wt_qualify as Q  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures" / "approval_trial_a"
FAIL = 0
CLEAN_REPO = None      # 判定器の版の記録に使う、きれいな一時の git リポジトリ (この試験は作業ツリーの状態に依らない; 判定専用の環境は run_gate_env_tests.py)


def clean_repo(td: Path) -> Path:
    import subprocess
    r = td / "clean_repo"
    r.mkdir()
    for args in (["init", "-q"], ["config", "user.email", "t@t"], ["config", "user.name", "t"]):
        subprocess.run(["git", "-C", str(r), *args], check=True)
    (r / "x.txt").write_text("x")
    subprocess.run(["git", "-C", str(r), "add", "x.txt"], check=True)
    subprocess.run(["git", "-C", str(r), "commit", "-q", "-m", "x"], check=True)
    return r


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def raises(fn, needle=None):
    try:
        fn()
    except ValueError as e:
        return needle is None or needle in str(e)
    return False


USER_REF = "ユーザ「この表でよい」(2026-10-17、試験の仮の承認)"
SOURCES = {"gas": {"kind": "user", "ref": "依頼の原文の組成"}, "conditions": {"kind": "user", "ref": "依頼の原文の入口条件"},
           "M_exit": {"kind": "user", "ref": "依頼の原文"}, "throat": {"kind": "user", "ref": "依頼の原文"},
           "inlet": {"kind": "user", "ref": "依頼の原文"}, "divergent_length_m": {"kind": "user", "ref": "依頼の原文"},
           "wall": {"kind": "user", "ref": "依頼の原文の壁温"}, "p_ambient": {"kind": "recipe_default"},
           "tolerances": {"kind": "recipe_default"}, "budget": {"kind": "user", "ref": "見積もりの承認 (GPU 4 h・50 GB)"}}


def v2(card_name, **edit) -> dict:
    d = yaml.safe_load((FIX / card_name).read_text())
    d["format"] = W.CARD_FORMAT_V2
    d["budget"] = {"gpu_hours": 4.0, "disk_gb": 50.0}
    d["sources"] = copy.deepcopy(SOURCES)
    for k, v in edit.items():
        if v is None:
            d.pop(k, None)
        else:
            d[k] = v
    return d


def put(td: Path, name: str, doc) -> Path:
    p = td / name
    if isinstance(doc, str):
        p.write_text(doc)
    else:
        p.write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True))
    return p


def prob(name="problem_p0.yaml", **edit):
    d = yaml.safe_load((FIX / name).read_text())
    for path, v in edit.items():
        cur = d
        ks = path.split("__")
        for k in ks[:-1]:
            cur = cur.setdefault(k, {})
        cur[ks[-1]] = v
    return d


def approve(td: Path, card_doc, prob_doc=None, tag="a") -> tuple[Path, Path, Path]:
    c = put(td, f"card_{tag}.yaml", card_doc)
    p = put(td, f"problem_{tag}.yaml", prob_doc if prob_doc is not None else prob())
    pv = AP.preview(c, p, repo_dir=CLEAN_REPO, judge_env=False)
    pvp = td / f"preview_{tag}.json"
    pvp.write_text(json.dumps(pv, ensure_ascii=False))
    if pv["problems"]:
        raise AssertionError(f"preview の問題 {pv['problems']}")
    ap = td / f"card_approval_{tag}.json"
    AP.write(pvp, USER_REF, ap)
    return c, p, ap


def preview_problems(td: Path, card_doc, prob_doc=None, tag="n") -> list:
    c = put(td, f"card_{tag}.yaml", card_doc)
    p = put(td, f"problem_{tag}.yaml", prob_doc if prob_doc is not None else prob())
    return AP.preview(c, p, repo_dir=CLEAN_REPO, judge_env=False)["problems"]


def trial(td: Path, ap: Path, name="card_p1.yaml", tag="t", **edit) -> Path:
    d = v2(name, **edit)
    d["approval"] = {"file": str(ap), "sha256": AP.sha256_file(ap)}
    return put(td, f"card_{tag}.yaml", d)


def copy_run(td: Path, tag: str, regenerate: bool = True) -> Path:
    """固定の入力の run を写す。regenerate なら停止規則の基準を今のコードで作り直す (承認の版で作った run の代わり; 実物の基準は
    10-17 に足した residual_finite.required の前の生成器のもの)。"""
    from forge_design.evaluate import settle_presets as SP
    r = td / f"run_{tag}"
    shutil.copytree(FIX / "run_0003_ns_p0", r)
    if regenerate:
        old = yaml.safe_load((r / "settle_criteria.yaml").read_text())
        es = (yaml.safe_load((FIX / "problem_p0.yaml").read_text())["evaluate"]["settle"])
        crit = SP.resolve(es, json.loads((r / "prepare_info.json").read_text()), yaml.safe_load((r / "solverConfig.yaml").read_text()),
                          yaml.safe_load((r / "bcondConfig.yaml").read_text()), segment=(old.get("grid") or {}).get("segment"))
        (r / "settle_criteria.yaml").write_text(yaml.safe_dump(crit, sort_keys=False, allow_unicode=True))
    return r


with tempfile.TemporaryDirectory() as tdn:
    td = Path(tdn)
    CLEAN_REPO = clean_repo(td)
    print("== 正例")
    c0, p0, ap0 = approve(td, v2("card_p0.yaml"))
    apd = json.loads(ap0.read_text())
    check("(P1) preview に問題が無く write で承認の記録ができる", apd["format"] == AP.FORMAT and apd["user_ref"] == USER_REF,
          f"protected sha256 {apd['protected_sha256'][:16]}…")
    check("(P1) 承認の記録に既定の値が出所つきで入る (tolerances・line_n・run_policy・mesh_limits)",
          apd["protected"]["qualification"]["line_n"] == 2001 and apd["protected"]["run_policy"] == {"max_iterations": 3, "max_recalibrations": 1}
          and apd["sources"]["qualification.line_n"]["kind"] == "recipe_default"
          and apd["sources"]["tolerances.length_m"]["kind"] == "recipe_default",
          f"{apd['protected']['qualification']}")
    r = AP.check(ap0, c0, p0)
    check("(P2) 承認の原本の check が OK (original)", r["status"] == "OK" and r["kind"] == "original", str(r["diffs"] + r["problems"]))
    p1 = put(td, "problem_p1.yaml", prob("problem_p1.yaml"))
    t1 = trial(td, ap0)
    r = AP.check(ap0, t1, p1)
    check("(P3) 較正の係数だけを変えた試行の入力の check が OK (trial、再承認なし)", r["status"] == "OK" and r["kind"] == "trial",
          str(r["diffs"] + r["problems"]))
    run_a = copy_run(td, "A")
    rr = AP.check_run(ap0, run_a)
    check("(P4) run の停止規則の基準 (今のコードで作り直した) の check_run が OK", rr["status"] == "OK", str(rr["problems"] + rr["diffs"])[:300])
    run_old = copy_run(td, "old", regenerate=False)
    ro = AP.check_run(ap0, run_old)
    check("(N15) 実物の run_0003 の基準 (古い生成器) は MISMATCH、違いは residual_finite.required だけ",
          ro["status"] == "MISMATCH" and ro["diffs"] == ["quantities.residual_finite.required (今の値に無い)"] and not ro["problems"],
          str(ro["problems"] + ro["diffs"])[:300])
    info = json.loads((run_a / "prepare_info.json").read_text())
    tol = Q.approved_tolerances(run_a, info, AP.load_approval(ap0), ap0)
    crit = yaml.safe_load((run_a / "settle_criteria.yaml").read_text())
    lm = next(q for q in crit["quantities"] if q["name"] == "line_mach")
    check("(P4) 承認の値の公差・窓が基準の線と同じ (η・点数・x の窓)",
          tol["line_n"]["value"] == lm["n"] == 2001 and tol["line_eta"]["value"] == lm["eta"]
          and all(abs(a - b) <= 1e-12 for a, b in zip(tol["line_x_range_m"]["value"], lm["x_range_m"]))
          and tol["exit_core_M_tol_rel"]["source"].startswith("承認の記録"),
          f"{tol['line_x_range_m']['value']} / {lm['x_range_m']}")

    print("== 負例 (preview)")
    d = v2("card_p0.yaml")
    del d["sources"]["conditions"]
    pr = preview_problems(td, d, tag="n1")
    check("(N1) 出所の無いユーザ欄で止まる", any("conditions.Pt: 出所が無い" in x for x in pr), str(pr)[:200])
    d = v2("card_p0.yaml")
    d["sources"]["conditions.Px"] = {"kind": "user", "ref": "x"}
    pr = preview_problems(td, d, tag="n2")
    check("(N2) 出所の未知のパスで止まる", any("conditions.Px: 未知のパス" in x for x in pr), str(pr)[:200])
    d = v2("card_p0.yaml", tolerances={"length_m": 2.0e-3, "diameter_m": 2.0e-5})
    pr = preview_problems(td, d, tag="n3")
    check("(N3) recipe_default なのに既定と違う値で止まる", any("tolerances.length_m: recipe_default だが値" in x for x in pr), str(pr)[:200])
    d = v2("card_p0.yaml", budget=None)
    del d["sources"]["budget"]
    pr = preview_problems(td, d, tag="n4")
    check("(N4) 予算の欠落で止まる", any("budget" in x and "必須" in x for x in pr), str(pr)[:200])
    d = v2("card_p0.yaml")
    d["sources"]["wall"] = {"kind": "proposal", "ref": "AI の提案"}
    pr = preview_problems(td, d, tag="n5")
    check("(N5) 出所 proposal で止まる", any("wall.Tw: 出所 'proposal'" in x for x in pr), str(pr)[:200])
    d = v2("card_p0.yaml", approval={"file": "x.json", "sha256": "0" * 64})
    pr = preview_problems(td, d, tag="n12")
    check("(N12) 承認の原本に approval の参照があると止まる", any("approval" in x and "原本" in x for x in pr), str(pr)[:200])
    pr = preview_problems(td, v2("card_p0.yaml"), prob(qualification={"exit_core_M_tol_rel": 0.1}), tag="n11p")
    check("(N11) 問題 YAML の qualification で止まる (preview)", any("qualification" in x for x in pr), str(pr)[:200])

    print("== 負例 (承認の後の変更; check)")
    rq = [{"name": "test_section_flow", "method": "未定", "tolerance": "未定", "evidence": "未定", "judge": "人"}]
    d6 = v2("card_p0.yaml", requirements_without_gate=rq)
    d6["sources"]["requirements_without_gate"] = {"kind": "user", "ref": "依頼の原文"}
    c6, p6, ap6 = approve(td, d6, tag="n6")
    t6 = trial(td, ap6, tag="t6")
    r = AP.check(ap6, t6, p6)
    check("(N6) requirements_without_gate の削除で MISMATCH", r["status"] == "MISMATCH" and any("requirements_without_gate" in x for x in r["diffs"]),
          str(r["diffs"])[:200])
    t8 = trial(td, ap0, tag="t8", qualification={"line_n": 2})
    r = AP.check(ap0, t8, p1)
    check("(N8) 承認の後の line_n の変更で MISMATCH", r["status"] == "MISMATCH" and "qualification.line_n" in r["diffs"], str(r["diffs"])[:200])
    p9 = put(td, "problem_n9.yaml", prob("problem_p1.yaml", evaluate__settle={"preset": "wt_nozzle_v3",
                                                                               "thresholds": {"line_mach": {"tol_maxabs": 0.05}}}))
    r = AP.check(ap0, t1, p9)
    check("(N9) 停止規則の閾値の変更で MISMATCH", r["status"] == "MISMATCH" and any(x.startswith("evaluation.settle") for x in r["diffs"]),
          str(r["diffs"])[:200])
    g = prob("problem_p1.yaml")["gas"]
    g["transport"]["H2O"] = "cea"
    p10 = put(td, "problem_n10.yaml", prob("problem_p1.yaml", gas=g))
    r = AP.check(ap0, t1, p10)
    check("(N10) 物性 (H2O の輸送) の変更で MISMATCH", r["status"] == "MISMATCH" and "physics.gas_transport.H2O" in r["diffs"], str(r["diffs"])[:200])
    p11 = put(td, "problem_n11.yaml", prob("problem_p1.yaml", qualification={"exit_core_M_tol_rel": 0.1}))
    r = AP.check(ap0, t1, p11)
    check("(N11) 問題 YAML の qualification で MISMATCH (check)", r["status"] == "MISMATCH" and any("qualification" in x for x in r["problems"]),
          str(r["problems"])[:200])
    d13 = v2("card_p1.yaml")
    t13 = put(td, "card_t13.yaml", d13)
    r = AP.check(ap0, t13, p1)
    check("(N13) 試行の入力に approval の参照が無いと MISMATCH", r["status"] == "MISMATCH" and any("参照) が無い" in x for x in r["problems"]),
          str(r["problems"])[:200])
    d13["approval"] = {"file": str(ap0), "sha256": "f" * 64}
    t13b = put(td, "card_t13b.yaml", d13)
    r = AP.check(ap0, t13b, p1)
    check("(N13) 参照の sha256 が違うと MISMATCH", r["status"] == "MISMATCH" and any("sha256" in x for x in r["problems"]), str(r["problems"])[:200])
    t13c = trial(td, ap0, tag="t13c", conditions={"Pt": 3.2e6, "Tt": 1187.5})
    r = AP.check(ap0, t13c, put(td, "problem_n13c.yaml", prob("problem_p1.yaml", spec__Pt=3.2e6)))
    check("(N13) 試行の入力で入口の全圧を変えると MISMATCH", r["status"] == "MISMATCH" and "requirements.conditions.Pt" in r["diffs"],
          str(r["diffs"])[:200])

    print("== 負例 (write)")
    check("(N14) 承認の記録は上書きしない", raises(lambda: AP.write(td / "preview_a.json", USER_REF, ap0), "既にある"))
    bad = AP.preview(put(td, "card_n14.yaml", v2("card_p0.yaml", budget=None)), p0, repo_dir=CLEAN_REPO, judge_env=False)
    (td / "preview_n14.json").write_text(json.dumps(bad, ensure_ascii=False))
    check("(N14) 問題のある preview は write しない", raises(lambda: AP.write(td / "preview_n14.json", USER_REF, td / "ap_n14.json"), "承認できない"))
    pv = json.loads((td / "preview_a.json").read_text())
    pv["protected"]["qualification"]["exit_core_M_tol_rel"] = 0.1
    (td / "preview_tamper.json").write_text(json.dumps(pv, ensure_ascii=False))
    check("(N14) 書き換えた preview は write しない", raises(lambda: AP.write(td / "preview_tamper.json", USER_REF, td / "ap_t.json"), "合わない"))
    check("(N14) user_ref が空なら write しない", raises(lambda: AP.write(td / "preview_a.json", " ", td / "ap_e.json"), "user_ref"))

    print("== 負例 (判定の正式な入口 wt_qualify.approval_mode)")
    none_cli = {"exit_core_M_tol_rel": None, "line_eta": None, "line_x_range_m": None, "line_tol_maxabs_pct": None, "line_tol_halfamp_pct": None}
    ok = Q.approval_mode(ap0, t1, p1, none_cli, (), 1, None, None)
    check("正例: CLI の指定なし・--iteration ありなら通り、上限は承認の値", ok[1:] == (3, 1, ()), str(ok[1:]))
    check("(N7) CLI の公差 (承認と同じ値 0.001) も拒否",
          raises(lambda: Q.approval_mode(ap0, t1, p1, {**none_cli, "exit_core_M_tol_rel": 1e-3}, (), 1, None, None), "CLI の公差"))
    check("(N7) --optional-parts を拒否", raises(lambda: Q.approval_mode(ap0, t1, p1, none_cli, ("excel",), 1, None, None), "optional-parts"))
    check("(N7) --max-iterations を拒否", raises(lambda: Q.approval_mode(ap0, t1, p1, none_cli, (), 1, 3, None), "max-iterations"))
    check("(N7) --iteration の欠落を拒否", raises(lambda: Q.approval_mode(ap0, t1, p1, none_cli, (), None, None, None), "--iteration"))
    check("(N8) 承認と違う試行の入力 (line_n 2) を判定の入口で拒否",
          raises(lambda: Q.approval_mode(ap0, t8, p1, none_cli, (), 1, None, None), "承認の記録と違う"))

    print("== 判別 A/B")
    cli01 = {**none_cli, "exit_core_M_tol_rel": 0.1}
    check("(i) 承認値 0.001 に CLI 0.1 → 拒否", raises(lambda: Q.approval_mode(ap0, t1, p1, cli01, (), 1, None, None), "CLI の公差"))
    base = Q.resolve_tolerances(run_a, info, yaml.safe_load(c0.read_text()), None, cli01)
    print(f"     (参考) 承認なしの今の経路は CLI を採用する: {base['exit_core_M_tol_rel']}")
    run_b = copy_run(td, "B")
    cb = yaml.safe_load((run_b / "settle_criteria.yaml").read_text())
    for q in cb["quantities"]:
        if q["name"] == "line_mach":
            q["n"] = 2
    (run_b / "settle_criteria.yaml").write_text(yaml.safe_dump(cb, sort_keys=False, allow_unicode=True))
    ra, rb = AP.check_run(ap0, run_a), AP.check_run(ap0, run_b)
    check("(ii) A (n 2001) は OK、B (n 2) は MISMATCH",
          ra["status"] == "OK" and rb["status"] == "MISMATCH" and any("line_mach.n" in x for x in rb["diffs"]),
          f"A {ra['status']}, B {rb['status']} {rb['diffs'][:3]}")
    if ra["status"] != "OK":
        print("     試験の不成立: A も止まった")

    print("== run の物理の照合 (§5.1 #3)")

    def phys_with(tag, edit):
        r_ = copy_run(td, tag)
        for fname, fn in edit.items():
            doc = yaml.safe_load((r_ / fname).read_text())
            fn(doc)
            (r_ / fname).write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True))
        return AP.check_run_physics(AP.load_approval(ap0), r_, p0)

    ph = AP.check_run_physics(AP.load_approval(ap0), run_a, p0)
    check("(R0) 正例: 試走 A の NS p0 の本段の設定が承認の物理と一致", ph["status"] == "OK", str(ph["problems"])[:300])
    check("(R0) 照合した項目 (入口 Pt・Tt・組成・熱力学・輸送・乱流・軸対称・壁・出口)",
          {"inlet.Pt", "inlet.Tt", "physProp.transport", "turbulence", "mesh.isAxisymmetric", "wall", "outlet.Ps"} <= set(ph["checked"]),
          str(ph["checked"]))
    neg = {
        "入口の全温 1200 K": {"bcondConfig.yaml": lambda d: d["inlet"]["floats"].update(Tt=1200.0)},
        "入口の全圧 3.2 MPa": {"bcondConfig.yaml": lambda d: d["inlet"]["floats"].update(Pt=3.2e6)},
        "入口の組成 (H2O 0.09)": {"bcondConfig.yaml": lambda d: d["inlet"]["floats"].update(Y0=0.91, Y1=0.09)},
        "入口の分布 (inletProfile)": {"bcondConfig.yaml": lambda d: d["inlet"].update(ints={"inletProfile": 1})},
        "壁温 400 K": {"bcondConfig.yaml": lambda d: d["wall"]["floats"].update(Ts=400.0)},
        "断熱壁": {"bcondConfig.yaml": lambda d: d["wall"].update(kind="wall", floats=None)},
        "出口の静圧": {"bcondConfig.yaml": lambda d: d["outlet"]["floats"].update(Ps=3000.0)},
        "平面 (isAxisymmetric 0)": {"solverConfig.yaml": lambda d: d["mesh"].update(isAxisymmetric=0)},
        "H2O の輸送 cea": {"solverConfig.yaml": lambda d: d["physProp"]["transport"].update(H2O="cea")},
        "粘性のモデル (viscMethod 1)": {"solverConfig.yaml": lambda d: d["physProp"].update(viscMethod=1)},
        "熱力学 (thermoHrefTemp 0)": {"solverConfig.yaml": lambda d: d["physProp"].update(thermoHrefTemp=0.0)},
        "壁関数 (wallTreatmentSST 1)": {"solverConfig.yaml": lambda d: d["turbulence"].update(wallTreatmentSST=1)},
        "圧縮性補正 (dilatationCorrection 0)": {"solverConfig.yaml": lambda d: d["turbulence"].update(dilatationCorrection=0)},
        "エネルギーに k (sstEnergyIncludesK 1)": {"solverConfig.yaml": lambda d: d["turbulence"].update(sstEnergyIncludesK=1)},
    }
    for i, (nm, ed) in enumerate(neg.items()):
        ph = phys_with(f"R{i + 1}", ed)
        check(f"(R{i + 1}) {nm} → MISMATCH", ph["status"] == "MISMATCH", str(ph["problems"])[:200])
    r_num = copy_run(td, "Rnum")
    cfg = yaml.safe_load((r_num / "solverConfig.yaml").read_text())
    cfg["time"]["deltaT"].update(cfl=2.0, cfl_pseudo=2.0)
    cfg["space"]["limiter"] = 1
    (r_num / "solverConfig.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))
    ph = AP.check_run_physics(AP.load_approval(ap0), r_num, p0)
    check("(R+) 数値の設定 (CFL・リミッタ) だけを変えた run は止めない (AI が直してよい範囲)", ph["status"] == "OK", str(ph["problems"])[:200])
    rr = AP.check_run(ap0, run_a, p0)
    check("(R+) check_run (停止規則の基準 + 物理) が正例で OK", rr["status"] == "OK", str(rr["problems"] + rr["diffs"])[:200])

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
