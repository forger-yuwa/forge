#!/usr/bin/env python3
"""実行の計画 (`forge_design.campaign.plan`、`runner_axismach.staged_ns_launch_configs`・`run_staged_ns(plan=…)`) の試験
(plan tooling-design-orchestrator §4.1.1 の 2・3、§5.1 #4b)。GPU・forge なし: 起動を記録するだけのダミーの forge を使う。

手元の case/71 run_0014_chain_v2 の写し (設定・準備の記録・停止規則の基準・節点の表; 格子 nozzle.h5 は symlink、場 res_* は持たない) で回す。

  (P1) 準備の後の計画の確定: 段の列 S1_soft・S2_mid・main、前の段と本段の差は recipe の許す差だけ、本段の nStepOuter は停止規則の上限
  (P2) 計画の各段の設定の FNV が run_0014 の実際の段の記録 (stage_manifest.json の cfg_fnv) と一致
  判別の A/B (codex 2026-10-17): 計画を確定した後、
  (A)  変更なし → 計画の照合を通り forge (ダミー) を 1 回起動する (場を書かないので段は失敗する — 試験の対象は起動まで)
  (B)  S1_soft の設定の CFL だけを計画の値から変える (段の設定を作る関数を差し替える) → 起動 0 回で止まる
       A も止まれば試験の不成立
  (N1) 確定した計画を書き換える → 読み込みで止まる
  (N2) 計画をもう一度確定する → 上書きしない
  (N3) 前の段の物理 (乱流の dilatationCorrection) が本段と違う → 確定しない
  (N4) 計画の後に境界条件を変える → 起動 0 回で止まる
  (N6) 計画の確定の後に停止規則の基準を変える → 起動 0 回で止まる (旧の台本の「準備の後に変わっていない」の検査)
  (N5) 走った後の照合: 段の記録が計画と違う (S1 だけで止まった run) → MISMATCH。計画どおりの段の記録・起動の記録なら OK

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_execution_plan_tests.py
所要: 数十秒。run_0014 が無ければ飛ばす (rc 0 と「飛ばした」)。
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from forge_design.campaign import plan as EP  # noqa: E402
from forge_design.evaluate import runner as RN  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402

R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", HERE.parents[1] / "case" / "71.m55_combustion_wt" / "run_0014_chain_v2"))
C71 = R14.parent
FAIL = 0
SKIP_FILES = {"STOP_REQUEST", "STOPPED_AT", "SNAPSHOT_REQUEST", "MAIN_STAGE_STARTED", "STOP_RULE_STATE.json", "execution_plan.json"}


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


if not (R14 / "solverConfig.yaml").is_file() or not (R14 / RA.PRESTAGE_CONFIG).is_file():
    print(f"飛ばした: {R14} に run が無い")
    sys.exit(0)

td = Path(tempfile.mkdtemp(prefix="exec_plan_tests_"))
launch_log = td / "launches.txt"
dummy = td / "dummy_forge.sh"
dummy.write_text(f"#!/bin/bash\necho \"$PWD\" >> {launch_log}\nexit 0\n")
dummy.chmod(0o755)
RN._ENV["FORGE_BIN"] = str(dummy)          # runner の起動の環境 (import の時点の写し) の forge をダミーに
os.environ["FORGE_BIN"] = str(dummy)
appr = td / "card_approval.json"
appr.write_text("{}")


def copy_run(tag: str) -> Path:
    r = td / f"run_{tag}"
    r.mkdir()
    for f in R14.iterdir():
        if f.name in SKIP_FILES or (f.name.startswith("res_") and f.suffix in (".h5", ".xmf")) or f.is_dir():
            continue
        if f.suffix == ".h5":
            (r / f.name).symlink_to(f)
        else:
            shutil.copy2(f, r / f.name)
    return r


def launches() -> int:
    return len(launch_log.read_text().splitlines()) if launch_log.exists() else 0


def inputs_for(run: Path) -> dict:
    return EP.fixed_inputs(card=C71 / "card_m55_3coef.yaml", problem=C71 / "problem_m55_p5_settle_v2.yaml", approval=appr,
                           mesh_msh=run / "nozzle.msh", forge_bin=dummy)


try:
    print("== 計画の確定")
    ra = copy_run("A")
    doc = EP.confirm(ra, inputs_for(ra), attempt_id="A")
    st = doc["confirmed"]["stages"]
    check("(P1) 段の列 S1_soft・S2_mid・main、本段の nStepOuter は停止規則の上限",
          [s["tag"] for s in st] == ["S1_soft", "S2_mid", "main"] and st[-1]["nsteps"] == doc["confirmed"]["stop_rule_cap"],
          f"{[(s['tag'], s['nsteps']) for s in st]}")
    orig_sm = json.loads((R14 / "stage_manifest.json").read_text())
    check("(P2) 計画の各段の設定 (FNV) が、run_0014 が実際に走った段の記録の cfg_fnv と一致 (段の設定を作り直すと同じ)",
          [(s["tag"], s["fnv"]) for s in st] == [(x["tag"], x["cfg_fnv"]) for x in orig_sm["stages"]],
          f"{[(s['tag'], s['fnv']) for s in st]}")
    diffs = sorted({d for s in st[:-1] for d in EP._diff_paths(s["config"], st[-1]["config"])})
    check("(P1) 前の段と本段の差は recipe の許す差だけ (物理のパスを含まない)",
          not EP.stage_diff_problems([{"tag": s["tag"], "config": s["config"]} for s in st], EP.recipe())
          and not any(d.split(".")[0] in ("physProp", "turbulence", "mesh") for d in diffs), str(diffs))

    print("== 判別の A/B")
    n0 = launches()
    try:
        RA.run_staged_ns(ra, plan=ra)
        resA = "完走"
    except (RuntimeError, ValueError) as e:
        resA = f"{type(e).__name__}: {str(e)[:120]}"
    nA = launches() - n0
    check("(A) 変更なし → 計画の照合を通り forge (ダミー) を 1 回起動する", nA == 1 and "計画と違う" not in resA, f"起動 {nA} 回、{resA}")
    rb = copy_run("B")
    EP.confirm(rb, inputs_for(rb), attempt_id="B")
    orig = RA._ns_prestage

    def patched(cfg, cfl):
        return orig(cfg, 0.6 if cfl == 0.5 else cfl)        # S1_soft の CFL だけを計画の値 (0.5) から変える
    RA._ns_prestage = patched
    n0 = launches()
    try:
        RA.run_staged_ns(rb, plan=rb)
        resB = "完走"
    except (RuntimeError, ValueError) as e:
        resB = f"{type(e).__name__}: {str(e)[:200]}"
    finally:
        RA._ns_prestage = orig
    nB = launches() - n0
    check("(B) S1_soft の CFL だけを変える → 起動 0 回で止まる (計画と違う)", nB == 0 and "計画と違う" in resB, f"起動 {nB} 回、{resB}")
    if nA != 1:
        print("     試験の不成立: A も止まった")

    print("== 負例")
    bad = json.loads((rb / EP.PLAN_FILE).read_text())
    bad["confirmed"]["stages"][0]["sha256"] = "0" * 64
    (td / "plan_bad.json").write_text(json.dumps(bad, ensure_ascii=False))
    try:
        EP.load_confirmed(td / "plan_bad.json")
        ok = False
    except ValueError as e:
        ok = "合わない" in str(e)
    check("(N1) 確定した計画を書き換える → 読み込みで止まる", ok)
    try:
        EP.confirm(rb, inputs_for(rb))
        ok = False
    except ValueError as e:
        ok = "上書きしない" in str(e)
    check("(N2) 計画をもう一度確定する → 上書きしない", ok)
    rc3 = copy_run("N3")
    pre = (rc3 / RA.PRESTAGE_CONFIG).read_text()
    (rc3 / RA.PRESTAGE_CONFIG).write_text(pre.replace("dilatationCorrection: 2", "dilatationCorrection: 0", 1))
    try:
        EP.confirm(rc3, inputs_for(rc3))
        ok, why = False, "確定した"
    except ValueError as e:
        ok, why = "turbulence" in str(e), str(e)[:200]
    check("(N3) 前の段の乱流 (dilatationCorrection) が本段と違う → 確定しない (段の設定を作る関数か計画の検査のどちらかで)", ok, why)
    st3 = [{"tag": "S1_soft", "config": st[0]["config"].replace("dilatationCorrection: 2", "dilatationCorrection: 0", 1)},
           {"tag": "main", "config": st[-1]["config"]}]
    pr3 = EP.stage_diff_problems(st3, EP.recipe())
    check("(N3) 計画の検査そのもの: 前の段と本段の差に turbulence があれば問題にする", any("turbulence" in x for x in pr3), str(pr3)[:200])
    rc4 = copy_run("N4")
    EP.confirm(rc4, inputs_for(rc4))
    bc = (rc4 / "bcondConfig.yaml").read_text()
    (rc4 / "bcondConfig.yaml").write_text(bc.replace("Ts: 350.0", "Ts: 351.0", 1))
    n0 = launches()
    try:
        RA.run_staged_ns(rc4, plan=rc4)
        res4 = "完走"
    except (RuntimeError, ValueError) as e:
        res4 = str(e)[:200]
    check("(N4) 計画の後に境界条件 (壁温) を変える → 起動 0 回で止まる", launches() - n0 == 0 and "境界条件" in res4, res4)
    rc6 = copy_run("N6")
    EP.confirm(rc6, inputs_for(rc6))
    (rc6 / "settle_criteria.yaml").write_text((rc6 / "settle_criteria.yaml").read_text() + "\n# 書き足し\n")
    n0 = launches()
    try:
        RA.run_staged_ns(rc6, plan=rc6)
        res6 = "完走"
    except (RuntimeError, ValueError) as e:
        res6 = str(e)[:200]
    check("(N6) 計画の確定の後に停止規則の基準を変える → 起動 0 回で止まる (旧の台本 p5_ns.py の検査)", launches() - n0 == 0
          and "settle_criteria.yaml が計画の確定の後に変わった" in res6, res6)
    v = EP.verify_after(ra)
    check("(N5) S1 だけで止まった run の走った後の照合 → MISMATCH (段の記録が計画と違う)", v["status"] == "MISMATCH"
          and any("段の記録" in p for p in v["problems"]), str(v["problems"])[:200])
    rc5 = copy_run("N5")
    d5 = EP.confirm(rc5, inputs_for(rc5))
    (rc5 / "stage_manifest.json").write_text(json.dumps({"manifest_version": 1, "stages": [
        {"order": s["order"], "tag": s["tag"], "cfg_fnv": s["fnv"]} for s in d5["confirmed"]["stages"]]}))
    (rc5 / "RUN_PROVENANCE.txt").write_text(f"forge_sha256: {EP.sha256_file(dummy)}\n")
    v = EP.verify_after(rc5)
    check("(N5) 正例: 計画どおりの段の記録・起動の記録なら OK", v["status"] == "OK", str(v["problems"])[:200])
finally:
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
