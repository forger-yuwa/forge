#!/usr/bin/env python3
"""段の見積もりの判別の A/B (plan tooling-design-orchestrator §6 V3 の 13a の事前登録、codex 諮問 v3、2026-10-18)。GPU・forge なし。

変えるのは Euler の段の見積もりの gpu_hours だけ: A = 無し (今の recipe と同じ)、B = 4/60。同じ「格子の段 (mesh_p0) が済んだ」台帳・承認の予算
GPU 0.01 時間。本物の driver・recipe (axismach.wt_card_chain/v1) で起動の判定まで各 1 回回し、起動 (`Driver._start_unit`) は到達の記録だけの
代わり (記録して止める)。見る量: Euler の試行の追加数・起動への到達数・止まった理由・先の試行の保持。

期待 (事前登録): A は追加 1・到達 1、B は追加 0・到達 0・予算の理由で止まる。A も止まれば対照の不成立。B の値は試験の入力で、生産の見積もりの
妥当性は示さない。

実行: /home/sano/work/forge/design/.venv-opt/bin/python notes/investigations/2026-10-18-estimate-ab/estimate_ab.py
結果は同じディレクトリの result.json。
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "design"))
from forge_design.campaign import driver as DR  # noqa: E402
from forge_design.campaign.ledger import Ledger  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

C71 = REPO / "case" / "71.m55_combustion_wt"


class Reached(Exception):
    pass


def make_case(td: Path, tag: str, estimates) -> Path:
    case = td / tag
    case.mkdir()
    d = yaml.safe_load((C71 / "card_m55_3coef.yaml").read_text())
    d.update(format=W.CARD_FORMAT_V2, budget={"gpu_hours": 0.01, "disk_gb": 50.0}, mesh_limits={"max_nodes_ns": 443000, "max_nodes_euler": 467000},
             sources={**{k: {"kind": "user", "ref": "A/B"} for k in ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall",
                                                                     "tolerances", "p_ambient", "budget")}, "mesh_limits": {"kind": "user", "ref": "A/B"}})
    card = case / "card_p0.yaml"
    card.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    prob = case / "problem_p0_card.yaml"
    shutil.copy(C71 / "problem_m55_p5_settle_v2.yaml", prob)
    clean = td / f"clean_{tag}"
    clean.mkdir()
    for c in (["init", "-q"], ["-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "--allow-empty", "-m", "t"]):
        subprocess.run(["git", "-C", str(clean), *c], check=True, capture_output=True)
    pv = AP.preview(card, prob, repo_dir=clean, judge_env=False)
    assert not pv["problems"], pv["problems"]
    (case / "pv.json").write_text(json.dumps(pv, ensure_ascii=False))
    AP.write(case / "pv.json", "A/B の仮の承認", case / "card_approval.json")
    camp = {"format": "forge.campaign/v1", "name": f"ab_{tag}", "role": "fixer",
            "approval": {"file": "card_approval.json", "sha256": AP.sha256_file(case / "card_approval.json")},
            "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger"}
    if estimates is not None:
        camp["estimates"] = estimates
    (case / "campaign.yaml").write_text(yaml.safe_dump(camp, sort_keys=False, allow_unicode=True))
    # 格子の段 (mesh_p0) が済んだ台帳: 本物の Ledger で作り、試行 1 つを取り込み済みにして結果のファイルを置く
    from forge_design.campaign import campaign as CP
    c = CP.load(case / "campaign.yaml")
    with Ledger(c["ledger"]) as L:
        L.init(c, {"gpu_hours": 0.01, "disk_gb": 50.0})
        L.doc["attempts"].append({"id": "0001-mesh_p0", "step": "mesh_p0", "state": "ingested", "run_dir": None,
                                  "interpretation": {"ok": True, "branch": None, "detail": "A/B の前提 (格子の段は済んだ)"}})
        L.doc["ingested"].append("0001-mesh_p0")
        L.save()
    ad = Path(c["ledger"]) / "attempts" / "0001-mesh_p0"
    ad.mkdir(parents=True)
    (ad / "result.json").write_text(json.dumps({"step": "mesh_p0", "ok": True, "branch": None, "detail": "A/B の前提", "outputs": {}}))
    return case


def run_arm(case: Path) -> dict:
    from forge_design.campaign import campaign as CP
    c = CP.load(case / "campaign.yaml")
    reached = []
    orig = DR.Driver._start_unit

    def stub(self, a):
        reached.append(a["id"])
        raise Reached()
    DR.Driver._start_unit = stub
    try:
        d = DR.Driver(case / "campaign.yaml", DR.recipe_for(c, REPO), systemd=False, design_dir=REPO / "design")
        try:
            out = d.run()
        except Reached:
            out = {"status": "reached (起動の代わりで止めた)"}
    finally:
        DR.Driver._start_unit = orig
    led = json.loads((Path(c["ledger"]) / "ledger.json").read_text())
    eu = [a["id"] for a in led["attempts"] if a["step"].startswith("euler")]
    return {"out": {k: v for k, v in out.items() if k != "ledger"}, "euler_attempts_added": len(eu), "start_reached": len(reached),
            "attempts": [a["id"] for a in led["attempts"]], "stops": led.get("stops") or [], "estimates_in_ledger": led.get("estimates")}


td = Path(tempfile.mkdtemp(prefix="estimate_ab_"))
try:
    A = run_arm(make_case(td, "A", None))
    B = run_arm(make_case(td, "B", {"euler": {"gpu_hours": 4 / 60, "disk_gb": 0.5, "basis": "A/B の試験の入力 (生産の見積もりではない)"}}))
    ok_A = A["euler_attempts_added"] == 1 and A["start_reached"] == 1
    ok_B = (B["euler_attempts_added"] == 0 and B["start_reached"] == 0 and B["out"].get("status") == "stopped"
            and "予算" in str(B["out"].get("why")) and "0001-mesh_p0" in B["attempts"])
    verdict = ("対照の不成立 (A が到達しない)" if not ok_A else ("期待どおり (見積もりの未接続の仮説を支持)" if ok_B else "B も到達・または別の理由で止まった"))
    res = {"A": A, "B": B, "verdict": verdict, "ok_A": ok_A, "ok_B": ok_B}
    (HERE / "result.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str) + "\n")
    print(json.dumps({"verdict": verdict, "A": {k: A[k] for k in ("out", "euler_attempts_added", "start_reached")},
                      "B": {k: B[k] for k in ("out", "euler_attempts_added", "start_reached")}}, ensure_ascii=False, indent=1))
finally:
    shutil.rmtree(td, ignore_errors=True)
