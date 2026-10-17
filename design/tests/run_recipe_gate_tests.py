#!/usr/bin/env python3
"""recipe の判定と成果物の段を、判定専用の環境の本物の道具で回す試験 (plan tooling-design-orchestrator §4.5.1、§5.1 #6c の判別の A/B)。

codex (diagnose 2026-10-17 executor) の判別の A/B: 同じ run の写しで、成果物が無い (B) と有る (A) とで判定と recipe の遷移を見る。
  (B) 成果物が無い run の qualify_p1 → incomplete・make_deliverables → recipe は deliverables_p1 へ進む (止まらない)
  (A) deliverables_p1 (判定専用の環境で成果物 → 最後の判定) → accepted → recipe は終わる
  A も失敗すれば試験の不成立。

case/71 の run_0014_chain_v2 (手元) の写し (h5 は symlink) と、case/71 のカードを版 2 にした承認の記録 (このリポジトリの HEAD を承認の commit) を使う。
**HEAD を承認の commit に使う**ので、追跡ファイルに変更の無い状態で回す。

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_recipe_gate_tests.py
所要: 約 5 分 (判定 2 回・成果物 1 回)。run・道具が無ければ飛ばす。
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
from forge_design.campaign.recipes import wt_card_chain as WC  # noqa: E402
from forge_design.feedback import gate_env as GE  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

C71 = REPO / "case" / "71.m55_combustion_wt"
R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", C71 / "run_0014_chain_v2"))
CARD, PROB = C71 / "card_m55_3coef.yaml", C71 / "problem_m55_p5_settle_v2.yaml"
CONV = os.environ.get("FORGE_CONVERTER") or str(REPO / "solver_density_cuda" / ".build-native" / "card4a" / "convertGmshToForge")
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True)


if any(not p.exists() for p in (R14 / "nozzle.h5", CARD, PROB, Path(CONV))):
    print("飛ばした: run_0014・カード・変換器が無い")
    sys.exit(0)
if git("status", "--porcelain", "--untracked-files=no").stdout.strip():
    print("止めた: 追跡ファイルに変更がある (HEAD を承認の commit に使う)")
    sys.exit(2)

td = Path(tempfile.mkdtemp(prefix="recipe_gate_tests_"))
os.environ[GE.GATE_ROOT_ENV] = str(td / "gates")
os.environ["FORGE_CONVERTER"] = CONV
head = git("rev-parse", "HEAD").stdout.strip()
try:
    case = td / "case"
    led = case / "_ledger"
    (led / "attempts").mkdir(parents=True)
    d = yaml.safe_load(CARD.read_text())
    src = {k: {"kind": "user", "ref": "case/71"} for k in ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall",
                                                         "tolerances", "p_ambient", "budget")}
    d.update(format=W.CARD_FORMAT_V2, budget={"gpu_hours": 8.0, "disk_gb": 50.0}, sources=src)
    card = case / "card_p1.yaml"
    card.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    prob = case / "problem_p1.yaml"
    shutil.copy(PROB, prob)
    pv = AP.preview(card, prob)
    check("前提: 承認の preview に問題が無い", not pv["problems"], str(pv["problems"])[:300])
    (td / "pv.json").write_text(json.dumps(pv, ensure_ascii=False))
    apath = case / "card_approval.json"
    AP.write(td / "pv.json", "試験の仮の承認", apath)
    run = case / "run_0014_ns_p1"
    run.mkdir()
    for f in R14.iterdir():
        if f.is_file():
            if f.suffix == ".h5":
                (run / f.name).symlink_to(f)
            else:
                shutil.copy2(f, run / f.name)
    camp = {"case_dir": str(case), "ledger": str(led), "approval": {"file": str(apath), "sha256": AP.sha256_file(apath)}, "name": "t"}
    base = {"k": 1, "case_dir": str(case), "approval": str(apath), "card": str(card), "problem": str(prob), "problem_card": str(prob),
            "problem_euler": str(prob), "repo": str(REPO), "iteration": 1, "recalibrations_done": 0, "ns_run": str(run), "prev_run": None}

    def doc_with(steps):
        doc = {"attempts": [], "ingested": [], "requests": [], "run_env": {}}
        for i, (st, r) in enumerate(steps, 1):
            aid = f"{i:04d}-{st}"
            (led / "attempts" / aid).mkdir(parents=True, exist_ok=True)
            (led / "attempts" / aid / "result.json").write_text(json.dumps({"step": st, **r}))
            doc["attempts"].append({"id": aid, "step": st, "run_dir": r.get("run_dir")})
            doc["ingested"].append(aid)
        return doc

    # 前の段 (p0・p1 の格子・Euler・NS) は合格の結果のファイルだけを置く (この試験は判定と成果物の段だけを本物で回す)
    eu = {"ok": True, "euler_ref": {"run": str(case / "run_e"), "res": str(case / "run_e" / "res_1.h5"), "sha256": None}}
    ns = {"ok": True, "run_dir": str(run), "final_field": {"path": str(run / "res_46171.h5"), "sha256": None}}
    prior = [("mesh_p0", {"ok": True}), ("euler_p0", eu), ("ns_p0", ns), ("calibrate_p0", {"ok": True}), ("problem_p1", {"ok": True}),
             ("mesh_p1", {"ok": True}), ("euler_p1", eu), ("ns_p1", ns)]
    R = WC.Recipe(repo_dir=REPO)
    qa = WC.exec_step({**base, "step": "qualify_p1", "out": str(td / "q.json")})
    nxt = R.next_step(doc_with(prior + [("qualify_p1", qa)]), camp)
    check("(B) 成果物が無い run の qualify_p1 → incomplete・make_deliverables、recipe は deliverables_p1 へ進む (止まらない)",
          qa["ok"] and qa["branch"] == "make_deliverables" and nxt is not None and nxt.get("step") == "deliverables_p1", f"{qa['detail']} → {nxt and nxt.get('step')}")
    dl = WC.exec_step({**base, "step": "deliverables_p1", "out": str(td / "d.json")})
    shutil.rmtree(led / "attempts")
    (led / "attempts").mkdir()
    nxt2 = R.next_step(doc_with(prior + [("qualify_p1", qa), ("deliverables_p1", dl)]), camp)
    check("(A) deliverables_p1 (判定専用の環境で成果物 → 最後の判定) → accepted、recipe は終わる",
          dl["ok"] and dl["branch"] == "accepted" and nxt2 is None, f"{dl['detail']} → {nxt2}")
    man = json.loads((run / "deliverables" / "MANIFEST.json").read_text()) if (run / "deliverables" / "MANIFEST.json").is_file() else {}
    check("(A) 成果物の MANIFEST の gate は承認の commit", (man.get("gate") or {}).get("commit") == head, str(man.get("gate"))[:200])
finally:
    co = td / "gates" / head
    if co.exists():
        git("worktree", "remove", "--force", str(co))
    git("worktree", "prune")
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
