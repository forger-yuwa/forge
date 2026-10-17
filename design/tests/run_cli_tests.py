#!/usr/bin/env python3
"""設計の実行役の入口 (`forge_design.campaign.cli`) の試験 (plan tooling-design-orchestrator §4.8、§5.1 #6c)。driver は起動しない。

  (C1) 台帳が無いときの status は「まだ無い」
  (C2) note は変更の要求のファイルを置く (台帳には書かない; 入力の sha256 と理由)
  (C3) run の起動の前の検査: この PC で Linger=no なら止める (Linger=yes なら飛ばす)
  (C4) --role が campaign の role と違えば起動の前に止める (§4.0; 役割は起動のときに人が決める)
  (C5) 台帳があるのに --resume なし・台帳が無いのに --resume → 止める (付け忘れで別の campaign を上書きしない)
  (C6) status の集計: 試行の数 (状態ごと)・人に返した回数 (台帳の stops)・直した回数 (変更の要求 note)・予算

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_cli_tests.py
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

DESIGN = Path(__file__).resolve().parents[1]
PY = sys.executable
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def cli(*args):
    return subprocess.run([PY, "-m", "forge_design.campaign.cli", *map(str, args)], cwd=str(DESIGN), capture_output=True, text=True,
                          env=dict(os.environ, PYTHONPATH=str(DESIGN)))


with tempfile.TemporaryDirectory() as tdn:
    case = Path(tdn)
    (case / "card_approval.json").write_text(json.dumps({"protected": {"budget": {"gpu_hours": 1, "disk_gb": 1}}}))
    sha = hashlib.sha256((case / "card_approval.json").read_bytes()).hexdigest()
    (case / "campaign.yaml").write_text(yaml.safe_dump({"format": "forge.campaign/v1", "name": "t", "role": "fixer",
                                                        "approval": {"file": "card_approval.json", "sha256": sha},
                                                        "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger",
                                                        "estimates": {"euler": {"gpu_hours": 0.1, "disk_gb": 1.0, "basis": "試験"}, "ns": {"gpu_hours": 1.0, "disk_gb": 5.0, "basis": "試験"}}}))
    r = cli("status", case / "campaign.yaml")
    check("(C1) 台帳が無いときの status", r.returncode == 0 and "まだ無い" in r.stdout, r.stdout + r.stderr[-200:])
    (case / "problem_p1.yaml").write_text("x: 1\n")
    r = cli("note", case / "campaign.yaml", "--path", case / "problem_p1.yaml", "--reason", "CFL を下げた")
    reqs = list((case / "_ledger" / "requests").glob("*.json"))
    doc = json.loads(reqs[0].read_text()) if reqs else {}
    check("(C2) note は変更の要求のファイルを置く (sha256 と理由)", r.returncode == 0 and len(reqs) == 1 and doc.get("reason") == "CFL を下げた"
          and doc.get("input_sha256") == hashlib.sha256(b"x: 1\n").hexdigest() and not (case / "_ledger" / "ledger.json").exists(), str(doc))
    lg = subprocess.run(["loginctl", "show-user", os.environ.get("USER", ""), "-p", "Linger"], capture_output=True, text=True).stdout.strip()
    r = cli("run", case / "campaign.yaml", "--role", "fixed")
    check("(C4) --role fixed が campaign の role fixer と違う → 起動の前に止める (終了 2)",
          r.returncode == 2 and "campaign の role fixer と違う" in r.stderr and not (case / "_ledger" / "ledger.json").exists(), r.stderr[-200:])
    r = cli("run", case / "campaign.yaml", "--role", "fixer", "--resume")
    check("(C5) 台帳が無いのに --resume → 止める (終了 2)", r.returncode == 2 and "台帳" in r.stderr and "無い" in r.stderr, r.stderr[-200:])
    led = {"format": "forge.ledger/v1", "campaign_id": "x", "name": "t", "role": "fixer", "recipe": {"name": "axismach.wt_card_chain", "version": "v1"},
           "attempts": [{"id": "0001-mesh_p0", "state": "ingested"}, {"id": "0002-euler_p0", "state": "ingested"},
                        {"id": "0003-ns_p0", "state": "launched"}],
           "requests": [{"kind": "note", "reason": "CFL を下げた"}, {"kind": "note", "reason": "格子を細かく"}],
           "stops": [{"time": "t", "status": "stopped", "why": "関門で止まった"}],
           "budget": {"approved": {"gpu_hours": 1, "disk_gb": 1}, "consumed": {"gpu_hours": 0.25, "disk_gb": 0.5}}}
    (case / "_ledger").mkdir(exist_ok=True)
    (case / "_ledger" / "ledger.json").write_text(json.dumps(led))
    r = cli("run", case / "campaign.yaml", "--role", "fixer")
    check("(C5) 台帳があるのに --resume なし → 止める (終了 2)", r.returncode == 2 and "--resume" in r.stderr, r.stderr[-200:])
    r = cli("status", case / "campaign.yaml", "--json")
    sm = json.loads(r.stdout.strip().splitlines()[-1]) if r.returncode == 0 else {}
    check("(C6) status の集計: 試行 3 (ingested 2・launched 1)・人に返した 1・直した 2・消費した予算",
          sm.get("attempts") == 3 and sm.get("attempts_by_state") == {"ingested": 2, "launched": 1} and sm.get("returned_to_human") == 1
          and sm.get("fixes") == 2 and sm.get("budget_consumed") == {"gpu_hours": 0.25, "disk_gb": 0.5}, r.stdout[-300:] + r.stderr[-200:])
    (case / "_ledger" / "ledger.json").unlink()
    if lg == "Linger=no":
        r = cli("run", case / "campaign.yaml", "--role", "fixer")
        check("(C3) Linger=no なら run は起動の前に止める", r.returncode == 2 and "Linger=no" in r.stderr and not (case / "_ledger" / "ledger.json").exists(),
              r.stderr[-200:])
    else:
        print(f"飛ばした (C3): {lg}")
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
