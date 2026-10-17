#!/usr/bin/env python3
"""設計の実行役の入口 (`forge_design.campaign.cli`) の試験 (plan tooling-design-orchestrator §4.8、§5.1 #6c)。driver は起動しない。

  (C1) 台帳が無いときの status は「まだ無い」
  (C2) note は変更の要求のファイルを置く (台帳には書かない; 入力の sha256 と理由)
  (C3) run の起動の前の検査: この PC で Linger=no なら止める (Linger=yes なら飛ばす)

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
                                                        "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger"}))
    r = cli("status", case / "campaign.yaml")
    check("(C1) 台帳が無いときの status", r.returncode == 0 and "まだ無い" in r.stdout, r.stdout + r.stderr[-200:])
    (case / "problem_p1.yaml").write_text("x: 1\n")
    r = cli("note", case / "campaign.yaml", "--path", case / "problem_p1.yaml", "--reason", "CFL を下げた")
    reqs = list((case / "_ledger" / "requests").glob("*.json"))
    doc = json.loads(reqs[0].read_text()) if reqs else {}
    check("(C2) note は変更の要求のファイルを置く (sha256 と理由)", r.returncode == 0 and len(reqs) == 1 and doc.get("reason") == "CFL を下げた"
          and doc.get("input_sha256") == hashlib.sha256(b"x: 1\n").hexdigest() and not (case / "_ledger" / "ledger.json").exists(), str(doc))
    lg = subprocess.run(["loginctl", "show-user", os.environ.get("USER", ""), "-p", "Linger"], capture_output=True, text=True).stdout.strip()
    if lg == "Linger=no":
        r = cli("run", case / "campaign.yaml")
        check("(C3) Linger=no なら run は起動の前に止める", r.returncode == 2 and "Linger=no" in r.stderr and not (case / "_ledger" / "ledger.json").exists(),
              r.stderr[-200:])
    else:
        print(f"飛ばした (C3): {lg}")
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
