#!/usr/bin/env python3
"""search の無い campaign (`forge_design.campaign.campaign`) の試験 (plan tooling-design-orchestrator §5.1 #4c)。正例 1・負例 7。

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_campaign_schema_tests.py
"""
import hashlib
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.campaign import campaign as CP  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def bad(doc, needle, base):
    try:
        CP.validate(doc, base)
    except ValueError as e:
        return needle in str(e), str(e)[:160]
    return False, "通った"


with tempfile.TemporaryDirectory() as tdn:
    td = Path(tdn)
    (td / "card_approval.json").write_text('{"x": 1}')
    sha = hashlib.sha256((td / "card_approval.json").read_bytes()).hexdigest()
    good = {"format": CP.FORMAT, "name": "m55_trial", "approval": {"file": "card_approval.json", "sha256": sha},
            "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger"}
    (td / "campaign.yaml").write_text(yaml.safe_dump(good))
    c = CP.load(td / "campaign.yaml")
    check("(P) 正例: 読めて、承認の記録・case・台帳のパスを解決する",
          c["approval"]["file"] == str((td / "card_approval.json").resolve()) and c["ledger"] == str((td / "_ledger").resolve()))
    for nm, edit, needle in (
            ("search", {"search": {"dv": ["L_c"]}}, "探索は未対応"),
            ("evaluations", {"evaluations": [{}]}, "複数の評価"),
            ("知らないキー", {"foo": 1}, "知らないキー"),
            ("知らない recipe", {"recipe": {"name": "axismach.other", "version": "v1"}}, "知らない recipe"),
            ("知らない版", {"recipe": {"name": "axismach.wt_card_chain", "version": "v9"}}, "知らない recipe"),
            ("承認の記録の sha256 の不一致", {"approval": {"file": "card_approval.json", "sha256": "0" * 64}}, "sha256 が campaign の値と違う"),
            ("必須のキーの欠け (ledger)", None, "キー ['ledger'] が無い")):
        doc = dict(good)
        if edit is None:
            doc.pop("ledger")
        else:
            doc.update(edit)
        ok, why = bad(doc, needle, td)
        check(f"(N) {nm} → 拒否", ok, why)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
