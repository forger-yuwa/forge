#!/usr/bin/env python3
"""search の無い campaign (`forge_design.campaign.campaign`) の試験 (plan tooling-design-orchestrator §5.1 #4c)。正例 1・負例 8。

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
    good = {"format": CP.FORMAT, "name": "m55_trial", "role": "fixed", "approval": {"file": "card_approval.json", "sha256": sha},
            "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger",
            "estimates": {"euler": {"gpu_hours": 0.1, "disk_gb": 1.0, "basis": "試験"}, "ns": {"gpu_hours": 1.0, "disk_gb": 5.0, "basis": "試験"}}}
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
            ("必須のキーの欠け (ledger)", None, "キー ['ledger'] が無い"),
            ("知らない役割", {"role": "admin"}, "role は"),
            # 段の見積もりの契約 (plan §4.5.1、2026-10-18): GPU の段 (euler・ns) の見積もりが無い・形が違う・負の値
            ("GPU の段の見積もりの欠け (estimates なし)", {"estimates": None}, "GPU の段 ['euler', 'ns']"),
            ("NS の見積もりの欠け", {"estimates": {"euler": {"gpu_hours": 0.1, "disk_gb": 1.0, "basis": "x"}}}, "GPU の段 ['ns']"),
            ("見積もりの根拠の欠け", {"estimates": {"euler": {"gpu_hours": 0.1, "disk_gb": 1.0}, "ns": {"gpu_hours": 1, "disk_gb": 1, "basis": "x"}}},
             "estimates.euler は"),
            ("負の見積もり", {"estimates": {"euler": {"gpu_hours": -1, "disk_gb": 1.0, "basis": "x"}, "ns": {"gpu_hours": 1, "disk_gb": 1, "basis": "x"}}},
             "0 以上の数"),
            ("知らない段の種類", {"estimates": {"euler": {"gpu_hours": 1, "disk_gb": 1, "basis": "x"}, "ns": {"gpu_hours": 1, "disk_gb": 1, "basis": "x"},
                                           "cfd": {"gpu_hours": 1, "disk_gb": 1, "basis": "x"}}}, "知らない段の種類")):
        doc = dict(good)
        if edit is None:
            doc.pop("ledger")
        else:
            doc.update(edit)
            if doc.get("estimates") is None:
                doc.pop("estimates", None)
        ok, why = bad(doc, needle, td)
        check(f"(N) {nm} → 拒否", ok, why)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
