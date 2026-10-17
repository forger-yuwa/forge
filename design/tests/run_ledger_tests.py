#!/usr/bin/env python3
"""台帳・予算の勘定・停止の記録 (`forge_design.campaign.ledger`・`budget`・`stop`) の試験 (plan tooling-design-orchestrator §4.5.1、§5.1 #6a)。GPU なし。

  (L1) 二つ目の書き手 (別のプロセス) はロックを取れない (LedgerBusy)、一つ目が離すと取れる
  (L2) 取り込みを 2 回呼んでも精算は 1 回 (GPU 時間・予約の解放)
  (L3) 試行は起動の前に記録され (state recorded・unit 名・予約)、同じ名前の別の campaign (別の台帳) は campaign の ID が違う
  (L4) 台帳を作った後に役割・承認・recipe を変えた campaign は拒否
  (L5) 変更の要求は CLI が置き、driver が取り込む (台帳の書き手は 1 つ)
  (B1) 終わりの時刻の無い起動は 0 にしない (根拠の時刻まで)、根拠も無ければ例外
  (B2) ディスクは重ならない所有の範囲 (入れ子は外側だけ、symlink の先は数えない)
  (B3) 見積もりが残り (承認 − 消費 − 予約) を超えれば始めない、内なら始める
  (S1) 停止の記録の欄がそろい、上書きしない。driver の補いは原因 unknown (根拠があればその分類)

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_ledger_tests.py
"""
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.campaign import budget as BG  # noqa: E402
from forge_design.campaign import ledger as LG  # noqa: E402
from forge_design.campaign import stop as ST  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


CAMP = {"name": "m55_trial", "role": "fixed", "approval": {"file": "/x/card_approval.json", "sha256": "a" * 64},
        "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": "/x"}
BUD = {"gpu_hours": 4.0, "disk_gb": 50.0}
TRY_LOCK = r'''
import sys
sys.path.insert(0, sys.argv[2])
from forge_design.campaign import ledger as LG
try:
    LG.Ledger(sys.argv[1]).acquire()
    print("ACQUIRED")
except LG.LedgerBusy:
    print("BUSY")
'''

with tempfile.TemporaryDirectory() as tdn:
    td = Path(tdn)
    led = td / "_ledger"
    L = LG.Ledger(led).acquire()
    doc = L.init(CAMP, BUD)
    other = subprocess.run([sys.executable, "-c", TRY_LOCK, str(led), str(Path(__file__).resolve().parents[1])], capture_output=True, text=True)
    check("(L1) 二つ目の書き手はロックを取れない", other.stdout.strip() == "BUSY", other.stdout + other.stderr[-200:])
    a = L.new_attempt("euler_p0", td / "run_0001_euler_p0", {"gpu_hours": 4.0})
    check("(L3) 試行は起動の前に記録 (recorded・unit 名・予約)", a["state"] == "recorded" and a["unit"].startswith(f"forge-step-{doc['campaign_id'][:8]}")
          and LG.read_ledger(led)["budget"]["reserved"]["gpu_hours"] == 4.0, str(a))
    st = {"launches": [{"start": 1000.0, "end": 1000.0 + 1800.0, "seconds": 1800.0, "end_basis": "終わりの時刻"}], "disk_gb": 1.5, "result": "ok"}
    r1 = L.ingest(a["id"], st, owned=[td / "run_0001_euler_p0"])
    r2 = L.ingest(a["id"], st, owned=[td / "run_0001_euler_p0"])
    b = LG.read_ledger(led)["budget"]
    check("(L2) 取り込みを 2 回呼んでも精算は 1 回", r1 is True and r2 is False and abs(b["consumed"]["gpu_hours"] - 0.5) < 1e-12
          and b["reserved"]["gpu_hours"] == 0.0 and len(b["launches"]) == 1, str(b))
    LG.put_request(led, "note", reason="CFL を 1.0 から 0.8 に下げた", input_sha256="b" * 64)
    taken = L.take_requests()
    check("(L5) 変更の要求は CLI が置き、driver が取り込む", len(taken) == 1 and LG.read_ledger(led)["requests"][0]["kind"] == "note"
          and not list((led / "requests").glob("*.json")))
    L.release()
    other = subprocess.run([sys.executable, "-c", TRY_LOCK, str(led), str(Path(__file__).resolve().parents[1])], capture_output=True, text=True)
    check("(L1) 一つ目が離すと取れる", other.stdout.strip() == "ACQUIRED", other.stdout)
    led2 = td / "other_case" / "_ledger"
    with LG.Ledger(led2) as L2:
        d2 = L2.init(dict(CAMP), BUD)
    check("(L3) 同じ名前の別の campaign は campaign の ID が違う", d2["campaign_id"] != doc["campaign_id"] and d2["name"] == doc["name"])
    with LG.Ledger(led) as L3:
        try:
            L3.init({**CAMP, "role": "fixer"}, BUD)
            ok = False
        except ValueError as e:
            ok = "role" in str(e)
    check("(L4) 台帳を作った後に役割を変えた campaign は拒否", ok)

    print("== 予算")
    r = BG.launch_seconds({"start": 1000.0, "end": None}, {"time": 1000.0 + 600.0, "basis": "unit が止まったと分かった時刻"})
    try:
        BG.launch_seconds({"start": 1000.0, "end": None}, None)
        ok = False
    except ValueError:
        ok = True
    check("(B1) 終わりの時刻の無い起動は根拠の時刻まで (600 s)、根拠も無ければ例外 (0 にしない)", r["seconds"] == 600.0 and "根拠" in r["end_basis"] and ok)
    own = td / "own"
    (own / "run_a" / "sub").mkdir(parents=True)
    (own / "run_a" / "f1").write_bytes(b"x" * 1000)
    (own / "run_a" / "sub" / "f2").write_bytes(b"y" * 500)
    big = td / "big_outside"
    big.write_bytes(b"z" * 100000)
    (own / "run_a" / "link").symlink_to(big)
    n = BG.owned_disk_bytes([own / "run_a", own / "run_a" / "sub", own / "run_a" / "f1"])
    check("(B2) 入れ子は外側だけ、symlink の先は数えない", 1500 <= n < 1500 + 200, f"{n} バイト")
    bud = {"approved": {"gpu_hours": 4.0, "disk_gb": 50.0}, "consumed": {"gpu_hours": 3.0, "disk_gb": 10.0}, "reserved": {"gpu_hours": 0.5, "disk_gb": 0.0}}
    c1, c2 = BG.check_start(bud, {"gpu_hours": 0.4, "disk_gb": 5.0}), BG.check_start(bud, {"gpu_hours": 0.6, "disk_gb": 5.0})
    check("(B3) 見積もりが残りの内なら始める、超えれば始めない", c1["ok"] and not c2["ok"] and "始めない" in c2["why"], str(c2["remaining"]))
    check("(B3) 試行の監視: 消費 + この試行 > 承認 を返す", BG.over_budget(bud, {"gpu_hours": 1.1}) == ["gpu_hours"] and BG.over_budget(bud, {"gpu_hours": 0.9}) == [])

    print("== 停止の記録")
    ad = td / "attempts" / "0002-ns_p0"
    ad.mkdir(parents=True)
    (ad / "evid.txt").write_text("VERDICT: NOT SETTLED")
    rec = ST.write_stop(ad, attempt="0002-ns_p0", step="ns_p0", verdicts=[{"tool": "check_settled", "verdict": "NOT SETTLED"}],
                        evidence_paths=[ad / "evid.txt"], budget={"gpu_hours": 1.2}, process={"unit": "u", "state": "stopped"},
                        cause="budget", detail="GPU 時間が承認の値に届いた")
    try:
        ST.write_stop(ad, attempt="x", step="x", verdicts=[], evidence_paths=[], budget={}, process={}, cause="budget", detail="x")
        ok = False
    except FileExistsError:
        ok = True
    check("(S1) 停止の記録の欄がそろい、上書きしない", not ST.missing_fields(rec) and rec["evidence"][0]["sha256"] and ok)
    ad2 = td / "attempts" / "0003-ns_p0"
    ad2.mkdir()
    sup = ST.supplement(ad2, attempt="0003-ns_p0", step="ns_p0", process={"unit": "u2", "state": "inactive"}, budget={})
    ad3 = td / "attempts" / "0004-ns_p0"
    ad3.mkdir()
    sup3 = ST.supplement(ad3, attempt="0004-ns_p0", step="ns_p0", process={}, budget={}, known_cause="resource", basis="forge_run.log に Write failed")
    check("(S1) driver の補い: 根拠が無ければ unknown、根拠があればその分類、既にあれば補わない",
          sup["cause"]["class"] == "unknown" and sup["recorded_by"] == "driver" and sup3["cause"]["class"] == "resource"
          and ST.supplement(ad, attempt="x", step="x", process={}, budget={}) is None)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
