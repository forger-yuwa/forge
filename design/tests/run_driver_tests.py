#!/usr/bin/env python3
"""driver と試行の service の障害の試験 (plan tooling-design-orchestrator §4.5.1、§5.1 #6b)。GPU・forge なし: ダミーの段 (bash の sleep) を
systemd の user service の試行で回す (試験の recipe: design/tests/fixtures/dummy_recipe.py)。

  (T1)  ふつう: 2 段を回して終わる、各段は 1 回だけ走る、GPU の段の時間を数える
  (T2)  台帳の記録の後・起動の前に driver が落ちる → 再開で起動する (段は 1 回)
  (T3)  起動の後・起動の印の前に driver が落ちる → 再開で動いている試行を見つけて待つ (段は 1 回、二つ目を起動しない)
  (T4)  driver のいない間に試行が終わる → 再開で 1 回だけ取り込む
  (T5)  取り込みの前に落ちる → 再開で 1 回取り込む。取り込みの後に落ちる → 再開で 2 回目を取り込まない (GPU 時間を二重に数えない)
  (T6)  二つ目の driver (二重の --resume) → busy を返し、一つ目を止めない
  (T7)  予算の超過 → 試行が先に停止の記録 (budget) を書いてから止まる、子のプロセスが残らない
  (T8)  driver (CLI) が死んでも試行の予算の監視は続く (driver なしで停止の記録 budget)
  (T9)  試行の実行担当を SIGKILL → 子のプロセスが残らない、driver が停止の記録を補う (原因 unknown)、GPU 時間を 0 にしない
  (T10) 停止の途中の再開 → 待ってから取り込む
  (T11) 空きの容量が下限より小さい → 起動の前に止まる (試行を記録しない)
  (T12) 出力先の run が既にある → 起動の前に止まる

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_driver_tests.py
所要: 約 2〜3 分。systemd の user service が要る (無ければ飛ばす)。
"""
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
DESIGN = HERE.parent
sys.path.insert(0, str(DESIGN))
from forge_design.campaign import driver as DR  # noqa: E402
from forge_design.campaign import ledger as LG  # noqa: E402
from forge_design.campaign import step as SP  # noqa: E402

PY = sys.executable
FAIL = 0
UNITS = []


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


if shutil.which("systemd-run") is None or subprocess.run(["systemctl", "--user", "is-system-running"], capture_output=True).returncode not in (0, 1):
    print("飛ばした: systemd の user service が使えない")
    sys.exit(0)

TD = Path(tempfile.mkdtemp(prefix="driver_tests_"))


def make_case(name, steps, gpu_hours=1.0, disk_gb=10.0):
    case = TD / name
    case.mkdir()
    (case / "card_approval.json").write_text(json.dumps({"protected": {"budget": {"gpu_hours": gpu_hours, "disk_gb": disk_gb}}}))
    sha = hashlib.sha256((case / "card_approval.json").read_bytes()).hexdigest()
    (case / "campaign.yaml").write_text(yaml.safe_dump({"format": "forge.campaign/v1", "name": name, "role": "fixed",
                                                        "approval": {"file": "card_approval.json", "sha256": sha},
                                                        "recipe": {"name": "axismach.wt_card_chain", "version": "v1"},
                                                        "case_dir": ".", "ledger": "_ledger"}))
    (case / "steps.json").write_text(json.dumps(steps))
    return case


def step(name, body, gpu=False):
    return {"step": name, "run_dir": f"run_{name}", "gpu": gpu,
            "command": ["bash", "-c", f"mkdir -p {{RUN}} && echo x >> {{RUN}}/../count_{name}.txt && {body} && echo ok > {{RUN}}/out"]}


def drv(case, crash=None, background=False, extra=(), timeout=180):
    env = dict(os.environ, FORGE_DUMMY_STEPS=str(case / "steps.json"), PYTHONPATH=str(DESIGN))
    env.pop(DR.CRASH_ENV, None)
    if crash:
        env[DR.CRASH_ENV] = crash
    cmd = [PY, "-m", "forge_design.campaign.driver", str(case / "campaign.yaml"), "--recipe-file", str(HERE / "fixtures" / "dummy_recipe.py"),
           "--poll", "0.2", "--step-poll", "0.2", *extra]
    if background:
        return subprocess.Popen(cmd, env=env, cwd=str(DESIGN), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return subprocess.run(cmd, env=env, cwd=str(DESIGN), capture_output=True, text=True, timeout=timeout)


def out_status(r):
    for ln in reversed((r.stdout or "").splitlines()):
        if ln.startswith("{"):
            return json.loads(ln)
    return {}


def ledger(case):
    d = LG.read_ledger(case / "_ledger")
    if d:
        for a in d["attempts"]:
            if a["unit"] not in UNITS:
                UNITS.append(a["unit"])
    return d


def count(case, name):
    f = case / f"count_{name}.txt"
    return len(f.read_text().splitlines()) if f.exists() else 0


def wait_for(fn, timeout=60, poll=0.2):
    t_end = time.time() + timeout
    while time.time() < t_end:
        v = fn()
        if v:
            return v
        time.sleep(poll)
    return None


def alive(pid):
    try:
        st = Path(f"/proc/{int(pid)}/stat").read_text()
        return st[st.rindex(")") + 2:].split()[0] not in ("Z", "X")
    except (OSError, ValueError, TypeError):
        return False


try:
    c = make_case("t1", [step("s1", "sleep 1", gpu=True), step("s2", "sleep 0.5")])
    r = drv(c)
    d = ledger(c)
    check("(T1) ふつう: 終わる・各段 1 回・GPU の段の時間を数える",
          out_status(r).get("status") == "done" and len(d["ingested"]) == 2 and count(c, "s1") == 1 and count(c, "s2") == 1
          and 0.0 < d["budget"]["consumed"]["gpu_hours"] < 0.01, f"{out_status(r)} {r.stdout[-200:]}")

    c = make_case("t2", [step("s1", "sleep 1")])
    r1 = drv(c, crash="after_record")
    d1 = ledger(c)
    r2 = drv(c)
    d2 = ledger(c)
    check("(T2) 記録の後・起動の前に落ちる → 再開で起動する (段は 1 回)",
          r1.returncode == 99 and d1["attempts"][0]["state"] == "recorded" and out_status(r2).get("status") == "done" and count(c, "s1") == 1
          and d2["attempts"][0].get("launched_by") == "reconcile", f"{r1.returncode} {out_status(r2)}")

    c = make_case("t3", [step("s1", "sleep 3")])
    r1 = drv(c, crash="after_start")
    r2 = drv(c)
    d2 = ledger(c)
    check("(T3) 起動の後・印の前に落ちる → 再開で待つ (段は 1 回、二つ目を起動しない)",
          r1.returncode == 99 and out_status(r2).get("status") == "done" and count(c, "s1") == 1 and len(d2["attempts"]) == 1,
          f"{r1.returncode} {out_status(r2)} count {count(c, 's1')}")

    c = make_case("t4", [step("s1", "sleep 3")])
    p = drv(c, background=True)
    wait_for(lambda: (ledger(c) or {}).get("attempts") and ledger(c)["attempts"][0]["state"] == "launched")
    p.kill()
    p.wait()
    ad = c / "_ledger" / "attempts" / ledger(c)["attempts"][0]["id"]
    fin = wait_for(lambda: (SP.read_status(ad) or {}).get("state") in SP.FINAL_STATES, timeout=30)
    r2 = drv(c)
    d2 = ledger(c)
    check("(T4) driver のいない間に試行が終わる → 再開で 1 回だけ取り込む",
          fin and out_status(r2).get("status") == "done" and d2["ingested"] == [d2["attempts"][0]["id"]] and count(c, "s1") == 1,
          str(out_status(r2)))

    c = make_case("t5a", [step("s1", "sleep 1", gpu=True)])
    r1 = drv(c, crash="before_ingest")
    d1 = ledger(c)
    r2 = drv(c)
    d2 = ledger(c)
    g = d2["budget"]["consumed"]["gpu_hours"]
    check("(T5) 取り込みの前に落ちる → 再開で 1 回取り込む",
          r1.returncode == 99 and not d1["ingested"] and len(d2["ingested"]) == 1 and len(d2["budget"]["launches"]) == 1 and g > 0,
          f"{r1.returncode} {d2['budget']}")
    c = make_case("t5b", [step("s1", "sleep 1", gpu=True)])
    r1 = drv(c, crash="after_ingest")
    g1 = ledger(c)["budget"]["consumed"]["gpu_hours"]
    r2 = drv(c)
    d2 = ledger(c)
    check("(T5) 取り込みの後に落ちる → 再開で 2 回目を取り込まない (GPU 時間を二重に数えない)",
          r1.returncode == 99 and out_status(r2).get("status") == "done" and d2["budget"]["consumed"]["gpu_hours"] == g1
          and len(d2["budget"]["launches"]) == 1, f"{g1} → {d2['budget']['consumed']['gpu_hours']}")

    c = make_case("t6", [step("s1", "sleep 5")])
    p1 = drv(c, background=True)
    wait_for(lambda: (ledger(c) or {}).get("attempts") and ledger(c)["attempts"][0]["state"] == "launched")
    r2 = drv(c)
    out1, _ = p1.communicate(timeout=60)
    s1 = next((json.loads(ln) for ln in reversed(out1.splitlines()) if ln.startswith("{")), {})
    check("(T6) 二つ目の driver → busy を返し、一つ目を止めない (一つ目は終わる)",
          out_status(r2).get("status") == "busy" and s1.get("status") == "done" and count(c, "s1") == 1, f"{out_status(r2)} / {s1}")

    c = make_case("t7", [step("s1", "sleep 30", gpu=True)], gpu_hours=0.0005)
    t0 = time.time()
    r = drv(c)
    el = time.time() - t0
    d = ledger(c)
    ad = c / "_ledger" / "attempts" / d["attempts"][0]["id"]
    st = SP.read_status(ad) or {}
    sr = json.loads((ad / "STOP_REASON.json").read_text()) if (ad / "STOP_REASON.json").exists() else {}
    check("(T7) 予算の超過 → 先に停止の記録 (budget) を書いてから止まる、子のプロセスが残らない、30 s を待たない",
          out_status(r).get("status") == "stopped" and sr.get("cause", {}).get("class") == "budget" and st.get("state") == "stopped"
          and st.get("stop_requested", {}).get("at") and sr.get("recorded_at") <= st.get("finished") and not alive(st.get("child_pid")) and el < 25,
          f"{out_status(r)} {el:.1f} s")

    c = make_case("t8", [step("s1", "sleep 30", gpu=True)], gpu_hours=0.0008)
    p = drv(c, background=True)
    wait_for(lambda: (ledger(c) or {}).get("attempts") and ledger(c)["attempts"][0]["state"] == "launched")
    p.kill()
    p.wait()
    ad = c / "_ledger" / "attempts" / ledger(c)["attempts"][0]["id"]
    got = wait_for(lambda: (ad / "STOP_REASON.json").exists(), timeout=40)
    st = wait_for(lambda: (SP.read_status(ad) or {}).get("state") in SP.FINAL_STATES and SP.read_status(ad), timeout=40) or {}
    r2 = drv(c)
    check("(T8) driver が死んでも試行の予算の監視は続く (driver なしで停止の記録 budget)、再開で取り込んで止まる",
          got and json.loads((ad / "STOP_REASON.json").read_text())["cause"]["class"] == "budget" and st.get("state") == "stopped"
          and out_status(r2).get("status") == "stopped", f"{out_status(r2)}")

    c = make_case("t9", [step("s1", "sleep 60", gpu=True)])
    p = drv(c, background=True)
    wait_for(lambda: (ledger(c) or {}).get("attempts") and ledger(c)["attempts"][0]["state"] == "launched")
    a = ledger(c)["attempts"][0]
    ad = c / "_ledger" / "attempts" / a["id"]
    st = wait_for(lambda: (SP.read_status(ad) or {}).get("child_pid") and SP.read_status(ad), timeout=20) or {}
    time.sleep(1.0)
    us = DR.unit_state(a["unit"])
    os.kill(int(us["main_pid"]), signal.SIGKILL)
    out, _ = p.communicate(timeout=90)
    d = ledger(c)
    sr = json.loads((ad / "STOP_REASON.json").read_text()) if (ad / "STOP_REASON.json").exists() else {}
    ln = (d["budget"]["launches"] or [{}])[0]
    check("(T9) 実行担当を SIGKILL → 子が残らない、driver が停止の記録を補う (unknown)、GPU 時間を 0 にしない",
          not alive(st.get("child_pid")) and sr.get("recorded_by") == "driver" and sr.get("cause", {}).get("class") == "unknown"
          and (ln.get("seconds") or 0) > 0.5 and "根拠" in str(ln.get("end_basis")), f"{sr.get('cause')} {ln}")

    c = make_case("t10", [step("s1", "trap '' TERM; sleep 30", gpu=True)], gpu_hours=0.0003)
    p = drv(c, background=True, extra=("--stop-timeout", "6"))
    wait_for(lambda: (ledger(c) or {}).get("attempts") and ledger(c)["attempts"][0]["state"] == "launched")
    ad = c / "_ledger" / "attempts" / ledger(c)["attempts"][0]["id"]
    stopping = wait_for(lambda: (SP.read_status(ad) or {}).get("state") == "stopping", timeout=30)
    p.kill()
    p.wait()
    r2 = drv(c, extra=("--stop-timeout", "6"))
    st = SP.read_status(ad) or {}
    check("(T10) 停止の途中の再開 → 待ってから取り込む",
          stopping and out_status(r2).get("status") == "stopped" and st.get("state") == "stopped" and ledger(c)["ingested"],
          f"{out_status(r2)} {st.get('state')}")

    c = make_case("t11", [step("s1", "sleep 1")])
    r = drv(c, extra=("--min-free-bytes", "1e18"))
    check("(T11) 空きの容量が下限より小さい → 起動の前に止まる (試行を記録しない)",
          out_status(r).get("status") == "stopped" and "空きの容量" in out_status(r).get("why", "") and not ledger(c)["attempts"])

    c = make_case("t12", [step("s1", "sleep 1")])
    (c / "run_s1").mkdir()
    r = drv(c)
    check("(T12) 出力先の run が既にある → 起動の前に止まる", out_status(r).get("status") == "stopped" and "既にある" in out_status(r).get("why", "")
          and not ledger(c)["attempts"])
finally:
    for u in UNITS:
        subprocess.run(["systemctl", "--user", "stop", u], capture_output=True)
    shutil.rmtree(TD, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
