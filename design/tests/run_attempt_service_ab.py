#!/usr/bin/env python3
"""試行ごとの service の契約の事前検査と、実行担当への `SIGKILL` の判別の A/B
(plan tooling-design-orchestrator §4.5「停止の責任は実行担当の外に置く」・§5.1 #15)。GPU・forge は使わない。

契約: 試行ごとに systemd の (user) service を 1 つ作り、実行担当とその全子孫 (`run_case.sh`・forge) を収める。
実行担当を main process とし、ExitType=main・KillMode=control-group・TimeoutStopSec (有限)・Restart=no。
停止規則 v2 は `run_case.sh` を start_new_session=True で起動する (stop_rule.py:1901) ので、ダミーも同じ形で起動する。

  事前検査: systemd-run --user が使えるか・systemd の版 (ExitType は 250 以降)・Linger (無いとログインのセッションが
           無くなったとき user service が止まる)。使えなければ UNAVAILABLE (rc 2)。契約では、そのとき nohup に黙って戻さず起動を止める。
  A (介入なし): 10 秒、100 ms ごとに 3 つ (実行担当・ダミーの run_case・ダミーの forge) の生存を観測。全観測で 3 つとも生きていること。
  B (実行担当だけを SIGKILL): 3 つがそろってから 2 秒の時点で実行担当だけを SIGKILL し、停止期限 5 秒の内に全子孫が消えること。
  判定: A が続き B の全子孫が期限の内に消えれば契約の成立 (PASS)。B が残れば不成立 (FAIL)。A も止まったら試験の不成立 (INVALID)。
  対照 (記録だけ、合否に入れない): service に入れず、同じ形の実行担当を直接起動して SIGKILL する。子孫が残れば、この試験は
  「残る」を検出できる (codex の観測 2026-10-17 の再現)。残った子孫は試験の最後に PID を指定して止める。

自分が起動した PID だけを止める (pgrep・pkill は使わない)。止める前に /proc/PID/cgroup がこの試験の unit に属することを確かめる。

実行: python3 design/tests/run_attempt_service_ab.py [--out result.json] [--no-control]
所要: 約 30 秒。rc: 0 PASS / 1 FAIL・INVALID / 2 UNAVAILABLE
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path

POLL_S = 0.1
A_SECONDS = 10.0
B_KILL_AT_S = 2.0
STOP_TIMEOUT_S = 5
B_OBSERVE_S = 10.0
START_WAIT_S = 10.0
EXIT_TYPE_MIN_VERSION = 250

EXECUTOR_PY = r'''
import subprocess, sys, time
from pathlib import Path
d = Path(sys.argv[1])
(d / "executor.pid").write_text(str(__import__("os").getpid()))
# 停止規則 v2 と同じく run_case を別のセッションで起動する (stop_rule.py:1901)
p = subprocess.Popen(["bash", str(d / "run_case_dummy.sh"), str(d)], stdin=subprocess.DEVNULL,
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
deadline = time.time() + 120
while time.time() < deadline:
    time.sleep(0.2)
'''

RUN_CASE_SH = r'''#!/bin/bash
d="$1"
echo $$ > "$d/run_case.pid"
# run_case.sh と同じく forge を前景の子として起動して待つ
python3 "$d/forge_dummy.py" "$d"
'''

FORGE_PY = r'''
import os, sys, time
from pathlib import Path
d = Path(sys.argv[1])
(d / "forge.pid").write_text(str(os.getpid()))
deadline = time.time() + 120
while time.time() < deadline:
    time.sleep(0.2)
'''

NAMES = ("executor", "run_case", "forge")


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def precheck() -> dict:
    rec = {"systemd_run": shutil.which("systemd-run"), "systemctl": shutil.which("systemctl")}
    if not rec["systemd_run"] or not rec["systemctl"]:
        rec["status"] = "UNAVAILABLE"
        rec["why"] = "systemd-run / systemctl が無い"
        return rec
    v = sh(["systemctl", "--version"])
    try:
        rec["systemd_version"] = int(v.stdout.split()[1])
    except (IndexError, ValueError):
        rec["systemd_version"] = None
    st = sh(["systemctl", "--user", "is-system-running"])
    rec["user_manager"] = (st.stdout or st.stderr).strip()
    lg = sh(["loginctl", "show-user", os.environ.get("USER", ""), "-p", "Linger"])
    rec["linger"] = (lg.stdout or lg.stderr).strip()
    probs = []
    if rec["systemd_version"] is None or rec["systemd_version"] < EXIT_TYPE_MIN_VERSION:
        probs.append(f"systemd の版 {rec['systemd_version']} < {EXIT_TYPE_MIN_VERSION} (ExitType が使えない)")
    if rec["user_manager"] not in ("running", "degraded"):
        probs.append(f"user の systemd が動いていない ({rec['user_manager']!r})")
    rec["problems"] = probs
    # Linger は合否に入れず警告にする (この試験は数十秒で、ログインのセッションの中で回す)。長い run の前の検査では必須
    rec["warnings"] = [] if rec["linger"] == "Linger=yes" else [
        f"{rec['linger']}: ログインのセッションが無くなると user service が止まる。長い run の前に loginctl enable-linger が要る"]
    rec["status"] = "UNAVAILABLE" if probs else "OK"
    return rec


def write_files(d: Path):
    (d / "executor.py").write_text(EXECUTOR_PY)
    (d / "run_case_dummy.sh").write_text(RUN_CASE_SH)
    (d / "forge_dummy.py").write_text(FORGE_PY)


def read_pids(d: Path) -> dict:
    out = {}
    for n in NAMES:
        f = d / f"{n}.pid"
        if f.exists():
            t = f.read_text().strip()
            if t.isdigit():
                out[n] = int(t)
    return out


def alive(pid: int) -> bool:
    try:
        st = Path(f"/proc/{pid}/stat").read_text()
    except (FileNotFoundError, ProcessLookupError):
        return False
    state = st[st.rindex(")") + 2:].split()[0]
    return state not in ("Z", "X")


def cgroup_of(pid: int) -> str:
    try:
        return Path(f"/proc/{pid}/cgroup").read_text().strip()
    except FileNotFoundError:
        return ""


def wait_pids(d: Path) -> dict:
    t_end = time.time() + START_WAIT_S
    while time.time() < t_end:
        p = read_pids(d)
        if len(p) == len(NAMES) and all(alive(x) for x in p.values()):
            return p
        time.sleep(POLL_S)
    return read_pids(d)


def observe(pids: dict, seconds: float, kill_at: float | None, kill_pid: int | None) -> dict:
    t0 = time.time()
    samples, killed_t = [], None
    while True:
        t = time.time() - t0
        if kill_at is not None and killed_t is None and t >= kill_at:
            os.kill(kill_pid, signal.SIGKILL)
            killed_t = time.time() - t0
        samples.append({"t": round(t, 3), **{n: alive(p) for n, p in pids.items()}})
        if t >= seconds:
            break
        time.sleep(POLL_S)
    return {"samples": samples, "killed_t": killed_t}


def stop_unit(unit: str):
    sh(["systemctl", "--user", "stop", unit])
    sh(["systemctl", "--user", "reset-failed", unit])


def kill_own(pids: dict):
    for p in pids.values():
        if alive(p):
            try:
                os.kill(p, signal.SIGKILL)
            except ProcessLookupError:
                pass


def arm_service(name: str, kill: bool) -> dict:
    d = Path(tempfile.mkdtemp(prefix=f"attempt_ab_{name}_"))
    write_files(d)
    unit = f"forge-attempt-ab-{name}-{uuid.uuid4().hex[:8]}"
    cmd = ["systemd-run", "--user", "--quiet", f"--unit={unit}", "--collect",
           "-p", "ExitType=main", "-p", "KillMode=control-group", "-p", f"TimeoutStopSec={STOP_TIMEOUT_S}", "-p", "Restart=no",
           sys.executable, str(d / "executor.py"), str(d)]
    rec = {"arm": name, "unit": unit, "dir": str(d), "cmd": cmd}
    r = sh(cmd)
    rec["systemd_run_rc"], rec["systemd_run_err"] = r.returncode, r.stderr.strip()
    pids = {}
    try:
        if r.returncode != 0:
            rec["verdict"] = "INVALID"
            rec["why"] = "systemd-run が失敗した"
            return rec
        pids = wait_pids(d)
        rec["pids"] = pids
        if len(pids) != len(NAMES):
            rec["verdict"] = "INVALID"
            rec["why"] = f"3 つのプロセスがそろわない ({pids})"
            return rec
        cg = {n: cgroup_of(p) for n, p in pids.items()}
        rec["cgroup"] = cg
        rec["in_unit"] = {n: (unit in c) for n, c in cg.items()}
        rec["sessions_differ"] = os.getsid(pids["executor"]) != os.getsid(pids["run_case"])
        if not all(rec["in_unit"].values()):
            rec["verdict"] = "INVALID"
            rec["why"] = "子孫がこの unit の cgroup に入っていない"
            return rec
        if kill:
            ob = observe(pids, B_OBSERVE_S, B_KILL_AT_S, pids["executor"])
            rec.update(ob)
            after = [s for s in ob["samples"] if s["t"] >= ob["killed_t"]]
            gone = [s for s in after if not any(s[n] for n in NAMES)]
            t_gone = None if not gone else gone[0]["t"] - ob["killed_t"]
            rec["all_gone_after_kill_s"] = None if t_gone is None else round(t_gone, 3)
            before_ok = all(all(s[n] for n in NAMES) for s in ob["samples"] if s["t"] < ob["killed_t"])
            if not before_ok:
                rec["verdict"], rec["why"] = "INVALID", "SIGKILL の前に止まったプロセスがある"
            elif t_gone is not None and t_gone <= STOP_TIMEOUT_S and not any(any(s[n] for n in NAMES) for s in after if s["t"] - ob["killed_t"] > t_gone):
                rec["verdict"], rec["why"] = "PASS", f"SIGKILL の {t_gone:.3f} 秒後に全子孫が消えた (期限 {STOP_TIMEOUT_S} 秒)"
            else:
                rec["verdict"], rec["why"] = "FAIL", f"期限 {STOP_TIMEOUT_S} 秒の内に全子孫が消えなかった"
        else:
            ob = observe(pids, A_SECONDS, None, None)
            rec.update(ob)
            ok = all(all(s[n] for n in NAMES) for s in ob["samples"])
            rec["verdict"] = "PASS" if ok else "INVALID"
            rec["why"] = (f"{A_SECONDS:g} 秒・{len(ob['samples'])} 回の観測で 3 つとも生きていた" if ok
                          else "介入なしで止まったプロセスがある (試験の不成立)")
        return rec
    finally:
        stop_unit(unit)
        time.sleep(0.5)
        rec["left_after_cleanup"] = {n: alive(p) for n, p in pids.items()}
        kill_own(pids)
        shutil.rmtree(d, ignore_errors=True)


def arm_control() -> dict:
    """service に入れずに同じ形で起動し、実行担当を SIGKILL する (記録だけ)。"""
    d = Path(tempfile.mkdtemp(prefix="attempt_ab_control_"))
    write_files(d)
    rec = {"arm": "control_no_service", "dir": str(d)}
    p = subprocess.Popen([sys.executable, str(d / "executor.py"), str(d)], stdin=subprocess.DEVNULL,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    pids = {}
    try:
        pids = wait_pids(d)
        rec["pids"] = pids
        if len(pids) != len(NAMES):
            rec["result"] = "INVALID"
            return rec
        ob = observe(pids, 4.0, 1.0, pids["executor"])
        p.wait(timeout=5)
        last = ob["samples"][-1]
        rec["killed_t"] = ob["killed_t"]
        rec["alive_at_end"] = {n: last[n] for n in NAMES}
        rec["result"] = ("子孫が残った (service が無いと残ることを検出できる)" if (last["run_case"] or last["forge"])
                         else "子孫も消えた (この環境では対照が残らない)")
        return rec
    finally:
        kill_own(pids)
        time.sleep(0.3)
        rec["left_after_cleanup"] = {n: alive(x) for n, x in pids.items()}
        shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", help="結果の JSON の書き出し先")
    ap.add_argument("--no-control", action="store_true", help="対照 (service なし) を回さない")
    a = ap.parse_args()
    res = {"host": os.uname().nodename, "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "precheck": precheck()}
    print("事前検査:", json.dumps({k: v for k, v in res["precheck"].items()}, ensure_ascii=False))
    if res["precheck"]["status"] != "OK":
        res["verdict"] = "UNAVAILABLE"
    else:
        res["A"] = arm_service("A", kill=False)
        print(f"A: {res['A']['verdict']} — {res['A'].get('why')}")
        res["B"] = arm_service("B", kill=True)
        print(f"B: {res['B']['verdict']} — {res['B'].get('why')}")
        if not a.no_control:
            res["control"] = arm_control()
            print(f"対照: {res['control'].get('result')}")
        va, vb = res["A"]["verdict"], res["B"]["verdict"]
        res["verdict"] = ("INVALID" if "INVALID" in (va, vb) else "PASS" if (va, vb) == ("PASS", "PASS") else "FAIL")
    left = {k: res[k].get("left_after_cleanup") for k in ("A", "B", "control") if k in res}
    print("後始末の後に残ったもの:", left)
    print(f"VERDICT: {res['verdict']}")
    if a.out:
        for k in ("A", "B", "control"):
            if k in res and "samples" in res[k]:
                res[k]["n_samples"] = len(res[k]["samples"])
        Path(a.out).write_text(json.dumps(res, ensure_ascii=False, indent=1))
    sys.exit({"PASS": 0, "UNAVAILABLE": 2}.get(res["verdict"], 1))


if __name__ == "__main__":
    main()
