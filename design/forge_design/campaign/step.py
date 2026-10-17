"""試行の実行担当 (plan tooling-design-orchestrator §4.5・§4.5.1、§5.1 #6b)。systemd の user service の main process として動く。

    python -m forge_design.campaign.step <試行のディレクトリ>

試行のディレクトリの `attempt.json` (driver が起動の前に書く) を読み、段の処理 (`command`) を子のプロセスとして回し、予算を監視する。
書くのは自分の状態 (`status.json`)・停止の記録 (`STOP_REASON.json`)・ログ (`step.log`、上書きしない) だけ (台帳は driver が書く)。

- GPU 時間: GPU の段 (`gpu: true`) は、段の処理の経過時間を GPU 時間として数える (保守的: 準備の時間を含む)。起動の時刻は始めに書き、
  終わりの時刻は終わったときに書く (この担当が死ねば終わりの時刻は残らず、driver が根拠の時刻で補う)。
- ディスク: 所有の範囲 (`owned_paths`) の大きさ (`budget.owned_disk_bytes`)。監視の間隔ごと。
- 予算の超過 (台帳の消費済み (起動の前の写し) + この試行 > 承認): **先に**停止の記録 (cause budget・測った値) と状態 (stop_requested) を書き、
  その後で子のプロセスのグループを止める (SIGTERM → 期限 → SIGKILL)。この担当が死んでも service が子孫ごと止める (§5.1 #15)。
- 段の処理の終了コードと時刻を状態に書いて終わる。結果の解釈 (次の段) は driver が段の出力を読んで決める。
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from . import budget as BG
from . import stop as ST
from .ledger import _now, atomic_write_json

STATUS_FILE, ATTEMPT_FILE, LOG_FILE = "status.json", "attempt.json", "step.log"
FINAL_STATES = ("finished", "stopped", "failed")


def read_status(attempt_dir) -> dict | None:
    p = Path(attempt_dir) / STATUS_FILE
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def _write(attempt_dir, st: dict) -> None:
    st["updated"] = _now()
    atomic_write_json(Path(attempt_dir) / STATUS_FILE, st)


def _kill_group(proc: subprocess.Popen, term_wait: float) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    t_end = time.time() + term_wait
    while time.time() < t_end and proc.poll() is None:
        time.sleep(0.1)
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def run(attempt_dir) -> int:
    ad = Path(attempt_dir)
    at = json.loads((ad / ATTEMPT_FILE).read_text())
    poll = float(at.get("poll_s", 5.0))
    st = {"attempt": at["attempt"], "step": at["step"], "state": "running", "pid": os.getpid(), "started": _now(),
          "launches": [], "disk_gb": None}
    _write(ad, st)
    log = open(ad / LOG_FILE, "a")
    t0 = time.time()
    launch = {"start": t0, "end": None, "seconds": None, "end_basis": None}
    if at.get("gpu"):
        st["launches"] = [launch]
        _write(ad, st)
    env = dict(os.environ, **{k: str(v) for k, v in (at.get("env") or {}).items()})
    proc = subprocess.Popen([str(x) for x in at["command"]], cwd=at.get("cwd") or str(ad), stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, env=env, start_new_session=True)
    st["child_pid"] = proc.pid
    _write(ad, st)
    stopped = None
    while proc.poll() is None:
        time.sleep(poll)
        usage = {"gpu_hours": (time.time() - t0) / 3600.0 if at.get("gpu") else 0.0,
                 "disk_gb": BG.owned_disk_bytes(at.get("owned_paths") or []) / 1e9}
        st["usage"] = usage
        over = BG.over_budget(at["budget"], usage)
        if over and stopped is None:
            # 先に記録してから止める (codex 2026-10-17 Major 3)
            ST.write_stop(ad, attempt=at["attempt"], step=at["step"], verdicts=[], evidence_paths=[],
                          budget={"approved": at["budget"]["approved"], "consumed_before": at["budget"]["consumed"], "attempt_usage": usage},
                          process={"unit": at.get("unit"), "child_pid": proc.pid, "state": "stopping"}, cause="budget",
                          detail=f"予算に届いた {over} (承認 {at['budget']['approved']}、消費済み {at['budget']['consumed']}、この試行 {usage})")
            stopped = {"at": _now(), "over": over, "usage": usage}
            st.update(state="stopping", stop_requested=stopped)
            _write(ad, st)
            _kill_group(proc, float(at.get("term_wait_s", 30.0)))
            break
        _write(ad, st)
    rc = proc.wait()
    t1 = time.time()
    if at.get("gpu"):
        launch.update(end=t1, seconds=t1 - t0, end_basis="終わりの時刻 (実行担当が書いた)")
        st["launches"] = [launch]
    st["disk_gb"] = BG.owned_disk_bytes(at.get("owned_paths") or []) / 1e9
    st.update(rc=rc, finished=_now(), state=("stopped" if stopped else ("finished" if rc == 0 else "failed")))
    _write(ad, st)
    log.close()
    return 0 if rc == 0 and not stopped else 1


if __name__ == "__main__":
    sys.exit(run(sys.argv[1]))
