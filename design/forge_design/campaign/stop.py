"""停止の記録 `STOP_REASON.json` (plan tooling-design-orchestrator §4.3・§4.5.1、§5.1 #6a)。

欄: 試行 ID・段、判定器と VERDICT (対象量・値・閾値・判定区間)、証拠のパスと sha256、消費した予算、プロセスの状態、原因の分類 (一過性・資源・数値・照合の不一致・
予算・不明)。上書きしない (一度だけ書く)。実行担当が書けずに死んだときは、driver が根拠つきで補い (`recorded_by: driver`)、原因が確かめられなければ「不明」。
予算で止めるときは、この記録を**先に**書いてから試行を止める。
"""
from __future__ import annotations

import datetime as _dt
import hashlib
from pathlib import Path

from .ledger import atomic_write_json

STOP_FILE = "STOP_REASON.json"
CAUSES = ("transient", "resource", "numeric", "mismatch", "budget", "unknown")
REQUIRED = ("attempt", "step", "verdicts", "evidence", "budget", "process", "cause", "recorded_at", "recorded_by")


def _sha(p) -> str | None:
    p = Path(p)
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def evidence(paths) -> list:
    return [{"path": str(p), "sha256": _sha(p)} for p in paths]


def write_stop(attempt_dir, *, attempt: str, step: str, verdicts: list, evidence_paths: list, budget: dict, process: dict,
               cause: str, detail: str, recorded_by: str = "supervisor") -> dict:
    if cause not in CAUSES:
        raise ValueError(f"原因の分類 {cause!r} は {CAUSES} のどれか (分からなければ unknown)")
    p = Path(attempt_dir) / STOP_FILE
    if p.exists():
        raise FileExistsError(f"{p} が既にある — 停止の記録は上書きしない")
    rec = {"format": "forge.stop_reason/v1", "attempt": attempt, "step": step, "verdicts": list(verdicts),
           "evidence": evidence(evidence_paths), "budget": dict(budget), "process": dict(process),
           "cause": {"class": cause, "detail": detail}, "recorded_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "recorded_by": recorded_by}
    atomic_write_json(p, rec)
    return rec


def supplement(attempt_dir, *, attempt: str, step: str, process: dict, budget: dict, known_cause: str | None = None,
               basis: str | None = None) -> dict | None:
    """実行担当が停止の記録を書けずに死んだ試行に、driver が補う。既にあれば何もしない (None)。根拠が無ければ原因は unknown。"""
    if (Path(attempt_dir) / STOP_FILE).exists():
        return None
    cause = known_cause if (known_cause in CAUSES and basis) else "unknown"
    return write_stop(attempt_dir, attempt=attempt, step=step, verdicts=[], evidence_paths=[], budget=budget, process=process,
                      cause=cause, detail=(basis or "実行担当が停止の記録を書かずに終わった (原因を確かめられない)"), recorded_by="driver")


def missing_fields(rec: dict) -> list:
    return [k for k in REQUIRED if k not in rec]
