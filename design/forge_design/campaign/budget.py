"""予算の勘定 (plan tooling-design-orchestrator §4.5.1、§5.1 #6a)。承認の記録の `budget` {gpu_hours, disk_gb}。一律の上限は置かない (ユーザ決定)。

- **GPU 時間**: GPU の段の forge の起動ごとの経過時間 (GPU 1 台)。終わりの時刻が無い起動 (実行担当が死んだ) は 0 にしない:
  根拠の時刻 (止まったと分かった時刻・出力の最後の更新の時刻) までで数え、根拠を書く (`launch_seconds`)。AWS のインスタンスの稼働時間とは別。
- **ディスク**: 依頼が新たに占めた出力の、重ならない所有の範囲の大きさ (`owned_disk_bytes`; 入れ子のパスは外側だけを数える、symlink は先を数えない)。
- **予約の検査** (`check_start`): 段の見積もりが残り (承認 − 消費 − 予約) の内か。見積もりとの差は警告に使い、固定の倍率では止めない。
"""
from __future__ import annotations

import datetime as _dt
import os
from pathlib import Path


def _ts(x) -> float:
    if isinstance(x, (int, float)):
        return float(x)
    return _dt.datetime.fromisoformat(str(x)).timestamp()


def launch_seconds(launch: dict, evidence_end=None) -> dict:
    """1 回の起動の経過時間。end が無ければ evidence_end (根拠の時刻) まで。どちらも無ければ ValueError (0 にしない)。"""
    st = _ts(launch["start"])
    if launch.get("end") is not None:
        return {"seconds": max(0.0, _ts(launch["end"]) - st), "end_basis": "終わりの時刻"}
    if evidence_end is None:
        raise ValueError("終わりの時刻も根拠の時刻も無い起動 (0 として数えない)")
    return {"seconds": max(0.0, _ts(evidence_end["time"]) - st), "end_basis": f"根拠: {evidence_end['basis']}"}


def last_output_mtime(run_dir) -> dict | None:
    """run の出力の最後の更新の時刻 (終わりの時刻の無い起動の根拠)。"""
    ts = [p.stat().st_mtime for p in Path(run_dir).rglob("*") if p.is_file() and not p.is_symlink()] if Path(run_dir).is_dir() else []
    return None if not ts else {"time": max(ts), "basis": f"{run_dir} の出力の最後の更新の時刻"}


def owned_disk_bytes(paths) -> int:
    """重ならない所有の範囲の大きさ。入れ子のパスは外側だけを数え、symlink は先を数えない (リンク自体の大きさだけ)。"""
    ps = sorted({Path(os.path.abspath(p)) for p in paths}, key=lambda p: len(p.parts))
    roots = []
    for p in ps:
        if not any(r == p or r in p.parents for r in roots):
            roots.append(p)
    total = 0
    for r in roots:
        if r.is_symlink() or r.is_file():
            total += r.lstat().st_size
        elif r.is_dir():
            for dp, dns, fns in os.walk(r, followlinks=False):
                for fn in fns:
                    try:
                        total += os.lstat(os.path.join(dp, fn)).st_size
                    except FileNotFoundError:
                        pass
    return total


def remaining(budget: dict) -> dict:
    a, c, r = budget["approved"], budget["consumed"], budget["reserved"]
    return {k: float(a.get(k, 0.0)) - float(c.get(k, 0.0)) - float(r.get(k, 0.0)) for k in ("gpu_hours", "disk_gb")}


def check_start(budget: dict, estimate: dict) -> dict:
    """段を始めてよいか (見積もりが残りの内か)。戻り {ok, remaining, estimate, why}。"""
    rem = remaining(budget)
    over = [k for k in ("gpu_hours", "disk_gb") if float(estimate.get(k, 0.0)) > rem[k] + 1e-12]
    return {"ok": not over, "remaining": rem, "estimate": dict(estimate),
            "why": None if not over else f"見積もり {estimate} が残りの予算 {rem} を超える ({over}) — 次の段を始めない"}


def over_budget(budget: dict, attempt_usage: dict) -> list:
    """消費済み + この試行の使用 > 承認 になった量の列 (試行の監視が使う)。"""
    a, c = budget["approved"], budget["consumed"]
    return [k for k in ("gpu_hours", "disk_gb") if float(c.get(k, 0.0)) + float(attempt_usage.get(k, 0.0)) > float(a.get(k, 0.0)) + 1e-12]
