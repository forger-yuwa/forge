"""campaign の台帳 (plan tooling-design-orchestrator §4.5.1、§5.1 #6a)。

- driver が**唯一の書き手**。生きている間ずっと `<ledger>/ledger.lock` の flock (排他・待たない) を持ち、取れなければ `LedgerBusy`
  (二つ目の `--resume` は今の driver を止めずに状態を返す)。更新は一時ファイル → fsync → rename。
- campaign の ID は初めて台帳を作るときの UUID (名前は表示だけ; 同じ名前の別の case の campaign と衝突しない)。
- 試行は「記録 → 起動 → 取り込み」の順。`new_attempt` が試行 ID・unit 名・run の出力先・予算の予約を**起動の前に**書き、`ingest` が終わった試行の結果を
  試行 ID ごとに**一度だけ**取り込んで精算する (2 回呼んでも精算は 1 回)。
- 変更の要求 (CLI の `note` など) は `<ledger>/requests/*.json` に置かれ、driver が `take_requests` で取り込む。

読むだけ (`read_only`) なら lock は要らない (CLI の `status`)。
"""
from __future__ import annotations

import datetime as _dt
import fcntl
import json
import os
import uuid
from pathlib import Path

FORMAT = "forge.campaign_ledger/v1"
LEDGER_FILE, LOCK_FILE = "ledger.json", "ledger.lock"


class LedgerBusy(RuntimeError):
    """ほかの driver が台帳のロックを持っている。"""


def _now() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def atomic_write_json(path: Path, doc) -> None:
    path = Path(path)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(tmp, "w") as f:
        json.dump(doc, f, indent=1, ensure_ascii=False, default=str)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def read_ledger(ledger_dir) -> dict | None:
    p = Path(ledger_dir) / LEDGER_FILE
    return json.loads(p.read_text()) if p.is_file() else None


class Ledger:
    def __init__(self, ledger_dir):
        self.dir = Path(ledger_dir)
        self.path = self.dir / LEDGER_FILE
        self._fd = None
        self.doc = None

    # ---- ロック
    def acquire(self):
        self.dir.mkdir(parents=True, exist_ok=True)
        fd = os.open(self.dir / LOCK_FILE, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd)
            raise LedgerBusy(f"台帳 {self.dir} はほかの driver が書いている (ロックを取れない)") from None
        os.ftruncate(fd, 0)
        os.write(fd, f"{os.getpid()} {_now()}\n".encode())
        self._fd = fd
        self.doc = read_ledger(self.dir)
        return self

    def release(self):
        if self._fd is not None:
            fcntl.flock(self._fd, fcntl.LOCK_UN)
            os.close(self._fd)
            self._fd = None

    def __enter__(self):
        return self.acquire()

    def __exit__(self, *exc):
        self.release()

    def _need_lock(self):
        if self._fd is None:
            raise RuntimeError("台帳のロックを持たずに書こうとした")

    def save(self):
        self._need_lock()
        self.doc["updated"] = _now()
        atomic_write_json(self.path, self.doc)

    # ---- 初期化
    def init(self, campaign: dict, approved_budget: dict) -> dict:
        """台帳が無ければ作る (campaign の ID を決める)。あれば campaign の承認の記録・recipe・役割が同じかを確かめる (役割は変えられない)。"""
        self._need_lock()
        ap = campaign["approval"]
        if self.doc is None:
            self.doc = {"format": FORMAT, "campaign_id": str(uuid.uuid4()), "name": campaign["name"], "role": campaign["role"],
                        "approval": {"file": ap["file"], "sha256": ap["sha256"]}, "recipe": campaign["recipe"], "case_dir": campaign["case_dir"],
                        "created": _now(), "attempts": [], "ingested": [], "owned_paths": [], "requests": [], "notes": [],
                        "budget": {"approved": dict(approved_budget), "reserved": {"gpu_hours": 0.0, "disk_gb": 0.0},
                                   "consumed": {"gpu_hours": 0.0, "disk_gb": 0.0}, "launches": []},
                        "drivers": [],
                        # 段の見積もり (plan §4.5.1「段の見積もりの契約」): 台帳を作るときに campaign から写して固定する (再開で違えば止まる)
                        "estimates": dict(campaign.get("estimates") or {})}
        else:
            if self.doc.get("format") != FORMAT:
                raise ValueError(f"{self.path} は台帳でない")
            diffs = [k for k, a, b in (("approval", self.doc["approval"]["sha256"], ap["sha256"]),
                                       ("recipe", self.doc["recipe"], campaign["recipe"]), ("role", self.doc["role"], campaign["role"]),
                                       ("estimates", self.doc.get("estimates") or {}, campaign.get("estimates") or {}))
                     if a != b]
            if diffs:
                raise ValueError(f"台帳と campaign が違う {diffs} (役割・承認・recipe・見積もりは台帳を作った後に変えられない; 新しい campaign にする)")
        self.doc["drivers"].append({"pid": os.getpid(), "started": _now()})
        self.save()
        return self.doc

    # ---- 試行
    def new_attempt(self, step: str, run_dir, reservation: dict, extra: dict | None = None) -> dict:
        """試行を起動の**前**に記録する (状態 recorded)。試行 ID は通し番号 + 段の名前、unit 名は campaign の ID と試行 ID から。"""
        self._need_lock()
        n = len(self.doc["attempts"]) + 1
        aid = f"{n:04d}-{step}"
        cid = self.doc["campaign_id"]
        rec = {"id": aid, "step": step, "unit": f"forge-step-{cid[:8]}-{aid}", "run_dir": None if run_dir is None else str(run_dir),
               "reservation": dict(reservation), "state": "recorded", "recorded": _now(), **(extra or {})}
        for k, v in reservation.items():
            self.doc["budget"]["reserved"][k] = float(self.doc["budget"]["reserved"].get(k, 0.0)) + float(v)
        self.doc["attempts"].append(rec)
        self.save()
        return rec

    def attempt(self, aid: str) -> dict:
        a = next((x for x in self.doc["attempts"] if x["id"] == aid), None)
        if a is None:
            raise KeyError(f"試行 {aid} が台帳に無い")
        return a

    def mark(self, aid: str, **fields) -> dict:
        self._need_lock()
        a = self.attempt(aid)
        a.update(fields)
        self.save()
        return a

    def ingest(self, aid: str, status: dict, owned: list | None = None) -> bool:
        """終わった試行の結果を一度だけ取り込み、予約を外して消費を精算する。既に取り込んだ試行なら何もせず False。
        status: 試行の状態 (`launches` = forge の起動ごとの {start, end, seconds, end_basis}、`disk_gb` = 所有の範囲の大きさ (測った値))。"""
        self._need_lock()
        if aid in self.doc["ingested"]:
            return False
        a = self.attempt(aid)
        b = self.doc["budget"]
        for k, v in (a.get("reservation") or {}).items():
            b["reserved"][k] = max(0.0, float(b["reserved"].get(k, 0.0)) - float(v))
        gpu_s = sum(float(x.get("seconds") or 0.0) for x in status.get("launches") or [])
        b["consumed"]["gpu_hours"] = float(b["consumed"]["gpu_hours"]) + gpu_s / 3600.0
        b["launches"] += [{"attempt": aid, **x} for x in status.get("launches") or []]
        for pth in owned or []:
            if str(pth) not in self.doc["owned_paths"]:
                self.doc["owned_paths"].append(str(pth))
        if status.get("disk_gb") is not None:
            b["consumed"]["disk_gb"] = float(status["disk_gb"])     # 所有の範囲全体を測った値 (足し合わせない; 二重に数えない)
        a.update(state="ingested", ingested=_now(), result=status.get("result"), stop=status.get("stop"))
        self.doc["ingested"].append(aid)
        self.save()
        return True

    # ---- 変更の要求
    def take_requests(self) -> list:
        self._need_lock()
        rq = self.dir / "requests"
        got = []
        for f in sorted(rq.glob("*.json")) if rq.is_dir() else []:
            doc = json.loads(f.read_text())
            doc["file"] = f.name
            doc["taken"] = _now()
            self.doc["requests"].append(doc)
            got.append(doc)
            f.rename(f.with_suffix(".taken"))
        if got:
            self.save()
        return got


def put_request(ledger_dir, kind: str, **fields) -> Path:
    """CLI が置く変更の要求 (台帳には書かない; driver が取り込む)。"""
    rq = Path(ledger_dir) / "requests"
    rq.mkdir(parents=True, exist_ok=True)
    p = rq / f"{_dt.datetime.now().strftime('%Y%m%dT%H%M%S%f')}-{uuid.uuid4().hex[:6]}.json"
    atomic_write_json(p, {"kind": kind, "made": _now(), **fields})
    return p
