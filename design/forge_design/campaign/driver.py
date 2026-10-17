"""campaign の driver (plan tooling-design-orchestrator §4.5.1、§5.1 #6b の最小の形; recipe の段は #6c)。

    python -m forge_design.campaign.driver <campaign.yaml> [--recipe-file F.py] [--no-systemd]

- 台帳の唯一の書き手。台帳のロックを取れなければ (ほかの driver が動いている) 何も変えずに状態を返す (`busy`)。
- 試行は「台帳に記録 → service を起動 → 起動の印 → 終わりを待つ → 一度だけ取り込み・精算」の順。
- 起動のたびに、取り込んでいない試行を照合する (`reconcile`): unit が動いている → 待つ。記録だけで起動していない (unit も状態も無い) → 今起動する
  (段の処理はまだ走っていない)。止まっていて状態が終わり → 取り込む。担当が死んで状態が終わっていない → 停止の記録を根拠つきで補い、
  終わりの時刻の無い起動は根拠の時刻まで数えて取り込む。unit の照会に失敗 → 人に返す。同じ試行の段の処理は呼び直さない。
- 起動の前の検査: 出力先の run が既にあれば止める、空きの容量 (`min_free_bytes`) が足りなければ止める、段の見積もりが残りの予算を超えれば止める。
- recipe は `next_step(ledger_doc, campaign) -> {step, command, gpu, run_dir, estimate, owned_paths, cwd, env} | {"stop": 理由} | None` と
  `interpret(step_spec, attempt, status) -> {ok, cause, detail}` を持つ物 (#6c で `axismach.wt_card_chain/v1`)。
- 試験用の障害: 環境変数 `FORGE_DRIVER_CRASH_AT` (after_record・after_start・before_ingest・after_ingest) のところで driver を落とす。
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import budget as BG
from . import campaign as CP
from . import stop as ST
from .ledger import Ledger, LedgerBusy, _now, atomic_write_json, read_ledger
from .step import ATTEMPT_FILE, FINAL_STATES, read_status

DESIGN_DIR = Path(__file__).resolve().parents[2]
CRASH_ENV = "FORGE_DRIVER_CRASH_AT"
ACTIVE = ("active", "activating", "deactivating", "reloading")


def _crash(point: str):
    if os.environ.get(CRASH_ENV) == point:
        print(f"[driver] 試験の障害: {point} で落ちる", flush=True)
        os._exit(99)


def unit_state(unit: str) -> dict:
    r = subprocess.run(["systemctl", "--user", "show", unit, "-p", "ActiveState", "-p", "LoadState", "-p", "MainPID"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return {"query_ok": False, "why": r.stderr.strip()[-200:]}
    kv = dict(ln.split("=", 1) for ln in r.stdout.splitlines() if "=" in ln)
    exists = kv.get("LoadState") not in (None, "not-found")
    return {"query_ok": True, "exists": exists, "active": kv.get("ActiveState") if exists else "missing", "main_pid": kv.get("MainPID")}


class Driver:
    def __init__(self, campaign_path, recipe, *, systemd: bool = True, poll_s: float = 1.0, step_poll_s: float = 1.0,
                 min_free_bytes: float = 200e6, stop_timeout_s: int = 30, python=None, design_dir=None):
        self.campaign_path = Path(campaign_path).resolve()
        self.campaign = CP.load(self.campaign_path)
        self.recipe = recipe
        self.systemd = systemd
        self.poll_s, self.step_poll_s = float(poll_s), float(step_poll_s)
        self.min_free_bytes = float(min_free_bytes)
        self.stop_timeout_s = int(stop_timeout_s)
        self.python = python or sys.executable
        self.design_dir = Path(design_dir or DESIGN_DIR)
        self.ledger_dir = Path(self.campaign["ledger"])
        self.L = None
        self.log = None

    # ---- 起動と待ち
    def _attempt_dir(self, aid: str) -> Path:
        return self.ledger_dir / "attempts" / aid

    def _write_attempt_file(self, a: dict, spec: dict) -> Path:
        ad = self._attempt_dir(a["id"])
        ad.mkdir(parents=True, exist_ok=True)
        b = self.L.doc["budget"]
        atomic_write_json(ad / ATTEMPT_FILE, {
            "attempt": a["id"], "step": a["step"], "unit": a["unit"], "command": spec["command"], "gpu": bool(spec.get("gpu")),
            "cwd": spec.get("cwd"), "env": spec.get("env") or {}, "owned_paths": [str(p) for p in spec.get("owned_paths") or []],
            "budget": {"approved": b["approved"], "consumed": b["consumed"]}, "poll_s": self.step_poll_s, "term_wait_s": self.stop_timeout_s})
        return ad

    def _start_unit(self, a: dict) -> None:
        ad = self._attempt_dir(a["id"])
        cmd = [self.python, "-m", "forge_design.campaign.step", str(ad)]
        if self.systemd:
            r = subprocess.run(["systemd-run", "--user", "--quiet", f"--unit={a['unit']}", "--collect",
                                "-p", "ExitType=main", "-p", "KillMode=control-group", "-p", f"TimeoutStopSec={self.stop_timeout_s}",
                                "-p", "Restart=no", f"--working-directory={self.design_dir}",
                                f"--setenv=PYTHONPATH={self.design_dir}", *cmd], capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"試行の service を起動できない ({a['unit']}): {r.stderr.strip()[-300:]}")
        else:
            subprocess.Popen(cmd, cwd=str(self.design_dir), env=dict(os.environ, PYTHONPATH=str(self.design_dir)),
                             stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)

    def _wait_inactive(self, a: dict, deadline_s: float | None = None) -> dict:
        """試行の unit が動いていない (止まった・消えた) まで待つ。systemd を使わない試験では試行の状態が終わりになるまで。"""
        t_end = None if deadline_s is None else time.time() + deadline_s
        while True:
            if self.systemd:
                us = unit_state(a["unit"])
                if not us["query_ok"] or us["active"] not in ACTIVE:
                    return us
            else:
                st = read_status(self._attempt_dir(a["id"]))
                if st and st.get("state") in FINAL_STATES:
                    return {"query_ok": True, "exists": False, "active": "missing"}
                us = {"query_ok": True}
            if t_end is not None and time.time() > t_end:
                return {**us, "timeout": True}
            time.sleep(self.poll_s)

    # ---- 取り込み
    def _ingest(self, a: dict, spec: dict | None = None) -> dict:
        ad = self._attempt_dir(a["id"])
        st = read_status(ad) or {}
        launches = []
        for ln in st.get("launches") or []:
            if ln.get("end") is None:
                ev = BG.last_output_mtime(ad) or {"time": time.time(), "basis": "driver が止まったと分かった時刻"}
                ln = {**ln, **BG.launch_seconds(ln, ev)}
            launches.append(ln)
        if st.get("state") not in FINAL_STATES:
            ST.supplement(ad, attempt=a["id"], step=a["step"], process={"unit": a["unit"], "state": "inactive", "last_status": st.get("state")},
                          budget={"launches": launches})
            st = {**st, "state": "failed", "rc": None, "supplemented": True}
        owned = [str(p) for p in (a.get("owned_paths") or [])]
        all_owned = list(dict.fromkeys(self.L.doc["owned_paths"] + owned))
        st_in = {"launches": launches, "disk_gb": BG.owned_disk_bytes(all_owned) / 1e9, "result": st.get("state"),
                 "stop": (json.loads((ad / ST.STOP_FILE).read_text()) if (ad / ST.STOP_FILE).exists() else None)}
        _crash("before_ingest")
        self.L.ingest(a["id"], st_in, owned=owned)
        _crash("after_ingest")
        return {**st, "launches": launches}

    # ---- 照合 (再開)
    def reconcile(self) -> dict | None:
        """取り込んでいない試行を照合する。人に返すときは理由の dict、それ以外は None。"""
        for a in list(self.L.doc["attempts"]):
            if a["id"] in self.L.doc["ingested"]:
                continue
            ad = self._attempt_dir(a["id"])
            us = unit_state(a["unit"]) if self.systemd else {"query_ok": True, "exists": False, "active": "missing"}
            if not us["query_ok"]:
                return {"status": "human", "why": f"試行 {a['id']} の unit の照会に失敗 ({us.get('why')})"}
            st = read_status(ad)
            if us["active"] in ACTIVE:
                if a["state"] == "recorded":
                    self.L.mark(a["id"], state="launched", launched=_now(), launched_by="reconcile")
                self._log(f"再開: 試行 {a['id']} は動いている — 待つ")
                us = self._wait_inactive(a)
                if not us["query_ok"]:
                    return {"status": "human", "why": f"試行 {a['id']} の待ちで照会に失敗"}
                self._ingest(a)
            elif a["state"] == "recorded" and st is None and not us["exists"]:
                # 記録の後・起動の前に落ちた: 段の処理はまだ走っていない → 今起動する (呼び直しではない)
                self._log(f"再開: 試行 {a['id']} は記録だけで起動していない — 起動する")
                if not (ad / ATTEMPT_FILE).exists():
                    return {"status": "human", "why": f"試行 {a['id']} の attempt.json が無い (起動の内容が分からない)"}
                self._start_unit(a)
                self.L.mark(a["id"], state="launched", launched=_now(), launched_by="reconcile")
                self._wait_inactive(a)
                self._ingest(a)
            elif st is not None:
                self._log(f"再開: 試行 {a['id']} は止まっている (状態 {st.get('state')}) — 取り込む")
                self._ingest(a)
            else:
                return {"status": "human", "why": f"試行 {a['id']} の unit も状態も無い (記録 {a['state']}) — 完了の証拠が無い"}
        return None

    def _check_run_env(self) -> dict | None:
        """走らせる環境 (forge の実行ファイル・変換器など; recipe の run_env) を初回に台帳へ写し、再開では同じかを確かめる。"""
        fn = getattr(self.recipe, "run_env", None)
        if fn is None:
            return None
        now = fn()
        rec = self.L.doc.get("run_env")
        if rec is None:
            self.L.doc["run_env"] = now
            self.L.save()
            return None
        diff = [k for k in set(rec) | set(now) if (rec.get(k) or {}).get("value") != (now.get(k) or {}).get("value")
                or (rec.get(k) or {}).get("sha256") != (now.get(k) or {}).get("sha256")]
        if diff:
            return {"status": "human", "why": f"走らせる環境が台帳を作ったときと違う {sorted(diff)} (同じ環境で再開する)"}
        return None

    def _sync_readme(self) -> None:
        try:
            from .readme import sync_runs
            sync_runs(Path(self.campaign["case_dir"]), self.L.doc)
        except Exception as e:  # noqa: BLE001 — README の同期の失敗は記録して続ける (run の判定には効かない)
            self._log(f"README の同期に失敗: {type(e).__name__}: {e}")

    def _log(self, msg: str):
        line = f"{_now()} {msg}"
        print(line, flush=True)
        if self.log:
            self.log.write(line + "\n")
            self.log.flush()

    # ---- 本体
    def run(self) -> dict:
        try:
            self.L = Ledger(self.ledger_dir).acquire()
        except LedgerBusy as e:
            return {"status": "busy", "why": str(e), "ledger": read_ledger(self.ledger_dir)}
        try:
            (self.ledger_dir / "logs").mkdir(parents=True, exist_ok=True)
            self.log = open(self.ledger_dir / "logs" / f"driver-{time.strftime('%Y%m%dT%H%M%S')}-{os.getpid()}.log", "a")
            apr = json.loads(Path(self.campaign["approval"]["file"]).read_text())
            bud = (apr.get("protected") or {}).get("budget") or {}
            self.L.init(self.campaign, {"gpu_hours": float(bud.get("gpu_hours", 0.0)), "disk_gb": float(bud.get("disk_gb", 0.0))})
            h = self._check_run_env()
            if h:
                return h
            h = self.reconcile()
            if h:
                return h
            while True:
                self.L.take_requests()
                spec = self.recipe.next_step(self.L.doc, self.campaign)
                if spec is None:
                    self._log("recipe の段が終わった")
                    return {"status": "done"}
                if "stop" in spec:
                    self._log(f"止まる: {spec['stop']}")
                    return {"status": "stopped", "why": spec["stop"]}
                free = shutil.disk_usage(self.ledger_dir).free
                if free < self.min_free_bytes:
                    return {"status": "stopped", "why": f"空きの容量 {free / 1e9:.2f} GB が下限 {self.min_free_bytes / 1e9:.2f} GB より小さい (ログ・記録を書けない)"}
                if spec.get("run_dir") and Path(spec["run_dir"]).exists():
                    return {"status": "stopped", "why": f"出力先 {spec['run_dir']} が既にある (既存の出力先に書かない)"}
                chk = BG.check_start(self.L.doc["budget"], spec.get("estimate") or {})
                if not chk["ok"]:
                    return {"status": "stopped", "why": chk["why"]}
                a = self.L.new_attempt(spec["step"], spec.get("run_dir"), BG.remaining(self.L.doc["budget"]),
                                       extra={"owned_paths": [str(p) for p in spec.get("owned_paths") or []], "estimate": spec.get("estimate")})
                self._write_attempt_file(a, spec)
                _crash("after_record")
                self._start_unit(a)
                _crash("after_start")
                self.L.mark(a["id"], state="launched", launched=_now())
                self._log(f"試行 {a['id']} を起動 ({a['unit']})")
                us = self._wait_inactive(a)
                if not us["query_ok"]:
                    return {"status": "human", "why": f"試行 {a['id']} の待ちで照会に失敗"}
                st = self._ingest(a, spec)
                res = self.recipe.interpret(spec, self.L.attempt(a["id"]), st)
                self.L.mark(a["id"], interpretation=res)
                self._sync_readme()
                if not res.get("ok"):
                    return {"status": "stopped", "why": res.get("detail"), "attempt": a["id"]}
        finally:
            if self.log:
                self.log.close()
            self.L.release()


def recipe_for(campaign: dict, repo_dir):
    """campaign の recipe の名前と版から recipe の物を作る (repo_dir: 段の処理を回す checkout の root)。"""
    rc = campaign["recipe"]
    if (rc["name"], rc["version"]) == ("axismach.wt_card_chain", "v1"):
        from .recipes.wt_card_chain import Recipe
        return Recipe(repo_dir=repo_dir)
    raise ValueError(f"知らない recipe {rc}")


def load_recipe_file(path):
    spec = importlib.util.spec_from_file_location("forge_recipe_file", str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.RECIPE


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m forge_design.campaign.driver", description=__doc__.split("\n\n")[0])
    ap.add_argument("campaign")
    ap.add_argument("--recipe-file", default=None, help="recipe の物 (RECIPE) を持つ Python のファイル (試験用; 無ければ campaign の recipe; #6c)")
    ap.add_argument("--no-systemd", action="store_true", help="試行を service にしない (試験用)")
    ap.add_argument("--poll", type=float, default=1.0)
    ap.add_argument("--step-poll", type=float, default=1.0)
    ap.add_argument("--min-free-bytes", type=float, default=200e6)
    ap.add_argument("--stop-timeout", type=int, default=30, help="試行を止めるときの期限 [s] (SIGTERM の後 SIGKILL まで・unit の停止期限)")
    a = ap.parse_args(argv)
    recipe = load_recipe_file(a.recipe_file) if a.recipe_file else recipe_for(CP.load(a.campaign), DESIGN_DIR.parent)
    d = Driver(a.campaign, recipe, systemd=not a.no_systemd, poll_s=a.poll, step_poll_s=a.step_poll,
               min_free_bytes=a.min_free_bytes, stop_timeout_s=a.stop_timeout)
    r = d.run()
    print(json.dumps({k: v for k, v in r.items() if k != "ledger"}, ensure_ascii=False))
    return {"done": 0, "busy": 0}.get(r["status"], 1)


if __name__ == "__main__":
    sys.exit(main())
