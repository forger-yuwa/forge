"""driver の試験用の recipe (design/tests/run_driver_tests.py)。段の列は環境変数 FORGE_DUMMY_STEPS の JSON ファイル。

各段 {step, command, gpu, run_dir}: command の中の "{RUN}" は run のディレクトリ (case_dir/run_dir) に置き換える。
取り込んだ試行の状態が finished の段は済み。取り込んだが finished でない段があれば止まる (自動の再試行はしない)。
"""
import json
import os
from pathlib import Path


class DummyRecipe:
    def __init__(self, spec_file):
        self.steps = json.loads(Path(spec_file).read_text())

    def next_step(self, doc, campaign):
        done = {a["step"] for a in doc["attempts"] if a.get("result") == "finished"}
        failed = {a["step"] for a in doc["attempts"] if a["id"] in doc["ingested"] and a.get("result") != "finished"}
        for s in self.steps:
            if s["step"] in done:
                continue
            if s["step"] in failed:
                return {"stop": f"段 {s['step']} が失敗した"}
            run = str(Path(campaign["case_dir"]) / s["run_dir"]) if s.get("run_dir") else None
            cmd = [c.replace("{RUN}", run or "") for c in s["command"]]
            return {"step": s["step"], "command": cmd, "gpu": bool(s.get("gpu")), "run_dir": run,
                    "owned_paths": [run] if run else [], "estimate": s.get("estimate") or {}, "cwd": campaign["case_dir"]}
        return None

    def run_env(self):
        """走らせる環境の固定の試験用 (FORGE_DUMMY_ENV の値; 無ければ固定しない)。"""
        v = os.environ.get("FORGE_DUMMY_ENV")
        return {} if v is None else {"FORGE_DUMMY_ENV": {"value": v, "sha256": None}}

    def interpret(self, spec, attempt, status):
        return {"ok": status.get("state") == "finished", "detail": f"段 {spec['step']}: {status.get('state')}"}


RECIPE = DummyRecipe(os.environ["FORGE_DUMMY_STEPS"])
