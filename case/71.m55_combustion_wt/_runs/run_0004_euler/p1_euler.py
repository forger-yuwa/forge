#!/usr/bin/env python3
"""case/71 δ* の補正の 1 回目の Euler の参照 `run_0004_euler` を回す (plan tooling-nozzle-card-local-mesh §5.1 #3e (4))。
`run_0001_euler` の台本 (A `~/m55_final/m55_euler_driver.py`) と同じ: `runner_axismach.run_staged(stages="full")` → collect。
準備は先に `python -m forge_design.evaluate.runner_axismach <problem_m55_p1_euler.yaml> <run_0004_euler> --prepare-only`。"""
import json
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "design"))
from forge_design.evaluate.runner_axismach import run_staged, collect
from forge_design.probdef import load_problem
P, R = sys.argv[1], sys.argv[2]
pp = load_problem(P)
t0 = time.time()
rc = run_staged(R, cfl_main=pp.evaluate.get("cfl_main"), mid_stage=bool(pp.evaluate.get("mid_stage", pp.is_semiperfect)), stages="full")
print(f"forge exit={rc} (wall {time.time() - t0:.0f} s)", flush=True)
m = collect(P, R)
print(json.dumps(m, indent=1))
sys.exit(0 if int(rc) == 0 else 1)
