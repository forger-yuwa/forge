"""plan verification-m6-axis-wave-mesh-su2 §4.1 腕 A の run 準備 (solver は回さない)。

run_0038 と同じ壁 (δ_r = run_0025 の delta_r_next.csv, ω 1.0)・同じ warm レシピ (stages none / cfl 5 / relax 0.7 / 12000 step) で、
メッシュの軸側最大間隔だけを変える。IC は run_0038 最終場から interp_field (cross-mesh)。
usage: design/.venv-opt/bin/python prep_axis_wave.py <run_dir> <problem.yaml> [--no-ic] [--delta-csv CSV] [--ic RUN]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import prepare_ns  # noqa: E402

C = Path(__file__).resolve().parent
run_dir, problem = Path(sys.argv[1]), Path(sys.argv[2])
no_ic = "--no-ic" in sys.argv
arg = lambda k, dflt: (sys.argv[sys.argv.index(k) + 1] if k in sys.argv else dflt)
delta_csv = C / arg("--delta-csv", "run_0025_ns_ib_pass2_q/delta_r_next.csv")
ic_run = C / arg("--ic", "run_0038_ns_final_rt77p02")
info = prepare_ns(problem, run_dir, nsteps=12000, ic_from=None if no_ic else ic_run,
                  delta_r_csv=delta_csv, offset="radial",
                  euler_ref=C / "run_0037_euler_rt77p02", omega=1.0, prev_run=C / "run_0025_ns_ib_pass2_q",
                  cfl_main=5.0, implicit_relax=0.7)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}
(run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps(info["mesh"]), (run_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
