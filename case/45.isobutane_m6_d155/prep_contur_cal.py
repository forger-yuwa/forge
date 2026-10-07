"""CONTUR 較正壁の run 準備 (plan verification-m6-axis-wave-mesh-su2 §5.1 #8f)。solver は回さない。
usage: design/.venv-opt/bin/python prep_contur_cal.py <run_dir> <k_f> <k_N> <a> [--ic RUN] [--problem YAML] [--cfl CFL] [--relax R] [--no-ic]
"""
import json, sys
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import prepare_ns
C = Path(__file__).resolve().parent
run_dir = Path(sys.argv[1]); kf, kN, a = (float(v) for v in sys.argv[2:5])
arg = lambda k, dflt: (sys.argv[sys.argv.index(k) + 1] if k in sys.argv else dflt)
ic = None if "--no-ic" in sys.argv else C / arg("--ic", "run_0046_ns_band_edge")
relax = arg("--relax", "0.7")
info = prepare_ns(C / arg("--problem", "problem_d155_ns_rt77p02.yaml"), run_dir, nsteps=12000, ic_from=ic,
                  initializer={"model": "contur", "a_crocco": a, "cf_scale": kf, "n_scale": kN},
                  cfl_main=float(arg("--cfl", "5.0")), implicit_relax=(None if relax == "none" else float(relax)))
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}
(run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps({k: info[k] for k in ("mesh", "dstar_source")}, default=str), (run_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
