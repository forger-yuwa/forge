"""P5 (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #7): CFD ピン設計の NS run 準備 (solver は回さない)。
usage: python3 prep_c2pin.py <run_dir> <k_f> [--problem YAML] [--ic RUN] [--stages full|none]
δ は CONTUR 積分法 (cf_scale = k_f、出口 δ に較正; C2 方式)。IC 無しなら段階起動 full (soft → mid → 本段)。"""
import json, sys
from pathlib import Path
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import prepare_ns  # noqa: E402
run_dir, kf = Path(sys.argv[1]), float(sys.argv[2])
arg = lambda k, dflt: (sys.argv[sys.argv.index(k) + 1] if k in sys.argv else dflt)
ic = arg("--ic", None)
info = prepare_ns(C / arg("--problem", "problem_d155_ns_c2pin.yaml"), run_dir, nsteps=12000, ic_from=(C / ic if ic else None),
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=5.0, implicit_relax=0.7)
info["stages"] = {"stages": arg("--stages", "full" if ic is None else "none"), "ramp": None, "ramp_steps": 1000}
(run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps({k: info.get(k) for k in ("mesh", "dstar_source", "initial_line", "Md_moc_offset", "wall_repr", "scale_m")}, default=str))
print((run_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
