#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11c (事前登録 2026-10-05): 第三水準の格子 (細分格子の全方向 1/1.5) で δ_E の格子ゲート。
# IC = run_0103 res_60000 (interp_field) → 段階起動 full → 本段 2 次 cfl 1・60000 step・5000 ごと → δ_E 時系列 → run_0103 との比較。② へは自動で進まない。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
R=run_0106_ns_finemesh3_pass_cfl1; SRC=run_0103_ns_finemesh_pass_cfl1
if [ ! -f $R/res_60000.h5 ]; then
  KF0=$(python3 -c "import json;print(json.load(open('c2pin_solve_pass2.json'))['k_f'])")
  python3 prep_c2pin.py $R $KF0 --problem problem_d155_ns_finemesh3_pin.yaml --ic $SRC --stages full
  sed -i -E "s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 1.0, cfl_pseudo: 1.0/; s/nStepOuter: [0-9]+/nStepOuter: 60000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" $R/solverConfig.yaml
  grep -E "cfl_pseudo|nStepOuter|outStepInterval" $R/solverConfig.yaml
  python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="full"); print("forge exit", rc); sys.exit(rc)' $R
fi
echo "done $R"
python3 - <<'PY'
import sys, os, shutil, glob, re, csv, json, subprocess, numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, "../../design")
from pathlib import Path
R, SRC = "run_0106_ns_finemesh3_pass_cfl1", "run_0103_ns_finemesh_pass_cfl1"
EU = Path("run_0086_euler_wallfit_pincal_r1_ext6k")
def one(step):
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    from forge_design.evaluate.runner_axismach import design_chain, load_problem
    xF = float(design_chain(load_problem(Path("problem_d155_ns_finemesh_pin.yaml")))["wall_inv"][-1, 0])
    src = Path(R); d = Path(f"/tmp/dE3_{step}"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
    for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"): shutil.copy(src / f, d / f)
    os.symlink((src / "nozzle.h5").resolve(), d / "nozzle.h5"); os.symlink((src / f"res_{step}.h5").resolve(), d / f"res_{step}.h5")
    sm = extract_and_merge(d, EU, band_select="edge"); nx = read_delta_r_next(d / "delta_r_next.csv"); shutil.rmtree(d)
    return dict(step=step, delta_E=float(np.interp(xF, nx["x_rt"], nx["delta_E"])), lam=sm["smooth"]["lam"], held=int(np.sum(nx["held"])))
steps = sorted(int(re.findall(r"res_(\d+)\.h5", f)[0]) for f in glob.glob(R + "/res_[0-9]*.h5")); steps = [s for s in steps if s > 0]
with ProcessPoolExecutor(3) as ex: rows = list(ex.map(one, steps))
with open(f"{R}/delta_E_series.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
S = {}
for r in (R, SRC):
    import pandas as pd
    t = pd.read_csv(f"{r}/delta_E_series.csv").delta_E.values[-5:]; S[r] = (t.mean(), (t.max() - t.min()) / 2)
D = S[R][0] / S[SRC][0] - 1; u = S[R][1] / S[R][0] + S[SRC][1] / S[SRC][0]
print(json.dumps({k: [float(v[0]), float(v[1])] for k, v in S.items()}))
print(f"GRID3 {R[:8]}/{SRC[:8]} - 1 = {100*D:+.4f} %  評価幅 [{100*(D-u):+.4f}, {100*(D+u):+.4f}] %")
PY
python3 ../../solver_density_cuda/tools/check_quasisteady.py $R --series-csv $R/delta_E_series.csv --series-cols delta_E 2>&1 | grep -E "^===.*csv|delta_E"
python3 ../../solver_density_cuda/tools/check_convergence.py $R 2>&1 | grep -E "^===|<--"
python3 ../../solver_density_cuda/tools/check_wall_resolution.py $R --mesh $R/nozzle.h5 --groups wall --over-frac 5 2>&1 | tail -3
echo ALLDONE
