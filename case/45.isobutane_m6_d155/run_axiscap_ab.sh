#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11h (事前登録 2026-10-06): 最終壁 (run_0117) を固定し ni 2000・nj 257 で axis_cap_frac だけを変える NS A/B。
#  腕 A run_0124_ns_axiscap020 (cap 0.02) / 腕 B run_0125_ns_axiscap0133 (cap 0.0133)。IC = run_0117 res_60000 を interp_field → 段階起動 full → 本段 2 次 cfl 1・60000 step・5000 ごと。
set -e
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_recal.json'))['k_f'])")
RUN='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="full"); print("forge exit", rc); sys.exit(rc)'
for arm in "run_0124_ns_axiscap020 cap020" "run_0125_ns_axiscap0133 cap0133"; do
  set -- $arm
  python3 prep_c2pin.py $1 $KF --problem problem_d155_ns_recal_final_$2.yaml --ic run_0117_ns_recal_final_ext --stages full
  sed -i -E "s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 1.0, cfl_pseudo: 1.0/; s/nStepOuter: [0-9]+/nStepOuter: 60000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" $1/solverConfig.yaml
  python3 -c "$RUN" $1 || echo "FAILED $1"
done
for R in run_0124_ns_axiscap020 run_0125_ns_axiscap0133; do
  python3 axisgrid_metrics.py $R || true
  EULER_REF=run_0114_euler_pin_G1_recal_ext6k SOLVE_JSON=c2pin_solve_recal.json FINAL_PROBLEM=problem_d155_ns_finemesh_recal_final.yaml python3 exitM_sampling_ab.py $R run_0117_ns_recal_final_ext || true
done
python3 axisgrid_metrics.py run_0117_ns_recal_final_ext || true
echo ALLDONE
