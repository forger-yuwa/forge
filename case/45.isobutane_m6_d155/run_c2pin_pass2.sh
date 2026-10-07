#!/bin/bash
# P5 手 1 (plan tooling-nozzle-cfd-pinned-initial-line §9 2026-10-05 事前登録): C2 の 2 pass 目。
# run_0090 の場から E で δ_E を取り直し k_f・r_t を解き直して、IC = run_0090、同設定 12000 step (run_0092)。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
EU=run_0086_euler_wallfit_pincal_r1_ext6k; FN=run_0090_ns_c2pin_final; P2=run_0092_ns_c2pin_pass2
cp c2pin_solve.json c2pin_solve_pass1.json
# c2pin_solve は基準 problem (c2pin, r_t 76.8075) から solve_rt するが、prev_run の scale_m を r_t,prev に使うので run_0090 (r_t 76.773) からでも正しく換算される
python3 c2pin_solve.py $FN $EU
cp c2pin_solve.json c2pin_solve_pass2.json
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve.json'))['k_f'])")
python3 prep_c2pin.py $P2 $KF --problem problem_d155_ns_c2pin_final.yaml --ic $FN
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)' $P2
echo "done $P2"; echo ALLDONE
