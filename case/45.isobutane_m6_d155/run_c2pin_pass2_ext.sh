#!/bin/bash
# 最終採用の条件 ① 未達 (オーバーシュート η0.1 DRIFTING・波 TRANSIENT-UNSETTLED、3 snapshot) → 登録どおり 6000 step 延長 1 回 (run_0094、1000 step ごと出力)。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
SRC=run_0092_ns_c2pin_pass2; RD=run_0094_ns_c2pin_pass2_ext6k
mkdir "$RD"
cp "$SRC"/{nozzle.h5,nozzle.xmf,bcondConfig.yaml,solverConfig.yaml,probe.yaml,prepare_info.json,target_axis_M.csv,wall_design.csv,wall_physical.csv,delta_r_initial.csv} "$RD"/
cp "$SRC"/*.yaml "$RD"/ 2>/dev/null || true
python3 ../../solver_density_cuda/tools/restart_field.py "$SRC/res_12000.h5" "$RD/nozzle.h5" --dst-run "$RD" > "$RD/restart_field.log" 2>&1
tail -1 "$RD/restart_field.log"
sed -i -E 's/nStepOuter: [0-9]+/nStepOuter: 6000/; s/outStepInterval: [0-9]+/outStepInterval: 1000/' "$RD/solverConfig.yaml"
echo "{\"extends\": \"$SRC\", \"steps\": 6000}" > "$RD/EXTENDS.json"
../../solver_density_cuda/tools/run_case.sh "$RD" > "$RD/run_case_stdout.log" 2>&1
echo "done $RD"; echo ALLDONE
