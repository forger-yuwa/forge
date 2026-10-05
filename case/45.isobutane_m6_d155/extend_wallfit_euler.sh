#!/bin/bash
# plan verification-m6-axis-wave-mesh-su2 §5.1 #15: 事前登録の「準定常未達なら両腕を揃えて 6000 step 延長を 1 回」。
# run_0053〜0058 の最終場 (res_12000) を同一メッシュの restart_field で run_0059〜0064 に引き継ぎ、本段と同じ設定で 6000 step。
# usage: bash extend_wallfit_euler.sh
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
TOOLS=../../solver_density_cuda/tools
i=59
for src in run_0053_euler_wallfit_interp_r1 run_0054_euler_wallfit_interp_r2 run_0055_euler_wallfit_interp_r3 \
           run_0056_euler_wallfit_fit_r1 run_0057_euler_wallfit_fit_r2 run_0058_euler_wallfit_fit_r3; do
  rd=run_00${i}_${src#run_00??_}_ext6k
  mkdir "$rd"
  cp "$src"/{nozzle.h5,nozzle.xmf,bcondConfig.yaml,solverConfig.yaml,probe.yaml,prepare_info.json,target_axis_M.csv,wall_design.csv} "$rd"/
  cp "$src"/*.yaml "$rd"/ 2>/dev/null || true
  python3 $TOOLS/restart_field.py "$src/res_12000.h5" "$rd/nozzle.h5" --dst-run "$rd" > "$rd/restart_field.log" 2>&1
  sed -i -E 's/nStepOuter: [0-9]+/nStepOuter: 6000/; s/outStepInterval: [0-9]+/outStepInterval: 1000/' "$rd/solverConfig.yaml"
  echo "{\"extends\": \"$src\", \"steps\": 6000}" > "$rd/EXTENDS.json"
  $TOOLS/run_case.sh "$rd" > "$rd/run_case_stdout.log" 2>&1 || echo "FAIL $rd"
  echo "done $rd $(grep -m1 -E 'VERDICT|-> ' "$rd/restart_field.log" | head -1)"
  i=$((i+1))
done
