#!/bin/bash
# plan tooling-nozzle-mesh-size-guard §6 V1b: Euler の解の ni 依存 (ni 1800・900・600、soft + 本段 cfl 1・48,000 step、run_0013 と同じ段取り)
set -e
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../.." && pwd)
export FORGE_BIN=$root/solver_density_cuda/.build-native/relwithdebinfo/forge
export FORGE_CONVERTER=$root/solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge
PY=/home/sano/work/forge/design/.venv-opt/bin/python
cd "$root/design"
for spec in "run_0015_euler_c2b_ni1800:problem_m5_air_c2b_euler.yaml" "run_0016_euler_c2b_ni900:problem_m5_air_c2b_euler_ni900.yaml" "run_0017_euler_c2b_ni600:problem_m5_air_c2b_euler_ni600.yaml"; do
  run=${spec%%:*}; prob=${spec##*:}
  echo "== $run $(date)"
  $PY -m forge_design.evaluate.runner_axismach "$here/$prob" "$here/$run" --steps 48000 --stages soft > "$here/_$run.log" 2>&1
  echo "rc=$? $run $(date)"
done
echo ALLDONE
