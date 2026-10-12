#!/bin/bash
# case/70: CFD 前の物理壁 (CONTUR 積分法の δ_r) で NS を回す (等温壁 350 K)。壁ダンプに界面の熱量 (interfaceDiag) を出し、
# 成果物の CFD 壁面量表の熱量の単位 (plan tooling-nozzle-wt-deliverables §6 V5) を半径の違う壁で確かめるための run。
# 使い方: run_ns_ifacediag.sh EULER_RUN NS_RUN [STEPS] [IC] [PROBLEM] [STAGES]
#   IC = euler (Euler の場を移植) | isentropic (既定: prepare_ns の 1 次元等エントロピー場) | run:<NS run> (前の NS の場を移植; 壁が違えば cross-mesh)
#   PROBLEM = 問題 YAML (既定 problem_m5_air_final.yaml)、STAGES = full (既定) | none (前の NS からの warm start)
#   2026-10-13: Euler → NS の場の移植は interp_field が化学種の照合で拒否する (Euler は輸送の記録が無く NS は transport schema を持つので
#   species_hash が違う)。生産の chain は NS の場から NS を始めるので当たらない。新しい設計の最初の NS は既定の isentropic で始める。
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
eu=$1; ns=$2; steps=${3:-12000}; ic=${4:-isentropic}; prob=${5:-problem_m5_air_final.yaml}; stages=${6:-full}
icarg=(); [ "$ic" = euler ] && icarg=(--ic-from "$here/$eu")
case "$ic" in run:*) icarg=(--ic-from "$here/${ic#run:}");; esac
export FORGE_CONVERTER=$root/solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge
export FORGE_BIN=$root/solver_density_cuda/.build-native/relwithdebinfo/forge
PY=/home/sano/work/forge/design/.venv-opt/bin/python
cd "$root/design"
$PY -m forge_design.feedback.deltastar_loop --problem "$here/$prob" --euler-ref "$here/$eu" \
    --init-integral "${icarg[@]}" --run-dir "$here/$ns" --steps "$steps" --prepare-only
cfg="$here/$ns/solverConfig.yaml"
if grep -q "^output:" "$cfg"; then echo "solverConfig に output がある — 手で interfaceDiag を足す" >&2; exit 2; fi
echo "output: {interfaceDiag: 1}   # 界面の熱量 (iface_*) を壁ダンプに出す (plan tooling-nozzle-wt-deliverables §6 V5)" >> "$cfg"
$PY -c "from forge_design.evaluate.runner_axismach import run_staged_ns; import sys; sys.exit(run_staged_ns('$here/$ns', '$stages'))"
