#!/bin/bash
# case/71 の 2 係数の壁: p2_chain.sh (NS の準備まで) の DONE を待って NS run_0007_ns を float で起動する
# (plan tooling-nozzle-card-local-mesh §5.1 #3o)。成果物は手元で作る (A に FreeCAD・gmsh の venv が無い)。AWS A。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p2
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
for i in $(seq 1 360); do grep -q "STOP:\|DONE" "$W/chain.log" 2>/dev/null && break; sleep 20; done
grep -q "DONE" "$W/chain.log" || { log "follow: 準備が DONE でない — NS を起動しない"; exit 1; }
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
log "(6) NS run_0007_ns を起動 (float、停止規則・上限 200,000)"
( cd "$D" && $PY "$C/_runs/run_0007_ns/p2_ns.py" run --run-dir "$C/run_0007_ns" ) > "$C/run_0007_ns/p2_ns.driver.log" 2>&1
log "(6) NS の台本 終了 rc $? (停止の状態は run_0007_ns/STOP_RULE_STATE.json)"
log "FOLLOW DONE"
