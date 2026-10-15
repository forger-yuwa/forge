#!/bin/bash
# case/71 δ* の補正の 1 回目: p1_chain.sh (NS の準備まで) の DONE を待って NS を float で起動し、終わったら成果物を作る
# (plan tooling-nozzle-card-local-mesh §5.1 #3e (5)・(8)、ユーザ「20万ステップ上限のルールで、とりあえず作ってもらえたら。フロー回して。」)。AWS A。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p1
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
for i in $(seq 1 360); do grep -q "STOP:\|DONE" "$W/chain.log" 2>/dev/null && break; sleep 20; done
grep -q "DONE" "$W/chain.log" || { log "follow: 準備が DONE でない — NS を起動しない"; exit 1; }
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
log "(5) NS run_0005_ns を起動 (float、停止規則・上限 200,000)"
( cd "$D" && $PY "$C/_runs/run_0005_ns/p1_ns.py" run --run-dir "$C/run_0005_ns" ) > "$C/run_0005_ns/p1_ns.driver.log" 2>&1
log "(5) NS の台本 終了 rc $? (停止の状態は run_0005_ns/STOP_RULE_STATE.json)"
log "(8) 成果物 (export.deliverables)"
( cd "$D" && $PY -m forge_design.export.deliverables "$C/run_0005_ns" --card "$C/card_m55.yaml" --problem "$C/problem_m55_p1.yaml" ) \
    > "$W/deliverables.log" 2>&1
log "(8) 成果物 終了 rc $? ($W/deliverables.log)"
log "FOLLOW DONE"
