#!/bin/bash
# case/71 run_0008_ns_b0 (本段の陰解法 B0、plan tooling-nozzle-card-local-mesh §5.1 #3p・§6 V6) を AWS B で準備して回す。
# 入力 (run_0005_ns の最後の場・run_0006_euler・バイナリ) は A から写した (sha256 を照合)。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p3
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD)"
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
[ "$FBSHA" = 3a9e347bfac51d003c4bfc9d384f8740b85f3058df9aa486aff4806a089231f8 ] || die "forge の sha256 が A と違う ($FBSHA)"
EUSHA=cbd236ec362eb40f9ae14da4e28422e01fe6752b74db84e71b58198a3b796245
[ "$(sha256sum "$C/run_0006_euler/res_12000.h5" | cut -d" " -f1)" = "$EUSHA" ] || die "run_0006_euler/res_12000.h5 の sha256 が違う"
test -f "$C/run_0005_ns/res_200000.h5" || die "run_0005_ns/res_200000.h5 が無い"
log "(1) NS prepare run_0008_ns_b0"
( cd "$D" && $PY "$C/_runs/run_0008_ns_b0/p3_ns.py" prepare --problem "$C/problem_m55_p2_b0.yaml" \
    --ic-from "$C/run_0005_ns" --run-dir "$C/run_0008_ns_b0" \
    --euler-run "$C/run_0006_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" ) > "$W/ns_prepare.log" 2>&1 || die "NS prepare rc $? (ns_prepare.log)"
tail -1 "$C/run_0008_ns_b0/MESH_QUALITY.txt" | tee -a "$W/chain.log"
log "(2) NS run_0008_ns_b0 を起動 (float、本段 B0、停止規則・上限 200,000)"
( cd "$D" && $PY "$C/_runs/run_0008_ns_b0/p3_ns.py" run --run-dir "$C/run_0008_ns_b0" ) > "$C/run_0008_ns_b0/p3_ns.driver.log" 2>&1
log "(2) NS の台本 終了 rc $? (停止の状態は run_0008_ns_b0/STOP_RULE_STATE.json)"
log "FOLLOW DONE"
