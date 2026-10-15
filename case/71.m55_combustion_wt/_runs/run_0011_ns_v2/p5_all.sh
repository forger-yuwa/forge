#!/bin/bash
# 収束の判定の新方式の検証 (plan tooling-settle-convergence-and-monitor §6 V3(b)・V4) を AWS B で:
# run_0011_ns_v2 (停止規則 v2) と run_0012_ns_free (停止規則なし・本段 20 万 step) を準備して起動する。
# codex diagnose 2026-10-16 settle-v3b-v4 (Minor 1) で 1 本ずつ回す: 引数の run を準備して起動する (run_0012 を先に、run_0011 はその結果 A の後)。
# 使い方: bash p5_all.sh run_0012_ns_free   /   bash p5_all.sh run_0011_ns_v2
set -uo pipefail
R=${1:?run_0012_ns_free か run_0011_ns_v2}
case "$R" in run_0011_ns_v2) M=v2;; run_0012_ns_free) M=off;; *) echo "知らない run: $R"; exit 1;; esac
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/venv-cardmesh/bin/python
S=$C/_runs/run_0011_ns_v2/p5_ns.py
W=$HOME/m55_p5
FP=$HOME/forge-cardmesh-fp64/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD) run $R (stop-rule $M)"
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
log "forge $FORGE_BIN sha256 $FBSHA"
EUSHA=$(sha256sum "$C/run_0009_euler/res_12000.h5" | cut -d" " -f1)
case "$EUSHA" in 2b2d02ce*) ;; *) die "run_0009_euler/res_12000.h5 の sha256 が記録 (2b2d02ce…) と違う ($EUSHA)";; esac
test -f "$C/run_0008_ns_b0/res_90000.h5" || die "run_0008_ns_b0/res_90000.h5 (初期場) が無い"
test -e "$C/$R" && die "$R が既にある"
log "prepare $R (stop-rule $M)"
( cd "$D" && $PY "$S" prepare --problem "$C/problem_m55_p4_settle.yaml" --ic-from "$C/run_0008_ns_b0" --run-dir "$C/$R" \
    --euler-run "$C/run_0009_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" --stop-rule "$M" --nsteps 200000 ) > "$W/prepare_$R.log" 2>&1 \
    || die "prepare $R rc $? (prepare_$R.log)"
tail -1 "$C/$R/MESH_QUALITY.txt" | tee -a "$W/chain.log"
log "run $R"
( cd "$D" && $PY "$S" run --run-dir "$C/$R" ) > "$C/$R/p5_ns.driver.log" 2>&1
log "run $R 終了 rc $?"
log "ALL DONE"
