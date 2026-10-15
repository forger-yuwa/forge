#!/bin/bash
# plan tooling-settle-convergence-and-monitor §4.9 の 6 (V6、事前登録 2026-10-16): case/71 を設計チェーンの関数で 1 pass 通す。
# deltastar_loop.run_pass0_integral (problem_m55_p5_settle_v2.yaml = 3 係数の壁・判定 wt_nozzle_v2、Euler の参照 run_0009_euler res_12000、
# 初期場 run_0008_ns_b0) → NS (停止規則 v2) → 新方式の受け渡しの関門 → delta_r_next.csv → 成果物 (deliverables)。AWS B で。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/venv-cardmesh/bin/python
R=$C/run_0014_chain_v2
W=$HOME/m55_p6; mkdir -p $W
FP=$HOME/forge-cardmesh-fp64/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
die() { log "STOP: $*"; exit 1; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD) forge sha256 $(sha256sum $FORGE_BIN | cut -c1-16)"
test -e "$R" && die "$R が既にある"
EUSHA=$(sha256sum "$C/run_0009_euler/res_12000.h5" | cut -d" " -f1)
case "$EUSHA" in 2b2d02ce*) ;; *) die "Euler の参照の sha256 が違う";; esac
log "run_pass0_integral → NS (停止規則 v2) → 関門 → delta_r_next.csv"
( cd "$D" && $PY - "$C/problem_m55_p5_settle_v2.yaml" "$C/run_0009_euler" "$EUSHA" "$R" "$C/run_0008_ns_b0" <<'PYC'
import sys, json
from forge_design.feedback import deltastar_loop as DL
prob, erun, esha, run, ic = sys.argv[1:6]
eref = {"run": erun, "res": "res_12000.h5", "sha256": esha}
out = DL.run_pass0_integral(prob, eref, run, ic_from=ic, stop_rule=True)
print("RESULT", json.dumps(out if isinstance(out, (dict, list, str, int, float, type(None))) else str(out), ensure_ascii=False, default=str)[:3000])
PYC
) > "$W/pass0.log" 2>&1
log "pass0 rc $?"
tail -5 "$W/pass0.log" | cut -c1-400 | tee -a "$W/chain.log"
log "deliverables"
( cd "$D" && $PY -m forge_design.export.deliverables "$R" --card "$C/card_m55_3coef.yaml" --problem "$C/problem_m55_p5_settle_v2.yaml" --skip cad ) > "$W/deliverables.log" 2>&1
log "deliverables rc $?"
tail -5 "$W/deliverables.log" | cut -c1-400 | tee -a "$W/chain.log"
log "ALL DONE"
