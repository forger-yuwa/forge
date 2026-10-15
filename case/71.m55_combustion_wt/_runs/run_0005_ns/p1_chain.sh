#!/bin/bash
# case/71 δ* の補正の 1 回目: 関門の後の段を順に回す (plan tooling-nozzle-card-local-mesh §5.1 #3e (2)〜(5) の準備まで)。AWS A。
# 各段の検査に落ちたら止める。NS の起動 (p1_ns.py run) はこの台本の外 (結果を見てから)。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p1
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD)"
test -f "$C/run_0003_ns_cont/delta_r_next.csv" || die "run_0003_ns_cont/delta_r_next.csv が無い (関門が通っていない)"
# (2) カード (r_t を物理スロート径 200 mm で解き直す)
log "(2) wt_card --ns-correction"
( cd "$D" && $PY -m forge_design.intake.wt_card "$C/card_m55.yaml" --out "$C/problem_m55_p1_card.yaml" \
    --ns-correction "$C/run_0003_ns_cont" ) > "$W/card_gen.log" 2>&1 || die "wt_card rc $? (card_gen.log)"
grep -q "throat_diameter_physical" "$C/problem_m55_p1_card.yaml" || die "カードの寸法の記録が無い"
grep "寸法" "$C/problem_m55_p1_card.yaml" | head -3 | tee -a "$W/chain.log"
test -f "$C/problem_m55_p1_card_delta_r.csv" || die "補正の表 problem_m55_p1_card_delta_r.csv が無い"
# (3) 格子 (#3b と同じ係数・ni 4,800・axis_cap_frac 0.03)
log "(3) hand_mesh build"
( cd "$C" && $PY _mesh/hand_mesh.py build --design "$D" --card card_m55.yaml --problem problem_m55_p1_card.yaml \
    --out-yaml problem_m55_p1.yaml --record "$W/build_record.json" --plan-item "§5.1 #3e" --axis-cap-frac 0.03 ) > "$W/build.log" 2>&1 \
    || die "hand_mesh build rc $? (build.log)"
test -f "$C/problem_m55_p1_euler.yaml" || die "problem_m55_p1_euler.yaml が無い"
log "(3) hand_mesh check"
rm -rf "$W/check"
( cd "$C" && $PY _mesh/hand_mesh.py check --design "$D" --card card_m55.yaml --problem problem_m55_p1.yaml \
    --record "$W/check_record.json" --run-dir "$W/check" ) > "$W/check.log" 2>&1 || die "hand_mesh check rc $? (check.log)"
grep -E "VERDICT|違反|超過|節点" "$W/check.log" | tee -a "$W/chain.log"
[ "$(grep -c 'VERDICT: PASS' "$W/check.log")" -ge 2 ] || die "格子の品質が NS・Euler とも PASS でない (check.log)"
# (4) Euler run_0004_euler
log "(4) Euler prepare + run"
test -e "$C/run_0004_euler" && die "run_0004_euler が既にある"
( cd "$D" && $PY -m forge_design.evaluate.runner_axismach "$C/problem_m55_p1_euler.yaml" "$C/run_0004_euler" --prepare-only ) \
    > "$W/euler_prepare.log" 2>&1 || die "Euler prepare rc $?"
tail -1 "$C/run_0004_euler/MESH_QUALITY.txt" | tee -a "$W/chain.log"
grep -q "VERDICT: PASS" "$C/run_0004_euler/MESH_QUALITY.txt" || die "Euler の格子の品質が PASS でない"
( cd "$D" && $PY "$C/_runs/run_0004_euler/p1_euler.py" "$C/problem_m55_p1_euler.yaml" "$C/run_0004_euler" ) \
    > "$C/run_0004_euler/euler_driver.log" 2>&1 || die "Euler rc $? (run_0004_euler/euler_driver.log)"
test -f "$C/run_0004_euler/res_12000.h5" || die "run_0004_euler/res_12000.h5 が無い"
python3 - "$C/run_0004_euler/res_12000.h5" <<'PYC' || die "Euler の最後の場に非有限"
import sys, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    bad = {k: int((~np.isfinite(np.asarray(h["VALUE/" + k][:], float))).sum()) for k in h["VALUE"]}
print("非有限", {k: v for k, v in bad.items() if v}); sys.exit(1 if any(bad.values()) else 0)
PYC
EUSHA=$(sha256sum "$C/run_0004_euler/res_12000.h5" | cut -d" " -f1)
log "Euler res_12000 sha256 $EUSHA"
# (5) NS の準備 (起動はしない)
log "(5) NS prepare"
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
( cd "$D" && $PY "$C/_runs/run_0005_ns/p1_ns.py" prepare --problem "$C/problem_m55_p1.yaml" \
    --table "$C/problem_m55_p1_card_delta_r.csv" --ic-from "$C/run_0003_ns_cont" --run-dir "$C/run_0005_ns" \
    --euler-run "$C/run_0004_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" ) > "$W/ns_prepare.log" 2>&1 || die "NS prepare rc $? (ns_prepare.log)"
tail -1 "$C/run_0005_ns/MESH_QUALITY.txt" | tee -a "$W/chain.log"
log "DONE (NS は未起動)"
