#!/bin/bash
# case/71 の 2 係数の壁 (CONTUR × c + δ₀、plan tooling-nozzle-card-local-mesh §4.3・§5.1 #3o、ユーザ決定 2026-10-15「1かええな」):
# カード → 格子 → Euler run_0006_euler → NS run_0007_ns の準備を順に回す (AWS A)。各段の検査に落ちたら止める。NS の起動は p2_follow.sh。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p2
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD)"
test -f "$C/run_0005_ns/res_200000.h5" || die "run_0005_ns/res_200000.h5 (初期場) が無い"
# (2) カード (2 係数の δ で r_t を物理スロート径 200 mm に解く)
log "(2) wt_card card_m55_2coef.yaml"
( cd "$D" && $PY -m forge_design.intake.wt_card "$C/card_m55_2coef.yaml" --out "$C/problem_m55_p2_card.yaml" ) > "$W/card_gen.log" 2>&1 \
    || die "wt_card rc $? (card_gen.log)"
$PY - "$C/problem_m55_p2_card.yaml" <<'PYC' || die "カードの YAML に 2 係数が無い"
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))["deltastar_initializer"]
print("deltastar_initializer", d); sys.exit(0 if d.get("model") == "contur" and "delta_scale" in d and "delta_offset_m" in d else 1)
PYC
grep "寸法" "$C/problem_m55_p2_card.yaml" | head -3 | tee -a "$W/chain.log"
# (3) 格子 (#3b・#3e と同じ係数・axis_cap_frac 0.03)
log "(3) hand_mesh build"
( cd "$C" && $PY _mesh/hand_mesh.py build --design "$D" --card card_m55_2coef.yaml --problem problem_m55_p2_card.yaml \
    --out-yaml problem_m55_p2.yaml --record "$W/build_record.json" --plan-item "§5.1 #3o" --axis-cap-frac 0.03 ) > "$W/build.log" 2>&1 \
    || die "hand_mesh build rc $? (build.log)"
test -f "$C/problem_m55_p2_euler.yaml" || die "problem_m55_p2_euler.yaml が無い"
log "(3) hand_mesh check"
rm -rf "$W/check"
( cd "$C" && $PY _mesh/hand_mesh.py check --design "$D" --card card_m55_2coef.yaml --problem problem_m55_p2.yaml \
    --record "$W/check_record.json" --run-dir "$W/check" ) > "$W/check.log" 2>&1 || die "hand_mesh check rc $? (check.log)"
grep -E "VERDICT|違反|超過|節点" "$W/check.log" | tee -a "$W/chain.log"
[ "$(grep -c 'VERDICT: PASS' "$W/check.log")" -ge 2 ] || die "格子の品質が NS・Euler とも PASS でない (check.log)"
# (4) Euler run_0006_euler (新しい理想の形; codex diagnose 2026-10-15 ④)
log "(4) Euler prepare + run"
test -e "$C/run_0006_euler" && die "run_0006_euler が既にある"
( cd "$D" && $PY -m forge_design.evaluate.runner_axismach "$C/problem_m55_p2_euler.yaml" "$C/run_0006_euler" --prepare-only ) \
    > "$W/euler_prepare.log" 2>&1 || die "Euler prepare rc $?"
tail -1 "$C/run_0006_euler/MESH_QUALITY.txt" | tee -a "$W/chain.log"
grep -q "VERDICT: PASS" "$C/run_0006_euler/MESH_QUALITY.txt" || die "Euler の格子の品質が PASS でない"
( cd "$D" && $PY "$C/_runs/run_0004_euler/p1_euler.py" "$C/problem_m55_p2_euler.yaml" "$C/run_0006_euler" ) \
    > "$C/run_0006_euler/euler_driver.log" 2>&1 || die "Euler rc $? (run_0006_euler/euler_driver.log)"
test -f "$C/run_0006_euler/res_12000.h5" || die "run_0006_euler/res_12000.h5 が無い"
python3 - "$C/run_0006_euler/res_12000.h5" <<'PYC' || die "Euler の最後の場に非有限"
import sys, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    bad = {k: int((~np.isfinite(np.asarray(h["VALUE/" + k][:], float))).sum()) for k in h["VALUE"]}
print("非有限", {k: v for k, v in bad.items() if v}); sys.exit(1 if any(bad.values()) else 0)
PYC
EUSHA=$(sha256sum "$C/run_0006_euler/res_12000.h5" | cut -d" " -f1)
log "Euler res_12000 sha256 $EUSHA"
# (5) NS の準備 (起動はしない; 初期場は run_0005_ns の最後の場、codex ⑥)
log "(5) NS prepare"
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
( cd "$D" && $PY "$C/_runs/run_0007_ns/p2_ns.py" prepare --problem "$C/problem_m55_p2.yaml" \
    --ic-from "$C/run_0005_ns" --run-dir "$C/run_0007_ns" \
    --euler-run "$C/run_0006_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" ) > "$W/ns_prepare.log" 2>&1 || die "NS prepare rc $? (ns_prepare.log)"
tail -1 "$C/run_0007_ns/MESH_QUALITY.txt" | tee -a "$W/chain.log"
log "DONE (NS は未起動)"
