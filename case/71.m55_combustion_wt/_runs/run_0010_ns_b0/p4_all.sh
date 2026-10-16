#!/bin/bash
# case/71 の 3 係数の壁 (plan tooling-nozzle-card-local-mesh §4.3 の末尾・§5.1 #3t、ユーザ決定 2026-10-15「bがよさそう」):
# カード → 格子 → B0 の問題 YAML → Euler run_0009_euler → NS run_0010_ns_b0 (本段 B0) を AWS B で順に回す。各段の検査に落ちたら止める。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$HOME/forge-cardmesh/design
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
W=$HOME/m55_p4
FP=$HOME/forge-integ-dflow2-fp64/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
log "git $(git -C $HOME/forge-cardmesh rev-parse HEAD)"
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
[ "$FBSHA" = 3a9e347bfac51d003c4bfc9d384f8740b85f3058df9aa486aff4806a089231f8 ] || die "forge の sha256 が A と違う ($FBSHA)"
test -f "$C/run_0008_ns_b0/res_90000.h5" || die "run_0008_ns_b0/res_90000.h5 (初期場) が無い"
log "(1) wt_card card_m55_3coef.yaml"
( cd "$D" && $PY -m forge_design.intake.wt_card "$C/card_m55_3coef.yaml" --out "$C/problem_m55_p3_card.yaml" ) > "$W/card_gen.log" 2>&1 \
    || die "wt_card rc $? (card_gen.log)"
$PY - "$C/problem_m55_p3_card.yaml" <<'PYC' || die "カードの YAML に 3 係数が無い"
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))["deltastar_initializer"]
print("deltastar_initializer", d)
sys.exit(0 if d.get("model") == "contur" and all(k in d for k in ("delta_scale", "delta_scale_slope_per_m", "delta_offset_m")) else 1)
PYC
grep "寸法\|r_throat:" "$C/problem_m55_p3_card.yaml" | head -3 | cut -c1-400 | tee -a "$W/chain.log"
log "(2) hand_mesh build + check"
( cd "$C" && $PY _mesh/hand_mesh.py build --design "$D" --card card_m55_3coef.yaml --problem problem_m55_p3_card.yaml \
    --out-yaml problem_m55_p3.yaml --record "$W/build_record.json" --plan-item "§5.1 #3t" --axis-cap-frac 0.03 ) > "$W/build.log" 2>&1 \
    || die "hand_mesh build rc $? (build.log)"
rm -rf "$W/check"
( cd "$C" && $PY _mesh/hand_mesh.py check --design "$D" --card card_m55_3coef.yaml --problem problem_m55_p3.yaml \
    --record "$W/check_record.json" --run-dir "$W/check" ) > "$W/check.log" 2>&1 || die "hand_mesh check rc $? (check.log)"
grep -E "VERDICT|違反|節点" "$W/check.log" | tee -a "$W/chain.log"
[ "$(grep -c 'VERDICT: PASS' "$W/check.log")" -ge 2 ] || die "格子の品質が NS・Euler とも PASS でない (check.log)"
log "(3) B0 の問題 YAML (cfl_main を外して main_line を足す; run_0008_ns_b0 と同じ値)"
$PY - "$C/problem_m55_p3.yaml" "$C/problem_m55_p3_b0.yaml" <<'PYC' || die "B0 の問題 YAML を作れない"
import sys
s = open(sys.argv[1]).read()
old = "  cfl_main: 1.0   # 雛形 axismach_final_v1\n"
if s.count(old) != 1:
    sys.exit("cfl_main の行がちょうど 1 つでない")
new = ("  # 本段だけの陰解法 B0 (plan §5.1 #3p・#3r; run_0008_ns_b0 と同じ値、リミッタの基準値も同じ固定)\n"
       "  main_line: {line: dir, cap: 50, itj: 5, cfl: 4, relax: 0.7, limiter_ref: {ro: 0.5810042584, p: 186732.713, a: 341.8295597}}\n")
open(sys.argv[2], "w").write("# run_0010_ns_b0 用: problem_m55_p3.yaml の写しで、evaluate の cfl_main を外して main_line (B0) を足しただけ。\n" + s.replace(old, new))
PYC
log "(4) Euler run_0009_euler"
test -e "$C/run_0009_euler" && die "run_0009_euler が既にある"
( cd "$D" && $PY -m forge_design.evaluate.runner_axismach "$C/problem_m55_p3_euler.yaml" "$C/run_0009_euler" --prepare-only ) \
    > "$W/euler_prepare.log" 2>&1 || die "Euler prepare rc $?"
grep -q "VERDICT: PASS" "$C/run_0009_euler/MESH_QUALITY.txt" || die "Euler の格子の品質が PASS でない"
( cd "$D" && $PY "$C/_runs/run_0004_euler/p1_euler.py" "$C/problem_m55_p3_euler.yaml" "$C/run_0009_euler" ) \
    > "$C/run_0009_euler/euler_driver.log" 2>&1 || die "Euler rc $? (run_0009_euler/euler_driver.log)"
test -f "$C/run_0009_euler/res_12000.h5" || die "run_0009_euler/res_12000.h5 が無い"
python3 - "$C/run_0009_euler/res_12000.h5" <<'PYC' || die "Euler の最後の場に非有限"
import sys, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    bad = {k: int((~np.isfinite(np.asarray(h["VALUE/" + k][:], float))).sum()) for k in h["VALUE"]}
print("非有限", {k: v for k, v in bad.items() if v}); sys.exit(1 if any(bad.values()) else 0)
PYC
EUSHA=$(sha256sum "$C/run_0009_euler/res_12000.h5" | cut -d" " -f1)
log "Euler res_12000 sha256 $EUSHA"
log "(5) NS prepare run_0010_ns_b0"
( cd "$D" && $PY "$C/_runs/run_0010_ns_b0/p4_ns.py" prepare --problem "$C/problem_m55_p3_b0.yaml" \
    --ic-from "$C/run_0008_ns_b0" --run-dir "$C/run_0010_ns_b0" \
    --euler-run "$C/run_0009_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" ) > "$W/ns_prepare.log" 2>&1 || die "NS prepare rc $? (ns_prepare.log)"
tail -1 "$C/run_0010_ns_b0/MESH_QUALITY.txt" | tee -a "$W/chain.log"
log "(6) NS run_0010_ns_b0 を起動 (float、本段 B0、停止規則・上限 200,000)"
( cd "$D" && $PY "$C/_runs/run_0010_ns_b0/p4_ns.py" run --run-dir "$C/run_0010_ns_b0" ) > "$C/run_0010_ns_b0/p4_ns.driver.log" 2>&1
log "(6) NS の台本 終了 rc $?"
log "FOLLOW DONE"
