#!/bin/bash
# plan tooling-settle-convergence-and-monitor §5.1 #14・§6 V3(c) (事前登録 2026-10-16): float の v3 の長い区間の検証の run を AWS B で準備して回す。
# run_0015 と同じ入力 (準備の design は run_0015 を準備した commit dfad0267 の worktree) で、停止規則なし・本段 130,000 step・途中 restart なし。
set -uo pipefail
R=run_0017_ns_float_free130
WT=$HOME/forge-cardmesh-r17
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
D=$WT/design
S=$WT/case/71.m55_combustion_wt/_runs/run_0011_ns_v2/p5_ns.py
PROB=$WT/case/71.m55_combustion_wt/problem_m55_p5_settle_v2.yaml
PY=$HOME/venv-cardmesh/bin/python
W=$HOME/m55_p17
FP=$HOME/forge-cardmesh-float/solver_density_cuda/build
export FORGE_BIN=$FP/forge FORGE_CONVERTER=$FP/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 PYTHONDONTWRITEBYTECODE=1
unset FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE
mkdir -p "$W"
die() { echo "STOP: $*" | tee -a "$W/chain.log"; exit 1; }
log() { echo "[$(date -Iseconds)] $*" | tee -a "$W/chain.log"; }
[ "$(git -C $WT rev-parse HEAD)" = dfad026762ebd7a346e66cd31e501e7936733bd1 ] || die "$WT が dfad0267 でない"
FBSHA=$(sha256sum "$FORGE_BIN" | cut -d" " -f1)
case "$FBSHA" in 12033eb3*) ;; *) die "float の forge の sha256 が run_0015 と違う ($FBSHA)";; esac
EUSHA=$(sha256sum "$C/run_0009_euler/res_12000.h5" | cut -d" " -f1)
case "$EUSHA" in 2b2d02ce*) ;; *) die "run_0009_euler/res_12000.h5 の sha256 が違う ($EUSHA)";; esac
test -f "$C/run_0008_ns_b0/res_90000.h5" || die "run_0008_ns_b0/res_90000.h5 (初期場) が無い"
test -e "$C/$R" && die "$R が既にある"
log "git $(git -C $WT rev-parse HEAD) forge $FBSHA prepare $R (stop-rule off, nsteps 130000, precision float)"
( cd "$D" && $PY "$S" prepare --problem "$PROB" --ic-from "$C/run_0008_ns_b0" --run-dir "$C/$R" \
    --euler-run "$C/run_0009_euler" --euler-res res_12000.h5 --euler-sha256 "$EUSHA" \
    --forge-bin "$FORGE_BIN" --forge-sha256 "$FBSHA" --stop-rule off --nsteps 130000 --precision float ) > "$W/prepare_$R.log" 2>&1 \
    || die "prepare $R rc $? (prepare_$R.log)"
tail -1 "$C/$R/MESH_QUALITY.txt" | tee -a "$W/chain.log"
# run_0015 と設定を比べる: 差は本段・前の段の nStepOuter と記録の出所だけのはず
$PY - "$C/run_0015_ns_float_a" "$C/$R" <<'PYC' 2>&1 | tee -a "$W/chain.log" || die "run_0015 と設定の差が要求と違う"
import sys, yaml, hashlib
from pathlib import Path
a, b = Path(sys.argv[1]), Path(sys.argv[2])
def flat(d, pre=()):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, pre + (k,)))
    else:
        out[pre] = d
    return out
bad = []
NS = {("time", "last", "nStepOuter")}
for name, allow in [("solverConfig.yaml", NS), ("solverConfig_prestage.yaml", NS), ("bcondConfig.yaml", set()), ("settle_criteria.yaml", None)]:
    fa, fb = flat(yaml.safe_load((a / name).read_text())), flat(yaml.safe_load((b / name).read_text()))
    diff = sorted(k for k in set(fa) | set(fb) if fa.get(k, "<無>") != fb.get(k, "<無>"))
    if allow is None:
        allow = {k for k in diff if k and k[0] == "source"}
    extra = [k for k in diff if k not in allow]
    print(name, "差:", [".".join(map(str, k)) for k in diff])
    if extra:
        bad.append((name, extra))
nodes = lambda d: [l for l in (d / "settle_nodes.txt").read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
if nodes(a) != nodes(b):
    bad.append(("settle_nodes.txt", "違う"))
sh = lambda f: hashlib.sha256(f.read_bytes()).hexdigest()
if sh(a / "delta_r_initial.csv") != sh(b / "delta_r_initial.csv"):
    bad.append(("delta_r_initial.csv", "sha256 が違う (壁が違う)"))
print("設定の比較:", "OK" if not bad else f"NG {bad}")
sys.exit(1 if bad else 0)
PYC
log "run $R"
( cd "$D" && $PY "$S" run --run-dir "$C/$R" ) > "$C/$R/p5_ns.driver.log" 2>&1
log "run $R 終了 rc $?"
log "ALL DONE"
