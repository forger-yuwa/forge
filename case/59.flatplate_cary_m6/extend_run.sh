#!/usr/bin/env bash
# 同一設定・同一メッシュで継続する (新しい run に restart_field.py で引き継ぐ)。完了後に中間の全場を消す。
#   bash extend_run.sh <src_run> <dst_run> <steps> <out_int> <forge_tools_dir>
set -eu
cd "$(dirname "$0")"
S=$1; D=$2; N=$3; OUT=$4; FT=$5
ROOT=$(git rev-parse --show-toplevel)
test ! -e "$D"
last=$(ls "$S" | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//;s/\.h5//' | sort -n | tail -1)
mkdir "$D"
cp "$S"/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json} "$D"/
[ -f "$S/species_db.yaml" ] && cp "$S/species_db.yaml" "$D"/   # 多成分 (燃焼ガス) の run に要る
sed -e "s/nStepOuter: [0-9]*/nStepOuter: $N/" -e "s/outStepInterval: [0-9]*/outStepInterval: $OUT/" "$S/solverConfig.yaml" > "$D/solverConfig.yaml"
python3 "$ROOT/solver_density_cuda/tools/restart_field.py" --keep-src-dtype "$S/res_$last.h5" "$D/mesh.h5" | tail -1   # FP64 の場は倍精度のまま渡す
echo "$S/res_$last.h5 (restart_field.py、同一設定 +$N step)" > "$D/CONTINUED_FROM"
python3 - "$D" "$S/res_$last.h5" <<PY
import sys; sys.path.insert(0, "$ROOT/solver_density_cuda/tools")
from stage_manifest import StageManifest
d, src = sys.argv[1], sys.argv[2]
sm = StageManifest(d); sm.add("main", open(d + "/solverConfig.yaml").read(), open(d + "/bcondConfig.yaml").read(),
                              history="residual_history.csv", restart_from=src); sm.write()
PY
FORGE_CUDA_BLOCKSIZE=128 "$FT/run_case.sh" "$(pwd)/$D" > "_ext_$D.log" 2>&1 || true
python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" "$D" > "$D/CONVERGENCE_VERDICT_check.txt" 2>&1 || true
lastd=$(ls "$D" | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//;s/\.h5//' | sort -n | tail -1)
for f in "$D"/res_[0-9]*.h5; do b=$(basename "$f" .h5); [ "$b" = "res_$lastd" ] || [ "$b" = "res_0" ] || rm -f "$f" "$D/$b.xmf"; done
echo "[$D] EXT 完了 $(grep -m1 '^===' "$D/CONVERGENCE_VERDICT_check.txt")"
