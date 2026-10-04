#!/usr/bin/env bash
# R1 (plan boundary-node-periodic-gradient-fix §6 R1) の 1 本分を回す: S1_spin → (成形 IC) → S2_dev。
# 旧 (1266aba1) / 新 (565959c7) で**同じスクリプト・同じメッシュ h5・同じ config** を使い、FORGE_BIN だけ変える。
#
#   r1_driver.sh RUN_DIR FORGE_BIN
#
# 前提: mesh/hill_r1_80x50x30.h5 (node 変換済み、品質 PASS) があること。
# 実行中にこのファイルを編集しないこと (run_case.sh と同じ理由: bash はバイト位置で逐次読む)。
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
RUN="$(cd "$1" && pwd)"; BIN="$2"
MESH=hill_r1_80x50x30.h5
cd "$RUN"
# S2_ONLY=1: S1 の結果 (S1_spin/) を残したまま S2 だけ回し直す (未収束で nStepOuter を伸ばしたとき。S2 は ic_S2_dev.h5 から最初から)
if [ "${S2_ONLY:-0}" = 1 ]; then
  [ -e residual_history_S1_spin.csv ] || { echo "S1 の跡が無いのに S2_ONLY: $RUN"; exit 1; }
  rm -f res_*.h5 res_*.xmf forge_run.log residual_history.csv residual_history.png residual_history_S2_dev.csv \
        residual_history_segment.csv stage_manifest.json CONVERGENCE_VERDICT.txt CONVERGENCE_VERDICT_segment.txt \
        RUN_PROVENANCE.txt r1_series.csv r1_bulk.csv
else
[ -e residual_history_S1_spin.csv ] && { echo "既に S1 を回した跡がある: $RUN"; exit 1; }
cp "$HERE/mesh/$MESH" "$MESH"

echo "[r1] S1_spin ($(date -Is))"
python3 "$HERE/r1_make_ic.py" static "$MESH" ic_S1_spin.h5
cp solverConfig_S1_spin.yaml solverConfig.yaml
FORGE_BIN="$BIN" "$ROOT/solver_density_cuda/tools/run_case.sh" "$RUN" || { echo "[r1] S1 forge 失敗"; exit 2; }
mkdir -p S1_spin
mv res_*.h5 res_*.xmf S1_spin/ 2>/dev/null || true
mv forge_run.log CONVERGENCE_VERDICT.txt residual_history.png RUN_PROVENANCE.txt S1_spin/ 2>/dev/null || true
mv residual_history.csv residual_history_S1_spin.csv
last=S1_spin/res_$(ls S1_spin/ | sed -n 's/^res_\([0-9]\+\)\.h5$/\1/p' | sort -n | tail -1).h5
echo "[r1] S1 last = $last"
fi
last=S1_spin/res_$(ls S1_spin/ | sed -n 's/^res_\([0-9]\+\)\.h5$/\1/p' | sort -n | tail -1).h5

echo "[r1] S2_dev ($(date -Is))"
python3 "$HERE/r1_make_ic.py" shaped "$last" "$MESH" ic_S2_dev.h5
cp solverConfig_S2_dev.yaml solverConfig.yaml
FORGE_BIN="$BIN" "$ROOT/solver_density_cuda/tools/run_case.sh" "$RUN" || { echo "[r1] S2 forge 失敗"; exit 3; }
cp residual_history.csv residual_history_S2_dev.csv
python3 "$HERE/r1_make_ic.py" manifest "$RUN"
python3 "$ROOT/solver_density_cuda/tools/stage_manifest.py" "$RUN" --segments
python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" "$RUN" --segment > CONVERGENCE_VERDICT_segment.txt 2>&1 || true
cat CONVERGENCE_VERDICT_segment.txt
echo "[r1] done ($(date -Is))"
