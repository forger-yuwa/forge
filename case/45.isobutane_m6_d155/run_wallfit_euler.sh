#!/bin/bash
# plan verification-m6-axis-wave-mesh-su2 §5.1 #15: 補間壁 (A) と同時当てはめ壁 (B) の Euler A/B を AWS で回す。
# 各腕の入力を 1 回作り (_prep_wallfit_{interp,fit})、そこから 3 回ずつ独立に再実行する。
# usage: bash run_wallfit_euler.sh [IC の run | none]   (none = 等エントロピー初期場 + soft 段, run_0037 と同じ起動)   (forge は runner → run_case.sh 経由で起動される)
# バイナリ: 既定は AWS の main 統合ビルド (~/forge-integ)。FORGE_BIN で上書き可。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
IC=${1:-none}
ARMS=${ARMS:-"interp fit"}   # v4 腕だけ追加するとき: ARMS=v4 START=65
ICARG=(); [ "$IC" != none ] && ICARG=(--ic "$IC")
for arm in $ARMS; do rm -rf _prep_wallfit_$arm; python3 prep_wallfit_euler.py _prep_wallfit_$arm $arm "${ICARG[@]}"; done
i=${START:-53}
for arm in $ARMS; do
  for k in 1 2 3; do
    rd=run_00${i}_euler_wallfit_${arm}_r${k}
    cp -r _prep_wallfit_${arm} "$rd"
    python3 prep_wallfit_euler.py --run "$rd" > "$rd/run_stdout.log" 2>&1 || echo "FAIL $rd"
    echo "done $rd $(tail -1 $rd/run_stdout.log)"
    i=$((i+1))
  done
done
