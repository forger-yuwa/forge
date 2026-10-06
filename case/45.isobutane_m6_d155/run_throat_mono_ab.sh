#!/bin/bash
# plan tooling-nozzle-throat-monotone-r2 §6 E1〜E4 (§5.1 #5・#6): 現行壁 (腕 A) と r″ 単調壁 (腕 B) の Euler A/B (生産 Euler 格子 G1、各腕 3 回)。
#  腕 A: run_0140〜0142 euler_wallfit_pinG1_r{1,2,3}   (problem_d155_euler_pin_G1_recal.yaml、IC = run_0114 最終場を restart_field)
#  腕 B: run_0143〜0145 euler_wallfit_monoG1_r{1,2,3}  (problem_d155_euler_pin_G1_recal_mono.yaml、IC = run_0114 最終場を interp_field + IC 写像の記録)
#  各腕の入力を 1 回作り (_prep_throat_mono_{pinG1,monoG1}; 壁の証拠 M4 を照合)、そこから 3 回ずつ独立に再実行する。
#  段: soft (1 次 cfl 0.5、3000 step) → 本段 2 次 cfl 2・implicitRelax 0.7・18000 step・1000 ごと出力 (forge は run_case.sh 経由)。
#  終了後: residual_history.png の確認 (無ければ plot_residual.py)、check_convergence --segment、評価器 (check_quasisteady を末尾 5 枚で内包)。
# usage: bash run_throat_mono_ab.sh [NPAR (同時実行本数, 既定 3)]
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み; run_case.sh の注記と同じ)。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
NPAR=${1:-3}
case "$NPAR" in ''|*[!0-9]*) echo "NPAR は正の整数: $NPAR"; exit 2;; esac
[ "$NPAR" -ge 1 ] || { echo "NPAR は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
IC=run_0114_euler_pin_G1_recal_ext6k
[ -f "$IC/res_6000.h5" ] || { echo "IC の run $IC/res_6000.h5 が無い — 止める"; exit 2; }
declare -A START=([pinG1]=140 [monoG1]=143)
RUNS=()
for arm in pinG1 monoG1; do
  for k in 1 2 3; do
    rd=run_0$((START[$arm] + k - 1))_euler_wallfit_${arm}_r${k}
    [ -e "$rd" ] && { echo "$rd が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; }
    RUNS+=("$rd")
  done
done
# 1. 準備 (腕ごとに 1 回)。_prep_throat_mono_* はこのスクリプト専用の作業 dir
for arm in pinG1 monoG1; do
  rm -rf "_prep_throat_mono_$arm"
  python3 throat_mono_ab.py prep "_prep_throat_mono_$arm" "$arm" 2>&1 | tee "_prep_throat_mono_$arm.log" | tail -3
done
# 2. run dir を作り、NPAR 本ずつ並列に回す
for rd in "${RUNS[@]}"; do
  arm=$(echo "$rd" | sed -E 's/.*_euler_wallfit_([A-Za-z0-9]+)_r[0-9]$/\1/')
  cp -r "_prep_throat_mono_$arm" "$rd"
done
# 終了コードは run dir の RUN_RC に書く (wait -n で回収済みの pid を後から wait すると状態が取れないため)
run_one() { local rc=0; python3 throat_mono_ab.py run "$1" > "$1/run_stdout.log" 2>&1 || rc=$?; echo "$rc" > "$1/RUN_RC"; }
for rd in "${RUNS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n || true; done
  run_one "$rd" &
  echo "started $rd (pid $!)"
done
wait || true
FAILED=()
for rd in "${RUNS[@]}"; do
  rc=$(cat "$rd/RUN_RC" 2>/dev/null || echo missing)
  if [ "$rc" = 0 ]; then echo "done $rd $(tail -1 "$rd/run_stdout.log")"; else echo "FAIL $rd (rc $rc)"; FAILED+=("$rd"); fi
done
# 3. 後処理: 残差図・収束判定 (本段の区間 = stage_manifest の最後の区間)・NaN の有無
for rd in "${RUNS[@]}"; do
  [ -f "$rd/residual_history.csv" ] || { echo "$rd: residual_history.csv が無い"; continue; }
  [ -f "$rd/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$rd/residual_history.csv" -o "$rd/residual_history.png" > /dev/null 2>&1 || echo "$rd: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$rd" --segment > "$rd/CONVERGENCE_VERDICT_segment.txt" 2>&1 || true
  echo "$rd: $(grep -m1 -E '^=== .*-> ' "$rd/CONVERGENCE_VERDICT_segment.txt" || tail -1 "$rd/CONVERGENCE_VERDICT_segment.txt")"
done
# 4. 評価 (§6 E2〜E4、固定係数の刻み感度、末尾 5 枚の check_quasisteady を内包)
mkdir -p _band_ab
python3 eval_wallfit_euler.py . --fixed-coef --e3 --pair=pinG1,monoG1 > _band_ab/throat_mono_euler_ab.log 2>&1 || echo "評価器が失敗 (rc $?) — _band_ab/throat_mono_euler_ab.log"
tail -14 _band_ab/throat_mono_euler_ab.log
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
