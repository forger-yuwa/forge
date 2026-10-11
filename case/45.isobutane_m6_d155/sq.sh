#!/bin/bash
# plan architecture-float-state-double-geometry §6.36 (事前登録、codex diagnose 2026-10-11): FP64 のビルドで step の末尾に状態を float32 に丸める十分性の試験。
# 同じ新しい FP64 のバイナリ (~/forge-fgeom10-fp64、A でビルドして B へ同じファイルを写す) で、A* = 量子化 OFF (run_0595〜0599)、G = 量子化 ON (run_0620〜0624)。
# 各 30,000 step・500 ごと、§6.26 の FP64 の腕 (run_0487) と同じ入力・設定・環境変数。G だけ FORGE_DIAG_QUANTIZE_STATE=1。
# 使い方: bash sq.sh <A|B> <run>:<0|1> ... (機械ごとに 5 本)。途中 (step 1〜19,500) の res は消す (判定は step 0 と 20,000〜30,000 だけ)。
# 系列の抽出と前提の記録は sq_series.py (各機械で)、判定は sq_an.py (2 機械の記録を合わせて)。
set -uo pipefail
M=$1; shift
case $M in A) IID=i-0b1a5e0b8dc152f00;; B) IID=i-0ba2b91ba659254c4;; *) echo "機械は A か B"; exit 1;; esac
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "$IID" ] || { echo "$M ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/sq_$M.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
F64=~/forge-fgeom10-fp64/solver_density_cuda/build/forge
cd $C45
echo "== 開始 $(date -Is) 機械 $M F64 $(sha256sum $F64) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 26843545600 ] || { echo "ディスクの空きが 25 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ -f run_0487_tr_fp64/solverConfig.yaml ] || { echo "run_0487_tr_fp64/solverConfig.yaml が無い (中止)" >> $LOG; exit 1; }
RUNS=()
for a in "$@"; do
  r=${a%%:*}; q=${a##*:}
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  diff $r/solverConfig.yaml run_0487_tr_fp64/solverConfig.yaml > $r/config_vs_run_0487.diff; nd=$(grep -c "^[<>]" $r/config_vs_run_0487.diff)
  echo "$r: run_0487 との設定の差 $nd 行" >> $LOG
  [ "$nd" = 0 ] || { echo "$r: 設定が run_0487 と違う (中止)" >> $LOG; exit 1; }
  if [ $q = 1 ]; then ( export FORGE_BIN=$F64 FORGE_DIAG_QUANTIZE_STATE=1 FORGE_DIAG_QUANTIZE_STATE_CHECK=1,30000; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$F64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) 量子化 $q" >> $LOG; RUNS+=($r)
done
sleep 240
for r in "${RUNS[@]}"; do echo "$r 起動の確認: $(grep -h "quantizeState\]" $r/forge_run.log | head -2 | cut -c1-160 | tr '\n' ' ')" >> $LOG; done
# 途中の res (step 1〜19,500) を消しながら待つ。判定は step 0 と 20,000〜30,000 だけを使う
while true; do
  left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done
  for r in "${RUNS[@]}"; do for f in $r/res_*.h5; do s=$(basename $f .h5); s=${s#res_}; [[ $s =~ ^[0-9]+$ ]] || continue
    if [ $s -gt 0 ] && [ $s -lt 20000 ]; then rm -f $f $r/res_$s.xmf; fi; done; done
  [ $left = 0 ] && break; sleep 120
done
for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 sq_series.py $M "${RUNS[@]}" >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch sq_$M.done
