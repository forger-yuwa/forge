#!/bin/bash
# plan architecture-float-state-double-geometry §6.34 (事前登録、codex diagnose 2026-10-11): 平均流の状態を FP64 で積む float の軌跡 (qAccumulatorFP64 1、軸対称対応)。
# §6.26 と同じ入力 (_tr/q32_init.h5)・環境変数・設定で、E = run_0592_qa_f32・E' = run_0593_qa_f32b (qAccumulatorFP64 1) と
# 対照 F = run_0594_qa_off (同じ新しいバイナリ、qAccumulatorFP64 0)。各 30,000 step・500 ごと。インスタンス A で 3 本同時 (B は他セッションの run で埋まっていたため。A と B は同じ型・同じ AMI)。判定は qa_an.py。
# バイナリ ~/forge-fgeom9-f32 (qAccumulatorFP64 の軸対称対応と採用の mask の診断を含む commit)。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/qa.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom9-f32/$B
cd $C45
echo "== 開始 $(date -Is) F32 $(sha256sum $F32) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 32212254720 ] || { echo "ディスクの空きが 30 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <qAccumulatorFP64 0|1>
  local r=$1 qa=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  [ $qa = 1 ] && sed -i "s/^  deltaT: {/  deltaT: {qAccumulatorFP64: 1, /" $r/solverConfig.yaml
  grep -q "nStepOuter: 30000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml && grep -q "omg_axisym" $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  if [ $qa = 1 ]; then grep -q "qAccumulatorFP64: 1," $r/solverConfig.yaml || { echo "$r: qAccumulatorFP64 の書き換えに失敗" >> $LOG; return 1; }; fi
  diff $r/solverConfig.yaml run_0488_tr_f32/solverConfig.yaml > $r/config_vs_run_0488.diff; echo "$r: run_0488 との設定の差 $(grep -c "^[<>]" $r/config_vs_run_0488.diff) 行" >> $LOG
  if [ $qa = 1 ]; then ( export FORGE_BIN=$F32 FORGE_QACC_AXISYM_DIAG=1 FORGE_DIAG_QACC_ADOPT=$C45/$r/qacc_adopt FORGE_DIAG_QACC_ADOPT_STEPS=1,30000; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$F32; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) qAccumulatorFP64 $qa" >> $LOG
}
arm run_0592_qa_f32 1 || exit 1
arm run_0593_qa_f32b 1 || exit 1
arm run_0594_qa_off 0 || exit 1
sleep 240
for r in run_0592_qa_f32 run_0593_qa_f32b run_0594_qa_off; do echo "$r 起動の確認: $(grep -h "\[qAccumulatorFP64\]" $r/forge_run.log | head -2 | cut -c1-120 | tr '\n' ' ')" >> $LOG; done
wait
for r in run_0592_qa_f32 run_0593_qa_f32b run_0594_qa_off; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 qa_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch qa.done
