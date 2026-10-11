#!/bin/bash
# plan architecture-float-state-double-geometry §6.32 (事前登録、codex diagnose 2026-10-11): block の系の組み立ての精度の感度試験。
# §6.26 の float の腕 (run_0488・0489) と同じバイナリ (~/forge-fgeom7-f32)・入力 (_tr/q32_init.h5)・環境変数・設定で、time.deltaT.implicitSolvePrecision だけを 0 → 1 にした C = run_0590_ds_f32・C' = run_0591_ds_f32b。
# 各 30,000 step・500 ごと。インスタンス B。判定は ds_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/ds.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
B=solver_density_cuda/build/forge; F64=~/forge-fgeom7-fp64/$B; F32=~/forge-fgeom7-f32/$B
SRC=~/tr_in/res_145000.h5
cd $C45
echo "== 開始 $(date -Is) F32 $(sha256sum $F32) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 32212254720 ] || { echo "ディスクの空きが 30 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  sed -i "s/^  deltaT: {/  deltaT: {implicitSolvePrecision: 1, /" $r/solverConfig.yaml
  grep -q "nStepOuter: 30000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml && grep -q "omg_axisym" $r/solverConfig.yaml && grep -q "implicitSolvePrecision: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  diff $r/solverConfig.yaml run_0488_tr_f32/solverConfig.yaml > $r/config_vs_run_0488.diff; echo "$r: run_0488 との設定の差 $(grep -c "^[<>]" $r/config_vs_run_0488.diff) 行" >> $LOG
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG
}
arm run_0590_ds_f32 $F32 || exit 1
arm run_0591_ds_f32b $F32 || exit 1
sleep 240
for r in run_0590_ds_f32 run_0591_ds_f32b; do echo "$r 起動の確認: $(grep -h "implicitSolvePrecision" $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG; done
wait
for r in run_0590_ds_f32 run_0591_ds_f32b; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 ds_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch ds.done
