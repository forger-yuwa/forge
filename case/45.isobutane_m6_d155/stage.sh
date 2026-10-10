#!/bin/bash
# plan architecture-float-state-double-geometry §6.30 (事前登録、codex diagnose 2026-10-11): 共通の Q32 の状態からの最初の 1 回の更新を段ごとに測る。インスタンス B。
# バイナリ ~/forge-fgeom8-{fp64,f32} (commit ba9e889a、FORGE_DIAG_STAGE_DUMP あり)。共通の入力 _tr/q32_init.h5 (§6.26 と同じ)、FORGE_FREEZE_TURB=1。
# 1. 診断の有無の確認 (_stage/): 各精度で、診断なし 2 本と診断あり 1 本 (本番の 1 本目) の 1 step の res_1 を比べる (stage_an.py)。
# 2. 本番: 各精度 3 回の 1 step (毎回同じ入力から)。FP64 = run_0494〜0496、float = run_0497〜0499。各 run の stage.h5 に記録。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_FREEZE_TURB=1
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/stage.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F64=~/forge-fgeom8-fp64/$B; F32=~/forge-fgeom8-f32/$B
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) F32 $(sha256sum $F32 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
one() {   # one <dir> <バイナリ> <診断 0|1>
  local r=$1 b=$2 dg=$3
  mkdir -p $(dirname $r); mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 1}/; s/outStepInterval: 2500/outStepInterval: 1/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe]}/' $r/solverConfig.yaml
  grep -q "nStepOuter: 1}" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  if [ $dg = 1 ]; then ( export FORGE_BIN=$b FORGE_DIAG_STAGE_DUMP=$C45/$r/stage.h5 FORGE_DIAG_STAGE_STEP=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1 )
  else ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1 ); fi
  echo "$r 診断 $dg rc=$? $(grep -h '\[stageDump' $r/forge_run.log 2>/dev/null | tail -1 | cut -c1-200)" >> $LOG
}
for r in run_0494_st_f64a run_0495_st_f64b run_0496_st_f64c; do one $r $F64 1 || exit 1; done
for r in run_0497_st_f32a run_0498_st_f32b run_0499_st_f32c; do one $r $F32 1 || exit 1; done
for r in _stage/off_f64_1 _stage/off_f64_2; do one $r $F64 0 || exit 1; done
for r in _stage/off_f32_1 _stage/off_f32_2; do one $r $F32 0 || exit 1; done
python3 stage_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch stage.done
