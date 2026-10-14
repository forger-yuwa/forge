#!/bin/bash
# plan architecture-float-state-double-geometry §6.42 (事前登録、codex diagnose 2026-10-11): FP64 のビルドで保存量 9 つを全部 step 末尾に float32 に丸める短い試験 (G9)。インスタンス A、1 本ずつ。
# G9 = run_0706〜0708。§6.40 の C と同じ FP64 のバイナリ・入力・設定・環境変数 + FORGE_DIAG_QUANTIZE_STATE=1 (ONLY なし)。各 2,000 step・500 ごと。判定は rg_an.py (C・F は既存の run)。
# arm() の pRef の引数は §6.42 では使わない (pRef の試験は別に登録する。run_0700〜0705 はそのために取っておく)。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
PREF=${PREF:-5499000.0}
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rp.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F64=~/forge-fgeom11-fp64/$B; F32=~/forge-fgeom7-f32/$B
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) F32 $(sha256sum $F32 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16) pRef $PREF" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 16106127360 ] || { echo "ディスクの空きが 15 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || { echo "float のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F64 | cut -c1-16)" = "74954e7be5ef38a8" ] || { echo "FP64 のバイナリが §6.40 と違う (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ> <pRef 0|1> <量子化 0|1>
  local r=$1 b=$2 p=$3 q=$4
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 2000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  [ $p = 1 ] && sed -i "s/^space: {/space: {pRef: $PREF, /" $r/solverConfig.yaml
  grep -q "nStepOuter: 2000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  if [ $p = 1 ]; then grep -q "^space: {pRef: $PREF," $r/solverConfig.yaml || { echo "$r: pRef の書き換えに失敗" >> $LOG; return 1; }; fi
  diff $r/solverConfig.yaml run_0595_rq_c1/solverConfig.yaml > $r/config_vs_run_0595.diff; echo "$r: run_0595 との設定の差 $(grep -c "^[<>]" $r/config_vs_run_0595.diff) 行" >> $LOG
  if [ $q = 1 ]; then ( export FORGE_BIN=$b FORGE_DIAG_QUANTIZE_STATE=1 FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2000; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) pRef $p 量子化 $q" >> $LOG
}
for spec in "run_0706_rp_g9a $F64 0 1" "run_0707_rp_g9b $F64 0 1" "run_0708_rp_g9c $F64 0 1"; do
  set -- $spec
  for w in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 4500 ] && break; echo "$(date -Is) $1: available ${av} MB < 4500 MB、待つ" >> $LOG; sleep 60; done
  [ "$av" -ge 4500 ] || { echo "$1: メモリが空かない (中止)" >> $LOG; exit 1; }
  echo "$1: available ${av} MB で投入" >> $LOG
  arm $1 $2 $3 $4 || exit 1
  wait
  echo "$1 終了 $(date -Is) RUN_RC=$(cat $1/RUN_RC) 最後の行: $(grep -h '^step' $1/forge_run.log | tail -1 | cut -c1-100)" >> $LOG
  [ $4 = 1 ] && grep -h "\[quantizeState\]" $1/forge_run.log | grep "step " | cut -c1-300 >> $LOG
  [ $3 = 1 ] && grep -h "pRef" $1/forge_run.log | head -2 | cut -c1-200 >> $LOG
done
python3 rg_an.py >> $LOG 2>&1
python3 rq_dT.py >> $LOG 2>&1
python3 qs_record.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rp.done
