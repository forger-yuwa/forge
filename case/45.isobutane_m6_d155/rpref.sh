#!/bin/bash
# plan architecture-float-state-double-geometry §6.44 (事前登録の予定、codex diagnose 2026-10-11): pRef (対流流束の圧力の基準値) の短い試験。インスタンス A、1 本ずつ。
# Cp = FP64 + space.pRef (run_0700〜0702)、Fp = float + space.pRef (run_0703〜0705)。§6.40 と同じ入力・設定・環境変数・バイナリ。各 2,000 step・500 ごと。判定は rpref_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
PREF=${PREF:-5499000.0}
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rpref.log
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
  elif [ $r = run_0700_rp_cp1 ]; then ( export FORGE_BIN=$b FORGE_DIAG_HOOP_CLOSURE=$C45/$r/hoop.h5; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) pRef $p 量子化 $q" >> $LOG
}
one() {   # one "<run> <バイナリ> <pRef> <量子化>"
  set -- $1
  for w in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 4500 ] && break; echo "$(date -Is) $1: available ${av} MB < 4500 MB、待つ" >> $LOG; sleep 60; done
  [ "$av" -ge 4500 ] || { echo "$1: メモリが空かない (中止)" >> $LOG; exit 1; }
  echo "$1: available ${av} MB で投入" >> $LOG
  arm $1 $2 $3 $4 || exit 1
  wait
  echo "$1 終了 $(date -Is) RUN_RC=$(cat $1/RUN_RC) 最後の行: $(grep -h '^step' $1/forge_run.log | tail -1 | cut -c1-100)" >> $LOG
  [ $4 = 1 ] && grep -h "\[quantizeState\]" $1/forge_run.log | grep "step " | cut -c1-300 >> $LOG
  [ $3 = 1 ] && grep -h "pRef" $1/forge_run.log | head -2 | cut -c1-200 >> $LOG
}
# Cp (FP64 + pRef) を先に 3 本。1 本目で閉性の診断 (FORGE_DIAG_HOOP_CLOSURE) も取る
for spec in "run_0700_rp_cp1 $F64 1 0" "run_0701_rp_cp2 $F64 1 0" "run_0702_rp_cp3 $F64 1 0"; do one "$spec" || exit 1; done
python3 cp_gate.py >> $LOG 2>&1 || { echo "前段の関門を通らない (前提不成立、Fp は回さない)" >> $LOG; python3 rpref_an.py >> $LOG 2>&1
python3 rq_dT.py >> $LOG 2>&1; touch rpref.done; exit 0; }
for spec in "run_0703_rp_fp1 $F32 1 0" "run_0704_rp_fp2 $F32 1 0" "run_0705_rp_fp3 $F32 1 0"; do one "$spec" || exit 1; done
python3 rpref_an.py >> $LOG 2>&1
python3 rq_dT.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rpref.done
