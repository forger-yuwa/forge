#!/bin/bash
# plan tooling-nozzle-core-grid §4.94 P-tf (事前登録、ユーザ決定・codex 諮問 2026-10-13 thermofloat): 熱物性の評価の経路 (physProp.thermoFloat) を 1 → 0 にして、
# 既存の 5 状態 (Q0・GI/GB の Q0 + δ★/2・GI/GB の Q0 + δ★) を各 1 step 再評価する。δ★ は計算し直さない。
#   Q0 = run_0693 の res_20000 を restart、半補正 = run_0702・run_0703 の書き込み済みの入力 h5 をそのまま複製、全補正 = run_0698・run_0701 の res_1 を restart。
#   他は §4.92 の B の run と同じ (nStepInner 5・nStepOuter 1・GI は precondEps 1.0、FORGE_DIAG_PRECOND_NO_DT_SCALE=1・FORGE_FREEZE_TURB=1)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_tf.sh <バイナリの sha256 の先頭 8 文字>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64-lmp2dg/solver_density_cuda/build/forge; SHA8=$1
BASE=run_0693_pdt_off; SEED=$BASE/res_20000.h5; MESHSRC=run_0608_cg_gref_cont_k1/nozzle.h5
LOG=$C45/cg_tf.log
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F $(sha256sum $F | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(sha256sum $F | cut -c1-8)" = "$SHA8" ] || { echo "バイナリの sha256 が $SHA8 でない (中止)" >> $LOG; exit 1; }
[ "$(free_gb)" -ge 5 ] || { echo "ディスクの空きが 5 GB 未満 (中止)" >> $LOG; exit 1; }
# prep <run> <eps: 1.0 または default> <nStepInner> <nStepOuter>
prep() {
  local r=$1 eps=$2 nin=$3 nout=$4 seed=${5:-$SEED}
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r
  for f in $BASE/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|nozzle.h5|restart.log|run_case_stdout.log|*.xmf|*.png) ;; *) cp -a $f $r/;; esac; done
  case "$seed" in
    copy:*) cp ${seed#copy:} $r/nozzle.h5; echo "入力 h5 を複製: ${seed#copy:} (sha256 $(sha256sum ${seed#copy:} | cut -c1-16))" > $r/restart.log ;;
    *) cp $MESHSRC $r/nozzle.h5
       python3 $TL/restart_field.py $seed $r/nozzle.h5 --dst-run $r --forge $F > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
       grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; } ;;
  esac
  sed -i "s/nStepOuter: 20000}/nStepOuter: $nout}/; s/outStepInterval: 1000/outStepInterval: 1/; s/nStepInner: 5$/nStepInner: $nin/" $r/solverConfig.yaml
  [ "$eps" = 1.0 ] && sed -i "s/lowMachPrecond: 3, /precondEps: 1.0, lowMachPrecond: 3, /" $r/solverConfig.yaml
  sed -i "s/physProp: {thermalMethod: 2, /physProp: {thermoFloat: 0, thermalMethod: 2, /" $r/solverConfig.yaml
  # 設定の検査: run_0693 との差は nStepOuter・outStepInterval・nStepInner と (GI なら) precondEps の行だけ
  local d=$(diff $BASE/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | grep -vc "nStepOuter: $nout}\|outStepInterval: 1$\|nStepInner: $nin$\|precondEps: 1.0, lowMachPrecond: 3, \|physProp: {thermoFloat: 0, thermalMethod: 2, ")
  [ "$d" = 0 ] && [ "$(grep -c "nStepInner: $nin$" $r/solverConfig.yaml)" = 1 ] && [ "$(grep -c "nStepOuter: $nout}" $r/solverConfig.yaml)" = 1 ] \
    && [ "$(grep -c 'outStepInterval: 1$' $r/solverConfig.yaml)" = 1 ] && [ "$(grep -o thermoFloat $r/solverConfig.yaml | wc -l)" = 1 ] && [ "$(grep -o precondEps $r/solverConfig.yaml | wc -l)" = "$([ "$eps" = 1.0 ] && echo 1 || echo 0)" ] \
    || { echo "$r: 設定の差が想定と違う (想定外の行 $d) (中止)" >> $LOG; return 1; }
  echo "$r: 準備 eps=$eps nStepInner=$nin nStepOuter=$nout 起点 $seed" >> $LOG
}
runit() {
  local r=$1
  ( export FORGE_BIN=$F FORGE_DIAG_PRECOND_NO_DT_SCALE=1 FORGE_FREEZE_TURB=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 出力 $(ls $r | grep -E '^res_[0-9]+\.h5$' | tr '\n' ' ') スイッチの行 $(cat $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | grep -c 'FORGE_DIAG_PRECOND_NO_DT_SCALE=1:')" >> $LOG
}
for spec in "run_0900_tf_q0 1.0 $SEED" "run_0901_tf_gih 1.0 copy:run_0702_pg_gih/nozzle.h5" "run_0902_tf_gbh default copy:run_0703_pg_gbh/nozzle.h5" \
            "run_0903_tf_gib 1.0 run_0698_pg_gi800/res_1.h5" "run_0904_tf_gbb default run_0701_pg_gb800/res_1.h5"; do
  set -- $spec
  prep $1 $2 5 1 $3 || exit 1
  runit $1
  echo "$1 熱物性の経路: $(grep -h "'thermoFloat'" $1/forge_run.log | head -1)" >> $LOG
done
echo "== 終了 $(date -Is)" >> $LOG
