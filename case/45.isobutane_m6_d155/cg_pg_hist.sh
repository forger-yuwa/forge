#!/bin/bash
# plan tooling-nozzle-core-grid §4.92 P-h (事前登録、codex 諮問 2026-10-13 pg-interp): 全補正の後の状態 (800 sweep の run の res_1) を
# 新しい run に restart_field で移し、1 step 回して「再起動直後」の残差 R_B を取る (継続中の評価 R_A は既存の res_2)。
# 設定・環境変数は §4.90 の半分の run と同じ (nStepInner 5・nStepOuter 1・GI は precondEps 1.0、FORGE_DIAG_PRECOND_NO_DT_SCALE=1・FORGE_FREEZE_TURB=1)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_pg_hist.sh <バイナリの sha256 の先頭 8 文字>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64-lmp2dg/solver_density_cuda/build/forge; SHA8=$1
BASE=run_0693_pdt_off; SEED=$BASE/res_20000.h5; MESHSRC=run_0608_cg_gref_cont_k1/nozzle.h5
LOG=$C45/cg_pg_hist.log
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
  cp $MESHSRC $r/nozzle.h5
  python3 $TL/restart_field.py $seed $r/nozzle.h5 --dst-run $r --forge $F > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  sed -i "s/nStepOuter: 20000}/nStepOuter: $nout}/; s/outStepInterval: 1000/outStepInterval: 1/; s/nStepInner: 5$/nStepInner: $nin/" $r/solverConfig.yaml
  [ "$eps" = 1.0 ] && sed -i "s/lowMachPrecond: 3, /precondEps: 1.0, lowMachPrecond: 3, /" $r/solverConfig.yaml
  # 設定の検査: run_0693 との差は nStepOuter・outStepInterval・nStepInner と (GI なら) precondEps の行だけ
  local d=$(diff $BASE/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | grep -vc "nStepOuter: $nout}\|outStepInterval: 1$\|nStepInner: $nin$\|precondEps: 1.0, lowMachPrecond: 3, ")
  [ "$d" = 0 ] && [ "$(grep -c "nStepInner: $nin$" $r/solverConfig.yaml)" = 1 ] && [ "$(grep -c "nStepOuter: $nout}" $r/solverConfig.yaml)" = 1 ] \
    && [ "$(grep -c 'outStepInterval: 1$' $r/solverConfig.yaml)" = 1 ] && [ "$(grep -o precondEps $r/solverConfig.yaml | wc -l)" = "$([ "$eps" = 1.0 ] && echo 1 || echo 0)" ] \
    || { echo "$r: 設定の差が想定と違う (想定外の行 $d) (中止)" >> $LOG; return 1; }
  echo "$r: 準備 eps=$eps nStepInner=$nin nStepOuter=$nout 起点 $seed" >> $LOG
}
runit() {
  local r=$1
  ( export FORGE_BIN=$F FORGE_DIAG_PRECOND_NO_DT_SCALE=1 FORGE_FREEZE_TURB=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 出力 $(ls $r | grep -E '^res_[0-9]+\.h5$' | tr '\n' ' ') スイッチの行 $(cat $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | grep -c 'FORGE_DIAG_PRECOND_NO_DT_SCALE=1:')" >> $LOG
}
for spec in "run_0704_pg_gib run_0698_pg_gi800 1.0" "run_0705_pg_gbb run_0701_pg_gb800 default"; do
  set -- $spec
  prep $1 $3 5 1 $2/res_1.h5 || exit 1
  runit $1
done
echo "== 終了 $(date -Is)" >> $LOG
