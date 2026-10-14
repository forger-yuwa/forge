#!/bin/bash
# plan tooling-nozzle-core-grid §4.48 M-o (事前登録、codex 諮問 2026-10-11 mismatch): 応答のずれ e のうち、右辺の再構成次数の寄与を測る。
#  base: 診断入りの FP64 バイナリ (FORGE_LINE_DUMP_FINAL) で M-r の上限 50 を回し直す (2 外反復・nStepInner 42・SST 凍結)。
#        確定する補正 δ に対する右辺 b(δ) = R₀ − Nδ を直接書き出す。Q₀ = res_0、Q₁ = res_1 の保存量、R₀ = res_1・R₁ = res_2 の残差。
#  eval: 状態 (Q₀ または base の Q₁) から 1 外反復だけ回し、res_1 の残差 R(Q) を取る。convMethod (1 = 現行の 2 次、0 = 1 次) だけ変える。
#        SST 凍結は base と同じ。nStepInner は 5 のまま (res_1 の残差は更新前の状態の値なので反復の設定に依らない)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155):
#   bash cg_mis.sh base <run>
#   bash cg_mis.sh eval <run> <seed の h5> <convMethod>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
FD=$HOME/forge-fgeom7-fp64-diag/solver_density_cuda/build/forge
SRC=run_0619_cg_gref_b0_ctl; LOG=$C45/cg_mis.log
mode=$1; r=$2
case $mode in
  base) SEED=$SRC/res_60000.h5; CONV=1 ;;
  eval) SEED=$3; CONV=$4 ;;
  *) echo "mode が不明: $mode"; exit 1 ;;
esac
echo "== 開始 $(date -Is) $mode $r FD $(sha256sum $FD | cut -c1-16) 起点 $SEED $(sha256sum $SEED | cut -c1-16) convMethod $CONV" >> $LOG
[ "$(df --output=avail -B1G . | tail -1 | tr -dc 0-9)" -ge 6 ] || { echo "ディスクの空きが 6 GB 未満 (中止)" >> $LOG; exit 1; }
[ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
mkdir $r; for f in $SRC/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|.done_*|*.png) ;; *) cp -a $f $r/;; esac; done
python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $FD > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
if [ $mode = base ]; then
  sed -i "s/nStepOuter: [0-9]*}/nStepOuter: 2}/; s/nStepInner: 5$/nStepInner: 42/; s/outStepInterval: [0-9]*/outStepInterval: 1/" $r/solverConfig.yaml
  want=3
else
  sed -i "s/nStepOuter: [0-9]*}/nStepOuter: 1}/; s/outStepInterval: [0-9]*/outStepInterval: 1/; s/convMethod: 1,/convMethod: $CONV,/" $r/solverConfig.yaml
  want=2; [ $CONV != 1 ] && want=3
fi
d=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
[ "$d" = "$want" ] && grep -q "convMethod: $CONV," $r/solverConfig.yaml && grep -q "lineDtDirectionalCap: 50.0," $r/solverConfig.yaml && diff $SRC/bcondConfig.yaml $r/bcondConfig.yaml > /dev/null \
  || { echo "$r: 元の run との設定の差が想定と違う ($d 行) (中止)" >> $LOG; exit 1; }
if [ $mode = base ]; then
  NODES=$(python3 cg_sweep_nodes.py $r) || { echo "$r: 書き出す節点を選べない (中止)" >> $LOG; exit 1; }
  mkdir $r/linedump
  echo "$r: $SRC と nStepOuter 2・nStepInner 42・outStepInterval 1 以外同一、FORGE_FREEZE_TURB=1・FORGE_LINE_DUMP_FINAL=1" >> $LOG
  ( export FORGE_BIN=$FD FORGE_FREEZE_TURB=1 FORGE_LINE_DUMP_DIR=$C45/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES FORGE_LINE_DUMP_FINAL=1
    bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
else
  echo "$r: $SRC と nStepOuter 1・outStepInterval 1・convMethod $CONV 以外同一、FORGE_FREEZE_TURB=1" >> $LOG
  ( export FORGE_BIN=$FD FORGE_FREEZE_TURB=1
    bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
fi
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
