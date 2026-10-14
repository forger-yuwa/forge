#!/bin/bash
# plan tooling-nozzle-core-grid §4.56 M-t (事前登録、codex 諮問 2026-10-12 linefaces): ライン面の熱伝導の近傍作用を固定 δ の予測に足す。
#  base: 診断入りの FP64 バイナリ diag5 (§4.56 の熱伝導の近傍作用を足した版) で M-o の base を再現し (2 外反復・nStepInner 42・SST 凍結)、FORGE_LINE_DUMP_FINAL_SPLIT で
#        左辺の対流の面の寄与を外した組み直し (ライン外・ライン面・両方・確認) も書き出す。
#  eval: 状態 (Q₀ または base の Q₁) から 1 外反復、convMethod 0 (1 次) で、FORGE_DIAG_CONV_FACES = 第 4 引数 (0 全面・1 ライン外だけ・2 ライン面だけ・3 なし)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155):
#   bash cg_conv.sh base <run>
#   bash cg_conv.sh eval <run> <seed の h5> <面の種類 0〜3>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
FD=$HOME/forge-fgeom7-fp64-diag5/solver_density_cuda/build/forge
SRC=run_0619_cg_gref_b0_ctl; LOG=$C45/cg_therm.log
mode=$1; r=$2
case $mode in
  base) SEED=$SRC/res_60000.h5; CONV=1 ;;
  eval) SEED=$3; CONV=0; FACES=$4 ;;
  *) echo "mode が不明: $mode"; exit 1 ;;
esac
echo "== 開始 $(date -Is) $mode $r FD $(sha256sum $FD | cut -c1-16) 起点 $SEED $(sha256sum $SEED | cut -c1-16) convMethod $CONV 面 ${FACES:-}" >> $LOG
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
  echo "$r: $SRC と nStepOuter 2・nStepInner 42・outStepInterval 1 以外同一、FORGE_FREEZE_TURB=1・FORGE_LINE_DUMP_FINAL=1・FORGE_LINE_DUMP_FINAL_SPLIT=1" >> $LOG
  ( export FORGE_BIN=$FD FORGE_FREEZE_TURB=1 FORGE_LINE_DUMP_DIR=$C45/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES FORGE_LINE_DUMP_FINAL=1 FORGE_LINE_DUMP_FINAL_SPLIT=1
    bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
else
  echo "$r: $SRC と nStepOuter 1・outStepInterval 1・convMethod $CONV 以外同一、FORGE_FREEZE_TURB=1・FORGE_DIAG_CONV_FACES=$FACES" >> $LOG
  ( export FORGE_BIN=$FD FORGE_FREEZE_TURB=1 FORGE_DIAG_CONV_FACES=$FACES
    bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
fi
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
