#!/bin/bash
# plan tooling-nozzle-core-grid §4.34 の D1 (事前登録、codex 諮問 2026-10-11 の 2 回目): Gref の B0 の場 (run_0608 の res_60000) から、
# 上限 50 と 1 の 2 本を 1 本ずつ (ディスクのため)、各 2032 step。2000〜2032 の場を毎 step 書き出す (level 0 + 残差・dq_block_old・SST の対角の材料)。
# FORGE_DIAG_COMMIT_LOSS=1 で commit の丸め・床の切り詰めを毎 step 数える。1 本終わるごとに cg_b0_diag_an.py で要約し、途中の場を消す。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_b0_diag.sh <上限 50 の run> <上限 1 の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRCRUN=run_0608_cg_gref_cont_k1; SEED=$SRCRUN/res_60000.h5; A=$1; B=$2
LOG=$C45/cg_b0_diag.log
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, dt_local, dq_block_old_0, dq_block_old_1, dq_block_old_2, dq_block_old_4, src_jac_k, src_jac_omega, transport_diag_k, transport_diag_omega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
for pair in "$A 50.0" "$B 1.0"; do
  set -- $pair; r=$1; cap=$2
  [ "$(free_gb)" -ge 8 ] || { echo "$r: ディスクの空きが 8 GB 未満 (中止)" >> $LOG; exit 1; }
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $SRCRUN/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: 2032}/; s/outStepStart: 0/outStepStart: 2000/; s/outStepInterval: 2500/outStepInterval: 1/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 0, extraFields: $EXTRA}/; s/deltaT: {lineDtDirectionalCap: 50.0, /deltaT: {lineDtDirectionalCap: $cap, /" $r/solverConfig.yaml
  grep -q "nStepOuter: 2032}" $r/solverConfig.yaml && grep -q "outStepStart: 2000" $r/solverConfig.yaml && grep -q "outStepInterval: 1$" $r/solverConfig.yaml \
    && grep -q "level: 0, extraFields: \[res_ro, res_roUx" $r/solverConfig.yaml && grep -q "lineDtDirectionalCap: $cap, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $r/solverConfig.yaml \
    || { echo "$r: 書き換えが想定と違う (中止)" >> $LOG; exit 1; }
  echo "$r: 上限 $cap、$(diff $SRCRUN/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>') 行を書き換え" >> $LOG
  ( export FORGE_BIN=$F7_64 FORGE_DIAG_COMMIT_LOSS=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 出力 $(ls $r/res_*.h5 | grep -c 'res_[0-9]*\.h5') 個 commit_loss $(wc -l < $r/commit_loss.csv 2>/dev/null) 行" >> $LOG
  grep -h "lineDtDirectionalCap\|commitLoss\] FORGE" $r/forge_run.log | head -2 | cut -c1-120 | sed "s/^/  $r: /" >> $LOG
  python3 cg_b0_diag_an.py $r >> $LOG 2>&1 || { echo "$r: 解析に失敗 (場は消さない)" >> $LOG; exit 1; }
  for k in $(seq 2001 2031); do rm -f $r/res_$k.h5 $r/res_$k.xmf; done
done
python3 cg_b0_diag_cmp.py $A $B >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
