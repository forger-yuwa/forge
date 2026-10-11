#!/bin/bash
# plan tooling-nozzle-core-grid §4.46 M-r (事前登録、codex 諮問 2026-10-11 response): 停滞した B0 の場から 2 外反復、nStepInner 42、SST の保存量の更新を凍結
# (FORGE_FREEZE_TURB=1、両腕共通)。最初の外反復を FORGE_LINE_DUMP で書き出し (cg_sweep.sh と同じライン)、res_1 に R(Q₀)、res_2 に R(Q₁) を出す。
# 腕: gref50 = Gref・B0 (上限 50) / gref1 = 同じ場で上限 1。起点は run_0619 の res_60000。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_resp.sh <run> <arm>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
r=$1; arm=$2; LOG=$C45/cg_resp.log; NIN=42
case $arm in
  gref50)  SRC=run_0619_cg_gref_b0_ctl; CAP=50.0 ;;
  gref1)   SRC=run_0619_cg_gref_b0_ctl; CAP=1.0 ;;
  *) echo "腕が不明: $arm"; exit 1 ;;
esac
SEED=$SRC/res_60000.h5
echo "== 開始 $(date -Is) $r $arm F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $SEED $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1G . | tail -1 | tr -dc 0-9)" -ge 6 ] || { echo "ディスクの空きが 6 GB 未満 (中止)" >> $LOG; exit 1; }
[ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
mkdir $r; for f in $SRC/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|.done_*|*.png) ;; *) cp -a $f $r/;; esac; done
python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
sed -i "s/nStepOuter: [0-9]*}/nStepOuter: 2}/; s/nStepInner: 5$/nStepInner: $NIN/; s/lineDtDirectionalCap: 50.0,/lineDtDirectionalCap: $CAP,/; s/outStepInterval: [0-9]*/outStepInterval: 1/" $r/solverConfig.yaml
# 元の run の設定との差は nStepOuter・nStepInner・outStepInterval と (gref1 なら) 上限の行だけであること (凍結は環境変数)
d=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
want=3; [ $arm = gref1 ] && want=4
[ "$d" = "$want" ] && grep -q "nStepOuter: 2}" $r/solverConfig.yaml && grep -q "nStepInner: $NIN$" $r/solverConfig.yaml && grep -q "lineDtDirectionalCap: $CAP," $r/solverConfig.yaml \
  && grep -q "implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $r/solverConfig.yaml && diff $SRC/bcondConfig.yaml $r/bcondConfig.yaml > /dev/null \
  || { echo "$r: 元の run との設定の差が想定と違う ($d 行) (中止)" >> $LOG; exit 1; }
NODES=$(python3 cg_sweep_nodes.py $r) || { echo "$r: 書き出す節点を選べない (中止)" >> $LOG; exit 1; }
echo "$r: $SRC と nStepOuter 2・nStepInner $NIN・outStepInterval 1・上限 $CAP 以外同一、FORGE_FREEZE_TURB=1。書き出す節点 $(echo $NODES | tr ',' '\n' | wc -l) 本ぶん" >> $LOG
mkdir $r/linedump
( export FORGE_BIN=$F7_64 FORGE_FREEZE_TURB=1 FORGE_LINE_DUMP_DIR=$C45/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
  bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC)、書き出し $(grep -c . $r/linedump/meta.txt) 行、$(du -sh $r/linedump | cut -f1)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
