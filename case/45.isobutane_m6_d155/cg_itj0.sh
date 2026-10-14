#!/bin/bash
# plan tooling-nozzle-core-grid §4.40 M-a (事前登録、codex 諮問 2026-10-11): Gref の B0 の場 (run_0608 の res_60000) から、B0 のまま
# implicitThermalJacobian だけを 5 → 0 にして 20,000 step。対照は run_0619 (同じ起点・同じ構成) の 0〜20k。
# 2500 step ごとの出力から縮流部の残差 (cg_region_res.py) と dt・μt の変化 (cg_freeze_an.py) を取り出し、前の出力の場を消す (res_0・最後は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_itj0.sh <run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; r=$1; CTL=run_0619_cg_gref_b0_ctl
LOG=$C45/cg_itj0.log; OUT=2500; STEPS=20000
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
echo "== 開始 $(date -Is) $r F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1G . | tail -1 | tr -dc 0-9)" -ge 8 ] || { echo "ディスクの空きが 8 GB 未満 (中止)" >> $LOG; exit 1; }
[ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/; s/implicitThermalJacobian: 5, /implicitThermalJacobian: 0, /" $r/solverConfig.yaml
# 対照 (run_0619) の設定との差は nStepOuter・extraFields (dt_local を足した)・ITJ の 3 行だけであること
d=$(diff $CTL/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
[ "$d" = 3 ] && grep -q "implicitThermalJacobian: 0, lineImplicit: 1, lineDtDirectional: 1," $r/solverConfig.yaml && grep -q "nStepOuter: $STEPS}" $r/solverConfig.yaml && grep -q "volume, dt_local\]" $r/solverConfig.yaml && diff $CTL/bcondConfig.yaml $r/bcondConfig.yaml > /dev/null \
  || { echo "$r: 対照との設定の差が想定と違う ($d 行) (中止)" >> $LOG; exit 1; }
echo "$r: 対照 $CTL と nStepOuter・extraFields・ITJ (5 → 0) 以外同一" >> $LOG
( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
PID=$!; echo "$r 起動 $(date -Is)" >> $LOG
sleep 120; echo "$r 起動の確認: $(grep -h "'implicitThermalJacobian'" $r/forge_run.log | head -1 | cut -c1-80)" >> $LOG
prev=""
while true; do
  for k in $(seq $OUT $OUT $STEPS); do
    [ -f $r/res_$k.h5 ] || continue
    [ -f $r/.done_$k ] && continue
    if [ -f $r/res_$((k + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
      python3 cg_region_res.py $r $k
      pk=$((k - OUT)); [ $pk -lt $OUT ] && pk=0
      [ -f $r/res_$pk.h5 ] && python3 cg_freeze_an.py $r $pk $k
      touch $r/.done_$k; echo "$r: $k を取り出した" >> $LOG
      [ $pk -gt 0 ] && rm -f $r/res_$pk.h5 $r/res_$pk.xmf
    fi
  done
  [ -f $r/RUN_RC ] && [ -f $r/.done_$STEPS ] && break
  sleep 30
done
wait $PID
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
