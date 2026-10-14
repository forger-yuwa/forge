#!/bin/bash
# plan tooling-nozzle-core-grid §4.26 (事前登録、codex 諮問 2026-10-11): 反復条件の感度 s。Gref (run_0608) と Grefx (run_0610) の res_60000 から
# 同じ格子に restart_field (9 量・ビット一致) し、陰解法の方式だけを point に替えて (lineImplicit 0・lineDtDirectional 0・implicitThermalJacobian 0、
# lineDtDirectionalCap を外す) 各 60,000 step。cfl 4・緩和 0.7・キー 1・格子・境界条件・バイナリは変えない。2 本同時。
# 20,000 step 未満の場は次の出力が出てから消す。ディスクの空きが 3 GB を切ったら自分の forge だけを止める。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_point.sh <Gref の point の run> <Grefx の point の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
LOG=$C45/cg_point.log; MIN_START_GB=12; MIN_RUN_GB=3; OUT=2500; KEEP_FROM=20000
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
A=$1; B=$2; declare -A SRC=([$A]=run_0608_cg_gref_cont_k1 [$B]=run_0610_cg_grefx_ext_k1)
echo "== 開始 $(date -Is) $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
for r in $A $B; do
  s=${SRC[$r]}; seed=$s/res_60000.h5
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  [ -f $seed ] || { echo "$r: 起点 $seed が無い (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $s/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $seed $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う: $(grep VERDICT $r/restart.log)" >> $LOG; exit 1; }
  sed -i 's/deltaT: {lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1, /deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, /' $r/solverConfig.yaml
  d=$(diff $s/solverConfig.yaml $r/solverConfig.yaml | grep '^[<>]' | wc -l)
  grep -q 'deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0,' $r/solverConfig.yaml && [ "$d" = 2 ] \
    && grep -q "nStepOuter: 60000}" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml && grep -q "implicitRelax: 0.7" $r/solverConfig.yaml \
    && diff $s/bcondConfig.yaml $r/bcondConfig.yaml > /dev/null || { echo "$r: point への書き換えが想定と違う (diff の行 $d)" >> $LOG; exit 1; }
  echo "$r: $s の res_60000 から、point に書き換え ($(diff $s/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | cut -c1-110))" >> $LOG
done
PIDS=()
for r in $A $B; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG
done
sleep 150
for r in $A $B; do echo "$r 起動の確認: $(grep -h 'lineImplicit\|implicitThermalJacobian' $r/forge_run.log | head -2 | cut -c1-100 | tr '\n' ' ')" >> $LOG; done
while [ ! -f $A/RUN_RC ] || [ ! -f $B/RUN_RC ]; do
  for r in $A $B; do
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue
      [ $n -gt 0 ] && [ $n -lt $KEEP_FROM ] && [ -f $r/res_$((n + OUT)).h5 ] && rm -f $f $r/res_$n.xmf
    done
  done
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $A $B; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
    break
  fi
  sleep 60
done
wait "${PIDS[@]}"
for r in $A $B; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
