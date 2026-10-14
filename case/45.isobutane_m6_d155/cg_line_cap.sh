#!/bin/bash
# plan tooling-nozzle-core-grid §4.32 (事前登録、codex 諮問 2026-10-11、ユーザ依頼「line で同じ結果を出すには」): Gref の B0 の場 (run_0608 の res_60000) から、
# ラインは有効のまま lineDtDirectionalCap だけを A = 50 (B0 のまま)・B = 1 にして各 20,000 step。1000 step ごとに出力し、dt_local も書く。
# 出力ごとに登録の 8 量の系列を取り出してから場を消す (res_1000 と最後の場は残す)。ディスクの空きが 3 GB を切ったら自分の forge だけを止める。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_line_cap.sh <A の run> <B の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRCRUN=run_0608_cg_gref_cont_k1; SEED=$SRCRUN/res_60000.h5; A=$1; B=$2
LOG=$C45/cg_line_cap.log; MIN_START_GB=8; MIN_RUN_GB=3; STEPS=20000; OUT=1000
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
for r in $A $B; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $SRCRUN/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う: $(grep VERDICT $r/restart.log)" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/outStepInterval: 2500/outStepInterval: $OUT/; s/extraFields: \[res_ro\]}/extraFields: [res_ro, dt_local]}/" $r/solverConfig.yaml
done
sed -i 's/deltaT: {lineDtDirectionalCap: 50.0, /deltaT: {lineDtDirectionalCap: 1.0, /' $B/solverConfig.yaml
# 書き換えの検査: A は 3 行、B は 4 行 (上限の行を含む) だけが元と違う。A と B の違いは上限の 1 行だけ
na=$(diff $SRCRUN/solverConfig.yaml $A/solverConfig.yaml | grep -c '^>'); nb=$(diff $SRCRUN/solverConfig.yaml $B/solverConfig.yaml | grep -c '^>'); nab=$(diff $A/solverConfig.yaml $B/solverConfig.yaml | grep -c '^>')
grep -q "lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $A/solverConfig.yaml && grep -q "lineDtDirectionalCap: 1.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $B/solverConfig.yaml \
  && [ "$na" = 3 ] && [ "$nb" = 4 ] && [ "$nab" = 1 ] && grep -q "dt_local" $B/solverConfig.yaml && grep -q "outStepInterval: $OUT" $B/solverConfig.yaml \
  || { echo "書き換えが想定と違う (A $na 行・B $nb 行・A と B $nab 行) (中止)" >> $LOG; exit 1; }
diff $A/solverConfig.yaml $B/solverConfig.yaml | grep '^[<>]' | cut -c1-90 | sed 's/^/  A と B の差: /' >> $LOG
PIDS=()
for r in $A $B; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG
done
sleep 120
for r in $A $B; do echo "$r 起動の確認: $(grep -h 'lineDtDirectionalCap' $r/forge_run.log | head -1 | cut -c1-100)" >> $LOG; done
extract() {   # extract <run> <step>: 1 出力の 8 量を取り出す (gref_series_<run>_<step>.csv)
  python3 cg_gref_judge.py series "$1_$2" $1 --from $2 --to $2 > /dev/null 2>&1
}
while true; do
  done_all=1
  for r in $A $B; do
    [ -f $r/RUN_RC ] || done_all=0
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && continue
      # 書き終わったことを、次の出力があるか run の終了で確かめる
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        extract $r $n && echo "$r: $n を取り出した" >> $LOG
        [ $n -ne $OUT ] && [ $n -ne $STEPS ] && [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && rm -f $f $r/res_$n.xmf
      fi
    done
  done
  [ $done_all = 1 ] && { left=0; for r in $A $B; do for f in $r/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && left=1; done; done; [ $left = 0 ] && break; }
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $A $B; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
  fi
  sleep 30
done
wait "${PIDS[@]}"
for r in $A $B; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
