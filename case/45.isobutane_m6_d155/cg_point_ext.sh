#!/bin/bash
# plan tooling-nozzle-core-grid §4.28 (事前登録、codex 諮問 2026-10-11): Gref の point (run_0611) だけを res_60000 から同じ設定で 6 万 step 延ばす。
# restart_field で 9 量・ビット一致を検査。20,000 step 未満の場は次の出力が出てから消す。ディスクの空きが 3 GB を切ったら自分の forge だけを止める。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_point_ext.sh <新しい run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRCRUN=run_0611_cg_gref_point; SEED=$SRCRUN/res_60000.h5; r=$1
LOG=$C45/cg_point_ext.log; MIN_START_GB=12; MIN_RUN_GB=3; OUT=2500; KEEP_FROM=20000
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $r F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
[ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
mkdir $r; for f in $SRCRUN/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う: $(grep VERDICT $r/restart.log)" >> $LOG; exit 1; }
diff $SRCRUN/solverConfig.yaml $r/solverConfig.yaml > /dev/null && diff $SRCRUN/bcondConfig.yaml $r/bcondConfig.yaml > /dev/null || { echo "$r: 設定が $SRCRUN と違う (中止)" >> $LOG; exit 1; }
echo "$r: 設定は $SRCRUN と同一 ($(grep -o 'nStepOuter: [0-9]*' $r/solverConfig.yaml)、$(grep -o 'axisSegmentRWeight: [01]' $r/solverConfig.yaml))" >> $LOG
( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
PID_RC=$!; echo "$r 起動 $(date -Is)" >> $LOG
sleep 150; echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG
while [ ! -f $r/RUN_RC ]; do
  for f in $r/res_*.h5; do
    n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue
    [ $n -gt 0 ] && [ $n -lt $KEEP_FROM ] && [ -f $r/res_$((n + OUT)).h5 ] && rm -f $f $r/res_$n.xmf
  done
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done
    break
  fi
  sleep 60
done
wait $PID_RC
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
