#!/bin/bash
# plan architecture-float-state-double-geometry §6.61 (事前登録、codex diagnose `2026-10-13-testg-result-diagnose.md`): 速さの比較 S の 1 本の腕。
# 上限 50 (生産の B0) に固定し、格子だけを変える: g0c = _g0c/q32_init.h5 (G0c)、g0 = _tr/q32_init.h5 (G0)。FP64 (~/forge-fgeom7-fp64)、停止規則 §6.11 (m9_watch.py)。
# 設定は試験 G (gg_arm.sh) と同じ: run_0717_c50_d1 の solverConfig.yaml から nStepOuter 200000・outStepInterval 2500 の 2 行だけ。ほかの yaml は run_0717_c50_d1 と _g0c で同一 (2026-10-13 に確認)。
# GPU の占有: 起動の前にほかの forge が無くなるまで待ち (最大 1 時間)、走行中は 30 秒ごとに時刻・step・forge の数を timing.tsv に書く (判定器が占有を確かめる)。
# 使い方: bash sp_arm.sh <新しい run> <g0|g0c>
set -uo pipefail
R=$1; G=$2
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/sp.log; TL=~/forge-wallfit/solver_density_cuda/tools
BIN=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45
die() { echo "$R: 中止: $* ($(date -Is))" >> $LOG; exit 1; }
[ "$(sha256sum $BIN | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || die "FP64 のバイナリが §6.51 と違う"
case $G in
  g0c) INP=_g0c/q32_init.h5; WANT=6f8fcda6d97482e3; YD=_g0c;;
  g0)  INP=_tr/q32_init.h5;  WANT=ba859e02f3ec236a; YD=run_0717_c50_d1;;
  *) die "格子は g0 か g0c";;
esac
[ "$(sha256sum $INP | cut -c1-16)" = "$WANT" ] || die "入力 $INP が登録と違う"
for f in bcondConfig.yaml probe.yaml species_meta.yaml resolved_species_fc8ba1a8828147e3.yaml resolved_species_fc8ba1a8828147e3_70cc58aa292c9504.yaml; do
  cmp -s run_0717_c50_d1/$f _g0c/$f || die "$f が G0 と G0c で違う"
done
mkdir $R || die "$R が既にある"
cp $YD/*.yaml $R/ && cp run_0717_c50_d1/solverConfig.yaml $R/solverConfig.yaml && cp $INP $R/nozzle.h5 || die "複製"
sed -i 's/nStepOuter: 30000}/nStepOuter: 200000}/; s/outStepInterval: 500$/outStepInterval: 2500/' $R/solverConfig.yaml
DL=$(diff run_0717_c50_d1/solverConfig.yaml $R/solverConfig.yaml | grep '^[<>]')
[ "$(echo "$DL" | wc -l)" = 4 ] && [ -z "$(echo "$DL" | grep -v -e 'nStepOuter: ' -e 'outStepInterval: ')" ] \
  && grep -q "nStepOuter: 200000}" $R/solverConfig.yaml && grep -q "outStepInterval: 2500$" $R/solverConfig.yaml && grep -q "lineDtDirectionalCap: 50.0, " $R/solverConfig.yaml || die "設定の書き換え"
for i in $(seq 1 120); do [ -z "$(pgrep -x forge)" ] && break; sleep 30; done
[ -z "$(pgrep -x forge)" ] || die "ほかの forge が 1 時間空かない (占有できない)"
echo "$R 起動 $(date -Is) 格子 $G バイナリ $(sha256sum $BIN | cut -c1-16) 入力 $(sha256sum $R/nozzle.h5 | cut -c1-16)" >> $LOG
date +%s.%N > $R/START_EPOCH
( export FORGE_BIN=$BIN; bash $TL/run_case.sh $C45/$R > $R/run_case_stdout.log 2>&1; echo $? > $R/RUN_RC ) &
( while [ ! -e $R/RUN_RC ]; do
    s=$(tail -1 $R/residual_history.csv 2>/dev/null | cut -d, -f1); n=$(pgrep -x forge | wc -l)
    o=$(for p in $(pgrep -x forge); do readlink /proc/$p/cwd; done | grep -vx "$C45/$R" | tr '\n' ' ')
    echo -e "$(date +%s.%N)\t${s:-}\t$n\t$o" >> $R/timing.tsv; sleep 30
  done ) &
python3 m9_watch.py $R --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $R/m9_watch.log 2>&1 &
wait
echo "$R 終了 $(date -Is) RUN_RC=$(cat $R/RUN_RC 2>/dev/null) 見張り $(python3 -c "import json; d=json.load(open('$R/m9_watch.json')); print(d['status'], d.get('reach_step'))" 2>&1)" >> $LOG
