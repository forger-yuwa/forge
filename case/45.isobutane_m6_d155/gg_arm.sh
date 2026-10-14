#!/bin/bash
# plan architecture-float-state-double-geometry §6.58 (事前登録、codex diagnose 2026-10-12): 試験 G の 1 本の腕。G0c・FP64 (~/forge-fgeom7-fp64)・G0c の Q32 の起点
# (_g0c/q32_init.h5) から、停止規則 (§6.11、m9_watch.py) まで回す。上限 (lineDtDirectionalCap) を 50 (生産の B0) か 1 で選ぶ。
# 使い方: bash gg_arm.sh <新しい run> <50|1>
#   設定は run_0717_c50_d1 (§6.51 の B0・上限 50・出力に dt_local) の solverConfig.yaml から nStepOuter 30000 → 200000・outStepInterval 500 → 2500 の 2 行、
#   上限 1 ではさらに lineDtDirectionalCap 50.0 → 1.0 の 1 行だけを変える。
set -uo pipefail
R=$1; CAP=$2
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/gg.log; TL=~/forge-wallfit/solver_density_cuda/tools
BIN=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45
die() { echo "$R: 中止: $* ($(date -Is))" >> $LOG; exit 1; }
[ "$(sha256sum $BIN | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || die "FP64 のバイナリが §6.51 と違う"
[ "$(sha256sum _g0c/q32_init.h5 | cut -c1-16)" = "6f8fcda6d97482e3" ] || die "G0c の入力が登録と違う"
case $CAP in 50|1) ;; *) die "上限は 50 か 1";; esac
mkdir $R || die "$R が既にある"
cp _g0c/*.yaml $R/ && cp run_0717_c50_d1/solverConfig.yaml $R/solverConfig.yaml && cp _g0c/q32_init.h5 $R/nozzle.h5 || die "複製"
sed -i 's/nStepOuter: 30000}/nStepOuter: 200000}/; s/outStepInterval: 500$/outStepInterval: 2500/' $R/solverConfig.yaml
[ $CAP = 1 ] && sed -i 's/lineDtDirectionalCap: 50.0, /lineDtDirectionalCap: 1.0, /' $R/solverConfig.yaml
DL=$(diff run_0717_c50_d1/solverConfig.yaml $R/solverConfig.yaml | grep '^[<>]')
NEXP=4; [ $CAP = 1 ] && NEXP=6
[ "$(echo "$DL" | wc -l)" = $NEXP ] && [ -z "$(echo "$DL" | grep -v -e 'nStepOuter: ' -e 'outStepInterval: ' -e 'lineDtDirectionalCap: ')" ] \
  && grep -q "nStepOuter: 200000}" $R/solverConfig.yaml && grep -q "outStepInterval: 2500$" $R/solverConfig.yaml && grep -q "lineDtDirectionalCap: $CAP.0, " $R/solverConfig.yaml || die "設定の書き換え"
echo "$R 起動 $(date -Is) 上限 $CAP バイナリ $(sha256sum $BIN | cut -c1-16) 入力 $(sha256sum $R/nozzle.h5 | cut -c1-16)" >> $LOG
( export FORGE_BIN=$BIN; bash $TL/run_case.sh $C45/$R > $R/run_case_stdout.log 2>&1; echo $? > $R/RUN_RC ) &
python3 m9_watch.py $R --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $R/m9_watch.log 2>&1 &
wait
echo "$R 終了 $(date -Is) RUN_RC=$(cat $R/RUN_RC 2>/dev/null) 見張り $(python3 -c "import json; d=json.load(open('$R/m9_watch.json')); print(d['status'], d.get('reach_step'))" 2>&1)" >> $LOG
