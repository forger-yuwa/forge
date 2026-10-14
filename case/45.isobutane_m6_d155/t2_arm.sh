#!/bin/bash
# plan architecture-float-state-double-geometry §6.56 (事前登録、codex diagnose 2026-10-12): 上限 1 の float と FP64 を停止規則まで比べる (T2)。
# 1 本の腕を、§6.51 の上限 1 の run の res_30000 から新しい run に restart_field.py で再開し、見張り m9_watch.py (§6.11 と同じ引数) を付けて回す。
# 使い方: bash t2_arm.sh <新しい run> <元の run のディレクトリ (res_30000.h5・*.yaml を含む)> <f32|fp64>
#   設定は元の run の solverConfig.yaml から nStepOuter 30000 → 170000・outStepInterval 500 → 2500 の 2 行だけを変える (通算 20 万 step が上限)。
set -uo pipefail
R=$1; SRC=$2; P=$3
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/t2.log; TL=~/forge-wallfit/solver_density_cuda/tools
B=solver_density_cuda/build/forge
case $P in f32) BIN=~/forge-fgeom7-f32/$B; WANT=793215917d156680;; fp64) BIN=~/forge-fgeom7-fp64/$B; WANT=7e95e5bbd8f50c9e;; *) echo "精度は f32 か fp64"; exit 1;; esac
cd $C45
die() { echo "$R: 中止: $* ($(date -Is))" >> $LOG; exit 1; }
[ "$(sha256sum $BIN | cut -c1-16)" = "$WANT" ] || die "バイナリが §6.51 と違う"
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || die "共通の入力が §6.26 と違う"
[ -e $SRC/res_30000.h5 ] && [ -e $SRC/solverConfig.yaml ] || die "元の run に res_30000.h5・solverConfig.yaml が無い"
grep -q "lineDtDirectionalCap: 1.0, " $SRC/solverConfig.yaml || die "元の run が上限 1 でない"
mkdir $R || die "$R が既にある"
cp $SRC/*.yaml $R/ && cp _tr/q32_init.h5 $R/nozzle.h5 || die "複製"
sed -i 's/nStepOuter: 30000}/nStepOuter: 170000}/; s/outStepInterval: 500/outStepInterval: 2500/' $R/solverConfig.yaml
# 差は 2 行の変更 = diff の < と > が 4 行で、どれも nStepOuter か outStepInterval の行 (2026-10-12: 最初は 2 と決め打ちして止まった)
DL=$(diff $SRC/solverConfig.yaml $R/solverConfig.yaml | grep '^[<>]')
[ "$(echo "$DL" | wc -l)" = 4 ] && [ -z "$(echo "$DL" | grep -v -e 'nStepOuter: ' -e 'outStepInterval: ')" ] && grep -q "nStepOuter: 170000}" $R/solverConfig.yaml && grep -q "outStepInterval: 2500$" $R/solverConfig.yaml || die "設定の書き換え (2 行だけのはず)"
KD=""; [ $P = fp64 ] && KD="--keep-src-dtype"
python3 $TL/restart_field.py $SRC/res_30000.h5 $R/nozzle.h5 $KD --dst-run $R --forge $BIN > $R/restart.log 2>&1 || die "restart_field.py"
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $R/restart.log || die "restart_field の VERDICT が 9 量のビット一致でない: $(grep VERDICT $R/restart.log)"
echo "$R 起動 $(date -Is) $P バイナリ $(sha256sum $BIN | cut -c1-16) 元 $SRC/res_30000 $(sha256sum $SRC/res_30000.h5 | cut -c1-16) 入力 $(sha256sum $R/nozzle.h5 | cut -c1-16)" >> $LOG
( export FORGE_BIN=$BIN; bash $TL/run_case.sh $C45/$R > $R/run_case_stdout.log 2>&1; echo $? > $R/RUN_RC ) &
python3 m9_watch.py $R --phase line --budget 170000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $R/m9_watch.log 2>&1 &
wait
echo "$R 終了 $(date -Is) RUN_RC=$(cat $R/RUN_RC 2>/dev/null) 見張り $(python3 -c "import json; d=json.load(open('$R/m9_watch.json')); print(d['status'], d.get('reach_step'))" 2>&1)" >> $LOG
