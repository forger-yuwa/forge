#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.6 の 6 (case/48) の取り直し (2026-10-10)。
# fg9.sh の run_0051〜0053 は一様な初期値 (run_0048 の入力) から始めたので、新旧のバイナリとも同じく step 9 で発散した (台本の誤り)。
# ここでは収束場 (主ワークツリーの case/48 run_0025_B_tw300_y3_fx05/res_48000.h5、sha256 489470e5…、境界条件は run_0048 と同一) から 20 step。
# 格子は run_0048 と同じ (_conv614/c48/new64/m.h5 = 新しい変換器、double の幾何)。差は判定と同じく 新 = fgeom7-f32 / 旧 = fgeom6-f32 ×2。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4
LOG=$C45/fg9.log; RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
F7_32=~/forge-fgeom7-f32/solver_density_cuda/build/forge; F6_32=~/forge-fgeom6-f32/solver_density_cuda/build/forge
SEED=$C48/_seed/run_0025_res_48000.h5
echo "== fg9c 開始 $(date -Is) seed $(sha256sum $SEED | cut -c1-16)" >> $LOG
run() { local r=$1 b=$2; ( export FORGE_BIN=$b; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
for r in run_0054_hp_reg_new run_0055_hp_reg_old run_0056_hp_reg_oldb; do
  mkdir $C48/$r || { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  cp $C48/run_0048_fg6_dump32/{solverConfig.yaml,bcondConfig.yaml,probe.yaml,mesh.h5} $C48/$r/
  sed -i 's/nStepOuter: 1}/nStepOuter: 20}/; s/outStepInterval: 1$/outStepInterval: 20/' $C48/$r/solverConfig.yaml
  grep -q 'nStepOuter: 20}' $C48/$r/solverConfig.yaml && grep -q 'outStepInterval: 20$' $C48/$r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗 (中止)" >> $LOG; exit 1; }
  python3 $TL/restart_field.py $SEED $C48/$r/mesh.h5 --dst-run $C48/$r > $C48/$r/restart.log 2>&1 || { echo "$r: restart_field に失敗 (中止)" >> $LOG; exit 1; }
  echo "$r: $(grep VERDICT $C48/$r/restart.log)" >> $LOG
done
run $C48/run_0054_hp_reg_new $F7_32; run $C48/run_0055_hp_reg_old $F6_32; run $C48/run_0056_hp_reg_oldb $F6_32
echo "== fg9c 終了 $(date -Is)" >> $LOG
touch $C45/fg9c.done
