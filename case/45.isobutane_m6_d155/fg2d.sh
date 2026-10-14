#!/bin/bash
# plan architecture-float-state-double-geometry §6.3 の 3' (2026-10-10 追加、実行前に登録): closure を使う経路で 3 を回し直す。
# 3 の run (run_0399〜0404) は case/45 の既定 (hoopAreaFromClosure 0・axisRFloor 0) で、closure は hoop ソースに使われていなかった
# (dump の closure_active = 0)。ここでは mesh に hoopAreaFromClosure: 1 を足し、段 ② の closure の変更が効く経路で同じ判定をする。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-fgeom-fp64/solver_density_cuda/build/forge;  NEW64=~/forge-fgeom2-fp64/solver_density_cuda/build/forge
OLD32=~/forge-fgeom-f32/solver_density_cuda/build/forge;   NEW32=~/forge-fgeom2-f32/solver_density_cuda/build/forge
LOG=$C45/fg2.log
echo "== 開始 (fg2d) $(date -Is)" >> $LOG
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用のバイナリ (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) && rm -f $1/nozzle.msh; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
STILL="run_0405_fg2_stillhc_old32 run_0406_fg2_stillhc_new32 run_0407_fg2_stillhc_old32b run_0408_fg2_stillhc_new64"
for r in $STILL; do
  prep45 $r 10 1 || { echo "$r: 準備に失敗 (中止)" >> $LOG; exit 1; }
  python3 still_field.py $r --state-from $SRC/res_100000.h5 >> $LOG 2>&1 || { echo "$r: still_field に失敗 (中止)" >> $LOG; exit 1; }
  if grep -q '^output\|hoopAreaFromClosure' $r/solverConfig.yaml; then echo "$r: output か hoopAreaFromClosure が既にある (中止)" >> $LOG; exit 1; fi
  grep -c '^mesh: {' $r/solverConfig.yaml | grep -qx 1 || { echo "$r: mesh がフロー形式 1 行でない (中止)" >> $LOG; exit 1; }
  sed -i 's/^mesh: {/mesh: {hoopAreaFromClosure: 1, /' $r/solverConfig.yaml
  echo 'output: {level: 1, extraFields: [res_roUy]}' >> $r/solverConfig.yaml
  grep -n '^mesh' $r/solverConfig.yaml >> $LOG
done
run $C45/run_0405_fg2_stillhc_old32 $OLD32; run $C45/run_0406_fg2_stillhc_new32 $NEW32
run $C45/run_0407_fg2_stillhc_old32b $OLD32; run $C45/run_0408_fg2_stillhc_new64 $NEW64
for r in $STILL; do for n in 2 3 4 5 6 7 8 9; do rm -f $r/res_$n.h5 $r/res_$n.xmf; done; rm -f $r/nozzle.h5; done
grep -h "hoopAreaFromClosure\|closure" run_0406_fg2_stillhc_new32/forge_run.log | head -5 >> $LOG
echo "== 終了 (fg2d) $(date -Is)" >> $LOG
touch fg2d.done
