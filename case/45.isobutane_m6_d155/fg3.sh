#!/bin/bash
# (1) plan time_integration-line-implicit-speed §6.23 (事前登録): B0 の保存場から hoopAreaFromClosure 0/1 (+0 の再実行) を 2000 step、
#     同じ出発点から 1 step (level 2) を 3 本。段 ② の FP64 のバイナリ。
# (2) plan architecture-float-state-double-geometry §6.5 の 4b: closure を使う経路のダンプ (float・FP64) と、
#     静止場 (3') の FP64 の旧 (段 ①)・旧の再実行 (新は run_0408)。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-fgeom-fp64/solver_density_cuda/build/forge;  NEW64=~/forge-fgeom2-fp64/solver_density_cuda/build/forge
NEW32=~/forge-fgeom2-f32/solver_density_cuda/build/forge
LOG=$C45/fg3.log
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD64 NEW64 NEW32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; B0=run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:4}" > /dev/null ) && rm -f $1/nozzle.msh || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
hc1() { grep -c '^mesh: {' $1/solverConfig.yaml | grep -qx 1 && ! grep -q 'hoopAreaFromClosure' $1/solverConfig.yaml || { echo "$1: mesh の形が想定外 (中止)" >> $LOG; exit 1; }; sed -i 's/^mesh: {/mesh: {hoopAreaFromClosure: 1, /' $1/solverConfig.yaml; }
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
# ---- (2) 4b ----
for r in run_0415_fg2_dumphc32 run_0416_fg2_dumphc64; do prep45 $r 1 1; hc1 $r; done
run $C45/run_0415_fg2_dumphc32 $NEW32 FORGE_DIAG_GEOM_STAGE2_DUMP=$C45/run_0415_fg2_dumphc32/stage2.h5
run $C45/run_0416_fg2_dumphc64 $NEW64 FORGE_DIAG_GEOM_STAGE2_DUMP=$C45/run_0416_fg2_dumphc64/stage2.h5
grep -h '^\[geomStage2\] axisym' run_0415_fg2_dumphc32/forge_run.log run_0416_fg2_dumphc64/forge_run.log >> $LOG
for r in run_0417_fg2_stillhc_old64 run_0418_fg2_stillhc_old64b; do
  prep45 $r 10 1; python3 still_field.py $r --state-from $SRC/res_100000.h5 >> $LOG 2>&1 || { echo "$r: still_field に失敗 (中止)" >> $LOG; exit 1; }
  hc1 $r; out $r '{level: 1, extraFields: [res_roUy]}'
done
run $C45/run_0417_fg2_stillhc_old64 $OLD64; run $C45/run_0418_fg2_stillhc_old64b $OLD64
# ---- (1) §6.23 ----
for r in run_0409_hc_A run_0410_hc_B run_0411_hc_A2; do prep45 $r 2000 100 --field-from $B0; out $r '{level: 1, extraFields: [res_ro, res_roUx, res_roUy]}'; done
hc1 run_0410_hc_B
for r in run_0412_hc1_A run_0413_hc1_B run_0414_hc1_A2; do prep45 $r 1 1 --field-from $B0; out $r '{level: 2}'; done
hc1 run_0413_hc1_B
for r in run_0409_hc_A run_0410_hc_B run_0411_hc_A2 run_0412_hc1_A run_0413_hc1_B run_0414_hc1_A2; do grep -n '^mesh\|^output' $r/solverConfig.yaml | sed "s#^#  $r: #" >> $LOG; done
run $C45/run_0412_hc1_A $NEW64; run $C45/run_0413_hc1_B $NEW64; run $C45/run_0414_hc1_A2 $NEW64
run $C45/run_0409_hc_A $NEW64; run $C45/run_0410_hc_B $NEW64; run $C45/run_0411_hc_A2 $NEW64
for r in run_0415_fg2_dumphc32 run_0416_fg2_dumphc64 run_0417_fg2_stillhc_old64 run_0418_fg2_stillhc_old64b run_0409_hc_A run_0410_hc_B run_0411_hc_A2 run_0412_hc1_A run_0413_hc1_B run_0414_hc1_A2; do rm -f $r/nozzle.h5; done
for r in run_0417_fg2_stillhc_old64 run_0418_fg2_stillhc_old64b; do for n in 2 3 4 5 6 7 8 9; do rm -f $r/res_$n.h5 $r/res_$n.xmf; done; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg3.done
