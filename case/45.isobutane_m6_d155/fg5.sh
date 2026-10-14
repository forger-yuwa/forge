#!/bin/bash
# plan architecture-float-state-double-geometry §6.6 (段 ④ の事前登録) の run (2026-10-10)。旧 = 段 ③、新 = 段 ④。
#   2  : FP64 での値の同一性 (FORGE_DIAG_GEOM_STAGE4_CHECK=1 の起動時の照合。case/45・case/39 (周期)・case/56)
#        float の case/56 (float32 の格子) の照合も記録する (内部面の e・r0・r1 が 0 のはず)
#   3  : FP64 isp 1 の B0 20 step (旧・新・旧の再実行) と、同じ状態からの 1 step (commit 前の残差と dq)
#   4  : float case/56 20 step (新、旧は段 ③ の run_0085、旧の再実行を足す)
#   5  : 記録: float case/45 B0 20 step (新、旧は段 ③ の run_0428)
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C56=~/forge-wallfit/case/56.gap_tp1187; C39=~/forge-wallfit/case/39.periodic_hills
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-fgeom3-fp64/solver_density_cuda/build/forge;  NEW64=~/forge-fgeom4-fp64/solver_density_cuda/build/forge
OLD32=~/forge-fgeom3-f32/solver_density_cuda/build/forge;   NEW32=~/forge-fgeom4-f32/solver_density_cuda/build/forge
LOG=$C45/fg5.log
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD64 NEW64 OLD32 NEW32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; B0=run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:4}" > /dev/null ) && rm -f $1/nozzle.msh || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
chk() { grep -h '^\[geomStage4\]' $1/forge_run.log $1/run_case_stdout.log 2>/dev/null | sort -u | sed "s#^#  $(basename $1): #" >> $LOG; }
# ---- 3 (と 2 の case/45: 新の FP64 の run に照合を付ける) ----
for r in run_0429_fg4_isp1_old run_0430_fg4_isp1_new run_0431_fg4_isp1_old2; do prep45 $r 20 20 --isp 1; done
run $C45/run_0429_fg4_isp1_old $OLD64; run $C45/run_0430_fg4_isp1_new $NEW64 FORGE_DIAG_GEOM_STAGE4_CHECK=1; chk $C45/run_0430_fg4_isp1_new
run $C45/run_0431_fg4_isp1_old2 $OLD64
for r in run_0432_fg4_1s_old run_0433_fg4_1s_new run_0434_fg4_1s_old2; do prep45 $r 1 1 --isp 1 --field-from $B0; out $r '{level: 2}'; done
run $C45/run_0432_fg4_1s_old $OLD64; run $C45/run_0433_fg4_1s_new $NEW64; run $C45/run_0434_fg4_1s_old2 $OLD64
# ---- 2 (case/39・case/56 の FP64 と、case/56 の float の照合。起動して 1 step) ----
mkdir $C39/run_0958_fg4_chk64 && cp $C39/run_0954_lay2cmp_kf0/{hill_xc_80x50x30.h5,bcondConfig.yaml,probe.yaml,solverConfig.yaml} $C39/run_0958_fg4_chk64/
sed -i 's/nStepOuter: 2/nStepOuter: 1/' $C39/run_0958_fg4_chk64/solverConfig.yaml
run $C39/run_0958_fg4_chk64 $NEW64 FORGE_DIAG_GEOM_STAGE4_CHECK=1; chk $C39/run_0958_fg4_chk64
mkdir $C56/run_0087_fg4_chk64 && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/run_0087_fg4_chk64/
sed -i 's/last: {nStepOuter: 20}/last: {nStepOuter: 1}/; s/outStepInterval: 20/outStepInterval: 1/' $C56/run_0087_fg4_chk64/solverConfig.yaml
run $C56/run_0087_fg4_chk64 $NEW64 FORGE_DIAG_GEOM_STAGE4_CHECK=1; chk $C56/run_0087_fg4_chk64
# ---- 4 (新の float の run に照合を付ける) ----
for r in run_0088_fg4_id32_new run_0089_fg4_id32_old2; do mkdir $C56/$r && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/$r/; done
run $C56/run_0088_fg4_id32_new $NEW32 FORGE_DIAG_GEOM_STAGE4_CHECK=1; chk $C56/run_0088_fg4_id32_new
run $C56/run_0089_fg4_id32_old2 $OLD32
# ---- 5 (記録) ----
prep45 run_0435_fg4_f32_new 20 20; run $C45/run_0435_fg4_f32_new $NEW32 FORGE_DIAG_GEOM_STAGE4_CHECK=1; chk $C45/run_0435_fg4_f32_new
for r in run_0429_fg4_isp1_old run_0430_fg4_isp1_new run_0431_fg4_isp1_old2 run_0432_fg4_1s_old run_0433_fg4_1s_new run_0434_fg4_1s_old2 run_0435_fg4_f32_new; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg5.done
