#!/bin/bash
# plan architecture-float-state-double-geometry §6.4 (段 ③ の事前登録) と §6.5 の 4d の run (2026-10-10)。
#   2  : FP64 の面流束の同一性 (段 ③ の FP64 で FORGE_DIAG_GEOMAB_DUMP、B0 の状態)
#   3  : FP64 isp 1 の B0 20 step (旧 = 段 ②・新 = 段 ③・旧の再実行)
#   4  : float case/56 20 step (新 = 段 ③、旧 = 段 ② の run_0084、旧の再実行を足す)
#   4d : FP64 isp 1 の B0 の状態から 1 step (level 2): commit 前の残差と dq (旧・新・旧の再実行)
#   5  : 記録: FP64 isp 0 の 20 step (新、旧は run_0398)、float case/45 B0 20 step (旧・新)
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C56=~/forge-wallfit/case/56.gap_tp1187
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-fgeom2-fp64/solver_density_cuda/build/forge;  NEW64=~/forge-fgeom3-fp64/solver_density_cuda/build/forge
OLD32=~/forge-fgeom2-f32/solver_density_cuda/build/forge;   NEW32=~/forge-fgeom3-f32/solver_density_cuda/build/forge
LOG=$C45/fg4.log
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD64 NEW64 OLD32 NEW32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; B0=run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:4}" > /dev/null ) && rm -f $1/nozzle.msh || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
# ---- 2 ----
prep45 run_0419_fg3_v0dump64 1 1 --field-from $B0
run $C45/run_0419_fg3_v0dump64 $NEW64 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0419_fg3_v0dump64/geomab.h5
grep -h "geomab" run_0419_fg3_v0dump64/forge_run.log | head -8 >> $LOG
# ---- 3 ----
for r in run_0420_fg3_isp1_old run_0421_fg3_isp1_new run_0422_fg3_isp1_old2; do prep45 $r 20 20 --isp 1; done
run $C45/run_0420_fg3_isp1_old $OLD64; run $C45/run_0421_fg3_isp1_new $NEW64; run $C45/run_0422_fg3_isp1_old2 $OLD64
# ---- 4d ----
for r in run_0423_fg3_1s_old run_0424_fg3_1s_new run_0425_fg3_1s_old2; do prep45 $r 1 1 --isp 1 --field-from $B0; out $r '{level: 2}'; done
run $C45/run_0423_fg3_1s_old $OLD64; run $C45/run_0424_fg3_1s_new $NEW64; run $C45/run_0425_fg3_1s_old2 $OLD64
# ---- 4 ----
for r in run_0085_fg3_id32_new run_0086_fg3_id32_old2; do mkdir $C56/$r && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/$r/; done
run $C56/run_0085_fg3_id32_new $NEW32; run $C56/run_0086_fg3_id32_old2 $OLD32
# ---- 5 (記録) ----
prep45 run_0426_fg3_isp0_new 20 20; run $C45/run_0426_fg3_isp0_new $NEW64
for r in run_0427_fg3_f32_old run_0428_fg3_f32_new; do prep45 $r 20 20; done
run $C45/run_0427_fg3_f32_old $OLD32; run $C45/run_0428_fg3_f32_new $NEW32
for r in run_0419_fg3_v0dump64 run_0420_fg3_isp1_old run_0421_fg3_isp1_new run_0422_fg3_isp1_old2 run_0423_fg3_1s_old run_0424_fg3_1s_new run_0425_fg3_1s_old2 run_0426_fg3_isp0_new run_0427_fg3_f32_old run_0428_fg3_f32_new; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg4.done
