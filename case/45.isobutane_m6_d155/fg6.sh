#!/bin/bash
# plan architecture-float-state-double-geometry §6.8 (V3、事前登録): 同じ状態 Q32 (B0 の状態を float に丸めたもの) から 1 step、
# commit 前の残差を float の旧 (段 ①)・新 (段 ④)・FP64 の参照 (段 ④) で比べる。各 2 本 (再実行の差)。
# 加えて、粘性・k/ω の面の流束を V0 の経路 (GEOMAB の DUMP/REF) で同じ Q32 から測る (記録)。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD32=~/forge-fgeom-f32/solver_density_cuda/build/forge;  NEW32=~/forge-fgeom4-f32/solver_density_cuda/build/forge
REF64=~/forge-fgeom4-fp64/solver_density_cuda/build/forge
LOG=$C45/fg6.log
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD32 NEW32 REF64; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; B0=run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:4}" > /dev/null ) && rm -f $1/nozzle.msh || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
ARMS="run_0440_v3_old32a run_0441_v3_old32b run_0442_v3_new32a run_0443_v3_new32b run_0444_v3_ref64a run_0445_v3_ref64b"
for r in $ARMS run_0446_v3_dump32 run_0447_v3_ref64dump; do
  prep45 $r 1 1 --isp 1 --field-from $B0
  python3 q32.py $r >> $LOG 2>&1 || { echo "$r: q32 に失敗 (中止)" >> $LOG; exit 1; }
done
for r in $ARMS; do out $r '{level: 2}'; done
# 丸めた場が腕どうしで同じこと (q32 の sha256 が 1 種類)
[ "$(grep -h '^\[q32\]' $LOG | tail -8 | awk '{print $NF}' | sort -u | wc -l)" = 1 ] || { echo "Q32 の場が腕で違う (中止)" >> $LOG; exit 1; }
run $C45/run_0440_v3_old32a $OLD32; run $C45/run_0441_v3_old32b $OLD32
run $C45/run_0442_v3_new32a $NEW32; run $C45/run_0443_v3_new32b $NEW32
run $C45/run_0444_v3_ref64a $REF64; run $C45/run_0445_v3_ref64b $REF64
run $C45/run_0446_v3_dump32 $NEW32 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0446_v3_dump32/geomab.h5
run $C45/run_0447_v3_ref64dump $REF64 FORGE_DIAG_GEOMAB_REF=$C45/run_0446_v3_dump32/geomab.h5
grep -h "geomab" run_0446_v3_dump32/forge_run.log run_0447_v3_ref64dump/forge_run.log | head -12 >> $LOG
for r in $ARMS run_0446_v3_dump32 run_0447_v3_ref64dump; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg6.done
