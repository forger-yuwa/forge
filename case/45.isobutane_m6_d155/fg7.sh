#!/bin/bash
# plan architecture-float-state-double-geometry §6.12 (事前登録): V3 と同じ Q32 から、段 ④ の FP64 で流束の直前の入力を
#   A = 自分で作る (FORGE_DIAG_GEOMAB_DUMP) / B = float の DUMP (run_0446) から移す (FORGE_DIAG_GEOMAB_REF) で 1 step、各 2 回。
# B の 1 本目は V3 の run_0447_v3_ref64dump (同じバイナリ・同じ入力)。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
REF64=~/forge-fgeom4-fp64/solver_density_cuda/build/forge
LOG=$C45/fg7.log
echo "== 開始 $(date -Is)" >> $LOG
echo "REF64 $(sha256sum $REF64 | cut -c1-16)" >> $LOG
cd $C45
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; B0=run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:4}" > /dev/null ) && rm -f $1/nozzle.msh || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
for r in run_0448_v3_self64a run_0449_v3_self64b run_0450_v3_ref64dump2; do
  prep45 $r 1 1 --isp 1 --field-from $B0
  python3 q32.py $r >> $LOG 2>&1 || { echo "$r: q32 に失敗 (中止)" >> $LOG; exit 1; }
done
[ "$(grep -h '^\[q32\]' $LOG | tail -3 | awk '{print $NF}' | sort -u)" = "66ecd4ba2a2bfd03" ] || { echo "Q32 の場が V3 と違う (中止)" >> $LOG; exit 1; }
run $C45/run_0448_v3_self64a $REF64 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0448_v3_self64a/geomab.h5
run $C45/run_0449_v3_self64b $REF64 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0449_v3_self64b/geomab.h5
cp run_0446_v3_dump32/geomab.h5 run_0450_v3_ref64dump2/geomab_src.h5
run $C45/run_0450_v3_ref64dump2 $REF64 FORGE_DIAG_GEOMAB_REF=$C45/run_0450_v3_ref64dump2/geomab_src.h5
grep -h "geomab" run_0448_v3_self64a/forge_run.log run_0450_v3_ref64dump2/forge_run.log | cut -c1-200 | head -12 >> $LOG
for r in run_0448_v3_self64a run_0449_v3_self64b run_0450_v3_ref64dump2; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg7.done
