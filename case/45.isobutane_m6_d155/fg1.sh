#!/bin/bash
# plan architecture-float-state-double-geometry §5.1 #3 (段 ①): 既定の経路の同一性 (T1) と V0 の評価 (2026-10-10)。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C56=~/forge-wallfit/case/56.gap_tp1187
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-linespeed-fp64/solver_density_cuda/build/forge; NEW64=~/forge-fgeom-fp64/solver_density_cuda/build/forge
OLD32=~/forge-linespeed-f32/solver_density_cuda/build/forge; NEW32=~/forge-fgeom-f32/solver_density_cuda/build/forge
LOG=$C45/fg1.log
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD64 NEW64 OLD32 NEW32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
prep45() { ( export FORGE_BIN=$OLD64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $2 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 "${@:3}" > /dev/null ) && rm -f $1/nozzle.msh; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$?" >> $LOG; }
# T1 FP64 case/45 (20 step)
for r in run_0390_fg_id64_old run_0391_fg_id64_new run_0392_fg_id64_old2; do prep45 $r 20; done
run run_0390_fg_id64_old $OLD64; run run_0391_fg_id64_new $NEW64; run run_0392_fg_id64_old2 $OLD64
# T1 float case/56 (20 step)
for r in run_0079_fg_id32_old run_0080_fg_id32_new run_0081_fg_id32_old2; do mkdir $C56/$r && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/$r/; done
run $C56/run_0079_fg_id32_old $OLD32; run $C56/run_0080_fg_id32_new $NEW32; run $C56/run_0081_fg_id32_old2 $OLD32
# V0 (B0 の最終状態 = run_0252 res_60000)
for r in run_0393_fg_v0_dump run_0394_fg_v0_dump2 run_0395_fg_v0_ref; do prep45 $r 1 --field-from run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2; done
run run_0393_fg_v0_dump $NEW32 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0393_fg_v0_dump/geomab.h5
run run_0394_fg_v0_dump2 $NEW32 FORGE_DIAG_GEOMAB_DUMP=$C45/run_0394_fg_v0_dump2/geomab.h5
run run_0395_fg_v0_ref $NEW64 FORGE_DIAG_GEOMAB_REF=$C45/run_0393_fg_v0_dump/geomab.h5
grep -h "geomab" run_0393_fg_v0_dump/forge_run.log run_0395_fg_v0_ref/forge_run.log | head -20 >> $LOG
for r in run_0390_fg_id64_old run_0391_fg_id64_new run_0392_fg_id64_old2 run_0393_fg_v0_dump run_0394_fg_v0_dump2 run_0395_fg_v0_ref; do rm -f $r/nozzle.h5 $r/res_0.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg1.done
