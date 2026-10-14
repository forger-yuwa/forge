#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.6 (事前登録) の 2〜6 の run (2026-10-10)。新 = ~/forge-fgeom7-* (区間ごとの r 重み、未 commit)、旧 = ~/forge-fgeom6-* (cb6a0d8b)。
#   2: 変換器 (case/45 の正本 run_0183 の nozzle.msh・case/48・case/39)
#   3・4: case/45 の新しい格子で一様な静止場 (P 1e5 Pa、hoopAreaFromClosure 0) を 10 step、FP64 と float、キー 0/1、FORGE_DIAG_HOOP_CLOSURE
#   5: 適用の分岐 (axisRFloor > 0・axisymMethod 1 でキー 1 も OFF、新キー 0 と旧のバイナリの 1 step)
#   6: 軸対称でないケースの回帰 (case/48・case/56 の 20 step、float)
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4; C56=~/forge-wallfit/case/56.gap_tp1187
W=$C45/_conv9; LOG=$C45/fg9.log; RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
CMP=~/forge-fgeom7-f32/solver_density_cuda/tools/compare_mesh_h5.py
V7_64=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge; V7_32=~/forge-fgeom7-f32/solver_density_cuda/build/convertGmshToForge
F7_64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; F7_32=~/forge-fgeom7-f32/solver_density_cuda/build/forge
F6_64=~/forge-fgeom6-fp64/solver_density_cuda/build/forge; F6_32=~/forge-fgeom6-f32/solver_density_cuda/build/forge
echo "== 開始 $(date -Is)" >> $LOG
for b in V7_64 V7_32 F7_64 F7_32 F6_64 F6_32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
conv() {   # conv <case> <arm> <変換器> <msh> <出力名>
  local d=$W/$1/$2; mkdir -p $d; cp $C45/_conv614/inp/$1/*.yaml $d/
  ( cd $d && $3 $C45/_conv614/inp/$1/$4 $5 > conv.log 2>&1 ); local rc=$?
  if [ $rc -ne 0 ] && grep -q "Write Input HDF5" $d/conv.log && grep -q "GPUassert: invalid argument" $d/conv.log; then rc=0; fi
  echo "conv $1 $2 rc=$rc" >> $LOG
}
# ---- 2 ----
conv c45 new64 $V7_64 nozzle.msh nozzle.h5; conv c45 new32 $V7_32 nozzle.msh nozzle.h5
conv c48 new64 $V7_64 fp_y1_3um.msh m.h5; conv c39 new64 $V7_64 hill_des_80x50x30.msh hill_des_80x50x30.h5
echo "== 2a case/45 の rSurfVect (FP64 の新しい変換器) の独立な作り直しとの照合" >> $LOG
python3 hoop_verify_conv.py check $W/c45/new64/nozzle.h5 --rw >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 2b case/45 新 float 対 新 FP64 (--geometry-only --require-float64)" >> $LOG
python3 $CMP $W/c45/new32/nozzle.h5 $W/c45/new64/nozzle.h5 --geometry-only --require-float64 >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 2c case/45 旧 (§6.14 の新 = fgeom6) 対 新 の FP64 の差 (変わるのは境界の半割面と rSurfVect の追加だけのはず)" >> $LOG
python3 hoop_verify_conv.py diff $C45/_conv614/c45/new64/nozzle.h5 $W/c45/new64/nozzle.h5 >> $LOG 2>&1
echo "== 2d case/48 (平面) 旧 対 新" >> $LOG; python3 hoop_verify_conv.py diff $C45/_conv614/c48/new64/m.h5 $W/c48/new64/m.h5 >> $LOG 2>&1
echo "== 2e case/39 (3D) 旧 対 新 (全件、ビット一致のはず)" >> $LOG; python3 $CMP $C45/_conv614/c39/new64/hill_des_80x50x30.h5 $W/c39/new64/hill_des_80x50x30.h5 >> $LOG 2>&1; echo "  rc=$?" >> $LOG
# ---- 3・4: 静止場 (新しい格子) ----
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
prepnew() {   # prepnew <run> <steps> <out>: B0 の構成で準備し、格子を新しい格子に替えて場を移す
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }
  rm -f $1/nozzle.msh; mv $1/nozzle.h5 $1/nozzle_oldmesh.h5; cp $W/c45/new64/nozzle.h5 $1/nozzle.h5
  python3 $TL/restart_field.py $1/nozzle_oldmesh.h5 $1/nozzle.h5 --dst-run $1 --forge $F7_64 > $1/restart_newmesh.log 2>&1 || { echo "$1: 新しい格子への場の移しに失敗 (中止)" >> $LOG; exit 1; }
  rm -f $1/nozzle_oldmesh.h5
}
key() {   # key <run> <0|1>
  grep -c '^mesh: {' $1/solverConfig.yaml | grep -qx 1 || { echo "$1: mesh の形が想定外 (中止)" >> $LOG; exit 1; }
  sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: $2, /" $1/solverConfig.yaml
}
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b)))) $(grep -h 'axisSegmentRWeight' $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | sort -u | head -1 | cut -c1-140)" >> $LOG; }
STILL="run_0463_hp_still64_k1 run_0464_hp_still64_k0 run_0465_hp_still32_k1 run_0466_hp_still32_k1b run_0467_hp_still32_k0 run_0468_hp_still32_k0b"
for r in $STILL; do
  prepnew $r 10 1; python3 still_field.py $r --state-from $SRC/res_100000.h5 >> $LOG 2>&1 || { echo "$r: still_field に失敗 (中止)" >> $LOG; exit 1; }
  case $r in *_k1*) key $r 1;; *) key $r 0;; esac; out $r '{level: 1, extraFields: [res_roUy, res_roUx]}'
done
run $C45/run_0463_hp_still64_k1 $F7_64 FORGE_DIAG_HOOP_CLOSURE=$C45/run_0463_hp_still64_k1/hoop.h5
run $C45/run_0464_hp_still64_k0 $F7_64 FORGE_DIAG_HOOP_CLOSURE=$C45/run_0464_hp_still64_k0/hoop.h5
for r in run_0465_hp_still32_k1 run_0466_hp_still32_k1b run_0467_hp_still32_k0 run_0468_hp_still32_k0b; do run $C45/$r $F7_32 FORGE_DIAG_HOOP_CLOSURE=$C45/$r/hoop.h5; done
for r in $STILL; do for n in 2 3 4 5 6 7 8 9; do rm -f $r/res_$n.h5 $r/res_$n.xmf; done; done
# ---- 5: 適用の分岐 (FP64、1 step、新しい格子) ----
for r in run_0469_hp_br_floor_k1 run_0470_hp_br_floor_k0 run_0471_hp_br_m1_k1 run_0472_hp_br_m1_k0 run_0473_hp_br_k0new run_0474_hp_br_old run_0475_hp_br_oldb; do prepnew $r 1 1; out $r '{level: 2}'; done
for r in run_0469_hp_br_floor_k1 run_0470_hp_br_floor_k0; do sed -i 's/^mesh: {/mesh: {axisRFloor: 1.0e-4, /' $r/solverConfig.yaml; done
for r in run_0471_hp_br_m1_k1 run_0472_hp_br_m1_k0; do sed -i 's/^mesh: {/mesh: {axisymMethod: 1, /' $r/solverConfig.yaml; done
key run_0469_hp_br_floor_k1 1; key run_0470_hp_br_floor_k0 0; key run_0471_hp_br_m1_k1 1; key run_0472_hp_br_m1_k0 0; key run_0473_hp_br_k0new 0
for r in run_0469_hp_br_floor_k1 run_0470_hp_br_floor_k0 run_0471_hp_br_m1_k1 run_0472_hp_br_m1_k0 run_0473_hp_br_k0new; do run $C45/$r $F7_64 FORGE_DIAG_HOOP_CLOSURE=$C45/$r/hoop.h5; done
run $C45/run_0474_hp_br_old $F6_64; run $C45/run_0475_hp_br_oldb $F6_64
# ---- 6: 軸対称でないケースの回帰 (float、20 step) ----
for r in run_0051_hp_reg_new run_0052_hp_reg_old run_0053_hp_reg_oldb; do mkdir $C48/$r && cp $C48/run_0048_fg6_dump32/{solverConfig.yaml,bcondConfig.yaml,probe.yaml,mesh.h5} $C48/$r/; sed -i 's/nStepOuter: 1}/nStepOuter: 20}/; s/outStepInterval: 1$/outStepInterval: 20/' $C48/$r/solverConfig.yaml; done
run $C48/run_0051_hp_reg_new $F7_32; run $C48/run_0052_hp_reg_old $F6_32; run $C48/run_0053_hp_reg_oldb $F6_32
for r in run_0093_hp_reg_new run_0094_hp_reg_old run_0095_hp_reg_oldb; do mkdir $C56/$r && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/$r/; done
run $C56/run_0093_hp_reg_new $F7_32; run $C56/run_0094_hp_reg_old $F6_32; run $C56/run_0095_hp_reg_oldb $F6_32
for r in $STILL run_0469_hp_br_floor_k1 run_0470_hp_br_floor_k0 run_0471_hp_br_m1_k1 run_0472_hp_br_m1_k0 run_0473_hp_br_k0new run_0474_hp_br_old run_0475_hp_br_oldb; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg9.done
