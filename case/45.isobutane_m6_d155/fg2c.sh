#!/bin/bash
# fg2b.sh の 3・4 の回し直し (2026-10-10): fg2b.sh は still_field.py を AWS に送っておらず、送った後も保存量だけの場で T が無く止まった (1・2 のダンプは済み)。
# plan architecture-float-state-double-geometry §5.1 #4 (段 ②): §6.3 の合格条件 1〜5 の run。
#   1・2: 読み込み時の量のダンプ (FORGE_DIAG_GEOM_STAGE2_DUMP) を新しい float / FP64 のビルドで case/45・56・39・48 に
#   3   : 一様な静止場 (still_field.py) で closure の非劣化。旧 = 段 ① のビルド、10 step・毎 step 出力 (res_1 と res_10 だけ残す)
#   4   : FP64 case/45 B0 20 step (新)。旧と再実行は段 ① の run_0390・run_0392 を使う
#   5   : float case/56 20 step (新)。旧と再実行は段 ① の run_0079・run_0081 を使う
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C56=~/forge-wallfit/case/56.gap_tp1187
C39=~/forge-wallfit/case/39.periodic_hills; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
OLD64=~/forge-fgeom-fp64/solver_density_cuda/build/forge;  NEW64=~/forge-fgeom2-fp64/solver_density_cuda/build/forge
OLD32=~/forge-fgeom-f32/solver_density_cuda/build/forge;   NEW32=~/forge-fgeom2-f32/solver_density_cuda/build/forge
CONV64=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
LOG=$C45/fg2.log   # fg2b: case/45 の部分の回し直し
echo "== 開始 $(date -Is)" >> $LOG
for b in OLD64 NEW64 OLD32 NEW32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
cd $C45
export REAL_CONVERTER=$CONV64 FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用のバイナリ (計算には使わない)
prep45() { ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) && rm -f $1/nozzle.msh; }
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
dump() { run $1 $2 FORGE_DIAG_GEOM_STAGE2_DUMP=$1/stage2.h5; grep -h '^\[geomStage2\]\|^\[wallRepPoint\]\|^\[setPeriodicPartner\]\|^\[lineImplicit\]' $1/forge_run.log $1/run_case_stdout.log 2>/dev/null | sort -u | sed "s#^#  $(basename $1): #" >> $LOG; }
rm -rf run_0399_fg2_still_old32   # fg2b が準備だけして止まった run (自分の run)
# ---- 3: 一様な静止場 (closure) ----
STILL="run_0399_fg2_still_old32 run_0400_fg2_still_new32 run_0401_fg2_still_old64 run_0402_fg2_still_new64 run_0403_fg2_still_old64b run_0404_fg2_still_old32b"
for r in $STILL; do
  prep45 $r 10 1 || { echo "$r: 準備に失敗 (中止)" >> $LOG; exit 1; }
  python3 still_field.py $r --state-from $SRC/res_100000.h5 >> $LOG 2>&1 || { echo "$r: still_field に失敗 (中止)" >> $LOG; exit 1; }
  if grep -q '^output' $r/solverConfig.yaml; then echo "$r: output キーが既にある (中止)" >> $LOG; exit 1; fi
  echo 'output: {level: 1, extraFields: [res_roUy]}' >> $r/solverConfig.yaml
done
run $C45/run_0399_fg2_still_old32 $OLD32; run $C45/run_0400_fg2_still_new32 $NEW32
run $C45/run_0401_fg2_still_old64 $OLD64; run $C45/run_0402_fg2_still_new64 $NEW64
run $C45/run_0403_fg2_still_old64b $OLD64; run $C45/run_0404_fg2_still_old32b $OLD32
for r in $STILL; do for n in 2 3 4 5 6 7 8 9; do rm -f $r/res_$n.h5 $r/res_$n.xmf; done; done
# ---- 4: FP64 case/45 B0 20 step ----
prep45 run_0398_fg2_id64_new 20 20 || { echo "run_0398: 準備に失敗 (中止)" >> $LOG; exit 1; }; run $C45/run_0398_fg2_id64_new $NEW64
for r in run_0396_fg2_dump32 run_0397_fg2_dump64 run_0398_fg2_id64_new $STILL; do rm -f $r/nozzle.h5; done
echo "== 終了 $(date -Is)" >> $LOG
touch fg2c.done
