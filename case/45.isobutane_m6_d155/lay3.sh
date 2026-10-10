#!/bin/bash
# plan time_integration-line-implicit-speed §6.19 (2026-10-10、§5.1 #13): LAYOUT2 を既定にする前の経路の確認。従来の並びとのビット一致 (FORGE_LINE_COMPARE=1)。
# (a) case/45 粘性入りのライン (lvc 3・マスク 5、FP64 の lineL) / (b) case/56 長さがばらつくライン・部分被覆 (float の lineL_f32) /
# (c) case/39 3D・周期・dual-time、lineKFreeze 0 と 1 (float の lineL_f32)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
LOG=$PWD/lay3.log
F32=$HOME/forge-linespeed-f32/solver_density_cuda/build/forge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
echo "== 開始 $(date -Is)" >> $LOG
[ -x $F32 ] || { echo "失敗: float のバイナリが無い" >> $LOG; touch lay3.done; exit 1; }
echo "float バイナリ sha256 $(sha256sum $F32 | cut -d' ' -f1)" >> $LOG
judge() {   # judge <run> <factors> <solves> <lines>
  python3 line_cmp_judge.py $1 --arm LAYOUT2 --eta-limit 1e-11 --bitwise --factors $2 --solves $3 --lines $4 >> $LOG 2>&1; local c=$?
  echo "$1 判定 rc=$c (0 PASS / 1 FAIL / 2 INDETERMINATE)" >> $LOG
}
# (a) case/45 L5
SRC=run_0183_ns_coldmesh_tw300_ext
r=run_0380_lay2cmp_L5
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
python3 cold_cfl.py prep $SRC $r --steps 20 --out 20 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --lvc 3 > /dev/null && rm -f $r/nozzle.msh \
  && ( export FORGE_LVC_TERMS=5 FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=2; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); echo "$r rc=$?" >> $LOG
judge $r 20 100 4719
rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_20.h5
unset FORGE_BIN COLD_ALT_BINARY REAL_CONVERTER FORGE_CONVERTER
# (b) case/56
r=../56.gap_tp1187/run_0073_lay2cmp_varlen
( export FORGE_BIN=$F32 FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=2; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$?" >> $LOG
grep -m1 "lineImplicit\] lines" $r/forge_run.log >> $LOG
judge $r 20 80 643
# (c) case/39
for k in 0 1; do
  r=../39.periodic_hills/run_095$((4 + k))_lay2cmp_kf$k
  ( export FORGE_BIN=$F32 FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=2; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$?" >> $LOG
  grep -m1 "lineImplicit\] lines" $r/forge_run.log >> $LOG
  if [ $k = 0 ]; then judge $r 40 200 2511; else judge $r 2 200 2511; fi
done
echo "== 終了 $(date -Is)" >> $LOG
touch lay3.done
