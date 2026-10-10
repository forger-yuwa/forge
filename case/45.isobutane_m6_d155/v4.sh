#!/bin/bash
# plan architecture-float-state-double-geometry §6.11 (V4、事前登録): 段 ④ + commit の診断 (56b1a6fe) の float と FP64 を、
# run_0183 の res_100000 から B0 の構成で、水準を 2 出力連続で満たした 2 つ目の出力まで (上限 20 万 step、2500 step ごと)。2 本を同時に投入する。
# 見張りは m9_watch.py (--interval 2500 --min-points 9 --keep-tail 9)。commit の診断 FORGE_DIAG_COMMIT_LOSS=2500 を両方に付ける。
set -uo pipefail
cd ~/forge-wallfit/case/45.isobutane_m6_d155
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B32=~/forge-fgeom5-f32/solver_density_cuda/build/forge; B64=~/forge-fgeom5-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$PWD/v4.log
echo "== 開始 $(date -Is)" >> $LOG
for b in B32 B64; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 10737418240 ] || { echo "ディスクの空きが 10 GiB 未満 (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps 200000 --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$r: 準備に失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  grep -n '^time\|nStepOuter\|outStepInterval\|lineImplicit\|implicitRelax\|extraFields' $r/solverConfig.yaml | sed "s#^#  $r: #" >> $LOG
  ( export FORGE_BIN=$b FORGE_DIAG_COMMIT_LOSS=2500; bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  python3 m9_watch.py $r --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $r/m9_watch.log 2>&1 &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG
}
arm run_0451_v4_f32 $B32
arm run_0452_v4_fp64 $B64
wait
for r in run_0451_v4_f32 run_0452_v4_fp64; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 見張り=$(python3 -c "import json; s=json.load(open('$r/m9_watch.json')); print(s['status'], s.get('reach_step'))")" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
touch v4.done
