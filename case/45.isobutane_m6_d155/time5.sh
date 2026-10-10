#!/bin/bash
# plan time_integration-line-implicit-speed §6.19 (既定化の後の速さ) と float ビルドの速さの下見 (2026-10-10、ユーザ「float ビルドできるの？早くなる？」)。
# B0 の構成 (値 0 + 方向別 + キー 5 + 上限 50) を、既定化した FP64 (D64) と既定化した float (F32) で交互に 1000 step、環境変数なし (既定の並び)。GPU 専有。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
LOG=$PWD/time5.log
declare -A BIN=([D64]=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge [F32]=$HOME/forge-linespeed-f32/solver_density_cuda/build/forge)
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
SRC=run_0183_ns_coldmesh_tw300_ext
gpucount() { local o; o=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null) || return 1; echo "$o" | grep -c '[0-9]'; true; }
echo "== 開始 $(date -Is)" >> $LOG
k=0
for v in D64 F32 D64 F32 D64 F32 D64; do
  k=$((k + 1)); r=run_0383_t5_$(printf %02d $k)_$v
  ( export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
    python3 cold_cfl.py prep $SRC $r --steps 1000 --out 1000 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) || { echo "prep $r 失敗" >> $LOG; continue; }
  rm -f $r/nozzle.msh
  n=$(gpucount) && echo "ok $n" > $r/gpu_pre.txt || echo fail > $r/gpu_pre.txt
  ( export FORGE_BIN=${BIN[$v]}; bash $RC $r > $r/run_case_stdout.log 2>&1 ) &
  job=$! mx=0 nf=0
  while kill -0 $job 2>/dev/null; do if c=$(gpucount); then [ "$c" -gt "$mx" ] && mx=$c; else nf=$((nf + 1)); fi; sleep 5; done
  wait $job; rc=$?; echo "done $nf $mx" > $r/gpu_watch.txt
  grep -q "Thomas の配列の並び: LAYOUT2 (既定" $r/forge_run.log && lay=既定LAYOUT2 || lay=違う
  echo "$r rc=$rc $v 並び $lay GPU 見張り (失敗 $nf・最大 $mx 本)" >> $LOG
  python3 scr_check.py $r 570999 121 5 >> $LOG 2>&1
  rm -f $r/nozzle.h5 $r/res_0.h5
done
echo "== 終了 $(date -Is)" >> $LOG
touch time5.done
