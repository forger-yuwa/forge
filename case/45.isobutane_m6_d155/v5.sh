#!/bin/bash
# plan architecture-float-state-double-geometry §6 の V5 (事前登録): 専有の GPU で、V4 と同じバイナリ (56b1a6fe、診断は off) の
# FP64 (D64) と float (F32) を交互に 1000 step (D64 F32 D64 F32 D64 F32 D64)。B0 の構成、run_0183 の res_100000 から。
# 判定: float の 1 step の時間の平均 ≤ FP64 の平均の 0.70 倍。各 run の GPU の最大の使用メモリも記録する (§5.1 #12)。
set -uo pipefail
cd ~/forge-wallfit/case/45.isobutane_m6_d155
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
declare -A BIN=([D64]=$HOME/forge-fgeom5-fp64/solver_density_cuda/build/forge [F32]=$HOME/forge-fgeom5-f32/solver_density_cuda/build/forge)
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$PWD/v5.log
gpucount() { local o; o=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null) || return 1; echo "$o" | grep -c '[0-9]'; true; }
echo "== 開始 $(date -Is)" >> $LOG
for v in D64 F32; do echo "$v $(sha256sum ${BIN[$v]} | cut -c1-16)" >> $LOG; done
k=0
for v in D64 F32 D64 F32 D64 F32 D64; do
  k=$((k + 1)); r=run_$(printf %04d $((452 + k)))_v5_$(printf %02d $k)_$v
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps 1000 --out 1000 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) || { echo "prep $r 失敗 (中止)" >> $LOG; exit 1; }
  rm -f $r/nozzle.msh
  w=0; until [ "$(gpucount)" = "0" ]; do sleep 30; w=$((w + 30)); [ $w -ge 3600 ] && { echo "$r: 1 時間待っても GPU が空かない (中止)" >> $LOG; exit 1; }; done
  ( export FORGE_BIN=${BIN[$v]}; bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1 ) &
  job=$! mx=0 mem=0 nf=0
  while kill -0 $job 2>/dev/null; do
    if c=$(gpucount); then [ "$c" -gt "$mx" ] && mx=$c; else nf=$((nf + 1)); fi
    m=$(nvidia-smi --query-compute-apps=used_memory --format=csv,noheader,nounits 2>/dev/null | awk '{s+=$1} END {print s+0}'); [ "${m:-0}" -gt "$mem" ] && mem=$m
    sleep 3
  done
  wait $job; rc=$?
  t=$(grep -h "ms/step" $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | tail -1)
  echo "$r rc=$rc $v GPU の本数の最大 $mx (照会の失敗 $nf)、GPU メモリの最大 ${mem} MiB | $t" >> $LOG
  rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1000.h5
done
echo "== 終了 $(date -Is)" >> $LOG
python3 - "$LOG" <<'PY' >> $LOG 2>&1
import re, sys, statistics as st
L = open(sys.argv[1]).read().split("== 開始")[-1]
rows = re.findall(r"run_\d+_v5_\d+_(D64|F32) rc=(\d+) .*?本数の最大 (\d+).*?\(wall, (\d+) steps, ([0-9.]+) ms/step\)", L)
d = {"D64": [], "F32": []}; bad = []
for v, rc, mx, n, ms in rows:
    if rc != "0" or mx != "1" or n != "1000": bad.append((v, rc, mx, n))
    d[v].append(float(ms))
if bad or len(d["D64"]) != 4 or len(d["F32"]) != 3:
    print(f"V5: 判定不能 (専有でない・失敗・本数の不足: {bad}、D64 {len(d['D64'])} 本・F32 {len(d['F32'])} 本)")
else:
    r = st.mean(d["F32"]) / st.mean(d["D64"])
    print(f"V5: F32 {d['F32']} ms/step、D64 {d['D64']} ms/step → 比 {r:.3f} (基準 ≤ 0.70) → {'PASS' if r <= 0.70 else 'FAIL'}")
PY
touch v5.done
