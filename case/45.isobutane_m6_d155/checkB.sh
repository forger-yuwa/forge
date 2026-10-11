#!/bin/bash
# plan time_integration-line-implicit-speed §6 案 B の判定と、plan time_integration-line-viscous-jacobian §6.2 のライン行列の書き出し (2026-10-09)。
# 新バイナリ = 64ed2cd6 + typedef double (lineB_fp64、~/forge-linevisc-fp64 を上書きしてビルド、B の並列 Thomas)、旧 = linespeed_fp64 (d8b06ebc、案 A まで)。run_0183 の res_100000 から、forge は cold_cfl.py run 経由。
#   (1) run_0275: 新、run_0223 の設定 (ライン + 方向別 + キー 5 + 上限 50) 20 step、FORGE_LINE_COMPARE=1 (同じ入力で v2 と v3 を両方解いて差を出す)
#   (2) run_0276〜0281: 短期の一致 (旧・新 20 step × 2、1 step × 1)
#   (3) run_0282: 新、値 2・上限なし 1 step、FORGE_LINE_DUMP_DIR (列 12・33・65・1640・2183 の壁の節点を含むライン)
#   (4) run_0283〜0285: 新の性能 (専有 GPU、1000 step × 3。案 A の旧 34.98 ms/step と比べる)
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
REF="--limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
C223="--line dir --itj 5 --cap 50"
use_new() { export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineB_fp64; }
use_old() { export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64; }
go() {  # $1 run, $2 steps, 以降 prep の追加引数
  local r=$1 st=$2; shift 2
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps $st --cfl 4 --out $st $REF "$@" > /dev/null && rm -f $r/nozzle.msh \
    && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  echo "$r rc=$? $(grep -o 'Time = .*' $r/forge_run.log | tail -1) $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> checkB.log
}
use_new; export FORGE_LINE_COMPARE=1; go run_0275_cmpB 20 $C223; unset FORGE_LINE_COMPARE
grep "\[lineCompare\]" run_0275_cmpB/forge_run.log > checkB_compare.txt
use_old; go run_0276_abB_old_s20a 20 $C223; go run_0277_abB_old_s20b 20 $C223
use_new; go run_0278_abB_new_s20a 20 $C223; go run_0279_abB_new_s20b 20 $C223
use_old; go run_0280_abB_old_s1 1 $C223
use_new; go run_0281_abB_new_s1 1 $C223
python3 ab_compare.py abB > /dev/null 2>&1
mkdir -p run_0282_dump_lvc2_linedump
use_new; export FORGE_LINE_DUMP_DIR=$PWD/run_0282_dump_lvc2_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$((12*121+120)),$((33*121+120)),$((65*121+120)),$((1640*121+120)),$((2183*121+120))
go run_0282_dump_lvc2 1 --line dir --itj 5 --lvc 2; unset FORGE_LINE_DUMP_DIR FORGE_LINE_DUMP_CALL FORGE_LINE_DUMP_NODES
python3 linedump_analyze.py run_0282_dump_lvc2_linedump > checkB_linedump.txt 2>&1
for k in a b c; do
  if [ -n "$(pgrep -x forge)" ]; then echo "性能: 他の forge が走っている — 止める" >> checkB.log; break; fi
  use_new; go run_028$((3 + $(printf '%d' "'$k") - 97))_timeB_new_$k 1000 $C223
done
touch checkB.done
