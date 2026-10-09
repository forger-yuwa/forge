#!/bin/bash
# (B2: 代入を行の分担に書き直した版、2026-10-09) plan time_integration-line-implicit-speed §6 案 B の判定と、plan time_integration-line-viscous-jacobian §6.2 のライン行列の書き出し (2026-10-09)。
# (1) run_0287: 新 (lineB2_fp64、874ae90a + typedef double) で run_0223 の設定 20 step、FORGE_LINE_COMPARE=1 (v2 と v3 を同じ入力で解いて差を出す)
# (2) run_0288: 値 2・上限なし 1 step、FORGE_LINE_DUMP_DIR (列 12・33・65・1640・2183 の壁の節点を含むライン)
# (3) run_0289〜0291: 性能 (専有 GPU、1000 step × 3、案 A の 34.98 ms/step と比べる)
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
REF="--limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
C223="--line dir --itj 5 --cap 50"
use_new() { export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineB2_fp64; }
use_old() { export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64; }
go() {  # $1 run, $2 steps, 以降 prep の追加引数
  local r=$1 st=$2; shift 2
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps $st --cfl 4 --out $st $REF "$@" > /dev/null && rm -f $r/nozzle.msh \
    && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  echo "$r rc=$? $(grep -o 'Time = .*' $r/forge_run.log | tail -1) $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> checkB2.log
}
use_new; export FORGE_LINE_COMPARE=1; go run_0287_cmpB2 20 $C223; unset FORGE_LINE_COMPARE
grep "\[lineCompare\]" run_0287_cmpB2/forge_run.log > checkB2_compare.txt
mkdir -p run_0288_dump_lvc2_linedump
use_new; export FORGE_LINE_DUMP_DIR=$PWD/run_0288_dump_lvc2_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$((12*121+120)),$((33*121+120)),$((65*121+120)),$((1640*121+120)),$((2183*121+120))
go run_0288_dump_lvc2 1 --line dir --itj 5 --lvc 2; unset FORGE_LINE_DUMP_DIR FORGE_LINE_DUMP_CALL FORGE_LINE_DUMP_NODES
python3 linedump_analyze.py run_0288_dump_lvc2_linedump > checkB2_linedump.txt 2>&1
for k in a b c; do
  if [ -n "$(pgrep -x forge)" ]; then echo "性能: 他の forge が走っている — 止める" >> checkB2.log; break; fi
  use_new; go run_02$((89 + $(printf '%d' "'$k") - 97))_timeB2_new_$k 1000 $C223
done
touch checkB2.done
