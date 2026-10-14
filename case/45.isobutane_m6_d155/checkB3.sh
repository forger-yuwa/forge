#!/bin/bash
# (B3: 1 ライン 1 スレッドの lu5 をレジスタに置いた版、2026-10-09) plan time_integration-line-implicit-speed §6 案 B の判定と、plan time_integration-line-viscous-jacobian §6.2 のライン行列の書き出し (2026-10-09)。
# (1) run_0287: 新 (lineB3_fp64、874ae90a + typedef double) で run_0223 の設定 20 step、FORGE_LINE_COMPARE=1 (v2 と v3 を同じ入力で解いて差を出す)
# (3) run_0293〜0295: 性能 (専有 GPU、1000 step × 3、案 A の 34.98 ms/step と比べる)
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
REF="--limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
C223="--line dir --itj 5 --cap 50"
use_new() { export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineB3_fp64; }
use_old() { export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64; }
go() {  # $1 run, $2 steps, 以降 prep の追加引数
  local r=$1 st=$2; shift 2
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps $st --cfl 4 --out $st $REF "$@" > /dev/null && rm -f $r/nozzle.msh \
    && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  echo "$r rc=$? $(grep -o 'Time = .*' $r/forge_run.log | tail -1) $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> checkB3.log
}
use_new; export FORGE_LINE_COMPARE=1; go run_0292_cmpB3 20 $C223; unset FORGE_LINE_COMPARE
grep "\[lineCompare\]" run_0292_cmpB3/forge_run.log > checkB3_compare.txt
for k in a b c; do
  if [ -n "$(pgrep -x forge)" ]; then echo "性能: 他の forge が走っている — 止める" >> checkB3.log; break; fi
  use_new; go run_02$((93 + $(printf '%d' "'$k") - 97))_timeB3_new_$k 1000 $C223
done
touch checkB3.done
