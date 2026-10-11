#!/bin/bash
# plan time_integration-line-implicit-speed §6 性能の判定 (案 A): 専有 GPU・FORGE_PROFILE なし・profiler なしで、run_0223 の設定 (ライン + 方向別 + キー 5 + 上限 50) を
# 旧 (thermjac_cap_fp64、35e498b1) と新 (linespeed_fp64、d8b06ebc = 案 A) で 3 回ずつ、point (新、ライン 0) を 1 回、run_0183 の res_100000 から 1000 step (出力は最後だけ)。
# 壁時計は forge_run.log の "Time = … (wall, 1000 steps, … ms/step)" (初期化は含まない)。forge は cold_cfl.py run 経由。同時に他の forge が無いことを確かめてから回す。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
n=268
one() {  # $1 = old|new|point, $2 = 接尾辞
  local r=$(printf "run_%04d_timeA_%s_%s" $n $1 $2); n=$((n+1))
  local c="--steps 1000 --cfl 4 --out 1000 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
  case $1 in
    old) export FORGE_BIN=$(cat prof_target.txt) COLD_ALT_BINARY=thermjac_cap_fp64; c="$c --line dir --itj 5 --cap 50";;
    new) export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64; c="$c --line dir --itj 5 --cap 50";;
    point) export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64;;
  esac
  if [ -n "$(pgrep -x forge)" ]; then echo "$r 他の forge が走っている — 計時を止める" >> time_lineA.log; return; fi
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $c > /dev/null && rm -f $r/nozzle.msh && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  echo "$r $(grep -o 'Time = .*' $r/forge_run.log | tail -1) $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> time_lineA.log
}
for k in a b c; do one old $k; one new $k; done
one point a
touch time_lineA.done
