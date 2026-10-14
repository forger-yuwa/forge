#!/bin/bash
# plan time_integration-line-implicit-speed §4.3 案 A の確認 (2026-10-09): ライン上の節点の対角を storeLU の sweep 以外で組まない変更 (4d394a71) が
# 数値を変えないこと。run_0223 と同じ設定 (ライン + 方向別 + キー 5 + 上限 50) で run_0183 の res_100000 から、
#   20 step: 旧 (thermjac_cap_fp64) 2 本・新 (linespeed_fp64) 2 本、1 step: 旧 1 本・新 1 本
# を回し、ab_compare.py で「旧 × 新」の差を「同じバイナリの再実行」の差と比べる。forge は cold_cfl.py run (= run_case.sh) 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
OLD=$(cat prof_target.txt)
NEW=$(cat prof_target_linespeed.txt)
n=254
one() {  # $1 = old|new, $2 = 接尾辞, $3 = step 数
  local r=$(printf "run_%04d_abA_%s_%s" $n $1 $2); n=$((n+1))
  if [ $1 = old ]; then export FORGE_BIN=$OLD COLD_ALT_BINARY=thermjac_cap_fp64; else export FORGE_BIN=$NEW COLD_ALT_BINARY=linespeed_fp64; fi
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps $3 --cfl 4 --out $3 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext \
      --line dir --itj 5 --cap 50 > /dev/null && rm -f $r/nozzle.msh && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  echo "$r rc=$? $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> ab_lineA.log
}
one old s20a 20; one old s20b 20; one new s20a 20; one new s20b 20; one old s1 1; one new s1 1
python3 ab_compare.py abA > /dev/null 2>&1
touch ab_lineA.done
