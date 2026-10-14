#!/bin/bash
# ライン陰解法の速度の内訳 (2026-10-09、ユーザ指示「LINE にすると計算速度がほぼ半減するのは課題、改善策を調べて」)。
# 同じ新バイナリ (thermjac_cap_fp64) で run_0183 の res_100000 から 500 step、GPU はこの 1 本だけで:
#   run_0249 point cfl 4 / run_0250 ライン + 方向別 + キー 5 + 上限 50 (run_0223 の設定) / run_0251 ラインだけ (方向別なし)
# run_case.sh 経由で、FORGE_BIN に prof_wrap.sh (nsys) を渡す。FORGE_PROFILE=1 で区間別の計時も出す。
cd "$(dirname "$0")"
export COLD_ALT_BINARY=thermjac_cap_fp64 FORGE_CUDA_BLOCKSIZE=128 FORGE_PROFILE=1
export FORGE_BIN=$(cat prof_target.txt)                      # prep のバイナリ照合用
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
common="--steps 500 --cfl 4 --out 500 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext run_0249_prof_point $common > /dev/null
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext run_0250_prof_linedir_tj5_cap50 $common --line dir --itj 5 --cap 50 > /dev/null
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext run_0251_prof_lineonly $common --line only > /dev/null
for r in run_0249_prof_point run_0250_prof_linedir_tj5_cap50 run_0251_prof_lineonly; do
  rm -f $r/nozzle.msh
  FORGE_BIN=$PWD/prof_wrap.sh bash $RC $r > $r/run_case_stdout.log 2>&1
  (cd $r && nsys stats --report cuda_gpu_kern_sum --format csv -o kern nsys_prof.nsys-rep > /dev/null 2>&1 \
     || nsys stats --report gpukernsum --format csv -o kern nsys_prof.nsys-rep > /dev/null 2>&1)
  echo "$r $(grep -o 'Time = .*' $r/forge_run.log | tail -1) $(ls $r/kern*.csv 2>/dev/null)" >> prof_runs.log
done
touch prof_runs.done
