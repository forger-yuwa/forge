#!/bin/bash
# ライン Thomas のカーネルの詳細 (ncu)。run_0250 と同じ設定で 20 step (ライン + 方向別 + キー 5 + 上限 50)。run_case.sh 経由。
cd "$(dirname "$0")"
export COLD_ALT_BINARY=thermjac_cap_fp64 FORGE_CUDA_BLOCKSIZE=128
export FORGE_BIN=$(cat prof_target.txt)
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES
r=run_0253_prof_ncu_linedir
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 20 --cfl 4 --out 20 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50 > /dev/null
rm -f $r/nozzle.msh
FORGE_BIN=$PWD/prof_wrap_ncu.sh bash $HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh $r > $r/run_case_stdout.log 2>&1
(cd $r && ncu --import ncu_prof.ncu-rep --page details --csv > ncu_details.csv 2>&1)
touch prof_ncu.done
