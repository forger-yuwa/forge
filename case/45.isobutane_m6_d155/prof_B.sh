#!/bin/bash
# 案 B の遅さの内訳 (nsys、2026-10-09): 新バイナリ (lineB_fp64) で run_0223 の設定を 300 step。run_case.sh 経由で FORGE_BIN に prof_wrap.sh。
cd "$(dirname "$0")"
export COLD_ALT_BINARY=lineB_fp64 FORGE_CUDA_BLOCKSIZE=128
export FORGE_BIN=$(cat prof_target_linevisc.txt)
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
r=run_0286_prof_lineB
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 300 --cfl 4 --out 300 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50 > /dev/null
rm -f $r/nozzle.msh; cp prof_target_linevisc.txt prof_target.txt.B
mv prof_target.txt prof_target.txt.orig; cp prof_target_linevisc.txt prof_target.txt
FORGE_BIN=$PWD/prof_wrap.sh bash $HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh $r > $r/run_case_stdout.log 2>&1
mv prof_target.txt.orig prof_target.txt; rm -f prof_target.txt.B
(cd $r && nsys stats --report cuda_gpu_kern_sum --format csv -o kern nsys_prof.nsys-rep > /dev/null 2>&1)
touch prof_B.done
