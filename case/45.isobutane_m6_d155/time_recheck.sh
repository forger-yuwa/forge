#!/bin/bash
# 計時の取り直し (2026-10-09): 案 A のバイナリ (linespeed_fp64) を 1 回 (run_0296) と、B3 (lineB3_fp64) を nsys で 300 step (run_0297)。専有 GPU。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
C="--cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50"
export FORGE_BIN=$(cat prof_target_linespeed.txt) COLD_ALT_BINARY=linespeed_fp64
r=run_0296_timeA_recheck
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 1000 --out 1000 $C > /dev/null && rm -f $r/nozzle.msh && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
echo "$r $(grep -o 'Time = .*' $r/forge_run.log | tail -1)" >> time_recheck.log
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineB3_fp64
r=run_0297_prof_lineB3
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 300 --out 300 $C > /dev/null && rm -f $r/nozzle.msh
mv prof_target.txt prof_target.txt.orig; cp prof_target_linevisc.txt prof_target.txt
FORGE_BIN=$PWD/prof_wrap.sh bash $HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh $r > $r/run_case_stdout.log 2>&1
mv prof_target.txt.orig prof_target.txt
(cd $r && nsys stats --report cuda_gpu_kern_sum --format csv -o kern nsys_prof.nsys-rep > /dev/null 2>&1)
echo "$r $(grep -o 'Time = .*' $r/forge_run.log | tail -1)" >> time_recheck.log
touch time_recheck.done
