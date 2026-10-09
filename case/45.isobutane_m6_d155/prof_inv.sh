#!/bin/bash
# §6.2 (3) で逆行列が遅かった (41.13 vs 34.63 ms/step) ので、カーネルの内訳を nsys で取る (事後の内訳、判定には使わない)。
# lineI で run_0183 の res_100000 から run_0223 の設定 300 step、従来 (run_0337_profLU) と逆行列 (run_0338_profINV)、GPU はこの 1 本だけで。
cd "$(dirname "$0")"
export COLD_ALT_BINARY=lineI_fp64 FORGE_CUDA_BLOCKSIZE=128
export FORGE_BIN=$(cat prof_target_linevisc.txt)
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_LINE_INV
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
for v in LU INV; do
  r=run_033$([ $v = LU ] && echo 7 || echo 8)_prof$v
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 300 --out 300 --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50 > /dev/null || { echo "$r prep 失敗" >> prof_inv.log; continue; }
  rm -f $r/nozzle.msh
  ( [ $v = INV ] && export FORGE_LINE_INV=1; FORGE_BIN=$PWD/prof_wrap_inv.sh bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$?" >> prof_inv.log
  (cd $r && nsys stats --report cuda_gpu_kern_sum --format csv -o kern nsys_prof.nsys-rep > /dev/null 2>&1)
  echo "$r $(grep -o 'Time = .*' $r/forge_run.log | tail -1)" >> prof_inv.log
  rm -f $r/nozzle.h5 $r/res_300.h5
done
touch prof_inv.done
