#!/bin/bash
# §5.1 #11: Thomas の律速の計測 (ncu)。lineJ (既定の double の LU の経路)、run_0223 の設定、run_0183 の res_100000 から 6 step。
cd "$(dirname "$0")"
export COLD_ALT_BINARY=lineJ_fp64 FORGE_CUDA_BLOCKSIZE=128
export FORGE_BIN=$(cat prof_target_linevisc.txt)
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_LINE_INV FORGE_LINE_F32 FORGE_LINE_PAR
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
r=run_0347_ncu_lu
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 6 --out 6 --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50 > /dev/null || { echo "prep 失敗" > ncu_run.log; touch ncu_run.done; exit 1; }
rm -f $r/nozzle.msh
FORGE_BIN=$PWD/ncu_wrap.sh bash $RC $r > $r/run_case_stdout.log 2>&1; echo "$r rc=$?" > ncu_run.log
sudo chown -R ubuntu:ubuntu $r
rm -f $r/nozzle.h5 $r/res_6.h5
touch ncu_run.done
