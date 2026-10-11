#!/bin/bash
# plan time_integration-line-viscous-jacobian §6 V-n1 (十分性の試験): run_0183 の res_100000 から、同じ新バイナリ (5ab83056 + typedef double) で
# lineViscCoupling 0 (A、run_0260) と 2 (B、run_0261) だけを変える。キー 5・方向別・上限なし・cfl 4、2000 step・200 ごと。forge は cold_cfl.py run 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=linevisc_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
for lvc in 0 2; do
  r=$( [ $lvc = 0 ] && echo run_0260_vn1_lvc0 || echo run_0261_vn1_lvc2 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 2000 --cfl 4 --out 200 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext \
      --line dir --itj 5 --lvc $lvc --extra res_ro,volume > /dev/null && rm -f $r/nozzle.msh
  (python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> vn1.log) &
done
wait
touch vn1.done
