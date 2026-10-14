#!/bin/bash
# V-n1 の B の発散の場所 (2026-10-09、計測のみ): 上限なし・方向別・キー 5 で lineViscCoupling 0 / 2 を 3 step・毎 step 出力 (run_0266・0267)。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=linevisc_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
for lvc in 0 2; do
  r=$( [ $lvc = 0 ] && echo run_0266_diag_s3_lvc0 || echo run_0267_diag_s3_lvc2 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 3 --cfl 4 --out 1 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext \
      --line dir --itj 5 --lvc $lvc --extra res_ro,volume > /dev/null && rm -f $r/nozzle.msh
  (python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> vn1_s3.log) &
done
wait; touch vn1_s3.done
