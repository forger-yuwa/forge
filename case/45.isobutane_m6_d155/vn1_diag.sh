#!/bin/bash
# V-n1 の B (lineViscCoupling 2) が 30 step で発散した切り分け (2026-10-09、計測のみ): 実装の誤りか、減衰を外したことによる不安定か。
#   run_0264: ライン + 値 2、方向別なし (Δτ は point と同じ)、キー 5、100 step・20 ごと
#   run_0265: ライン + 方向別 + 上限 50 + 値 2、キー 5 (run_0223 と同じ条件 + 値 2)、300 step・50 ごと
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=linevisc_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
c="--cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --itj 5 --lvc 2 --extra res_ro,volume"
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext run_0264_diag_lineonly_lvc2 --steps 100 --out 20 --line only $c > /dev/null
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext run_0265_diag_linedir_cap50_lvc2 --steps 300 --out 50 --line dir --cap 50 $c > /dev/null
for r in run_0264_diag_lineonly_lvc2 run_0265_diag_linedir_cap50_lvc2; do rm -f $r/nozzle.msh
  (python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> vn1_diag.log) &
done
wait; touch vn1_diag.done
