#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.6 (2026-10-09): lineViscCoupling 0 に固定し、implicitThermalJacobian 5 (run_0303) / 7 (run_0304) だけを変える A/B。
# 同じバイナリ (lineE_fp64)、run_0183 の res_100000 から、方向別・上限なし・cfl 4、2000 step・200 ごと、序盤 200 step は 112 節点の帳簿を毎 step。forge は cold_cfl.py run 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineE_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
C="--steps 2000 --cfl 4 --out 200 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --extra res_ro,volume"
for k in 5 7; do
  r=$( [ $k = 5 ] && echo run_0303_wallA_tj5 || echo run_0304_wallB_tj7 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $C --itj $k > /dev/null && rm -f $r/nozzle.msh
  (export FORGE_DUMP_LEDGER=$PWD/$r/ledger.csv FORGE_DUMP_LEDGER_NODES=1572,1571,1570,1569,1568,1567,1566,1565,1564,1563,1562,1561,1560,1559,1558,1557,4113,4112,4111,4110,4109,4108,4107,4106,4105,4104,4103,4102,4101,4100,4099,4098,4960,4959,4958,4957,4956,4955,4954,4953,4952,4951,4950,4949,4948,4947,4946,4945,6170,6169,6168,6167,6166,6165,6164,6163,6162,6161,6160,6159,6158,6157,6156,6155,7985,7984,7983,7982,7981,7980,7979,7978,7977,7976,7975,7974,7973,7972,7971,7970,198560,198559,198558,198557,198556,198555,198554,198553,198552,198551,198550,198549,198548,198547,198546,198545,264263,264262,264261,264260,264259,264258,264257,264256,264255,264254,264253,264252,264251,264250,264249,264248 FORGE_DUMP_LEDGER_CALLS=200
   python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> wallAB.log) &
done
wait
touch wallAB.done
