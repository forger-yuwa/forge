#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.8 改訂 E1 (2026-10-09): バイナリ lineF (4fdc2c0e + typedef double)。run_0183 の res_100000 から、キー 5・方向別・上限なし・cfl 4。
# (1) 1 step 目のライン行列の書き出し: 値 3 マスク 7 (run_0307)・5 (run_0308)・0 (run_0309)、値 0 + キー 7 (run_0310) → e1_compare.py
# (2) A = 値 3 マスク 7 (run_0305_e1_m7)、B = マスク 5 (run_0306_e1_m5) を逐次に 2000 step・200 ごと、序盤 200 step は 112 節点の帳簿。forge は cold_cfl.py run 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineF_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS
B="--cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --extra res_ro,volume"
dump() {  # $1 run, $2 マスク (空なら設定しない), 以降 prep の追加
  local r=$1 m=$2; shift 2
  mkdir -p ${r}_linedump
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 1 --out 1 $B "$@" > /dev/null && rm -f $r/nozzle.msh
  ( export FORGE_LINE_DUMP_DIR=$PWD/${r}_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=1572,4113,7985,198560,264263
    [ -n "$m" ] && export FORGE_LVC_TERMS=$m
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> e1.log )
}
dump run_0307_e1dump_m7 7 --itj 5 --lvc 3
dump run_0308_e1dump_m5 5 --itj 5 --lvc 3
dump run_0309_e1dump_m0 0 --itj 5 --lvc 3
dump run_0310_e1dump_v0k7 "" --itj 7
python3 e1_compare.py run_0307_e1dump_m7_linedump run_0308_e1dump_m5_linedump run_0309_e1dump_m0_linedump run_0310_e1dump_v0k7_linedump > e1_compare.txt 2>&1
for m in 7 5; do
  r=$( [ $m = 7 ] && echo run_0305_e1_m7 || echo run_0306_e1_m5 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 2000 --out 200 $B --itj 5 --lvc 3 > /dev/null && rm -f $r/nozzle.msh
  ( export FORGE_LVC_TERMS=$m FORGE_DUMP_LEDGER=$PWD/$r/ledger.csv FORGE_DUMP_LEDGER_NODES=1572,1571,1570,1569,1568,1567,1566,1565,1564,1563,1562,1561,1560,1559,1558,1557,4113,4112,4111,4110,4109,4108,4107,4106,4105,4104,4103,4102,4101,4100,4099,4098,4960,4959,4958,4957,4956,4955,4954,4953,4952,4951,4950,4949,4948,4947,4946,4945,6170,6169,6168,6167,6166,6165,6164,6163,6162,6161,6160,6159,6158,6157,6156,6155,7985,7984,7983,7982,7981,7980,7979,7978,7977,7976,7975,7974,7973,7972,7971,7970,198560,198559,198558,198557,198556,198555,198554,198553,198552,198551,198550,198549,198548,198547,198546,198545,264263,264262,264261,264260,264259,264258,264257,264256,264255,264254,264253,264252,264251,264250,264249,264248 FORGE_DUMP_LEDGER_CALLS=200
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> e1.log )
done
touch e1.done
