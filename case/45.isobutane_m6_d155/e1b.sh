#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.10 (2026-10-09): バイナリ lineG。run_0183 の res_100000 から値 3・キー 5・方向別・上限なし・cfl 4。
# (1) 1 step 目の書き出し: マスク 7 (run_0313)・15 (run_0314)・5 (run_0315) → e1b_gate.py (外れたら止める)
# (2) A = マスク 7 (run_0311_e1b_m7)、B = 15 (run_0312_e1b_m15) を逐次に 2000 step・200 ごと (場は監査まで消さない)、序盤 200 step は 112 節点の帳簿。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineG_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS
B="--cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --extra res_ro,volume --itj 5 --lvc 3"
dump() {
  local r=$1 m=$2
  mkdir -p ${r}_linedump
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 1 --out 1 $B > /dev/null && rm -f $r/nozzle.msh || return 1
  ( export FORGE_LVC_TERMS=$m FORGE_LINE_DUMP_DIR=$PWD/${r}_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=1572,4113,7985,198560,264263
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> e1b.log )
}
dump run_0313_e1bdump_m7 7 && dump run_0314_e1bdump_m15 15 && dump run_0315_e1bdump_m5 5 || { echo "書き出しの準備に失敗" >> e1b.log; touch e1b.done; exit 1; }
if ! python3 e1b_gate.py run_0313_e1bdump_m7_linedump run_0314_e1bdump_m15_linedump run_0315_e1bdump_m5_linedump > e1b_gate.txt 2>&1; then
  echo "ゲート不通過 — 比較の run を起動しない" >> e1b.log; touch e1b.done; exit 1
fi
for m in 7 15; do
  r=$( [ $m = 7 ] && echo run_0311_e1b_m7 || echo run_0312_e1b_m15 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 2000 --out 200 $B > /dev/null && rm -f $r/nozzle.msh || { echo "$r 準備に失敗" >> e1b.log; continue; }
  ( export FORGE_LVC_TERMS=$m FORGE_DUMP_LEDGER=$PWD/$r/ledger.csv FORGE_DUMP_LEDGER_NODES=1572,1571,1570,1569,1568,1567,1566,1565,1564,1563,1562,1561,1560,1559,1558,1557,4113,4112,4111,4110,4109,4108,4107,4106,4105,4104,4103,4102,4101,4100,4099,4098,4960,4959,4958,4957,4956,4955,4954,4953,4952,4951,4950,4949,4948,4947,4946,4945,6170,6169,6168,6167,6166,6165,6164,6163,6162,6161,6160,6159,6158,6157,6156,6155,7985,7984,7983,7982,7981,7980,7979,7978,7977,7976,7975,7974,7973,7972,7971,7970,198560,198559,198558,198557,198556,198555,198554,198553,198552,198551,198550,198549,198548,198547,198546,198545,264263,264262,264261,264260,264259,264258,264257,264256,264255,264254,264253,264252,264251,264250,264249,264248 FORGE_DUMP_LEDGER_CALLS=200
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> e1b.log )
done
touch e1b.done
