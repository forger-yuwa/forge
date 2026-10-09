#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.4 (2026-10-09): 同じ新バイナリ (lineE_fp64、6d49738e + typedef double) で、run_0183 の res_100000 から
# キー 5・方向別・上限なし・cfl 4 の A = lineViscCoupling 2 (run_0300) と B = 3 (run_0301) を 2000 step・200 ごと (並行)、
# その後 B の 1 step 目のライン行列の書き出し (run_0302、列 12・33・65・1640・2183 のライン)。forge は cold_cfl.py run 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineE_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR
C="--cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --extra res_ro,volume"
for v in 2 3; do
  r=$( [ $v = 2 ] && echo run_0300_vn1b_lvc2 || echo run_0301_vn1b_lvc3 )
  python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 2000 --out 200 $C --lvc $v > /dev/null && rm -f $r/nozzle.msh
  (python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> vn1b.log) &
done
wait
mkdir -p run_0302_dump_lvc3_linedump
export FORGE_LINE_DUMP_DIR=$PWD/run_0302_dump_lvc3_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$((12*121+120)),$((33*121+120)),$((65*121+120)),$((1640*121+120)),$((2183*121+120))
r=run_0302_dump_lvc3
python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r --steps 1 --out 1 $C --lvc 3 > /dev/null && rm -f $r/nozzle.msh && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
echo "$r rc=$?" >> vn1b.log
unset FORGE_LINE_DUMP_DIR FORGE_LINE_DUMP_CALL FORGE_LINE_DUMP_NODES
python3 linedump_analyze.py run_0302_dump_lvc3_linedump > vn1b_linedump.txt 2>&1
touch vn1b.done
