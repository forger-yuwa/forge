#!/bin/bash
# plan time_integration-line-implicit-speed §6.2 (2026-10-09): 逆行列の保存 (FORGE_LINE_INV=1) の照合・後退誤差・性能・同じ品質までの総時間。バイナリ lineI。
# 設定は run_0223 と同じ (値 0・方向別・キー 5・上限 50・cfl 4)、初期場は run_0183 の res_100000。各 forge の終了コードを確かめ、外れたら止める。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineI_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG FORGE_DIAG_FACE_H_DOUBLE FORGE_LINE_INV FORGE_LINE_PAR
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50"
LOG=$PWD/jinv.log
fail() { echo "失敗: $*" | tee -a $LOG; touch jinv.done; exit 1; }
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
run() {   # run <run> [env...]: env を付けて回し、終了コードを返す
  local r=$1; shift
  ( for kv in "$@"; do export "$kv"; done; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc $*" >> $LOG; return $rc
}
echo "== 開始 $(date -Is)" >> $LOG
# (1) 全ラインの照合
r=run_0324_invcmp
prep $r --steps 20 --out 20 $B || fail "prep $r"
run $r FORGE_LINE_COMPARE=1 FORGE_LINE_INV=1 || fail $r
python3 inv_cmp_judge.py $r >> $LOG 2>&1; c1=$?; echo "(1) rc=$c1" >> $LOG
rm -f $r/nozzle.h5
# (2) 後退誤差 (書き出し)
NODES=1572,4113,5321,7985,198560,264263
for spec in "run_0325_invdump_c1 1 1" "run_0326_invdump_c20 20 1" "run_0327_ludump_c1 1 0" "run_0328_ludump_c20 20 0"; do
  set -- $spec; r=$1; call=$2; inv=$3
  mkdir ${r}_linedump || fail "${r}_linedump が既にある"
  prep $r --steps $call --out $call $B || fail "prep $r"
  if [ $inv = 1 ]; then run $r FORGE_LINE_INV=1 FORGE_LINE_DUMP_DIR=$PWD/${r}_linedump FORGE_LINE_DUMP_CALL=$call FORGE_LINE_DUMP_NODES=$NODES || fail $r
  else run $r FORGE_LINE_DUMP_DIR=$PWD/${r}_linedump FORGE_LINE_DUMP_CALL=$call FORGE_LINE_DUMP_NODES=$NODES || fail $r; fi
  [ -f ${r}_linedump/dqnew_s4.f64 ] || fail "$r の書き出しが無い"
  rm -f $r/nozzle.h5
done
python3 line_eta.py run_0325_invdump_c1_linedump run_0326_invdump_c20_linedump >> $LOG 2>&1; echo "(2) 逆行列 (CPU の照合) rc=$?" >> $LOG
python3 line_eta.py run_0327_ludump_c1_linedump run_0328_ludump_c20_linedump >> $LOG 2>&1; echo "(2) 従来 (CPU の照合) rc=$?" >> $LOG
if [ $c1 -ne 0 ]; then echo "(1) が不合格または判定不能 — 性能と到達時間は回さない" >> $LOG; touch jinv.done; exit 0; fi
# (3) 性能 (交互、他の forge が無いことを記録)
for i in 1 2 3; do
  for v in LU INV; do
    n=$([ $v = LU ] && echo $((328 + i)) || echo $((331 + i)))
    r=run_0${n}_time${v}_$i
    prep $r --steps 1000 --out 1000 $B || fail "prep $r"
    echo "$r 投入前の GPU の計算プロセス: $(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader | wc -l) 本、forge $(pgrep -x forge | wc -l) 本" >> $LOG
    if [ $v = INV ]; then run $r FORGE_LINE_INV=1 || fail $r; else run $r || fail $r; fi
    echo "$r $(grep -h '^Time =' $r/forge_run.log | tail -1)" >> $LOG
    rm -f $r/nozzle.h5 $r/res_1000.h5
  done
done
python3 inv_time_judge.py run_0329_timeLU_1 run_0332_timeINV_1 run_0330_timeLU_2 run_0333_timeINV_2 run_0331_timeLU_3 run_0334_timeINV_3 >> $LOG 2>&1; echo "(3) rc=$?" >> $LOG
# (4) 同じ品質までの総時間
for v in LU INV; do
  r=run_033$([ $v = LU ] && echo 5 || echo 6)_qual$v
  prep $r --steps 10000 --out 500 $B --extra res_ro || fail "prep $r"
  if [ $v = INV ]; then run $r FORGE_LINE_INV=1 || fail $r; else run $r || fail $r; fi
  python3 cold_series.py $r >> $LOG 2>&1 || echo "$r cold_series rc=$?" >> $LOG
  python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$?" >> $LOG
  rm -f $r/nozzle.h5
done
python3 inv_qual_judge.py run_0335_qualLU run_0336_qualINV >> $LOG 2>&1; echo "(4) rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch jinv.done
