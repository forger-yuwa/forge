#!/bin/bash
# plan time_integration-line-implicit-speed §6.9 (2026-10-10): 仕上げまでの総時間。(1) 1 step の専有時間 (P・L0・L5 × 3)、(2) L5 のライン段 → (3) point の段、(4) 判定。バイナリ lineK。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineK_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG FORGE_DIAG_FACE_H_DOUBLE FORGE_LINE_INV FORGE_LINE_PAR FORGE_LINE_F32 FORGE_LINE_LAYOUT
SRC=run_0183_ns_coldmesh_tw300_ext
REF="--cfl 4 --limiter-ref-from $SRC"
OPT_P=""; OPT_L0="--line dir --itj 5 --cap 50"; OPT_L5="--line dir --itj 5 --lvc 3"
LOG=$PWD/m9.log
fail() { echo "失敗: $*" | tee -a $LOG; touch m9.done; exit 1; }
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
echo "== 開始 $(date -Is)" >> $LOG
# (1) 1 step の専有時間
if [ ! -f _band_ab/cold_pair/m9_unit.json ]; then
  for i in 1 2 3; do
    for m in P L0 L5; do
      r=run_035$((4 + i))_m9unit_${m}_$i
      eval opt=\$OPT_$m
      prep $r --steps 1000 --out 1000 $REF $opt || fail "prep $r"
      nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l > $r/gpu_procs.txt
      t0=$(python3 -c "import time;print(time.time())")
      ( [ $m = L5 ] && export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); rc=$?
      python3 -c "import time;print(time.time()-$t0)" > $r/wall.txt
      echo "$r rc=$rc GPU の他のプロセス $(cat $r/gpu_procs.txt)" >> $LOG
      [ $rc -eq 0 ] || fail $r
      rm -f $r/nozzle.h5 $r/res_1000.h5
    done
  done
  python3 m9_unit.py >> $LOG 2>&1 || fail "m9_unit (rc $?)"
fi
# (2) L5 のライン段
r=run_0353_m9_L5
if [ ! -f $r/m9_watch.json ] || ! grep -q '"REACHED"\|"DIVERGED"\|"ENDED"' $r/m9_watch.json; then
  [ -d $r ] || prep $r --steps 200000 --out 5000 $REF $OPT_L5 --extra res_ro || fail "prep $r"
  ( export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r forge の終了 rc=$?" >> $LOG ) &
  sleep 60
  python3 m9_watch.py $r >> $LOG 2>&1 || fail "見張り $r (rc $?)"
  wait
fi
grep -q '"REACHED"' $r/m9_watch.json || { echo "$r: ライン段は水準に届かず ($(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])"))" >> $LOG; python3 m9_judge.py >> $LOG 2>&1; touch m9.done; exit 0; }
python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$?" >> $LOG
# (3) point の段
c=run_0354_m9_L5cut
if [ ! -f $c/m9_watch.json ] || ! grep -q '"REACHED"\|"DIVERGED"\|"ENDED"' $c/m9_watch.json; then
  [ -d $c ] || prep $c --steps 200000 --out 5000 $REF --extra res_ro --field-from $r || fail "prep $c"
  echo "$c の出発点: $(python3 -c "import json;print(json.load(open('$c/COLD_PAIR.json'))['parent_res'])") ($r の水準の出力)" >> $LOG
  ( python3 cold_cfl.py run $c > $c/cold_pair_run_stdout.log 2>&1; echo "$c forge の終了 rc=$?" >> $LOG ) &
  sleep 60
  python3 m9_watch.py $c >> $LOG 2>&1 || fail "見張り $c (rc $?)"
  wait
fi
python3 ../../solver_density_cuda/tools/check_convergence.py $c > $c/CONVERGENCE_VERDICT.txt 2>&1; echo "$c check_convergence rc=$?" >> $LOG
python3 m9_judge.py >> $LOG 2>&1; echo "(4) 判定 rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch m9.done
