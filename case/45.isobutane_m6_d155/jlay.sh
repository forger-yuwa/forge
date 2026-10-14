#!/bin/bash
# plan time_integration-line-implicit-speed §6.7 (2026-10-10): Thomas の配列の並べ替え (FORGE_LINE_LAYOUT=1) の照合 (ビット一致)・性能・ncu の内訳。バイナリ lineK。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineK_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG FORGE_DIAG_FACE_H_DOUBLE FORGE_LINE_INV FORGE_LINE_PAR FORGE_LINE_F32 FORGE_LINE_LAYOUT
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50"
LOG=$PWD/jlay.log
fail() { echo "失敗: $*" | tee -a $LOG; touch jlay.done; exit 1; }
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
run() { local r=$1; shift; ( for kv in "$@"; do export "$kv"; done; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?; echo "$r rc=$rc $*" >> $LOG; return $rc; }
echo "== 開始 $(date -Is)" >> $LOG
rm -f _band_ab/cold_pair/lay_time_judge.json
# (1) 照合 (ビット一致)
r=run_0348_lay_cmp
prep $r --steps 20 --out 20 $B || fail "prep $r"
run $r FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=1 || fail $r
python3 line_cmp_judge.py $r --arm LAYOUT --eta-limit 1e-11 --bitwise --factors 20 --solves 100 --lines 4719 >> $LOG 2>&1; c1=$?
echo "(1) 判定 rc=$c1 (0 PASS / 1 FAIL / 2 INDETERMINATE)" >> $LOG
rm -f $r/nozzle.h5 $r/res_20.h5
[ $c1 -eq 0 ] || fail "(1) が PASS でない — 性能は測らない"
# (2) 性能 (3 組)
for i in 1 2 3; do
  for v in LU LAY; do
    r=run_03$((48 + i))_timeLAY_${v}_$i
    prep $r --steps 1000 --out 1000 $B || fail "prep $r"
    nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l > $r/gpu_procs.txt
    echo "$r 投入前の GPU の計算プロセス: $(cat $r/gpu_procs.txt) 本、forge $(pgrep -x forge | wc -l) 本" >> $LOG
    if [ $v = LAY ]; then run $r FORGE_LINE_LAYOUT=1 || fail $r; else run $r || fail $r; fi
    rm -f $r/nozzle.h5 $r/res_1000.h5
  done
done
python3 time_pairs_judge.py _band_ab/cold_pair/lay_time_judge.json LU LAYOUT run_0349_timeLAY_LU_1 run_0349_timeLAY_LAY_1 run_0350_timeLAY_LU_2 run_0350_timeLAY_LAY_2 run_0351_timeLAY_LU_3 run_0351_timeLAY_LAY_3 >> $LOG 2>&1 || fail "(2) の判定器 (rc $?)"
# (3) ncu の内訳 (並べ替え版、記録)
r=run_0352_ncu_lay
prep $r --steps 6 --out 6 $B || fail "prep $r"
( export FORGE_LINE_LAYOUT=1; FORGE_BIN=$PWD/ncu_wrap_lay.sh bash $HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh $r > $r/run_case_stdout.log 2>&1 ); rc=$?; echo "$r rc=$rc" >> $LOG
sudo chown -R ubuntu:ubuntu $r
[ $rc -eq 0 ] && [ -f $r/ncu_prof.ncu-rep ] || fail "$r: ncu が失敗したかレポートが無い"
grep -q "モード LAYOUT" $r/forge_run.log || fail "$r: 実効のモードが LAYOUT でない"
/usr/local/cuda/bin/ncu --import $r/ncu_prof.ncu-rep --page raw --csv --metrics gpu__time_duration.sum > $r/ncu_kernels.csv 2>/dev/null || fail "$r: レポートを読めない"
cnt=$(python3 -c "
import csv; rows=list(csv.DictReader(open('$r/ncu_kernels.csv')))[1:]
n=lambda s: sum(1 for x in rows if s in x['Kernel Name'])
print(n('lineThomasFactorL_d'), n('lineThomasSolveL_d'), n('implicit_defect_correction_block_d'))")
echo "$r ncu の対象: 分解・代入・block = $cnt" >> $LOG
[ "$cnt" = "1 5 5" ] || fail "$r: ncu の対象のカーネルの内訳が 1・5・5 でない ($cnt)"
rm -f $r/nozzle.h5 $r/res_6.h5
echo "== 終了 $(date -Is)" >> $LOG
touch jlay.done
