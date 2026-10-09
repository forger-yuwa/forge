#!/bin/bash
# plan time_integration-line-implicit-speed §6.4 (2026-10-09): Thomas を float に (FORGE_LINE_F32=1/2) の照合・性能・途中の到達。バイナリ lineJ。
# 設定は run_0223 と同じ (値 0・方向別・キー 5・上限 50・cfl 4)、初期場は run_0183 の res_100000。各 forge の終了コードを確かめ、外れたら止める。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineJ_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG FORGE_DIAG_FACE_H_DOUBLE FORGE_LINE_INV FORGE_LINE_PAR FORGE_LINE_F32
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50"
LOG=$PWD/jf32.log
fail() { echo "失敗: $*" | tee -a $LOG; touch jf32.done; exit 1; }
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
run() { local r=$1; shift; ( for kv in "$@"; do export "$kv"; done; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?; echo "$r rc=$rc $*" >> $LOG; return $rc; }
echo "== 開始 $(date -Is)" >> $LOG
# (1) 全ラインの照合
declare -A PASS
for m in 1 2; do
  r=$([ $m = 1 ] && echo run_0339_f32c_cmp || echo run_0340_f32cs_cmp)
  prep $r --steps 20 --out 20 $B || fail "prep $r"
  run $r FORGE_LINE_COMPARE=1 FORGE_LINE_F32=$m || fail $r
  python3 inv_cmp_judge.py $r --arm $arm --eta-limit 1e-5 --dq-limit -1 >> $LOG 2>&1; PASS[$arm]=$?; echo "(1) $arm rc=${PASS[$arm]}" >> $LOG
  rm -f $r/nozzle.h5 $r/res_20.h5
done
# (2) 性能 (3 組、従来・F32c・F32cs の順)
for i in 1 2 3; do
  for arm in LU F32c F32cs; do
    r=run_034$((i))_time3_${arm}_$i
    prep $r --steps 1000 --out 1000 $B || fail "prep $r"
    echo "$r 投入前の GPU の計算プロセス: $(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader | wc -l) 本、forge $(pgrep -x forge | wc -l) 本" >> $LOG
    case $arm in LU) run $r || fail $r;; F32c) run $r FORGE_LINE_F32=1 || fail $r;; F32cs) run $r FORGE_LINE_F32=2 || fail $r;; esac
    rm -f $r/nozzle.h5 $r/res_1000.h5
  done
done
python3 f32_time_judge.py >> $LOG 2>&1; echo "(2) rc=$?" >> $LOG
# (3) 途中の到達 ((1) 合格かつ (2) で速い腕だけ)
FAST=$(python3 -c "import json;d=json.load(open('_band_ab/cold_pair/f32_time_judge.json'));print(' '.join(a for a in ('F32c','F32cs') if d[a]['verdict']=='速い'))")
SEL=""; for arm in $FAST; do [ "${PASS[$arm]}" = 0 ] && SEL="$SEL $arm"; done
echo "(3) 対象の腕: ${SEL:-なし}" >> $LOG
if [ -n "$SEL" ]; then
  for arm in LU $SEL; do
    r=$(case $arm in LU) echo run_0344_qualLU2;; F32cs) echo run_0345_qualF32cs;; F32c) echo run_0346_qualF32c;; esac)
    prep $r --steps 10000 --out 1000 $B --extra res_ro || fail "prep $r"
    case $arm in LU) run $r || fail $r;; F32c) run $r FORGE_LINE_F32=1 || fail $r;; F32cs) run $r FORGE_LINE_F32=2 || fail $r;; esac
    python3 cold_series.py $r >> $LOG 2>&1 || echo "$r cold_series rc=$?" >> $LOG
    python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$?" >> $LOG
    rm -f $r/nozzle.h5
  done
  for arm in $SEL; do
    r=$([ $arm = F32cs ] && echo run_0345_qualF32cs || echo run_0346_qualF32c)
    python3 inv_qual_judge.py run_0344_qualLU2 $r >> $LOG 2>&1; cp _band_ab/cold_pair/inv_qual_judge.json _band_ab/cold_pair/f32_qual_judge_$arm.json; echo "(3) $arm rc=$?" >> $LOG
  done
fi
echo "== 終了 $(date -Is)" >> $LOG
touch jf32.done
