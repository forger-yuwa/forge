#!/bin/bash
# plan time_integration-line-implicit-speed §6.4 (2026-10-09、codex plan-3 の反映後): Thomas を float に (FORGE_LINE_F32=1/2) の照合・性能・途中の到達。バイナリ lineJ。
# 設定は run_0223 と同じ (値 0・方向別・キー 5・上限 50・cfl 4)、初期場は run_0183 の res_100000。forge・判定器の失敗では止める。古い判定の JSON は使わない。
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
rm -f _band_ab/cold_pair/f32_time_judge.json _band_ab/cold_pair/f32_qual_judge_*.json
# (1) 全ラインの照合
declare -A CMP
for m in 1 2; do
  if [ $m = 1 ]; then arm=F32c; r=run_0339_f32c_cmp; else arm=F32cs; r=run_0340_f32cs_cmp; fi
  prep $r --steps 20 --out 20 $B || fail "prep $r"
  run $r FORGE_LINE_COMPARE=1 FORGE_LINE_F32=$m || fail $r
  python3 line_cmp_judge.py $r --arm $arm --eta-limit 1e-5 --dq-limit -1 --factors 20 --solves 100 --lines 4719 >> $LOG 2>&1; CMP[$arm]=$?
  echo "(1) $arm 判定 rc=${CMP[$arm]} (0 PASS / 1 FAIL / 2 INDETERMINATE)" >> $LOG
  rm -f $r/nozzle.h5 $r/res_20.h5
done
# (2) 性能 (3 組、従来・F32c・F32cs の順)。律速の切り分けの記録のため、(1) の結果によらず 3 腕とも測る
for i in 1 2 3; do
  for arm in LU F32c F32cs; do
    r=run_034${i}_time3_${arm}_$i
    prep $r --steps 1000 --out 1000 $B || fail "prep $r"
    echo "$r 投入前の GPU の計算プロセス: $(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader | wc -l) 本、forge $(pgrep -x forge | wc -l) 本" >> $LOG
    case $arm in LU) run $r || fail $r;; F32c) run $r FORGE_LINE_F32=1 || fail $r;; F32cs) run $r FORGE_LINE_F32=2 || fail $r;; esac
    rm -f $r/nozzle.h5 $r/res_1000.h5
  done
done
python3 f32_time_judge.py >> $LOG 2>&1 || fail "(2) の判定器 (rc $?)"
# (3) 途中の到達 ((1) が PASS かつ (2) で速い腕だけ)
SEL=$(python3 -c "
import json;d=json.load(open('_band_ab/cold_pair/f32_time_judge.json'))
import sys; cmp={'F32c':${CMP[F32c]},'F32cs':${CMP[F32cs]}}
print(' '.join(a for a in ('F32c','F32cs') if d[a]['verdict']=='速い' and cmp[a]==0))") || fail "(3) の選別"
echo "(3) 対象の腕: ${SEL:-なし}" >> $LOG
if [ -n "$SEL" ]; then
  for arm in LU $SEL; do
    r=$(case $arm in LU) echo run_0344_qualLU2;; F32cs) echo run_0345_qualF32cs;; F32c) echo run_0346_qualF32c;; esac)
    prep $r --steps 10000 --out 1000 $B --extra res_ro || fail "prep $r"
    case $arm in LU) run $r || fail $r;; F32c) run $r FORGE_LINE_F32=1 || fail $r;; F32cs) run $r FORGE_LINE_F32=2 || fail $r;; esac
    python3 cold_series.py $r >> $LOG 2>&1 || fail "$r cold_series (rc $?)"
    python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$? $(grep -o 'NOT CONVERGED\|DIVERGED\|PASS' $r/CONVERGENCE_VERDICT.txt | head -1)" >> $LOG
    rm -f $r/nozzle.h5
  done
  for arm in $SEL; do
    r=$([ $arm = F32cs ] && echo run_0345_qualF32cs || echo run_0346_qualF32c)
    python3 inv_qual_judge.py run_0344_qualLU2 $r f32_qual_judge_$arm.json >> $LOG 2>&1 || fail "(3) $arm の判定器 (rc $?)"
    echo "(3) $arm 判定済み → _band_ab/cold_pair/f32_qual_judge_$arm.json" >> $LOG
  done
fi
echo "== 終了 $(date -Is)" >> $LOG
touch jf32.done
