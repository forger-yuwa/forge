#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.15 (2026-10-09): TP の面エンタルピーの精度だけを変える A/B。バイナリ lineH (5197e00e + typedef double)。
# A = 切替なし、B = FORGE_DIAG_FACE_H_DOUBLE=1。S0・方向 p7 (run_0313 の書き出し)・ε = 1e-6 (ε/2 = 5e-7)、各側 q0・q0b・±ε・±ε/2 の 6 本 (run_0323_jph_{a,b}_*)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineH_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG FORGE_DIAG_FACE_H_DOUBLE
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --lvc 3"
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,wi_eheat,wi_ework
LOG=$PWD/jprobe.log
fail() { echo "失敗: $*" | tee -a $LOG; touch jprobe_fh.done; exit 1; }
probe() {   # probe <run> <field_from|-> <fh 0|1>
  local r=$1 ff=$2 fh=$3
  if [ -f _jprobe/npz/$r.npz ] && [ -d $r ]; then echo "$r は抜き出し済み — 飛ばす" >> $LOG; return 0; fi
  local opt=""; [ "$ff" != "-" ] && opt="--field-from $ff"
  python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 $B --extra $EXTRA $opt > _jprobe/logs/$r.prep.log 2>&1 || return 1
  rm -f $r/nozzle.msh
  local t0; t0=$(date +%s)
  ( export FORGE_WI_FORCE_DIAG=1; [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc fh=$fh" >> $LOG
  [ $rc -eq 0 ] || return 1
  [ -f $r/res_1.h5 ] && [ $(stat -c %Y $r/res_1.h5) -ge $t0 ] || { echo "$r: res_1 が無いか古い" >> $LOG; return 1; }
  if [ "$fh" = 1 ]; then grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log || { echo "$r: 切替の表示が無い" >> $LOG; return 1; }
  else grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log && { echo "$r: 切替が入っている" >> $LOG; return 1; }; fi
  python3 jprobe.py extract $r >> $LOG 2>&1 || return 1
  rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1.h5
}
echo "== 面エンタルピーの A/B 開始 $(date -Is)" >> $LOG
[ -f _jprobe/s0p7h_fields.json ] || python3 jprobe.py fields $SRC/res_100000.h5 run_0313_e1bdump_m7_linedump s0p7h --eps 1e-6 --no-pp >> $LOG 2>&1 || fail "fields s0p7h"
for side in a b; do
  fh=$([ $side = a ] && echo 0 || echo 1)
  probe run_0323_jph_${side}_q0 - $fh || fail "${side} q0"
  probe run_0323_jph_${side}_q0b - $fh || fail "${side} q0b"
  for k in pe me ph mh; do probe run_0323_jph_${side}_$k _jprobe/s0p7h_$k $fh || fail "${side} $k"; done
done
rm -f _jprobe/s0p7h_??/res_0.h5
for side in a b; do
  cp _jprobe/s0p7h_fields.json _jprobe/s0p7h${side}_fields.json
  python3 jprobe.py compare run_0313_e1bdump_m7_linedump s0p7h${side} --eps 1e-6 --ops 7=run_0313_e1bdump_m7_linedump,5=run_0315_e1bdump_m5_linedump \
    --runs "q0=run_0323_jph_${side}_q0,q0b=run_0323_jph_${side}_q0b,pe=run_0323_jph_${side}_pe,me=run_0323_jph_${side}_me,ph=run_0323_jph_${side}_ph,mh=run_0323_jph_${side}_mh,pp=run_0317_jp_s0p7_pp" \
    > _jprobe/s0p7h${side}_compare.txt 2>&1 || fail "compare $side"
done
python3 jprobe.py fhjudge >> $LOG 2>&1; echo "fhjudge rc=$?" >> $LOG
echo "== 面エンタルピーの A/B 終了 $(date -Is)" >> $LOG
touch jprobe_fh.done
