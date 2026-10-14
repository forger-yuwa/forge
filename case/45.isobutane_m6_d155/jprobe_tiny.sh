#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.13 の事後の探索 (判定には使わない、2026-10-09): ε = 1e-2〜1e-1 では残差の応答に折れが多数入り差分が定まらなかったので、
# 折れの影響が小さい ε = 1e-6 (ε/2 = 5e-7) で S0 の p5、S1 の p7・p5 も測る (S0 の p7 は run_0322_jp_s0p7n_* で取得済み)。run_0322_jp_<s><p>n_*。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineG_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --lvc 3"
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,wi_eheat,wi_ework
LOG=$PWD/jprobe.log
fail() { echo "失敗: $*" | tee -a $LOG; touch jprobe_tiny.done; exit 1; }
probe() {
  local r=$1 ff=$2
  if [ -f _jprobe/npz/$r.npz ] && [ -d $r ]; then echo "$r は抜き出し済み — 飛ばす" >> $LOG; return 0; fi
  python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 $B --extra $EXTRA --field-from $ff > _jprobe/logs/$r.prep.log 2>&1 || return 1
  rm -f $r/nozzle.msh
  local t0; t0=$(date +%s)
  ( export FORGE_WI_FORCE_DIAG=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc" >> $LOG
  [ $rc -eq 0 ] || return 1
  [ -f $r/res_1.h5 ] && [ $(stat -c %Y $r/res_1.h5) -ge $t0 ] || { echo "$r: res_1 が無いか古い" >> $LOG; return 1; }
  python3 jprobe.py extract $r >> $LOG 2>&1 || return 1
  rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1.h5
}
family4() {   # family4 <run 接頭辞> <tag> <state> <dump> <eps>
  local pre=$1 tag=$2 st=$3 dump=$4 eps=$5
  [ -f _jprobe/${tag}_fields.json ] || python3 jprobe.py fields $st $dump $tag --eps $eps --no-pp >> $LOG 2>&1 || fail "fields $tag"
  for k in pe me ph mh; do
    probe ${pre}_${k} _jprobe/${tag}_${k} || fail "${pre}_${k}"
    rm -f _jprobe/${tag}_${k}/res_0.h5
  done
}
echo "== 探索 (ε = 1e-6) 開始 $(date -Is)" >> $LOG
S0=$SRC/res_100000.h5; S1=run_0316_jp_m7_s1/res_20.h5
D07=run_0313_e1bdump_m7_linedump; D05=run_0315_e1bdump_m5_linedump; D17=run_0318_jpdump_s1_m7_linedump; D15=run_0319_jpdump_s1_m5_linedump
family4 run_0322_jp_s0p5n s0p5n $S0 $D05 1e-6
family4 run_0322_jp_s1p7n s1p7n $S1 $D17 1e-6
family4 run_0322_jp_s1p5n s1p5n $S1 $D15 1e-6
OPS0="7=$D07,5=$D05"; OPS1="7=$D17,5=$D15"
for st in s0 s1; do for t in p7 p5; do
  if [ $st = s0 ]; then ops=$OPS0; d=$([ $t = p7 ] && echo $D07 || echo $D05); q="q0=run_0317_jp_s0_q0,q0b=run_0317_jp_s0_q0b"; pp=run_0317_jp_s0${t}_pp
  else ops=$OPS1; d=$([ $t = p7 ] && echo $D17 || echo $D15); q="q0=run_0320_jp_s1_q0,q0b=run_0320_jp_s1_q0b"; pp=run_0320_jp_s1${t}_pp; fi
  python3 jprobe.py compare $d ${st}${t}n --eps 1e-6 --ops $ops \
    --runs "$q,pe=run_0322_jp_${st}${t}n_pe,me=run_0322_jp_${st}${t}n_me,ph=run_0322_jp_${st}${t}n_ph,mh=run_0322_jp_${st}${t}n_mh,pp=$pp" > _jprobe/${st}${t}n_compare.txt 2>&1 || fail "compare ${st}${t}n"
done; done
for st in s0 s1; do for t in p7 p5; do      # ε = 1e-2・1e-1 の比較も実補正の予測差を足して出し直す
  if [ $st = s0 ]; then ops=$OPS0; d=$([ $t = p7 ] && echo $D07 || echo $D05); q="q0=run_0317_jp_s0_q0,q0b=run_0317_jp_s0_q0b"; pre=run_0317_jp_s0${t}; pp=run_0317_jp_s0${t}_pp
  else ops=$OPS1; d=$([ $t = p7 ] && echo $D17 || echo $D15); q="q0=run_0320_jp_s1_q0,q0b=run_0320_jp_s1_q0b"; pre=run_0320_jp_s1${t}; pp=run_0320_jp_s1${t}_pp; fi
  python3 jprobe.py compare $d ${st}${t} --eps 1e-2 --ops $ops --runs "$q,pe=${pre}_pe,me=${pre}_me,ph=${pre}_ph,mh=${pre}_mh,pp=$pp" > _jprobe/${st}${t}_compare.txt 2>&1 || fail "compare ${st}${t}"
  python3 jprobe.py compare $d ${st}${t}r --eps 1e-1 --ops $ops \
    --runs "$q,pe=run_0321_jp_${st}${t}r_pe,me=run_0321_jp_${st}${t}r_me,ph=run_0321_jp_${st}${t}r_ph,mh=run_0321_jp_${st}${t}r_mh,pp=$pp" > _jprobe/${st}${t}r_compare.txt 2>&1 || fail "compare ${st}${t}r"
done; done
echo "== 探索 (ε = 1e-6) 終了 $(date -Is)" >> $LOG
touch jprobe_tiny.done
