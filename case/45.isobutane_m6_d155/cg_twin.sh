#!/bin/bash
# plan tooling-nozzle-core-grid §4.62 (事前登録、codex 諮問 2026-10-12 twin、ユーザ決定「ほい」): 双子の試験。Gref・同じ起点 (run_0619 の res_60000) で、
# 上限 50 と 1 の各腕に基準 Q₀ と摂動 Q₀ + δQ₀ (cg_twin_perturb.py、島の ρE に max 1e−7) を置き、差の成長を測る。B0 の他の設定・SST 更新はそのまま。
#   1 step の確認 (各上限): base1・dup1 (base1 と同じ入力、丸めの雑音の床)・pert1。1 外反復・outStepInterval 1。
#   本番 (各上限): base・pert を同時に 6,000 step・250 ごと。出力ごとに cg_twin_diff.py で差を取り、1000 の倍数・0・6000 以外の場を消す。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_twin.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRC=run_0619_cg_gref_b0_ctl; SEED=$SRC/res_60000.h5; LOG=$C45/cg_twin.log; AMP=1e-7; STEPS=6000; OUT=250
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge 10 ] || { echo "ディスクの空きが 10 GB 未満 (中止)" >> $LOG; exit 1; }
mk() {  # mk <run> <上限> <step 数> <出力の間隔> <摂動 0/1>
  local r=$1 cap=$2 n=$3 out=$4 pert=$5
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r; for f in $SRC/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|.done_*|*.png) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  sed -i "s/nStepOuter: [0-9]*}/nStepOuter: $n}/; s/outStepInterval: [0-9]*/outStepInterval: $out/; s/lineDtDirectionalCap: 50.0,/lineDtDirectionalCap: $cap,/; s/output: {level: 1, extraFields: \[[^]]*\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
  local d=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>') want=3; [ "$cap" != 50.0 ] && want=4
  [ "$d" = "$want" ] && grep -q "lineDtDirectionalCap: $cap," $r/solverConfig.yaml && grep -q "nStepOuter: $n}" $r/solverConfig.yaml && grep -q "volume, dt_local\]" $r/solverConfig.yaml \
    || { echo "$r: 設定の差が想定と違う ($d 行) (中止)" >> $LOG; return 1; }
  if [ $pert = 1 ]; then python3 cg_twin_perturb.py $r $AMP >> $LOG 2>&1 || { echo "$r: 摂動に失敗" >> $LOG; return 1; }; fi
  echo "$r: 上限 $cap・$n step・出力 $out ごと・摂動 $pert" >> $LOG
}
run() { ( export FORGE_BIN=$F7_64; bash $RC $C45/$1 > $1/run_case_stdout.log 2>&1; echo $? > $1/RUN_RC ); }
for cap in 50.0 1.0; do
  c=${cap%.0}; case $c in 50) n0=664 ;; 1) n0=669 ;; esac
  B1=run_0$((n0))_tw${c}_base1; D1=run_0$((n0+1))_tw${c}_dup1; P1=run_0$((n0+2))_tw${c}_pert1; BM=run_0$((n0+3))_tw${c}_base; PM=run_0$((n0+4))_tw${c}_pert
  mk $B1 $cap 1 1 0 && mk $D1 $cap 1 1 0 && mk $P1 $cap 1 1 1 || exit 1
  for r in $B1 $D1 $P1; do run $r; echo "$r 終了 RUN_RC=$(cat $r/RUN_RC)" >> $LOG; done
  python3 cg_twin_diff.py $B1 $D1 1 noise_cap$c && python3 cg_twin_diff.py $B1 $P1 0 step01_cap$c && python3 cg_twin_diff.py $B1 $P1 1 step01_cap$c || { echo "1 step の差の取り出しに失敗 (中止)" >> $LOG; exit 1; }
  mk $BM $cap $STEPS $OUT 0 && mk $PM $cap $STEPS $OUT 1 || exit 1
  ( run $BM ) & PB=$!; sleep 15; ( run $PM ) & PP=$!; echo "$BM・$PM 起動 $(date -Is)" >> $LOG
  python3 cg_twin_diff.py $BM $PM 0 main_cap$c 2>> $LOG || true
  while true; do
    for f in $BM/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f $BM/.tw_$n ] && continue; [ -f $PM/res_$n.h5 ] || continue
      if { [ -f $BM/res_$((n + OUT)).h5 ] && [ -f $PM/res_$((n + OUT)).h5 ]; } || { [ -f $BM/RUN_RC ] && [ -f $PM/RUN_RC ]; }; then
        python3 cg_twin_diff.py $BM $PM $n main_cap$c && touch $BM/.tw_$n
        if [ $((n % 1000)) -ne 0 ] && [ $n -ne $STEPS ]; then rm -f $BM/res_$n.h5 $BM/res_$n.xmf $PM/res_$n.h5 $PM/res_$n.xmf; fi
      fi
    done
    if [ -f $BM/RUN_RC ] && [ -f $PM/RUN_RC ]; then left=0; for f in $BM/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f $BM/.tw_$n ] && left=1; done; [ $left = 0 ] && break; fi
    [ "$(free_gb)" -lt 3 ] && { echo "ディスク 3 GB 未満 (中止)" >> $LOG; kill $PB $PP 2>/dev/null; exit 1; }
    sleep 20
  done
  wait $PB $PP
  echo "$BM 終了 RUN_RC=$(cat $BM/RUN_RC)、$PM 終了 RUN_RC=$(cat $PM/RUN_RC) $(date -Is)" >> $LOG
done
echo "== 終了 $(date -Is)" >> $LOG
