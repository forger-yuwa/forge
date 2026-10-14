#!/bin/bash
# plan tooling-nozzle-core-grid §4.70 (事前登録、codex 諮問 2026-10-12 next-direction、ユーザ決定「susumete」): Gref の上限 50 の双子の回し直し。
# §4.62 と同じ起点 (run_0619 の res_60000)・バイナリ・設定・摂動 (列の番号で直接、a = 1e−7、sha256 05a38657… を確認) で、各 3,000 step・250 ごと。
# 列ごとの振幅と帯の擬似 dt を残す cg_twin_diff.py (§4.65・§4.70 版) で差を取る。1 step の確認はしない (雑音の床は §4.64 の値)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_twin_r.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRC=run_0619_cg_gref_b0_ctl; SEED=$SRC/res_60000.h5; LOG=$C45/cg_twin_r.log; AMP=1e-7; STEPS=3000; OUT=250
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge 8 ] || { echo "ディスクの空きが 8 GB 未満 (中止)" >> $LOG; exit 1; }
mk() {  # mk <run> <step 数> <出力の間隔> <摂動 0/1>
  local r=$1 n=$2 out=$3 pert=$4
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r; for f in $SRC/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|.done_*|*.png) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  sed -i "s/nStepOuter: [0-9]*}/nStepOuter: $n}/; s/outStepInterval: [0-9]*/outStepInterval: $out/; s/output: {level: 1, extraFields: \[[^]]*\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
  local d=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
  [ "$d" = 3 ] && grep -q "lineDtDirectionalCap: 50.0," $r/solverConfig.yaml && grep -q "nStepOuter: $n}" $r/solverConfig.yaml && grep -q "volume, dt_local\]" $r/solverConfig.yaml \
    || { echo "$r: 設定の差が想定と違う ($d 行) (中止)" >> $LOG; return 1; }
  if [ $pert = 1 ]; then python3 cg_twin_perturb.py $r $AMP $SEED >> $LOG 2>&1 && grep -q '"sha256_delta": "05a38657e40c38f72e178f0239fe610877049ca5a8edf1d55e14c4487e176e23"' $r/perturb.json || { echo "$r: 摂動に失敗 (または §4.62 の配列と違う)" >> $LOG; return 1; }; fi
  echo "$r: Gref・上限 50・$n step・出力 $out ごと・摂動 $pert" >> $LOG
}
run() { ( export FORGE_BIN=$F7_64; bash $RC $C45/$1 > $1/run_case_stdout.log 2>&1; echo $? > $1/RUN_RC ); }
BM=run_0679_tw50r_base; PM=run_0680_tw50r_pert
mk $BM $STEPS $OUT 0 && mk $PM $STEPS $OUT 1 || exit 1
( run $BM ) & PB=$!; sleep 15; ( run $PM ) & PP=$!; echo "$BM・$PM 起動 $(date -Is)" >> $LOG
while true; do
  for f in $BM/res_*.h5; do
    n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
    [ -f $BM/.tw_$n ] && continue; [ -f $PM/res_$n.h5 ] || continue
    if { [ -f $BM/res_$((n + OUT)).h5 ] && [ -f $PM/res_$((n + OUT)).h5 ]; } || { [ -f $BM/RUN_RC ] && [ -f $PM/RUN_RC ]; }; then
      python3 cg_twin_diff.py $BM $PM $n main_gref_r && touch $BM/.tw_$n
      if [ $((n % 1000)) -ne 0 ] && [ $n -ne $STEPS ]; then rm -f $BM/res_$n.h5 $BM/res_$n.xmf $PM/res_$n.h5 $PM/res_$n.xmf; fi
    fi
  done
  if [ -f $BM/RUN_RC ] && [ -f $PM/RUN_RC ]; then left=0; for f in $BM/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f $BM/.tw_$n ] && left=1; done; [ $left = 0 ] && break; fi
  [ "$(free_gb)" -lt 3 ] && { echo "ディスク 3 GB 未満 (中止)" >> $LOG; kill $PB $PP 2>/dev/null; exit 1; }
  sleep 20
done
wait $PB $PP
[ -f $BM/res_0.h5 ] && [ -f $PM/res_0.h5 ] && python3 cg_twin_diff.py $BM $PM 0 main_gref_r
echo "$BM 終了 RUN_RC=$(cat $BM/RUN_RC)、$PM 終了 RUN_RC=$(cat $PM/RUN_RC) $(date -Is)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
