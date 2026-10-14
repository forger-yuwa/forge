#!/bin/bash
# plan tooling-nozzle-core-grid §4.96 P-ts (事前登録、ユーザ決定・codex 諮問 2026-10-13 thermofloat-stop): 熱物性の評価の経路 (physProp.thermoFloat) が
# B0 の止まる先を変えるか。起点 run_0608 の res_60000 から入力を 1 回作り 4 本に複製、旧バイナリ (7e95e5bb) で B0 のまま各 60,000 step を 1 本ずつ:
#   T1a → T0a → T0b → T1b。違いは physProp.thermoFloat の 1 (明示) と 0 だけ。出力 1000 ごと、8 量・縮流部の残差・物理性を取り出して場を消す (res_0・60000 は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_tfs.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; STAGE=$C45/_tfs_input_nozzle.h5
LOG=$C45/cg_tfs.log; STEPS=60000; OUT=1000
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F $(sha256sum $F | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(sha256sum $F | cut -c1-8)" = 7e95e5bb ] || { echo "バイナリの sha256 が 7e95e5bb でない (中止)" >> $LOG; exit 1; }
prep() {   # prep <run> <thermoFloat 0|1>
  local r=$1 tf=$2
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|*.xmf|*.png) ;; *) cp -a $f $r/;; esac; done
  if [ ! -f $STAGE ]; then
    python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
    grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
    cp $r/nozzle.h5 $STAGE; echo "入力を作った: $STAGE sha $(sha256sum $STAGE | cut -c1-16)" >> $LOG
  else
    cp $STAGE $r/nozzle.h5; echo "入力 h5 を複製: $STAGE (sha256 $(sha256sum $STAGE | cut -c1-16))" > $r/restart.log
  fi
  sed -i "s/outStepInterval: 2500/outStepInterval: $OUT/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/; s/physProp: {thermalMethod: 2, /physProp: {thermoFloat: $tf, thermalMethod: 2, /" $r/solverConfig.yaml
  # 設定の検査: run_0608 との差は outStepInterval・output・physProp の 3 行だけ (pipefail の下なので数える形)
  local d=$(diff $GREF/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | grep -vc "outStepInterval: $OUT$\|extraFields: \[res_ro, res_roUx\|physProp: {thermoFloat: $tf, thermalMethod: 2, ")
  [ "$d" = 0 ] && [ "$(diff $GREF/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')" = 3 ] && [ "$(grep -o thermoFloat $r/solverConfig.yaml | wc -l)" = 1 ] \
    && [ "$(grep -c "nStepOuter: $STEPS}" $r/solverConfig.yaml)" = 1 ] && [ "$(grep -c 'lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1,' $r/solverConfig.yaml)" = 1 ] \
    && [ "$(sha256sum $r/nozzle.h5 | cut -c1-64)" = "$(sha256sum $STAGE | cut -c1-64)" ] || { echo "$r: 設定・入力が想定と違う (想定外の行 $d) (中止)" >> $LOG; return 1; }
  echo "$r: 準備 thermoFloat=$tf" >> $LOG
}
runone() {
  local r=$1
  [ "$(free_gb)" -ge 8 ] || { echo "$r: ディスクの空きが 8 GB 未満 (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$F; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  local pid=$!; echo "$r 起動 $(date -Is)" >> $LOG
  sleep 90; echo "$r 起動の確認: $(grep -h "'thermoFloat'\|'lineDtDirectionalCap'" $r/forge_run.log | head -2 | tr '\n' ' ' | cut -c1-200)" >> $LOG
  while true; do
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f $r/.done_$n ] && continue
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        if python3 cg_gref_judge.py series "${r}_$n" $r --from $n --to $n > /dev/null 2>&1 && python3 cg_region_res.py $r $n && python3 cg_phys.py $r $n; then
          [ $n -ne $STEPS ] && rm -f $f $r/res_$n.xmf
        else
          echo "$r: $n の取り出しに失敗 (場は残す)" >> $LOG
        fi
        touch $r/.done_$n
      fi
    done
    if [ -f $r/RUN_RC ]; then left=0; for f in $r/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f $r/.done_$n ] && left=1; done; [ $left = 0 ] && break; fi
    if [ "$(free_gb)" -lt 3 ]; then kill $pid 2>/dev/null; for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && kill $p; done; echo "$r: ディスク 3 GB 未満で停止 $(date -Is)" >> $LOG; fi
    sleep 20
  done
  wait $pid
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の出力 $(ls $r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//; s/.h5//' | sort -n | tail -1)" >> $LOG
}
for spec in "run_0905_tfs_t1a 1" "run_0906_tfs_t0a 0" "run_0907_tfs_t0b 0" "run_0908_tfs_t1b 1"; do
  set -- $spec
  prep $1 $2 || exit 1
  runone $1 || exit 1
done
echo "== 終了 $(date -Is)" >> $LOG
