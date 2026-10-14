#!/bin/bash
# plan tooling-nozzle-core-grid §4.104 P-hs (事前登録、ユーザ指示・codex 諮問 2026-10-14 hi-settle): 入口半径を半分にした形 H を、
# 共通の起点 run_0916_hi_h_cont/res_40000.h5 (本段の通算 10 万 step) から 2 本同時に回す。
#   腕 b0    run_0917_hs_b0    : 今の B0 のまま (cfl_pseudo 4・緩和 0.7・ライン陰解法・上限 50・ITJ 5)。上限 14 万 step (通算 24 万)。
#   腕 point run_0918_hs_point : point (lineImplicit 0・lineDtDirectional 0・ITJ 0・上限のキーを外す、§4.26 と同じ書き換え)。上限 20 万 step。
# 両腕とも出力 2,500 ごと。出力ごとに 8 量 (cg_gref_judge.py series)・領域の残差・物理性・停止規則の量 (cg_hs_stoprow.py) を取り出し、
# cg_hs_check.py が STOP を返したらその腕の forge を止める (止めた出力は残す)。それ以外の途中の場は消す (最後の出力は残す)。
# 停止規則の帯は共通の起点で 1 回だけ固定する (cg_hs_band.py → _band_ab/core_grid/hs_band.json)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_hs.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRC=run_0916_hi_h_cont; SEED=$SRC/res_40000.h5; HPREP=_band_ab/core_grid/prep_Hi
LOG=$C45/cg_hs.log; OUT=2500
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F $(sha256sum $F | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(sha256sum $F | cut -c1-8)" = 7e95e5bb ] || { echo "バイナリの sha256 が 7e95e5bb でない (中止)" >> $LOG; exit 1; }
[ -f $SEED ] && [ "$(cat $SRC/RUN_RC)" = 0 ] || { echo "起点が無い (中止)" >> $LOG; exit 1; }
if [ ! -f _band_ab/core_grid/hs_band.json ]; then
  python3 cg_hs_band.py >> $LOG 2>&1 || { echo "帯の固定に失敗 (中止)" >> $LOG; exit 1; }
fi
# mk <run> <arm b0|point> <nStepOuter>
mk() {
  local r=$1 arm=$2 n=$3
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r
  for f in $HPREP/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|*.xmf|*.png|nozzle.h5|nozzle.msh|solverConfig.yaml|bcondConfig.yaml) ;; *) cp -a $f $r/;; esac; done
  cp $SRC/solverConfig.yaml $SRC/bcondConfig.yaml $r/
  cp $HPREP/nozzle.h5 $r/nozzle.h5
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F --keep-src-dtype > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  sed -i "s/nStepOuter: 40000}/nStepOuter: $n}/; s/outStepInterval: 1000/outStepInterval: $OUT/" $r/solverConfig.yaml
  if [ $arm = point ]; then
    sed -i 's/deltaT: {lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1, /deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, /' $r/solverConfig.yaml
    [ "$(grep -c 'deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, ' $r/solverConfig.yaml)" = 1 ] || { echo "$r: point への書き換えが想定と違う (中止)" >> $LOG; return 1; }
  fi
  local nd=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
  local exp=2; [ $arm = point ] && exp=3
  [ "$nd" = "$exp" ] && cmp -s $SRC/bcondConfig.yaml $r/bcondConfig.yaml \
    && [ "$(grep -c "nStepOuter: $n}" $r/solverConfig.yaml)" = 1 ] && [ "$(grep -c "outStepInterval: $OUT" $r/solverConfig.yaml)" = 1 ] \
    || { echo "$r: 設定の差が想定と違う (追加行 $nd、想定 $exp) (中止)" >> $LOG; return 1; }
  echo "$r: 準備 ($arm nStepOuter $n 起点 $SEED)。$SRC との差: $(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | cut -c1-140 | tr '\n' ' ')" >> $LOG
}
killrun() {   # killrun <run>: cwd がその run の forge だけを止める
  local r=$1 p w
  for sig in TERM KILL; do
    for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && kill -$sig $p; done
    for w in $(seq 1 36); do local alive=0; for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && alive=1; done; [ $alive = 0 ] && return 0; sleep 5; done
  done
}
runone() {   # runone <run> <arm> <nStepOuter>
  local r=$1 arm=$2 n=$3
  [ "$(free_gb)" -ge 6 ] || { echo "$r: ディスクの空きが 6 GB 未満 (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$F; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  local pid=$!; echo "$r 起動 $(date -Is)" >> $LOG
  while true; do
    for f in $r/res_*.h5; do
      local k=$(basename $f .h5); k=${k#res_}; [[ $k =~ ^[0-9]+$ ]] || continue; [ $k -gt 0 ] || continue
      [ -f $r/.done_$k ] && continue
      if [ -f $r/res_$((k + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        if python3 cg_gref_judge.py series "${r}_$k" $r --from $k --to $k > /dev/null 2>&1 && python3 cg_region_res.py $r $k && python3 cg_phys.py $r $k && python3 cg_hs_stoprow.py $r $k; then
          local dec=$(python3 cg_hs_check.py $r $arm 2>>$r/hs_check.err | tail -1)
          if [ "$dec" = STOP ] && [ ! -f $r/STOP_RULE ]; then
            echo "$k" > $r/STOP_RULE; echo "$r: $k で止める条件を満たした (cg_hs_check) $(date -Is)" >> $LOG; killrun $r
          fi
          local keep=0; [ -f $r/STOP_RULE ] && [ "$(cat $r/STOP_RULE)" = $k ] && keep=1; [ $k = $n ] && keep=1
          [ $keep = 1 ] || rm -f $f $r/res_$k.xmf
        else
          echo "$r: $k の取り出しに失敗 (場は残す)" >> $LOG
        fi
        touch $r/.done_$k
      fi
    done
    if [ -f $r/RUN_RC ]; then local left=0; for f in $r/res_*.h5; do local k=$(basename $f .h5); k=${k#res_}; [[ $k =~ ^[0-9]+$ ]] && [ $k -gt 0 ] && [ ! -f $r/.done_$k ] && left=1; done; [ $left = 0 ] && break; fi
    if [ "$(free_gb)" -lt 3 ]; then killrun $r; echo "$r: ディスク 3 GB 未満で停止 $(date -Is)" >> $LOG; fi
    sleep 20
  done
  wait $pid
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) STOP_RULE=$(cat $r/STOP_RULE 2>/dev/null || echo -) 最後の出力 $(ls $r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//; s/.h5//' | sort -n | tail -1)" >> $LOG
}
mk run_0917_hs_b0 b0 140000 && mk run_0918_hs_point point 200000 || { echo "準備に失敗 (中止)" >> $LOG; exit 1; }
runone run_0917_hs_b0 b0 140000 &
P1=$!; sleep 60
runone run_0918_hs_point point 200000 &
P2=$!
wait $P1 $P2
echo "== 終了 $(date -Is)" >> $LOG
