#!/bin/bash
# plan tooling-nozzle-core-grid §4.82 P-m3 (事前登録、ユーザ決定・codex 諮問 2026-10-13 point-precond): point 陰解法で低マッハの時間項の前処理 (lowMachPrecond 3、左辺だけ・右辺不変) の A/B。
# 2 腕を同じ起点 (run_0608 の res_60000、restart_field 9 量・ビット一致) から各 20,000 step、同時に回す。バイナリは node の前処理の境界修正 (a60d8eb5) 入りの FP64。
#   共通: point (lineImplicit 0・lineDtDirectional 0・ITJ 0、lineDtDirectionalCap は消す)・cfl_pseudo 1・緩和 0.7・nStepInner 5。C = lowMachPrecond 0、P = lowMachPrecond 3 (この 1 行だけ)。
# 2500 step ごとの出力から 8 量の系列 (cg_gref_judge.py series)・縮流部と外の残差 (cg_region_res.py)・物理性 (cg_phys.py) を取り出し、
# res_0・2500・20000 以外の場は取り出した後に消す (2500 は資格の検査 cg_pp3_check.py、20000 は組成の確認に使う)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_pp3.sh <C の run> <P の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64-lmp2/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; R=$1; L=$2; PT=run_0611_cg_gref_point
LOG=$C45/cg_pp3.log; MIN_START_GB=5; MIN_RUN_GB=2; STEPS=20000; OUT=2500; KEEP=" 2500 20000 "
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $R $L F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F7_64 | cut -c1-8)" = d58ed787 ] || { echo "バイナリの sha256 が d58ed787 でない (中止)" >> $LOG; exit 1; }
for r in $R $L; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
done
for r in $R $L; do sed -i "s/deltaT: {lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0,/deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 1.0,/" $r/solverConfig.yaml; done
sed -i "s/lowMachPrecond: 0, /lowMachPrecond: 3, /" $L/solverConfig.yaml
# 設定の検査: C と P の差は lowMachPrecond の値の 1 行だけ (C 0、P 3、キーは各 1 個)。C と run_0611 (point の既存 run) の差は cfl_pseudo・extraFields・nStepOuter の行だけ
dRL=$(diff $R/solverConfig.yaml $L/solverConfig.yaml | grep -c '^>'); dPT=$(diff $PT/solverConfig.yaml $R/solverConfig.yaml | grep '^>' | grep -vc 'cfl_pseudo: 1.0\|extraFields\|nStepOuter')
[ "$dRL" = 1 ] && [ "$dPT" = 0 ] && [ "$(grep -o lowMachPrecond $R/solverConfig.yaml | wc -l)" = 1 ] && [ "$(grep -o lowMachPrecond $L/solverConfig.yaml | wc -l)" = 1 ] && [ "$(grep -c 'lowMachPrecond: 0, ' $R/solverConfig.yaml)" = 1 ] && [ "$(grep -c 'lowMachPrecond: 3, ' $L/solverConfig.yaml)" = 1 ] \
  && grep -q "deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 1.0," $R/solverConfig.yaml && [ "$(grep -c lineDtDirectionalCap $R/solverConfig.yaml)" = 0 ] \
  && grep -q "nStepOuter: $STEPS}" $R/solverConfig.yaml && grep -q "volume, dt_local\]" $R/solverConfig.yaml && diff $R/bcondConfig.yaml $L/bcondConfig.yaml > /dev/null \
  || { echo "設定の差が想定と違う (C-P $dRL 行、run_0611 との想定外の差 $dPT 行) (中止)" >> $LOG; exit 1; }
echo "C と P の差は lowMachPrecond の値の 1 行 (C 0、P 3)。C と run_0611 の差は cfl_pseudo (4 → 1)・extraFields・nStepOuter だけ" >> $LOG
PIDS=()
for r in $R $L; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG; sleep 20
done
sleep 120
for r in $R $L; do echo "$r 起動の確認: $(grep -h "'lowMachPrecond'\|'precondEps'\|'speciesImplicitCoupling'\|'speciesPrecondDt'\|'lineImplicit'\|'cfl_pseudo'" $r/forge_run.log | head -6 | tr '\n' ' ' | cut -c1-400)" >> $LOG; done
while true; do
  done_all=1
  for r in $R $L; do
    [ -f $r/RUN_RC ] || done_all=0
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f $r/.done_$n ] && continue
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        # 取り出しに失敗しても (発散した場など) 印を付けて先へ進む (待ち合わせが止まらないように)。失敗はログに残し、場は消さない
        if python3 cg_gref_judge.py series "${r}_$n" $r --from $n --to $n > /dev/null 2>&1 && python3 cg_region_res.py $r $n && python3 cg_phys.py $r $n; then
          echo "$r: $n を取り出した" >> $LOG
          case "$KEEP" in *" $n "*) echo "$r: $n の場を残す" >> $LOG ;; *) rm -f $f $r/res_$n.xmf ;; esac
        else
          echo "$r: $n の取り出しに失敗 (場は残す)" >> $LOG
        fi
        touch $r/.done_$n
      fi
    done
  done
  [ $done_all = 1 ] && { left=0; for r in $R $L; do for f in $r/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f $r/.done_$n ] && left=1; done; done; [ $left = 0 ] && break; }
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $R $L; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
  fi
  sleep 30
done
wait "${PIDS[@]}"
for r in $R $L; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
