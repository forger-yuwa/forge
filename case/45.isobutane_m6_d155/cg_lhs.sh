#!/bin/bash
# plan tooling-nozzle-core-grid §4.59 (事前登録、codex 諮問 2026-10-12 consistent-lhs、ユーザ決定「整合させる案でいいや」):
# 左辺の壁法線方向の整合 — ライン面の熱伝導の近傍結合だけの効果。3 腕を同じ起点 (run_0608 の res_60000、restart_field 9 量・ビット一致) から各 20,000 step、同時に回す。
#   R = B0 そのまま (lineViscCoupling 0)、A = lineViscCoupling 3・FORGE_LVC_TERMS=0 (スカラー + 熱伝導の D + 値 3 の壁拘束、熱伝導の近傍 K なし)、
#   B = lineViscCoupling 3・FORGE_LVC_TERMS=2 (A に熱伝導の近傍 K だけを足す)。他は B0 (ITJ 5・方向別 dt・上限 50・cfl 4・緩和 0.7・nStepInner 5・SST 更新あり)。
# A・B は最初の外反復のライン行列を FORGE_LINE_DUMP で書き出す (資格: D・dt・初期の右辺が同じで、K の差がエネルギー行だけ)。
# 2500 step ごとの出力から 8 量の系列 (cg_gref_judge.py series)・縮流部の残差 (cg_region_res.py)・物理性 (cg_phys.py) を取り出して場を消す (res_0 と最後は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_lhs.sh <R の run> <A の run> <B の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; R=$1; A=$2; B=$3
LOG=$C45/cg_lhs.log; MIN_START_GB=8; MIN_RUN_GB=3; STEPS=20000; OUT=2500
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $R $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
strings $F7_64 | grep -q FORGE_LVC_TERMS || { echo "バイナリに FORGE_LVC_TERMS が無い (中止)" >> $LOG; exit 1; }
for r in $R $A $B; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
done
for r in $A $B; do sed -i "s/deltaT: {lineDtDirectionalCap: 50.0, /deltaT: {lineDtDirectionalCap: 50.0, lineViscCoupling: 3, /" $r/solverConfig.yaml; done
# 設定の検査: R と A の差は lineViscCoupling の 1 行だけ、A と B は同一 (違いは環境変数のマスクだけ)
dRA=$(diff $R/solverConfig.yaml $A/solverConfig.yaml | grep -c '^>'); dAB=$(diff $A/solverConfig.yaml $B/solverConfig.yaml | grep -c '^>')
[ "$dRA" = 1 ] && [ "$dAB" = 0 ] && grep -q "lineViscCoupling: 3, " $A/solverConfig.yaml && grep -q "nStepOuter: $STEPS}" $R/solverConfig.yaml && grep -q "volume\]" $R/solverConfig.yaml \
  && grep -q "implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $R/solverConfig.yaml && diff $R/bcondConfig.yaml $B/bcondConfig.yaml > /dev/null \
  || { echo "設定の差が想定と違う (R-A $dRA 行、A-B $dAB 行) (中止)" >> $LOG; exit 1; }
NODES=$(python3 cg_sweep_nodes.py $A) || { echo "書き出す節点を選べない (中止)" >> $LOG; exit 1; }
echo "R と A の差は lineViscCoupling の 1 行、A と B は同一 (マスク 0 / 2 は環境変数)" >> $LOG
PIDS=()
for r in $R $A $B; do
  case $r in $A) TERMS=0 ;; $B) TERMS=2 ;; *) TERMS="" ;; esac
  if [ -n "$TERMS" ]; then mkdir $r/linedump
    ( export FORGE_BIN=$F7_64 FORGE_LVC_TERMS=$TERMS FORGE_LINE_DUMP_DIR=$C45/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
      bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else
    ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  PIDS+=($!); echo "$r 起動 $(date -Is) マスク ${TERMS:-なし}" >> $LOG; sleep 20
done
sleep 120
for r in $R $A $B; do echo "$r 起動の確認: $(grep -h "'lineViscCoupling'\|FORGE_LVC_TERMS=" $r/forge_run.log | head -2 | tr '\n' ' ' | cut -c1-200)" >> $LOG; done
while true; do
  done_all=1
  for r in $R $A $B; do
    [ -f $r/RUN_RC ] || done_all=0
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && continue
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        python3 cg_gref_judge.py series "${r}_$n" $r --from $n --to $n > /dev/null 2>&1 && python3 cg_region_res.py $r $n && python3 cg_phys.py $r $n && echo "$r: $n を取り出した" >> $LOG
        [ $n -ne $STEPS ] && [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && rm -f $f $r/res_$n.xmf
      fi
    done
  done
  [ $done_all = 1 ] && { left=0; for r in $R $A $B; do for f in $r/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && left=1; done; done; [ $left = 0 ] && break; }
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $R $A $B; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
  fi
  sleep 30
done
wait "${PIDS[@]}"
for r in $R $A $B; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
