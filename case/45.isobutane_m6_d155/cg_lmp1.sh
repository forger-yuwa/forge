#!/bin/bash
# plan tooling-nozzle-core-grid §4.80 L-m1 (事前登録、codex 諮問 2026-10-13 lowmach): 低マッハの散逸 (lowMachPrecond 1) だけを変えた A/B。
# 2 腕を同じ起点 (run_0608 の res_60000、restart_field 9 量・ビット一致) から各 20,000 step、同時に回す。
#   R = B0 そのまま (lowMachPrecond 0)、L = time.deltaT.lowMachPrecond: 1 の 1 行だけ (precondEps は既定 0.15)。他は B0 (ITJ 5・方向別 dt・上限 50・cfl 4・緩和 0.7・nStepInner 5)。
# 2500 step ごとの出力から 8 量の系列 (cg_gref_judge.py series)・縮流部と外の残差 (cg_region_res.py)・物理性 (cg_phys.py) を取り出し、
# 10k・15k・20k では P の格子間隔の上下 (cg_lmp1_saw.py) も取る。10k・15k・20k 以外の場は取り出した後に消す (res_0 は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_lmp1.sh <R の run> <L の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; R=$1; L=$2
LOG=$C45/cg_lmp1.log; MIN_START_GB=5; MIN_RUN_GB=2; STEPS=20000; OUT=2500; KEEP=" 10000 15000 20000 "
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $R $L F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
for r in $R $L; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
done
sed -i "s/deltaT: {lineDtDirectionalCap: 50.0, /deltaT: {lineDtDirectionalCap: 50.0, lowMachPrecond: 1, /" $L/solverConfig.yaml
# 設定の検査: R と L の差は lowMachPrecond の 1 行だけ、R にはこのキーが無い (既定 0)
dRL=$(diff $R/solverConfig.yaml $L/solverConfig.yaml | grep -c '^>')
[ "$dRL" = 1 ] && [ "$(grep -c lowMachPrecond $R/solverConfig.yaml)" = 0 ] && [ "$(grep -c 'lowMachPrecond: 1, ' $L/solverConfig.yaml)" = 1 ] \
  && grep -q "nStepOuter: $STEPS}" $R/solverConfig.yaml && grep -q "volume\]" $R/solverConfig.yaml \
  && grep -q "implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $R/solverConfig.yaml && diff $R/bcondConfig.yaml $L/bcondConfig.yaml > /dev/null \
  || { echo "設定の差が想定と違う (R-L $dRL 行) (中止)" >> $LOG; exit 1; }
echo "R と L の差は lowMachPrecond の 1 行 (R は既定 0)" >> $LOG
PIDS=()
for r in $R $L; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG; sleep 20
done
sleep 120
for r in $R $L; do echo "$r 起動の確認: $(grep -h "'lowMachPrecond'\|'precondEps'" $r/forge_run.log | head -2 | tr '\n' ' ' | cut -c1-200)" >> $LOG; done
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
          case "$KEEP" in *" $n "*) python3 cg_lmp1_saw.py $r $n && echo "$r: $n の格子間隔の上下を取った" >> $LOG ;; *) rm -f $f $r/res_$n.xmf ;; esac
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
