#!/bin/bash
# plan tooling-nozzle-core-grid §4.88 P-k (事前登録、codex 諮問 2026-10-13 pdt-interp): Γ = I (precondEps 1.0) で、前処理のカーネル (lowMachPrecond 3) と標準のカーネル (0) を比べる。
# 2 腕を同じ起点 (run_0608 の res_60000、restart_field 9 量・ビット一致) から各 20,000 step、1000 ごとに出力、同時に回す。
#   共通: point (lineImplicit 0・lineDtDirectional 0・ITJ 0)・cfl_pseudo 1・緩和 0.7・nStepInner 5・precondEps 1.0・両腕とも FORGE_DIAG_PRECOND_NO_DT_SCALE=1。
#   K0 = lowMachPrecond 0、K3 = lowMachPrecond 3 (設定の違いはこの 1 行だけ)。
# 出力ごとに 8 量の系列・領域の残差 (縮流部・外・throatwall)・物理性を取り出し、res_0・1000・20000 以外の場は消す。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_pk.sh <K0 の run> <K3 の run> <バイナリの sha256 の先頭 8 文字>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64-lmp2dg/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; R=$1; L=$2; SHA8=$3; 
LOG=$C45/cg_pk.log; MIN_START_GB=8; MIN_RUN_GB=3; STEPS=20000; OUT=1000; KEEP=" 1000 20000 "
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $R $L F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F7_64 | cut -c1-8)" = "$SHA8" ] || { echo "バイナリの sha256 が $SHA8 でない (中止)" >> $LOG; exit 1; }
[ "$(strings $F7_64 | grep -c FORGE_DIAG_PRECOND_NO_DT_SCALE)" -ge 1 ] || { echo "バイナリに FORGE_DIAG_PRECOND_NO_DT_SCALE が無い (中止)" >> $LOG; exit 1; }
for r in $R $L; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  mkdir $r; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; exit 1; }
  sed -i "s/nStepOuter: 60000}/nStepOuter: $STEPS}/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
done
for r in $R $L; do sed -i "s/deltaT: {lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0,/deltaT: {implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 1.0,/" $r/solverConfig.yaml; done
for r in $R $L; do sed -i "s/lowMachPrecond: 0, /precondEps: 1.0, lowMachPrecond: 0, /; s/outStepInterval: 2500/outStepInterval: 1000/" $r/solverConfig.yaml; done
sed -i "s/lowMachPrecond: 0, /lowMachPrecond: 3, /" $L/solverConfig.yaml
# 設定の検査: K0 と K3 の差は lowMachPrecond の値 (0 → 3) の 1 行だけ。precondEps と lowMachPrecond のキーは各 1 個 (pipefail の下なので数える形)
dRL=$(diff $R/solverConfig.yaml $L/solverConfig.yaml | grep -c '^>')
[ "$dRL" = 1 ] && [ "$(diff $R/solverConfig.yaml $L/solverConfig.yaml | grep '^>' | grep -c 'lowMachPrecond: 3, ')" = 1 ] \
  && [ "$(grep -o precondEps $R/solverConfig.yaml | wc -l)" = 1 ] && [ "$(grep -c 'precondEps: 1.0, lowMachPrecond: 0, ' $R/solverConfig.yaml)" = 1 ] && [ "$(grep -o lowMachPrecond $L/solverConfig.yaml | wc -l)" = 1 ] \
  && [ "$(grep -c 'outStepInterval: 1000' $R/solverConfig.yaml)" = 1 ] && [ "$(grep -c 'implicitThermalJacobian: 0, lineImplicit: 0, lineDtDirectional: 0,.*cfl_pseudo: 1.0,' $R/solverConfig.yaml)" = 1 ] \
  && diff $R/bcondConfig.yaml $L/bcondConfig.yaml > /dev/null || { echo "設定の差が想定と違う (K0-K3 $dRL 行) (中止)" >> $LOG; exit 1; }
echo "K0 と K3 の差は lowMachPrecond の値 (0 → 3) の 1 行。両腕とも precondEps 1.0・FORGE_DIAG_PRECOND_NO_DT_SCALE=1" >> $LOG
PIDS=()
for r in $R $L; do
  ( export FORGE_BIN=$F7_64 FORGE_DIAG_PRECOND_NO_DT_SCALE=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG; sleep 20
done
sleep 120
for r in $R $L; do echo "$r 起動の確認: $(grep -h "'lowMachPrecond'\|'precondEps'\|'lineImplicit'\|'cfl_pseudo'" $r/forge_run.log | head -4 | tr '\n' ' ' | cut -c1-300) 診断の切替の行 $(cat $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | grep -c 'FORGE_DIAG_PRECOND_NO_DT_SCALE=1:')" >> $LOG; done
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
