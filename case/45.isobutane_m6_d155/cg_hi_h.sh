#!/bin/bash
# plan tooling-nozzle-core-grid §4.98 P-hi の H の鎖の続き (段 2 → 本段)。cg_hi.sh の H の鎖は段 2 の restart_field で止まった
# (H の格子ファイルの roY0・roY1 の型が元の場と違い、相対 5e−8 の丸めが入って「ビット一致」の判定文にならなかった)。
# ここでは restart_field に --keep-src-dtype を付けて SRC の型のまま移す (ビット一致)。他は cg_hi.sh と同じ。C の鎖は cg_hi.sh のまま走らせる。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_hi_h.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; HPREP=_band_ab/core_grid/prep_Hi
LOG=$C45/cg_hi_h.log; OUT=1000
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F $(sha256sum $F | cut -c1-16) Gref の起点 $(sha256sum $GREF/res_60000.h5 | cut -c1-16) H の入力 $(sha256sum $HPREP/nozzle.h5 | cut -c1-16)" >> $LOG
[ "$(sha256sum $F | cut -c1-8)" = 7e95e5bb ] || { echo "バイナリの sha256 が 7e95e5bb でない (中止)" >> $LOG; exit 1; }
grep -q '"pass": true' $HPREP/hi_mesh_check.json || { echo "H の格子の資格が PASS でない (中止)" >> $LOG; exit 1; }
# mk <run> <arm C|H> <cfl_pseudo> <nStepOuter> <seed: res の h5 または "-"(H の準備の初期値をそのまま)>
mk() {
  local r=$1 arm=$2 cflp=$3 n=$4 seed=$5 base meshsrc
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r
  if [ $arm = C ]; then base=$GREF; meshsrc=$GREF/nozzle.h5; else base=$HPREP; meshsrc=$HPREP/nozzle.h5; fi
  for f in $base/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|*.xmf|*.png|nozzle.h5|nozzle.msh|solverConfig.yaml|bcondConfig.yaml) ;; *) cp -a $f $r/;; esac; done
  cp $GREF/solverConfig.yaml $GREF/bcondConfig.yaml $r/
  cp $meshsrc $r/nozzle.h5
  if [ "$seed" != "-" ]; then
    python3 $TL/restart_field.py $seed $r/nozzle.h5 --dst-run $r --forge $F --keep-src-dtype > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
    grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  else
    echo "H の準備の初期値をそのまま使う ($meshsrc sha256 $(sha256sum $meshsrc | cut -c1-16))" > $r/restart.log
  fi
  sed -i "s/cfl: 4.0, cfl_pseudo: 4.0,/cfl: 4.0, cfl_pseudo: $cflp,/; s/nStepOuter: 60000}/nStepOuter: $n}/; s/outStepInterval: 2500/outStepInterval: $OUT/; s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
  local d=$(diff $GREF/solverConfig.yaml $r/solverConfig.yaml | grep '^>' | grep -vc "cfl_pseudo: $cflp,\|nStepOuter: $n}\|outStepInterval: $OUT$\|extraFields: \[res_ro, res_roUx")
  [ "$d" = 0 ] && [ "$(grep -c "cfl: 4.0, cfl_pseudo: $cflp," $r/solverConfig.yaml)" = 1 ] && [ "$(grep -c "nStepOuter: $n}" $r/solverConfig.yaml)" = 1 ] \
    && [ "$(grep -c 'lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1,' $r/solverConfig.yaml)" = 1 ] \
    || { echo "$r: 設定が想定と違う (想定外の行 $d) (中止)" >> $LOG; return 1; }
  echo "$r: 準備 ($arm cfl_pseudo $cflp nStepOuter $n 起点 $seed)" >> $LOG
}
runone() {   # runone <run> <nStepOuter> <残す step (空白区切り)>
  local r=$1 n=$2 keep=" $3 "
  [ "$(free_gb)" -ge 6 ] || { echo "$r: ディスクの空きが 6 GB 未満 (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$F; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  local pid=$!; echo "$r 起動 $(date -Is)" >> $LOG
  while true; do
    for f in $r/res_*.h5; do
      local k=$(basename $f .h5); k=${k#res_}; [[ $k =~ ^[0-9]+$ ]] || continue; [ $k -gt 0 ] || continue
      [ -f $r/.done_$k ] && continue
      if [ -f $r/res_$((k + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        if python3 cg_gref_judge.py series "${r}_$k" $r --from $k --to $k > /dev/null 2>&1 && python3 cg_region_res.py $r $k && python3 cg_phys.py $r $k; then
          case "$keep" in *" $k "*) ;; *) rm -f $f $r/res_$k.xmf ;; esac
        else
          echo "$r: $k の取り出しに失敗 (場は残す)" >> $LOG
        fi
        touch $r/.done_$k
      fi
    done
    if [ -f $r/RUN_RC ]; then local left=0; for f in $r/res_*.h5; do local k=$(basename $f .h5); k=${k#res_}; [[ $k =~ ^[0-9]+$ ]] && [ $k -gt 0 ] && [ ! -f $r/.done_$k ] && left=1; done; [ $left = 0 ] && break; fi
    if [ "$(free_gb)" -lt 3 ]; then kill $pid 2>/dev/null; for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && kill $p; done; echo "$r: ディスク 3 GB 未満で停止 $(date -Is)" >> $LOG; fi
    sleep 20
  done
  wait $pid
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の出力 $(ls $r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//; s/.h5//' | sort -n | tail -1)" >> $LOG
  [ "$(cat $r/RUN_RC)" = 0 ] && [ -f $r/res_$n.h5 ] || { echo "$r: 完走していない (中止)" >> $LOG; return 1; }
}
chain() {   # chain <arm> <s1> <s2> <main> <最初の起点>
  local arm=$1 s1=$2 s2=$3 m=$4 seed0=$5
  mk $s1 $arm 0.5 2000 $seed0 && runone $s1 2000 "2000" \
  && mk $s2 $arm 1.0 2000 $s1/res_2000.h5 && runone $s2 2000 "2000" \
  && mk $m $arm 4.0 60000 $s2/res_2000.h5 && runone $m 60000 "60000" \
  && echo "鎖 $arm 終了 $(date -Is)" >> $LOG || echo "鎖 $arm 失敗 $(date -Is)" >> $LOG
}
[ -f run_0913_hi_h_s1/res_2000.h5 ] && [ "$(cat run_0913_hi_h_s1/RUN_RC)" = 0 ] || { echo "H の段 1 が無い (中止)" >> $LOG; exit 1; }
mk run_0914_hi_h_s2 H 1.0 2000 run_0913_hi_h_s1/res_2000.h5 && runone run_0914_hi_h_s2 2000 "2000" \
  && mk run_0915_hi_h_main H 4.0 60000 run_0914_hi_h_s2/res_2000.h5 && runone run_0915_hi_h_main 60000 "60000" \
  && echo "鎖 H 終了 $(date -Is)" >> $LOG || echo "鎖 H 失敗 $(date -Is)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
