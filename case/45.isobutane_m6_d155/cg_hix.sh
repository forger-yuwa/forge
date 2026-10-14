#!/bin/bash
# plan tooling-nozzle-core-grid §4.100 P-hix (事前登録、codex 諮問 2026-10-14 hi-next): 入口半径を半分にした形 H の本段 (run_0915_hi_h_main) を、
# 同じ格子・バイナリ・実効設定 (B0、cfl_pseudo 4.0) のまま 40,000 step 継続する (本段の通算 (60k, 100k])。C は延長しない (run_0912 の (40k, 60k] を固定の対照にする)。
# 起点は run_0915_hi_h_main/res_60000.h5 を H の格子 (prep_Hi/nozzle.h5) へ restart_field --keep-src-dtype (9 量・ビット一致)。
# 設定は run_0915 の solverConfig と nStepOuter (60000 → 40000) だけが違うこと、bcondConfig は同一であることを確かめてから回す。
# 出力は 1000 ごと。8 量・縮流部と外と throatwall の残差・物理性を取り出して場を消す (res_40000 は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_hix.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; HPREP=_band_ab/core_grid/prep_Hi
LOG=$C45/cg_hix.log; OUT=1000
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
M=run_0915_hi_h_main; R=run_0916_hi_h_cont
[ -f $M/res_60000.h5 ] && [ "$(cat $M/RUN_RC)" = 0 ] || { echo "H の本段の res_60000 が無い (中止)" >> $LOG; exit 1; }
mk $R H 4.0 40000 $M/res_60000.h5 || { echo "延長の準備に失敗 $(date -Is)" >> $LOG; exit 1; }
# 実効設定が H の本段と同じか (違いは nStepOuter だけ・bcond は同一)
dd=$(diff $M/solverConfig.yaml $R/solverConfig.yaml | grep '^[<>]' | grep -vc 'nStepOuter: [46]0000}')
nd=$(diff $M/solverConfig.yaml $R/solverConfig.yaml | grep -c '^>')
if [ "$dd" != 0 ] || [ "$nd" != 1 ] || ! cmp -s $M/bcondConfig.yaml $R/bcondConfig.yaml; then
  echo "$R: H の本段と設定が違う (想定外の行 $dd・追加行 $nd) (中止)" >> $LOG; exit 1
fi
echo "$R: 設定は H の本段と nStepOuter だけが違う (bcond 同一)" >> $LOG
runone $R 40000 "40000" && echo "延長 終了 $(date -Is)" >> $LOG || echo "延長 失敗 $(date -Is)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
