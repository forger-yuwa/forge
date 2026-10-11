#!/bin/bash
# plan time_integration-line-implicit-speed §6.11 (2026-10-10、v2): 仕上げまでの総時間。(1) 1 step・出力・起動の専有時間 (P・L0・L5 × 3)、
# (2) L5 のライン段 (水準で止める) → (3) point の段 (水準 + 参照から 0.1 % 以内で止める)、(4) 判定。バイナリ lineL。
set -uo pipefail
cd "$(dirname "$0")"
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE|FORGE_BIN|FORGE_CONVERTER) ;; *) unset $v;; esac; done   # 許した FORGE_* だけ
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES
SRC=run_0183_ns_coldmesh_tw300_ext
REF="--cfl 4 --limiter-ref-from $SRC"
OPT_P=""; OPT_L0="--line dir --itj 5 --cap 50"; OPT_L5="--line dir --itj 5 --lvc 3"
LOG=$PWD/m9.log
fail() { echo "失敗: $*" | tee -a $LOG; exit 1; }
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
gpu_others() { local n; n=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader) || return 1; echo "$n" | awk 'NF' | wc -l; return 0; }   # 取得の失敗 (rc 1) と件数を分ける
modes_ok() {   # modes_ok <run> <P|L0|L5>: ログから実効のモードを確かめる
  local r=$1 m=$2 f=$1/forge_run.log
  grep -q "比較あり" $f && return 1
  case $m in
    P)  ! grep -q "^\[line\] factor 1 回目" $f ;;
    L0) grep -q "^\[line\] factor 1 回目: モード LU$" $f && grep -q "^\[line\] solve 1 回目: モード LU$" $f && ! grep -q "FORGE_LVC_TERMS" $f ;;
    L5) grep -q "^\[line\] factor 1 回目: モード LU$" $f && grep -q "FORGE_LVC_TERMS=5" $f && grep -q "lineViscCoupling: 3" $1/solverConfig.yaml ;;
  esac
}
echo "== 開始 $(date -Is)" >> $LOG
# (1) 専有時間 (再開: 既にある run は、バイナリが lineL で、終了コード 0・モード OK なら使い直す。そうでなければ止める)
SHA=$(python3 -c "import cold_cfl as C;print(C.ALT_BINARIES['lineL_fp64'][0])")
for i in 1 2 3; do
  for m in P L0 L5; do
    r=run_035$((4 + i))_m9unit_${m}_$i
    if [ -e $r ]; then
      if [ "$(cat $r/RUN_RC 2>/dev/null)" = 0 ] && [ "$(cat $r/m9_modes.txt 2>/dev/null)" = OK ] && grep -q "$SHA" $r/COLD_PAIR.json; then echo "$r は使い直す" >> $LOG; continue; fi
      fail "$r が既にあり、使い直せない"
    fi
    eval opt=\$OPT_$m
    g=$(gpu_others) || fail "nvidia-smi が失敗"
    [ "$g" = 0 ] || fail "GPU に他の計算プロセスが $g 本ある — 専有でないので測らない"
    prep $r --steps 1000 --out 1000 $REF $opt || fail "prep $r"
    echo $g > $r/gpu_procs.txt
    t0=$(python3 -c "import time;print(time.time())")
    ( [ $m = L5 ] && export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); rc=$?
    python3 -c "import time;print(time.time()-$t0)" > $r/wall.txt
    if modes_ok $r $m; then echo OK > $r/m9_modes.txt; else echo NG > $r/m9_modes.txt; fi
    echo "$r rc=$rc GPU の他のプロセス $g モード $(cat $r/m9_modes.txt)" >> $LOG
    [ $rc -eq 0 ] || fail $r
    rm -f $r/nozzle.h5 $r/res_1000.h5
  done
done
python3 m9_unit.py >> $LOG 2>&1 || fail "m9_unit (rc $?)"
# (2)(3) L5 の 2 段
phase() {   # phase <run> <line|point> <prep の追加引数...>
  local r=$1 ph=$2; shift 2
  if [ -f $r/m9_watch.json ] && grep -q '"status": "\(REACHED\|CENSORED\|DIVERGED\|EXEC_ERROR\|DATA_ERROR\)"' $r/m9_watch.json; then return 0; fi
  [ -d $r ] || prep $r --steps 200000 --out 5000 $REF --extra res_ro "$@" || fail "prep $r"
  ( [ $ph = line ] && export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo $? > $r/forge_rc ) &
  local job=$!
  sleep 120
  if ! modes_ok $r $([ $ph = line ] && echo L5 || echo P); then
    echo "$r: 実効のモードが期待と違う — 止める" >> $LOG
    python3 m9_stop.py $r || echo "$r: forge を止められない" >> $LOG
    wait $job; fail "$r のモード"
  fi
  python3 m9_watch.py $r --phase $ph --budget 200000 >> $LOG 2>&1; local wrc=$?
  wait $job
  echo "$r 見張り rc=$wrc、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])")" >> $LOG
  python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$? $(grep -o 'NOT CONVERGED\|DIVERGED\|PASS' $r/CONVERGENCE_VERDICT.txt | head -1)" >> $LOG
}
r=run_0353_m9_L5
phase $r line $OPT_L5
if grep -q '"status": "REACHED"' $r/m9_watch.json; then
  n=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['reach_step'])")
  sha=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['reach_sha256'])")
  last=$(ls $r/res_[0-9]*.h5 | sed 's/.*res_//;s/.h5//' | sort -n | tail -1)
  [ "$last" = "$n" ] || fail "$r: 最新の出力 $last ≠ 到達 $n"
  [ "$(sha256sum $r/res_$n.h5 | cut -d' ' -f1)" = "$sha" ] || fail "$r: 到達の出力の sha256 が違う"
  c=run_0354_m9_L5cut
  phase $c point --field-from $r
  echo "$c の出発点: $(python3 -c "import json;print(json.load(open('$c/COLD_PAIR.json'))['parent_res'])")" >> $LOG
fi
python3 m9_judge.py >> $LOG 2>&1 || fail "m9_judge (rc $?)"
echo "== 終了 $(date -Is)" >> $LOG
touch m9.done
