#!/bin/bash
# plan time_integration-line-implicit-speed §6.13 (2026-10-10、ユーザ「0.85 と 1.0 を値 0 ライン・粘性ラインで回し総時間比較、並行してやれたら」):
# implicitRelax だけを変えた 4 腕 (値 0 + 方向別 + キー 5 + 上限 50 / 値 3・マスク 5 + 方向別・上限なし + キー 5) × (0.85 / 1.0) を並行で回す。
# 各腕: ライン段 (緩和 r) を水準まで → その出力から point cfl 4・キー 0・緩和 0.7 (P と同じ仕上げ) を E (水準 + 参照から 0.1 % 以内) まで。バイナリ lineL。
set -uo pipefail
cd "$(dirname "$0")"
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE|FORGE_BIN|FORGE_CONVERTER) ;; *) unset $v;; esac; done
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES
SRC=run_0183_ns_coldmesh_tw300_ext
REF="--cfl 4 --limiter-ref-from $SRC"
LOG=$PWD/m9r.log
modes_ok() {   # modes_ok <run> <P|L0|L5> <relax>
  local r=$1 m=$2 rx=$3 f=$1/forge_run.log
  grep -q "比較あり" $f && return 1
  grep -Eq "implicitRelax: *$rx([,} ]|$)" $r/solverConfig.yaml || return 1
  case $m in
    P)  ! grep -q "^\[line\] factor 1 回目" $f ;;
    L0) grep -q "^\[line\] factor 1 回目: モード LU$" $f && ! grep -q "FORGE_LVC_TERMS" $f ;;
    L5) grep -q "^\[line\] factor 1 回目: モード LU$" $f && grep -q "FORGE_LVC_TERMS=5" $f && grep -q "lineViscCoupling: 3" $1/solverConfig.yaml ;;
  esac
}
phase() {   # phase <run> <line|point> <L0|L5|P> <relax> <prep の追加引数...>  → 0 = 見張りが最終の状態まで進んだ
  local r=$1 ph=$2 m=$3 rx=$4; shift 4
  if [ -f $r/m9_watch.json ] && grep -q '"status": "\(REACHED\|CENSORED\|DIVERGED\|EXEC_ERROR\|DATA_ERROR\)"' $r/m9_watch.json; then return 0; fi
  if [ ! -d $r ]; then python3 cold_cfl.py prep $SRC $r --steps 200000 --out 5000 $REF --extra res_ro "$@" > /dev/null || { echo "$r prep 失敗" >> $LOG; return 1; }; rm -f $r/nozzle.msh; fi
  ( [ $m = L5 ] && export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo $? > $r/forge_rc ) &
  local job=$!
  sleep 180
  if ! modes_ok $r $m $rx; then echo "$r: 実効のモード・緩和が期待と違う — 止める" >> $LOG; python3 m9_stop.py $r; wait $job; return 1; fi
  python3 m9_watch.py $r --phase $ph --budget 200000 >> $LOG 2>&1; local wrc=$?
  [ $wrc -eq 0 ] || { python3 m9_stop.py $r; echo "$r: 見張りが rc=$wrc で終わった — forge を止めた" >> $LOG; }
  wait $job
  echo "$r 見張り rc=$wrc、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])" 2>/dev/null)" >> $LOG
  python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$? $(grep -o 'NOT CONVERGED\|DIVERGED\|PASS' $r/CONVERGENCE_VERDICT.txt | head -1)" >> $LOG
  return 0
}
arm() {   # arm <line の run> <point の run> <L0|L5> <relax> <line の構成...>
  local rl=$1 rp=$2 m=$3 rx=$4; shift 4
  phase $rl line $m $rx "$@" --relax $rx || { echo "$rl: 失敗" >> $LOG; return; }
  grep -q '"status": "REACHED"' $rl/m9_watch.json || { echo "$rl: ライン段は水準に届かず" >> $LOG; return; }
  local n sha last
  n=$(python3 -c "import json;print(json.load(open('$rl/m9_watch.json'))['reach_step'])")
  sha=$(python3 -c "import json;print(json.load(open('$rl/m9_watch.json'))['reach_sha256'])")
  last=$(ls $rl/res_[0-9]*.h5 | sed 's/.*res_//;s/.h5//' | sort -n | tail -1)
  [ "$last" = "$n" ] || { echo "$rl: 最新の出力 $last ≠ 到達 $n — point の段は回さない" >> $LOG; return; }
  [ "$(sha256sum $rl/res_$n.h5 | cut -d' ' -f1)" = "$sha" ] || { echo "$rl: 到達の出力の sha256 が違う" >> $LOG; return; }
  phase $rp point P 0.7 --field-from $rl || echo "$rp: 失敗" >> $LOG
}
echo "== 開始 $(date -Is)" >> $LOG
OPT_L0="--line dir --itj 5 --cap 50"; OPT_L5="--line dir --itj 5 --lvc 3"
arm run_0364_m9r_L0_r085 run_0365_m9r_L0_r085cut L0 0.85 $OPT_L0 &
sleep 240
arm run_0366_m9r_L0_r100 run_0367_m9r_L0_r100cut L0 1.0 $OPT_L0 &
sleep 240
arm run_0368_m9r_L5_r085 run_0369_m9r_L5_r085cut L5 0.85 $OPT_L5 &
sleep 240
arm run_0370_m9r_L5_r100 run_0371_m9r_L5_r100cut L5 1.0 $OPT_L5 &
wait
python3 m9r_judge.py >> $LOG 2>&1 || echo "m9r_judge rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch m9r.done
