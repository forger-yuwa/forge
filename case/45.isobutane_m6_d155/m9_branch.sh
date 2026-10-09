#!/bin/bash
# plan time_integration-line-implicit-speed §6.11 の追加の腕 L5L (2026-10-10、ユーザの問い「point のほうが粘性込み line より収束は早いと思ってるの？」、到達の前に登録):
# run_0353_m9_L5 が水準に到達したら、同じ出力から同じ構成 (値 3・マスク 5・方向別・上限なし・キー 5) のライン run_0363_m9_L5line を E (水準 + 参照から 0.1 % 以内) まで回す。
# m9.sh の point の仕上げ (run_0354_m9_L5cut) と同時に走る (総時間は専有の単価で換算するので同時でよい)。終わったら m9_judge.py をもう一度回す。
set -uo pipefail
cd "$(dirname "$0")"
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE|FORGE_BIN|FORGE_CONVERTER) ;; *) unset $v;; esac; done
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/m9_branch.log
fail() { echo "失敗: $*" | tee -a $LOG; exit 1; }
r=run_0353_m9_L5; b=run_0363_m9_L5line
echo "== 開始 $(date -Is) ($r の終わりを待つ)" >> $LOG
until [ -f $r/m9_watch.json ] && grep -q '"status": "\(REACHED\|CENSORED\|DIVERGED\|EXEC_ERROR\|DATA_ERROR\)"' $r/m9_watch.json; do sleep 30; done
grep -q '"status": "REACHED"' $r/m9_watch.json || { echo "$r が水準に届かなかった — 追加の腕は回さない" >> $LOG; exit 0; }
n=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['reach_step'])")
sha=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['reach_sha256'])")
sleep 10
last=$(ls $r/res_[0-9]*.h5 | sed 's/.*res_//;s/.h5//' | sort -n | tail -1)
[ "$last" = "$n" ] || fail "$r: 最新の出力 $last ≠ 到達 $n"
[ "$(sha256sum $r/res_$n.h5 | cut -d' ' -f1)" = "$sha" ] || fail "$r: 到達の出力の sha256 が違う"
[ -e $b ] && fail "$b が既にある"
python3 cold_cfl.py prep $SRC $b --steps 200000 --out 5000 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --lvc 3 --extra res_ro --field-from $r > /dev/null || fail "prep $b"
rm -f $b/nozzle.msh
echo "$b の出発点: $(python3 -c "import json;print(json.load(open('$b/COLD_PAIR.json'))['parent_res'])")" >> $LOG
( export FORGE_LVC_TERMS=5; python3 cold_cfl.py run $b > $b/cold_pair_run_stdout.log 2>&1; echo $? > $b/forge_rc ) &
job=$!
sleep 120
f=$b/forge_run.log
if grep -q "比較あり" $f || ! grep -q "^\[line\] factor 1 回目: モード LU$" $f || ! grep -q "FORGE_LVC_TERMS=5" $f || ! grep -q "lineViscCoupling: 3" $b/solverConfig.yaml; then
  python3 m9_stop.py $b; wait $job; fail "$b の実効のモードが期待と違う"
fi
python3 m9_watch.py $b --phase line_e2 --budget 200000 >> $LOG 2>&1; wrc=$?
wait $job
echo "$b 見張り rc=$wrc、forge rc=$(cat $b/forge_rc 2>/dev/null)、状態 $(python3 -c "import json;print(json.load(open('$b/m9_watch.json'))['status'])")" >> $LOG
python3 ../../solver_density_cuda/tools/check_convergence.py $b > $b/CONVERGENCE_VERDICT.txt 2>&1; echo "$b check_convergence rc=$? $(grep -o 'NOT CONVERGED\|DIVERGED\|PASS' $b/CONVERGENCE_VERDICT.txt | head -1)" >> $LOG
until [ -f m9.done ] || grep -q "^失敗" m9.log; do sleep 30; done
python3 m9_judge.py >> $LOG 2>&1 || fail "m9_judge (rc $?)"
echo "== 終了 $(date -Is)" >> $LOG
touch m9_branch.done
