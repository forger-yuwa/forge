#!/bin/bash
# plan time_integration-line-implicit-speed §6.18 (2026-10-10、diagnostician の推奨): 総時間の 3 腕 (S2・M64・M64+S2)。バイナリ lineL。
# (1) M64+S2 の単価を B0 で挟んで 3 本 (専有) → (2) 損益分岐の予算を決める → (3) 3 腕を同時に run_0183 の res_100000 から水準 (2 出力連続) まで → (4) 判定。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE|FORGE_BIN|FORGE_CONVERTER) ;; *) unset $v;; esac; done
unset FORGE_ALLOW_UNVERIFIED_SPECIES
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50"
LOG=$PWD/tt.log
fail() { echo "失敗: $*" | tee -a $LOG; touch tt.done; exit 1; }
diskok() { [ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 3221225472 ]; }
gpucount() { local o; o=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null) || return 1; echo "$o" | grep -c '[0-9]'; true; }
declare -A ENVS=([B0]="" [M64S2]="FORGE_LINE_MAXLEN=64" [S2]="" [M64]="FORGE_LINE_MAXLEN=64")
declare -A OPTS=([B0]="" [M64S2]="--inner 2" [S2]="--inner 2" [M64]="")
declare -A COV=([B0]=570999 [M64S2]=302016 [S2]=570999 [M64]=302016)
declare -A MLEN=([B0]=121 [M64S2]=64 [S2]=121 [M64]=64)
declare -A INNER=([B0]=5 [M64S2]=2 [S2]=2 [M64]=5)
echo "== 開始 $(date -Is)" >> $LOG
# (1) 単価 (§6.16 と同じ手順・確認)
unit() {   # unit <run> <案>
  local r=$1 v=$2 n
  diskok || fail "ディスクの空きが 3 GiB 未満"
  python3 cold_cfl.py prep $SRC $r --steps 1000 --out 1000 $B ${OPTS[$v]} > /dev/null || fail "prep $r"; rm -f $r/nozzle.msh
  n=$(gpucount) && echo "ok $n" > $r/gpu_pre.txt || echo "fail" > $r/gpu_pre.txt
  [ "$(cat $r/gpu_pre.txt)" = "ok 0" ] || fail "$r: 投入前の GPU が専有でないか照会に失敗 ($(cat $r/gpu_pre.txt))"
  ( export FORGE_LINE_LAYOUT=2; [ -n "${ENVS[$v]}" ] && export "${ENVS[$v]}"; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ) &
  local job=$! mx=0 nf=0 c
  while kill -0 $job 2>/dev/null; do if c=$(gpucount); then [ "$c" -gt "$mx" ] && mx=$c; else nf=$((nf + 1)); fi; sleep 5; done
  wait $job; local rc=$?; echo "done $nf $mx" > $r/gpu_watch.txt
  echo "$r rc=$rc $v GPU 見張り (失敗 $nf・最大 $mx 本)" >> $LOG
  python3 scr_check.py $r ${COV[$v]} ${MLEN[$v]} ${INNER[$v]} >> $LOG 2>&1
  rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1000.h5
}
k=0
for v in B0 M64S2 B0 M64S2 B0 M64S2 B0; do k=$((k + 1)); unit run_0376_ttu_$(printf %02d $k)_$v $v; done
# (2) 予算
python3 tt_judge.py budgets > tt_budgets.txt 2>> $LOG || fail "予算の計算 (rc $?): $(tail -1 tt_budgets.txt)"
cat tt_budgets.txt >> $LOG
declare -A BUD; while read -r k n; do BUD[$k]=$n; done < tt_budgets.txt
for v in S2 M64 M64S2; do [ -n "${BUD[$v]:-}" ] && [ "${BUD[$v]}" -gt 125000 ] || fail "$v の予算が読めない・B0 以下 (${BUD[$v]:-なし})"; done
# (3) 3 腕を同時に
PIDS=()
arm() {   # arm <run> <案>
  local r=$1 v=$2 CUR=""
  trap '[ -n "$CUR" ] && python3 m9_stop.py $CUR; exit 1' TERM INT HUP
  trap '[ -n "$CUR" ] && python3 m9_stop.py $CUR' EXIT
  diskok || { echo "$r: ディスクの空きが 3 GiB 未満 — 起動しない" >> $LOG; return 1; }
  python3 cold_cfl.py prep $SRC $r --steps ${BUD[$v]} --out 2500 $B ${OPTS[$v]} --extra res_ro > /dev/null || { echo "$r prep 失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  CUR=$r
  ( export FORGE_LINE_LAYOUT=2; [ -n "${ENVS[$v]}" ] && export "${ENVS[$v]}"; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo $? > $r/forge_rc ) &
  local job=$!
  python3 m9_watch.py $r --phase line --budget ${BUD[$v]} --consec 2 >> $LOG 2>&1 &
  local wjob=$!
  sleep 180
  local f=$r/forge_run.log ok=1
  grep -q "^\[line\] factor 1 回目: モード LAYOUT2$" $f || ok=0
  grep -q "lines=4719  covered CVs=${COV[$v]}/570999 .*maxLen=${MLEN[$v]}" $f || ok=0
  grep -Eq "nStepInner: *${INNER[$v]}([,} ]|$)" $r/solverConfig.yaml || ok=0
  grep -q "比較あり" $f && ok=0
  [ $ok = 1 ] || { echo "$r: 実効のモード・被覆・nStepInner が期待と違う — 止める" >> $LOG; python3 m9_stop.py $r; }
  wait $wjob; local wrc=$?
  [ $wrc -eq 0 ] || { python3 m9_stop.py $r; echo "$r: 見張りが rc=$wrc で終わった — forge を止めた" >> $LOG; }
  wait $job; CUR=""
  echo "$r 見張り rc=$wrc、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])" 2>/dev/null)" >> $LOG
  python3 ../../solver_density_cuda/tools/check_convergence.py $r > $r/CONVERGENCE_VERDICT.txt 2>&1; echo "$r check_convergence rc=$? $(grep -o 'NOT CONVERGED\|DIVERGED\|PASS' $r/CONVERGENCE_VERDICT.txt | head -1)" >> $LOG
  rm -f $r/nozzle.h5 $r/res_0.h5
}
diskguard() {   # 空きが 1 GB を切ったら自分の 3 腕の forge を止める。arm が全部終わったら自分で抜ける (m9r.sh の wait の欠陥を避ける)
  while :; do
    local alive=0; for p in "$@"; do kill -0 $p 2>/dev/null && alive=1; done
    [ $alive = 1 ] || return 0
    if [ "$(df --output=avail -BG . | tail -1 | tr -dc 0-9)" -lt 1 ]; then
      echo "ディスクの空き < 1 GB — 自分の run の forge を止める $(date -Is)" >> $LOG
      for r in run_0377_tt_S2 run_0378_tt_M64 run_0379_tt_M64S2; do [ -d $r ] && python3 m9_stop.py $r; done; return 0
    fi
    sleep 60
  done
}
arm run_0377_tt_S2 S2 & PIDS+=($!); sleep 240
arm run_0378_tt_M64 M64 & PIDS+=($!); sleep 240
arm run_0379_tt_M64S2 M64S2 & PIDS+=($!)
diskguard "${PIDS[@]}" &
DG=$!
for p in "${PIDS[@]}"; do wait $p; done
wait $DG
# (4) 判定
python3 tt_judge.py final >> $LOG 2>&1 || echo "tt_judge rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch tt.done
