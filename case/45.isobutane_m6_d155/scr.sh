#!/bin/bash
# plan time_integration-line-implicit-speed §6.16 (2026-10-10、codex plan-6 の反映): 1 step の時間のふるい (ライン長の上限 FORGE_LINE_MAXLEN・sweep の回数 nStepInner) と
# 部分被覆での LAYOUT2 のビット一致 (§5.1 #13 の一部)。バイナリ lineL、GPU 専有。基準 B0 = 値 0 + 方向別 + キー 5 + 上限 50 + LAYOUT2。
# 各巡で B0 と案を交互に並べ (B0 v B0 v … B0)、各案を前後の B0 で挟む。案の順は巡ごとに変える (順序効果をならす)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE|FORGE_BIN|FORGE_CONVERTER) ;; *) unset $v;; esac; done
unset FORGE_ALLOW_UNVERIFIED_SPECIES
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50"
LOG=$PWD/scr.log
fail() { echo "失敗: $*" | tee -a $LOG; touch scr.done; exit 1; }
diskok() { [ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 3221225472 ]; }   # 出力先の空き ≥ 3 GiB (バイト単位、codex plan-6 m1)
prep() { diskok || return 2; python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
gpucount() { local o; o=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null) || return 1; echo "$o" | grep -c '[0-9]'; true; }
echo "== 開始 $(date -Is)" >> $LOG
# (1) 部分被覆 (ライン長 64) での LAYOUT2 のビット一致
r=run_0372_ml64_cmp
prep $r --steps 20 --out 20 $B || fail "prep $r"
( export FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=2 FORGE_LINE_MAXLEN=64; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); echo "$r rc=$?" >> $LOG
python3 line_cmp_judge.py $r --arm LAYOUT2 --eta-limit 1e-11 --bitwise --factors 20 --solves 100 --lines 4719 >> $LOG 2>&1; c1=$?
grep -q "lines=4719  covered CVs=302016/570999 .*maxLen=64" $r/forge_run.log || echo "$r: ライン長 64 の被覆になっていない" >> $LOG
echo "(1) 判定 rc=$c1 (0 PASS / 1 FAIL / 2 INDETERMINATE)" >> $LOG
rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_20.h5
# (2) 1 step の時間
declare -A ENVS=([B0]="" [M100]="FORGE_LINE_MAXLEN=100" [M64]="FORGE_LINE_MAXLEN=64" [M48]="FORGE_LINE_MAXLEN=48" [S3]="" [S2]="")
declare -A OPTS=([B0]="" [M100]="" [M64]="" [M48]="" [S3]="--inner 3" [S2]="--inner 2")
declare -A COV=([B0]=570999 [M100]=471900 [M64]=302016 [M48]=226512 [S3]=570999 [S2]=570999)
declare -A MLEN=([B0]=121 [M100]=100 [M64]=64 [M48]=48 [S3]=121 [S2]=121)
declare -A INNER=([B0]=5 [M100]=5 [M64]=5 [M48]=5 [S3]=3 [S2]=2)
ORD=("M100 M64 M48 S3 S2" "S2 S3 M48 M64 M100" "M64 S2 M100 S3 M48")
one() {   # one <run> <案>
  local r=$1 v=$2 n
  prep $r --steps 1000 --out 1000 $B ${OPTS[$v]} || fail "prep $r (rc $?、2 = ディスク)"
  n=$(gpucount) && echo "ok $n" > $r/gpu_pre.txt || echo "fail" > $r/gpu_pre.txt
  [ "$(cat $r/gpu_pre.txt)" = "ok 0" ] || fail "$r: 投入前の GPU が専有でないか照会に失敗 ($(cat $r/gpu_pre.txt))"
  ( export FORGE_LINE_LAYOUT=2; [ -n "${ENVS[$v]}" ] && export "${ENVS[$v]}"; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ) &
  local job=$! mx=0 nf=0 c
  while kill -0 $job 2>/dev/null; do                       # 計測中の見張り (codex plan-6 M2)
    if c=$(gpucount); then [ "$c" -gt "$mx" ] && mx=$c; else nf=$((nf + 1)); fi
    sleep 5
  done
  wait $job; local rc=$?
  echo "done $nf $mx" > $r/gpu_watch.txt
  echo "$r rc=$rc $v GPU 見張り (失敗 $nf・最大 $mx 本)" >> $LOG
  python3 scr_check.py $r ${COV[$v]} ${MLEN[$v]} ${INNER[$v]} >> $LOG 2>&1
  rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1000.h5
}
for i in 1 2 3; do
  k=1; one run_037$((2 + i))_scr_r${i}_01_B0 B0
  for v in ${ORD[$((i - 1))]}; do
    k=$((k + 1)); one run_037$((2 + i))_scr_r${i}_$(printf %02d $k)_$v $v
    k=$((k + 1)); one run_037$((2 + i))_scr_r${i}_$(printf %02d $k)_B0 B0
  done
done
python3 scr_judge.py >> $LOG 2>&1 || echo "scr_judge rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch scr.done
