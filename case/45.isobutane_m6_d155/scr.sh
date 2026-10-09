#!/bin/bash
# plan time_integration-line-implicit-speed §6.16 (2026-10-10): 1 step の時間のふるい (ライン長の上限 FORGE_LINE_MAXLEN・sweep の回数 nStepInner) と
# 部分被覆での LAYOUT2 のビット一致 (§5.1 #13 の一部)。バイナリ lineL、GPU 専有。基準 B0 = 値 0 + 方向別 + キー 5 + 上限 50 + LAYOUT2。
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
prep() { python3 cold_cfl.py prep $SRC "$@" > /dev/null || return 1; rm -f $1/nozzle.msh; }
run() { local r=$1; shift; ( for kv in "$@"; do export "$kv"; done; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?; echo "$r rc=$rc $*" >> $LOG; return $rc; }
echo "== 開始 $(date -Is)" >> $LOG
[ "$(df --output=avail -BG ~ | tail -1 | tr -dc 0-9)" -ge 3 ] || fail "ディスクの空きが 3 GB 未満"
# (1) 部分被覆 (ライン長 64) での LAYOUT2 のビット一致
r=run_0372_ml64_cmp
prep $r --steps 20 --out 20 $B || fail "prep $r"
run $r FORGE_LINE_COMPARE=1 FORGE_LINE_LAYOUT=2 FORGE_LINE_MAXLEN=64 || fail $r
grep -q "lines=4719  covered CVs=302016/570999 .* maxLen=64" $r/forge_run.log || fail "$r: ライン長 64 の被覆になっていない ($(grep -m1 'lineImplicit\] lines' $r/forge_run.log))"
python3 line_cmp_judge.py $r --arm LAYOUT2 --eta-limit 1e-11 --bitwise --factors 20 --solves 100 --lines 4719 >> $LOG 2>&1; c1=$?
echo "(1) 判定 rc=$c1 (0 PASS / 1 FAIL / 2 INDETERMINATE)" >> $LOG
rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_20.h5
# (2) 1 step の時間 (3 巡、各巡で B0 と 5 案を交互に)
declare -A ENVS=([B0]="" [M100]="FORGE_LINE_MAXLEN=100" [M64]="FORGE_LINE_MAXLEN=64" [M48]="FORGE_LINE_MAXLEN=48" [S3]="" [S2]="")
declare -A OPTS=([B0]="" [M100]="" [M64]="" [M48]="" [S3]="--inner 3" [S2]="--inner 2")
declare -A COV=([B0]="570999/570999 .* maxLen=121" [M100]="471900/570999 .* maxLen=100" [M64]="302016/570999 .* maxLen=64" [M48]="226512/570999 .* maxLen=48" [S3]="570999/570999 .* maxLen=121" [S2]="570999/570999 .* maxLen=121")
declare -A INNER=([B0]=5 [M100]=5 [M64]=5 [M48]=5 [S3]=3 [S2]=2)
for i in 1 2 3; do
  for v in B0 M100 M64 M48 S3 S2; do
    r=run_037$((2 + i))_scr_${v}_$i
    prep $r --steps 1000 --out 1000 $B ${OPTS[$v]} || fail "prep $r"
    nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l > $r/gpu_procs.txt
    echo "$r 投入前の GPU の計算プロセス: $(cat $r/gpu_procs.txt) 本" >> $LOG
    run $r FORGE_LINE_LAYOUT=2 ${ENVS[$v]} || fail $r
    ok=OK
    grep -q "lines=4719  covered CVs=${COV[$v]}" $r/forge_run.log || ok="被覆が違う: $(grep -m1 'lineImplicit\] lines' $r/forge_run.log)"
    grep -Eq "nStepInner: *${INNER[$v]}([,} ]|$)" $r/solverConfig.yaml || ok="$ok nStepInner が ${INNER[$v]} でない"
    echo "$ok" > $r/scr_modes.txt
    rm -f $r/nozzle.h5 $r/res_0.h5 $r/res_1000.h5
  done
done
for v in M100 M64 M48 S3 S2; do
  python3 time_pairs_judge.py _band_ab/cold_pair/scr_time_$v.json LAYOUT2 LAYOUT2 run_0373_scr_B0_1 run_0373_scr_${v}_1 run_0374_scr_B0_2 run_0374_scr_${v}_2 run_0375_scr_B0_3 run_0375_scr_${v}_3 >> $LOG 2>&1 || echo "$v: 組の判定器 rc=$?" >> $LOG
done
python3 scr_judge.py >> $LOG 2>&1 || echo "scr_judge rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch scr.done
