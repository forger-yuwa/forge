#!/bin/bash
# plan architecture-float-state-double-geometry §6.20 (事前登録、§5.1 #15 (1)): hoop の修正 (キー 1) を入れた格子で、float を V4 (§6.11) と同じ停止規則で回す。
# 比べる FP64 は hp7.sh の run_0483_hp7_k1 (同じソース commit 0b4dff4e・同じ格子・同じ出発点・同じ B0 の実効設定)。格子の差し替えは hp7.sh と同じ。
# V4 との違い: 格子 (区間ごとの r 重み・半割面の端点)、キー 1、commit の診断 (FORGE_DIAG_COMMIT_LOSS) なし (FP64 の腕に合わせた)。判定は fl1_judge.py。
set -uo pipefail
cd ~/forge-wallfit/case/45.isobutane_m6_d155
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools; W=$PWD/_conv9
F7_64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; F7_32=~/forge-fgeom7-f32/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$PWD/fl1.log
echo "== 開始 $(date -Is) F7_32 $(sha256sum $F7_32 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 10737418240 ] || { echo "ディスクの空きが 10 GiB 未満 (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <キー>
  local r=$1 k=$2
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps 200000 --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$r: 準備に失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh; mv $r/nozzle.h5 $r/nozzle_oldmesh.h5; cp $W/c45/new64/nozzle.h5 $r/nozzle.h5
  python3 - $r/nozzle_oldmesh.h5 $r/nozzle.h5 <<'EOF' || { echo "$r: roY* のデータセットの作成に失敗" >> $LOG; return 1; }
import sys, re, h5py
with h5py.File(sys.argv[1], "r") as s, h5py.File(sys.argv[2], "r+") as d:
    ks = sorted(k for k in s["VALUE"] if re.fullmatch(r"roY\d+", k))
    if not ks: raise SystemExit("SRC に roY* が無い")
    for k in ks:
        if k not in d["VALUE"]: d["VALUE"].create_dataset(k, shape=s["VALUE/" + k].shape, dtype=s["VALUE/" + k].dtype)
EOF
  python3 $TL/restart_field.py $r/nozzle_oldmesh.h5 $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart_newmesh.log 2>&1 || { echo "$r: 場の移しに失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した" $r/restart_newmesh.log || { echo "$r: 移した量が 9 でない: $(grep VERDICT $r/restart_newmesh.log)" >> $LOG; return 1; }
  rm -f $r/nozzle_oldmesh.h5
  grep -c '^mesh: {' $r/solverConfig.yaml | grep -qx 1 || { echo "$r: mesh の形が想定外" >> $LOG; return 1; }
  sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: $k, /" $r/solverConfig.yaml
  grep -n '^mesh\|nStepOuter\|outStepInterval\|lineImplicit\|implicitRelax\|extraFields' $r/solverConfig.yaml | sed "s#^#  $r: #" >> $LOG
  ( export FORGE_BIN=$F7_32; bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  python3 m9_watch.py $r --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $r/m9_watch.log 2>&1 &
  echo "$r 起動 $(date -Is) キー $k" >> $LOG
}
arm run_0484_fl1_f32_k1 1 || exit 1
sleep 120
for r in run_0484_fl1_f32_k1; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG; done
wait
for r in run_0484_fl1_f32_k1; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 見張り=$(python3 -c "import json; s=json.load(open('$r/m9_watch.json')); print(s['status'], s.get('reach_step'))")" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
touch fl1.done
