#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.6 の 7 (事前登録、既定化の条件・非静止場): case/45 の B0 の構成 (FP64) で、
# 同じバイナリ (~/forge-fgeom7-fp64 = commit 0b4dff4e と同じソース)・同じ格子 (_conv9/c45/new64/nozzle.h5)・同じ保存量 (run_0183 の res_100000 を restart_field.py) から、
# キー mesh.axisSegmentRWeight 0 と 1 を同じ停止規則 (float 化 plan §6.11: 水準を 2 出力連続、上限 20 万 step、2500 step ごと) で回す。2 本を同時に投入する。
# 格子の差し替えは fg9d.sh と同じ (roY* を先に作り、移した量 9 を検査)。判定は hp7_judge.py。
set -uo pipefail
cd ~/forge-wallfit/case/45.isobutane_m6_d155
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools; W=$PWD/_conv9
F7_64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$PWD/hp7.log
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16)" >> $LOG
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
  ( export FORGE_BIN=$F7_64; bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  python3 m9_watch.py $r --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $r/m9_watch.log 2>&1 &
  echo "$r 起動 $(date -Is) キー $k" >> $LOG
}
arm run_0482_hp7_k0 0 || exit 1
arm run_0483_hp7_k1 1 || exit 1
sleep 120
for r in run_0482_hp7_k0 run_0483_hp7_k1; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG; done
wait
for r in run_0482_hp7_k0 run_0483_hp7_k1; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 見張り=$(python3 -c "import json; s=json.load(open('$r/m9_watch.json')); print(s['status'], s.get('reach_step'))")" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
touch hp7.done
