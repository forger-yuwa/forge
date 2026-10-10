#!/bin/bash
# plan tooling-nozzle-core-grid §5.1 #4 と plan axisymmetric-freestream-hoop-gauge §4.6 の 7 (ユーザ決定 2026-10-10: キー 0/1 の比較を core-grid の格子で行う)。
# hp7.sh と同じ構成 (FP64 の ~/forge-fgeom7-fp64、B0 の構成、段階起動なし、水準 2 出力連続・上限 20 万 step・2500 step ごと) で、
# 格子 G1 / G1x の上でキー mesh.axisSegmentRWeight 0 と 1 を回す。格子の比較の G0′ は別セッションの run_0483_hp7_k1 (G0・キー 1)。
# hp7.sh との違い: 格子が違うので、場は restart_field ではなく interp_field (run_0183 の res_100000 から、roY* は interp_field が作る。移した量 9 を検査)、
# prepare_info.json の格子の情報を新しい格子のものに替える (見張りの cold_series が ni を読む)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_runs.sh <G1|G1x> <キー 1 の run> <キー 0 の run> [<格子> <run_k1> <run_k0> ...]
# 格子ごとに移送は 1 回 (キー 1 の run で行う)。キー 0 の run はその複製で、初期 HDF5 の sha256 の一致を確かめる (codex 諮問 2026-10-10)。全腕を準備してから一斉に起動する。
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$C45/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; SRC_RES=$SRC/res_100000.h5
PREP64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$C45/cg_runs.log; MIN_START_GB=8; MIN_RUN_GB=3
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16)" >> $LOG
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
PIDS=()
prep_arm() {   # prep_arm <run> <格子> <キー>  (準備だけ、起動しない)
  local r=$1 g=$2 k=$3 mesh=$C45/_band_ab/core_grid/prep_$2/new64/nozzle.h5
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  [ -f $mesh ] || { echo "$r: 格子 $mesh が無い (中止)" >> $LOG; return 1; }
  python3 -c "import h5py,sys; assert 'PLANES/rSurfVect' in h5py.File(sys.argv[1],'r')" $mesh || { echo "$r: 格子に rSurfVect が無い (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps 200000 --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$r: 準備に失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh $r/nozzle.h5; cp $mesh $r/nozzle.h5
  python3 $TL/interp_field.py $SRC_RES $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/interp_newmesh.log 2>&1 || { echo "$r: 場の移しに失敗" >> $LOG; return 1; }
  python3 - $r/interp_newmesh.log <<'EOF' || { echo "$r: 移した量が想定と違う: $(grep -o 'moved.*' $r/interp_newmesh.log)" >> $LOG; return 1; }
import re, sys
t = open(sys.argv[1]).read(); m = re.search(r"moved \[(.*?)\]", t)
moved = sorted(x.strip().strip("'").replace("(new)", "") for x in m.group(1).split(",")) if m else []
want = sorted(["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"])
if moved != want: raise SystemExit(f"moved {moved} != {want}")
EOF
  python3 - $r/prepare_info.json $C45/_band_ab/core_grid/prep_$g/prepare_info.json $g <<'EOF' || { echo "$r: prepare_info の差し替えに失敗" >> $LOG; return 1; }
import json, sys
p, q, g = sys.argv[1:4]; a = json.load(open(p)); b = json.load(open(q))
if abs(a["scale_m"] - b["scale_m"]) > 0: raise SystemExit("scale_m が違う")
a["mesh"] = b["mesh"]; a["core_grid"] = {"grid": g, "mesh_from": f"_band_ab/core_grid/prep_{g}/new64/nozzle.h5", "ic": "interp_field run_0183/res_100000"}
json.dump(a, open(p, "w"), ensure_ascii=False, indent=1)
EOF
  grep -c '^mesh: {' $r/solverConfig.yaml | grep -qx 1 || { echo "$r: mesh の形が想定外" >> $LOG; return 1; }
  sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: $k, /" $r/solverConfig.yaml
  grep -n '^mesh\|nStepOuter\|outStepInterval\|lineImplicit\|implicitRelax\|extraFields' $r/solverConfig.yaml | sed "s#^#  $r: #" >> $LOG
}
copy_arm() {   # copy_arm <元の run> <新しい run> <キー>  (同じ初期 HDF5 でキーだけ替える)
  local a=$1 b=$2 k=$3
  [ -e $b ] && { echo "$b: 既にある (中止)" >> $LOG; return 1; }
  mkdir $b; for f in $a/*; do case $(basename $f) in *.log|RUN_*|forge_*|residual_*|m9_watch*) ;; *) cp -a $f $b/;; esac; done
  sed -i "s/^mesh: {axisSegmentRWeight: [01], /mesh: {axisSegmentRWeight: $k, /" $b/solverConfig.yaml
  grep -q "^mesh: {axisSegmentRWeight: $k, " $b/solverConfig.yaml || { echo "$b: キーの書き換えに失敗" >> $LOG; return 1; }
  [ "$(sha256sum < $a/nozzle.h5)" = "$(sha256sum < $b/nozzle.h5)" ] || { echo "$b: 初期 HDF5 が $a と違う" >> $LOG; return 1; }
  diff <(sed 's/axisSegmentRWeight: [01], //' $a/solverConfig.yaml) <(sed 's/axisSegmentRWeight: [01], //' $b/solverConfig.yaml) > /dev/null || { echo "$b: キー以外の設定が違う" >> $LOG; return 1; }
  echo "$b: $a の複製 (初期 HDF5 の sha256 一致、キー $k)" >> $LOG
}
launch() {   # launch <run>
  local r=$1
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!)
  python3 m9_watch.py $r --phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9 >> $r/m9_watch.log 2>&1 &
  echo "$r 起動 $(date -Is) $(grep -o 'axisSegmentRWeight: [01]' $r/solverConfig.yaml)" >> $LOG
}
RUNS=()
while [ $# -ge 3 ]; do
  prep_arm $2 $1 1 || exit 1
  copy_arm $2 $3 0 || exit 1
  RUNS+=($2 $3); shift 3
done
for r in "${RUNS[@]}"; do launch $r; done
sleep 150
for r in "${RUNS[@]}"; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-140)" >> $LOG; done
# ディスクの見張り: 空きが MIN_RUN_GB を切ったら、自分の run の forge だけを止める (cwd が自分の run のもの)
( while true; do
    alive=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || alive=1; done; [ $alive -eq 0 ] && break
    if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
      for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in "${RUNS[@]}"; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
      break
    fi
    sleep 60
  done ) &
wait "${PIDS[@]}"
for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 見張り=$(python3 -c "import json; s=json.load(open('$r/m9_watch.json')); print(s['status'], s.get('reach_step'))" 2>/dev/null)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
