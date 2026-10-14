#!/bin/bash
# plan tooling-nozzle-core-grid §4.21 (事前登録、codex 諮問 2026-10-11・ユーザ承認 2026-10-11): x 方向だけを粗くした Grefx の A/B。
#   A = Gref の継続: run_0607 の res_200000 から同じ格子へ restart_field (9 量・ビット一致を検査)
#   B = Grefx (3086 × 160): run_0607 の res_200000 から interp_field (9 量を検査)。準備は cg_gref.sh と同じ (cold_cfl.py prep → 格子の差し替え)
# 両腕とも ~/forge-fgeom7-fp64・B0 の構成・キー 1・60,000 step の固定長・2500 step ごと。2 本同時。
# 20,000 step 未満の場は次の出力が出てから消す (20k〜60k の 17 出力は残す)。ディスクの空きが 3 GB を切ったら自分の forge だけを止める。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_grefx.sh <A の run> <B の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$C45/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext                                   # cold_cfl.py prep の設定の元 (場は使わない)
PREP64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
GREF=run_0607_cg_gref_k1; SEED=$GREF/res_200000.h5
A=$1; B=$2; MESH_B=$C45/_band_ab/core_grid/prep_Grefx/new64/nozzle.h5
LOG=$C45/cg_grefx.log; MIN_START_GB=15; MIN_RUN_GB=3; STEPS=60000; OUT=2500; KEEP_FROM=20000
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16) 格子 B $(sha256sum $MESH_B | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
for r in $A $B; do [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }; done
[ -f $SEED ] && [ -f $MESH_B ] || { echo "起点か格子 B が無い (中止)" >> $LOG; exit 1; }
python3 -c "import h5py,sys; assert 'PLANES/rSurfVect' in h5py.File(sys.argv[1],'r')" $MESH_B || { echo "格子 B に rSurfVect が無い (中止)" >> $LOG; exit 1; }
# --- A: Gref の継続 (同じ格子 → restart_field)
mkdir $A; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $A/;; esac; done
python3 $TL/restart_field.py $SEED $A/nozzle.h5 --dst-run $A --forge $F7_64 > $A/restart.log 2>&1 || { echo "$A: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $A/restart.log || { echo "$A: 移した量が想定と違う: $(grep VERDICT $A/restart.log)" >> $LOG; exit 1; }
sed -i "s/nStepOuter: 200000}/nStepOuter: $STEPS}/" $A/solverConfig.yaml
grep -q "nStepOuter: $STEPS}" $A/solverConfig.yaml && grep -q "outStepInterval: $OUT" $A/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $A/solverConfig.yaml || { echo "$A: 設定の書き換えに失敗" >> $LOG; exit 1; }
# --- B: Grefx (別の格子 → interp_field)
( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $B --steps $STEPS --out $OUT --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$B: 準備に失敗" >> $LOG; exit 1; }
rm -f $B/nozzle.msh $B/nozzle.h5; cp $MESH_B $B/nozzle.h5
python3 $TL/interp_field.py $SEED $B/nozzle.h5 --dst-run $B --forge $F7_64 > $B/interp_newmesh.log 2>&1 || { echo "$B: 場の移しに失敗" >> $LOG; exit 1; }
python3 - $B/interp_newmesh.log <<'PY' || { echo "$B: 移した量が想定と違う: $(grep -o 'moved.*' $B/interp_newmesh.log)" >> $LOG; exit 1; }
import re, sys
t = open(sys.argv[1]).read(); m = re.search(r"moved \[(.*?)\]", t)
moved = sorted(x.strip().strip("'").replace("(new)", "") for x in m.group(1).split(",")) if m else []
want = sorted(["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"])
if moved != want: raise SystemExit(f"moved {moved} != {want}")
PY
python3 - $B/prepare_info.json $C45/_band_ab/core_grid/prep_Grefx/prepare_info.json <<'PY' || { echo "$B: prepare_info の差し替えに失敗" >> $LOG; exit 1; }
import json, sys
p, q = sys.argv[1:3]; a = json.load(open(p)); b = json.load(open(q))
if abs(a["scale_m"] - b["scale_m"]) > 0: raise SystemExit("scale_m が違う")
a["mesh"] = b["mesh"]; a["core_grid"] = {"grid": "Grefx", "mesh_from": "_band_ab/core_grid/prep_Grefx/new64/nozzle.h5", "ic": "interp_field run_0607/res_200000"}
json.dump(a, open(p, "w"), ensure_ascii=False, indent=1)
PY
grep -c '^mesh: {' $B/solverConfig.yaml | grep -qx 1 || { echo "$B: mesh の形が想定外" >> $LOG; exit 1; }
sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: 1, /" $B/solverConfig.yaml
# 両腕の設定が同じ (格子以外) ことを確かめる
diff $A/solverConfig.yaml $B/solverConfig.yaml > $C45/cg_grefx_config.diff; diff $A/bcondConfig.yaml $B/bcondConfig.yaml >> $C45/cg_grefx_config.diff
[ -s $C45/cg_grefx_config.diff ] && { echo "A と B の設定が違う (中止): $(head -5 $C45/cg_grefx_config.diff | tr '\n' ' ')" >> $LOG; exit 1; }
echo "A と B の solverConfig・bcondConfig は同一" >> $LOG
grep -n '^mesh\|nStepOuter\|outStepInterval\|lineImplicit\|implicitRelax\|extraFields\|cfl_pseudo' $A/solverConfig.yaml | sed "s#^#  $A: #" >> $LOG
PIDS=()
for r in $A $B; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG
done
sleep 150
for r in $A $B; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG; done
while [ ! -f $A/RUN_RC ] || [ ! -f $B/RUN_RC ]; do
  for r in $A $B; do
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue
      [ $n -gt 0 ] && [ $n -lt $KEEP_FROM ] && [ -f $r/res_$((n + OUT)).h5 ] && rm -f $f $r/res_$n.xmf
    done
  done
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $A $B; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
    break
  fi
  sleep 60
done
wait "${PIDS[@]}"
for r in $A $B; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
