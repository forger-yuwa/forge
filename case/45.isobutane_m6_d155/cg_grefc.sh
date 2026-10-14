#!/bin/bash
# plan tooling-nozzle-core-grid §4.36 (事前登録、codex 諮問 2026-10-11、ユーザ決定「すいしょうでいいよ」): 縮流部 x∈[−8, −2] r_t だけ流れ方向を粗くした Grefc と、
# Gref の B0 の継続 (対照) を、同じ起点 (run_0608 の res_60000) から各 60,000 step。B0 の構成 (ライン dir・ITJ 5・上限 50・cfl 4・緩和 0.7・キー 1)。
#   A = 対照: 同じ格子へ restart_field (9 量・ビット一致)。B = Grefc: prep_Grefc/new64 へ interp_field (9 量)。
# 2500 step ごとの出力から、登録の 8 量の系列 (cg_gref_judge.py series) と縮流部の残差 (cg_region_res.py) を取り出してから場を消す (res_0 と最後は残す)。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_grefc.sh <A の run> <B の run>
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$C45/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; PREP64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
GREF=run_0608_cg_gref_cont_k1; SEED=$GREF/res_60000.h5; A=$1; B=$2; MESH_B=$C45/_band_ab/core_grid/prep_Grefc/new64/nozzle.h5
LOG=$C45/cg_grefc.log; MIN_START_GB=10; MIN_RUN_GB=3; STEPS=60000; OUT=2500
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $A $B F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16) 格子 B $(sha256sum $MESH_B | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
for r in $A $B; do [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }; done
python3 -c "import h5py,sys; assert 'PLANES/rSurfVect' in h5py.File(sys.argv[1],'r')" $MESH_B || { echo "格子 B に rSurfVect が無い (中止)" >> $LOG; exit 1; }
# A: 対照
mkdir $A; for f in $GREF/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt) ;; *) cp -a $f $A/;; esac; done
python3 $TL/restart_field.py $SEED $A/nozzle.h5 --dst-run $A --forge $F7_64 > $A/restart.log 2>&1 || { echo "$A: restart_field に失敗" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $A/restart.log || { echo "$A: 移した量が想定と違う" >> $LOG; exit 1; }
# B: Grefc
( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $B --steps $STEPS --out $OUT --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$B: 準備に失敗" >> $LOG; exit 1; }
rm -f $B/nozzle.msh $B/nozzle.h5; cp $MESH_B $B/nozzle.h5
python3 $TL/interp_field.py $SEED $B/nozzle.h5 --dst-run $B --forge $F7_64 > $B/interp_newmesh.log 2>&1 || { echo "$B: 場の移しに失敗" >> $LOG; exit 1; }
python3 - $B/interp_newmesh.log <<'PY' || { echo "$B: 移した量が想定と違う" >> $LOG; exit 1; }
import re, sys
t = open(sys.argv[1]).read(); m = re.search(r"moved \[(.*?)\]", t)
moved = sorted(x.strip().strip("'").replace("(new)", "") for x in m.group(1).split(",")) if m else []
if moved != sorted(["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]): raise SystemExit(f"moved {moved}")
PY
python3 - $B/prepare_info.json $C45/_band_ab/core_grid/prep_Grefc/prepare_info.json <<'PY' || { echo "$B: prepare_info の差し替えに失敗" >> $LOG; exit 1; }
import json, sys
p, q = sys.argv[1:3]; a = json.load(open(p)); b = json.load(open(q))
if abs(a["scale_m"] - b["scale_m"]) > 0: raise SystemExit("scale_m が違う")
a["mesh"] = b["mesh"]; a["core_grid"] = {"grid": "Grefc", "mesh_from": "_band_ab/core_grid/prep_Grefc/new64/nozzle.h5", "ic": "interp_field run_0608/res_60000"}
json.dump(a, open(p, "w"), ensure_ascii=False, indent=1)
PY
grep -c '^mesh: {' $B/solverConfig.yaml | grep -qx 1 || { echo "$B: mesh の形が想定外" >> $LOG; exit 1; }
sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: 1, /" $B/solverConfig.yaml
for r in $A $B; do sed -i "s/output: {level: 1, extraFields: \[res_ro\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml; done
diff $A/solverConfig.yaml $B/solverConfig.yaml > $C45/cg_grefc_config.diff; diff $A/bcondConfig.yaml $B/bcondConfig.yaml >> $C45/cg_grefc_config.diff
[ -s $C45/cg_grefc_config.diff ] && { echo "A と B の設定が違う (中止): $(head -5 $C45/cg_grefc_config.diff | tr '\n' ' ')" >> $LOG; exit 1; }
grep -q "res_roe, res_roK, res_roOmega, volume" $A/solverConfig.yaml && grep -q "nStepOuter: $STEPS}" $A/solverConfig.yaml || { echo "設定の書き換えが想定と違う (中止)" >> $LOG; exit 1; }
echo "A と B の solverConfig・bcondConfig は同一 ($(grep -o 'lineDtDirectionalCap: [0-9.]*' $A/solverConfig.yaml)、$(grep -o 'axisSegmentRWeight: [01]' $A/solverConfig.yaml))" >> $LOG
PIDS=()
for r in $A $B; do
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  PIDS+=($!); echo "$r 起動 $(date -Is)" >> $LOG
done
sleep 150
for r in $A $B; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-110)" >> $LOG; done
while true; do
  done_all=1
  for r in $A $B; do
    [ -f $r/RUN_RC ] || done_all=0
    for f in $r/res_*.h5; do
      n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue; [ $n -gt 0 ] || continue
      [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && continue
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then
        python3 cg_gref_judge.py series "${r}_$n" $r --from $n --to $n > /dev/null 2>&1 && python3 cg_region_res.py $r $n && echo "$r: $n を取り出した" >> $LOG
        [ $n -ne $STEPS ] && [ -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && rm -f $f $r/res_$n.xmf
      fi
    done
  done
  [ $done_all = 1 ] && { left=0; for r in $A $B; do for f in $r/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ $n -gt 0 ] && [ ! -f _band_ab/core_grid/gref_series_${r}_$n.csv ] && left=1; done; done; [ $left = 0 ] && break; }
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in $A $B; do [ "$c" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done; done
  fi
  sleep 30
done
wait "${PIDS[@]}"
for r in $A $B; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
