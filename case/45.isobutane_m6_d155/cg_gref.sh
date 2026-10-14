#!/bin/bash
# plan tooling-nozzle-core-grid §4.18 (事前登録、codex 諮問 2026-10-11): 参照の格子 Gref (4719 × 160) を 1 本、キー 1・B0 の構成で 20 万 step の固定長で回す。
# 準備は cg_runs.sh の prep_arm と同じ (cold_cfl.py prep → 新しい変換の格子に差し替え → interp_field で 9 量 → prepare_info の mesh を差し替え → キー 1)。
# cg_runs.sh との違い: m9_watch を使わない (水準で止めない・出力を消さない)。160k 未満の場は、次の出力が出てから消す (残差・ログは残す)。
# 160k〜200k の 17 出力はすべて残す。ディスクの空きが 3 GB を切ったら自分の forge だけを止める。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_gref.sh <run>
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
G=Gref; r=$1; mesh=$C45/_band_ab/core_grid/prep_$G/new64/nozzle.h5
LOG=$C45/cg_gref.log; MIN_START_GB=20; MIN_RUN_GB=3; BUDGET=200000; OUT=2500; KEEP_FROM=160000
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) $r F7_64 $(sha256sum $F7_64 | cut -c1-16) mesh $(sha256sum $mesh | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge $MIN_START_GB ] || { echo "ディスクの空きが ${MIN_START_GB} GB 未満 (中止)" >> $LOG; exit 1; }
[ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
[ -f $mesh ] || { echo "$r: 格子 $mesh が無い (中止)" >> $LOG; exit 1; }
python3 -c "import h5py,sys; assert 'PLANES/rSurfVect' in h5py.File(sys.argv[1],'r')" $mesh || { echo "$r: 格子に rSurfVect が無い (中止)" >> $LOG; exit 1; }
( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps $BUDGET --out $OUT --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$r: 準備に失敗" >> $LOG; exit 1; }
rm -f $r/nozzle.msh $r/nozzle.h5; cp $mesh $r/nozzle.h5
python3 $TL/interp_field.py $SRC_RES $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/interp_newmesh.log 2>&1 || { echo "$r: 場の移しに失敗" >> $LOG; exit 1; }
python3 - $r/interp_newmesh.log <<'PY' || { echo "$r: 移した量が想定と違う: $(grep -o 'moved.*' $r/interp_newmesh.log)" >> $LOG; exit 1; }
import re, sys
t = open(sys.argv[1]).read(); m = re.search(r"moved \[(.*?)\]", t)
moved = sorted(x.strip().strip("'").replace("(new)", "") for x in m.group(1).split(",")) if m else []
want = sorted(["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"])
if moved != want: raise SystemExit(f"moved {moved} != {want}")
PY
python3 - $r/prepare_info.json $C45/_band_ab/core_grid/prep_$G/prepare_info.json $G <<'PY' || { echo "$r: prepare_info の差し替えに失敗" >> $LOG; exit 1; }
import json, sys
p, q, g = sys.argv[1:4]; a = json.load(open(p)); b = json.load(open(q))
if abs(a["scale_m"] - b["scale_m"]) > 0: raise SystemExit("scale_m が違う")
a["mesh"] = b["mesh"]; a["core_grid"] = {"grid": g, "mesh_from": f"_band_ab/core_grid/prep_{g}/new64/nozzle.h5", "ic": "interp_field run_0183/res_100000"}
json.dump(a, open(p, "w"), ensure_ascii=False, indent=1)
PY
grep -c '^mesh: {' $r/solverConfig.yaml | grep -qx 1 || { echo "$r: mesh の形が想定外" >> $LOG; exit 1; }
sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: 1, /" $r/solverConfig.yaml
grep -n '^mesh\|nStepOuter\|outStepInterval\|lineImplicit\|implicitRelax\|extraFields\|cfl_pseudo\|detectNaN' $r/solverConfig.yaml | sed "s#^#  $r: #" >> $LOG
( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
PID_RC=$!
echo "$r 起動 $(date -Is) $(grep -o 'axisSegmentRWeight: [01]' $r/solverConfig.yaml)" >> $LOG
sleep 150
echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-140)" >> $LOG
# 160k 未満の場の間引きとディスクの見張り (自分の run の forge だけを止める)
while [ ! -f $r/RUN_RC ]; do
  for f in $r/res_*.h5; do
    n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue
    [ $n -gt 0 ] && [ $n -lt $KEEP_FROM ] && [ -f $r/res_$((n + OUT)).h5 ] && rm -f $f $r/res_$n.xmf
  done
  if [ "$(free_gb)" -lt $MIN_RUN_GB ]; then
    for p in $(pgrep -x forge); do [ "$(readlink /proc/$p/cwd)" = "$C45/$r" ] && kill $p && echo "ディスク ${MIN_RUN_GB} GB 未満で停止: $r (pid $p) $(date -Is)" >> $LOG; done
    break
  fi
  sleep 60
done
wait $PID_RC
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の出力 $(ls $r/res_*.h5 | sed 's/.*res_//; s/.h5//' | sort -n | tail -1)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
