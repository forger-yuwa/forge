#!/bin/bash
# plan tooling-nozzle-core-grid §4.13 の起動の A/B (事前登録、codex 諮問 2026-10-10): G1x・キー 1・同じ初期 HDF5・同じバイナリで、
# time.deltaT.cfl_pseudo だけを A = 4 (本段)、B = 0.5 にして各 2000 step を順に回す。判定は cg_startab_judge.py。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_startab.sh
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
PREP64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
G=G1x; MESH=$C45/_band_ab/core_grid/prep_$G/new64/nozzle.h5
A=run_0601_cg_start_g1x_cfl4; B=run_0602_cg_start_g1x_cfl05
LOG=$C45/cg_startab.log
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1G . | tail -1 | tr -dc 0-9)" -ge 8 ] || { echo "ディスクの空きが 8 GB 未満 (中止)" >> $LOG; exit 1; }
for r in $A $B; do [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }; done
( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $A --steps 2000 --out 1000 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro > /dev/null ) || { echo "$A: 準備に失敗" >> $LOG; exit 1; }
rm -f $A/nozzle.msh $A/nozzle.h5; cp $MESH $A/nozzle.h5
python3 $TL/interp_field.py $SRC_RES $A/nozzle.h5 --dst-run $A --forge $F7_64 > $A/interp_newmesh.log 2>&1 || { echo "$A: 場の移しに失敗" >> $LOG; exit 1; }
python3 - $A/interp_newmesh.log <<'EOF' || { echo "$A: 移した量が想定と違う: $(grep -o 'moved.*' $A/interp_newmesh.log)" >> $LOG; exit 1; }
import re, sys
t = open(sys.argv[1]).read(); m = re.search(r"moved \[(.*?)\]", t)
moved = sorted(x.strip().strip("'").replace("(new)", "") for x in m.group(1).split(",")) if m else []
want = sorted(["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"])
if moved != want: raise SystemExit(f"moved {moved} != {want}")
EOF
python3 - $A/prepare_info.json $C45/_band_ab/core_grid/prep_$G/prepare_info.json $G <<'EOF' || { echo "$A: prepare_info の差し替えに失敗" >> $LOG; exit 1; }
import json, sys
p, q, g = sys.argv[1:4]; a = json.load(open(p)); b = json.load(open(q))
if a["scale_m"] != b["scale_m"]: raise SystemExit("scale_m が違う")
a["mesh"] = b["mesh"]; a["core_grid"] = {"grid": g, "mesh_from": f"_band_ab/core_grid/prep_{g}/new64/nozzle.h5", "ic": "interp_field run_0183/res_100000"}
json.dump(a, open(p, "w"), ensure_ascii=False, indent=1)
EOF
grep -c '^mesh: {' $A/solverConfig.yaml | grep -qx 1 || { echo "$A: mesh の形が想定外" >> $LOG; exit 1; }
sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: 1, /" $A/solverConfig.yaml
# B = A の複製で cfl_pseudo だけ 0.5 (同じ初期 HDF5)
mkdir $B; for f in $A/*; do case $(basename $f) in *.log|RUN_*|forge_*|residual_*) ;; *) cp -a $f $B/;; esac; done
python3 - $A $B <<'EOF' || { echo "$B: cfl_pseudo の書き換えに失敗" >> $LOG; exit 1; }
import sys, hashlib
from pathlib import Path
sys.path.insert(0, ".")
import ns_n012 as NS
a, b = Path(sys.argv[1]), Path(sys.argv[2])
ys = NS.yaml_strict(); t = (a / "solverConfig.yaml").read_text()
nt = ys.replace_scalars(t, {NS.CFLP: "0.5"})
d = sorted(NS.MK.diff_paths(ys.load(t), ys.load(nt)))
if d != [NS.CFLP]: raise SystemExit(f"差が cfl_pseudo だけでない: {d}")
(b / "solverConfig.yaml").write_text(nt)
h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
if h(a / "nozzle.h5") != h(b / "nozzle.h5"): raise SystemExit("初期 HDF5 が違う")
print("B: cfl_pseudo 0.5、初期 HDF5 は A と同じ", h(a / "nozzle.h5")[:16])
EOF
grep -n 'cfl\|^mesh' $A/solverConfig.yaml $B/solverConfig.yaml | sed 's#^#  #' >> $LOG
for r in $A $B; do
  echo "$r 起動 $(date -Is)" >> $LOG
  ( export FORGE_BIN=$F7_64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-140)" >> $LOG
done
python3 cg_startab_judge.py $A $B >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch cg_startab.done
