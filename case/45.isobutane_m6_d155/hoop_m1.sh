#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.11 (m1、事前登録): 小さな試験の格子 3 つを新しい変換器 (~/forge-fgeom7-fp64、0b4dff4e と同じソース) で
# isAxisymmetric 1・node として変換し、hoop_verify_conv.py check --rw で照合する。判定は hoop_m1_an.py。インスタンス B で回す。
# 使い方: hoop_m1.sh <格子 (.msh) のディレクトリ> <作業ディレクトリ>
set -u
M=$1; W=$2; mkdir -p $W
CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge
SC=~/forge-wallfit/case/48.flat_plate_cooled_m4/_smoke/solverConfig.yaml   # CPG の物性 (case/48)
HERE=$(cd $(dirname $0) && pwd)
echo "== $(date -Is) 変換器 $(sha256sum $CONV | cut -c1-16)" > $W/hoop_m1.log
for T in T1_tri_curved T2_mixed T3_step_multicorner; do
  d=$W/$T; rm -rf $d; mkdir -p $d; cp $M/$T.msh $d/
  python3 - $SC $d/solverConfig.yaml <<'PY' || { echo "$T: solverConfig の書き換えに失敗" >> $W/hoop_m1.log; continue; }
import re, sys
s = open(sys.argv[1]).read()   # case/48 の mesh: は複数行にまたがる (1 行目だけ替えると YAML が壊れる)
t = re.sub(r"^mesh:\s*\{.*?\}", 'mesh: {discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, meshFileName: "m.h5", valueFileName: "m.h5"}', s, count=1, flags=re.S | re.M)
assert t != s, "mesh: が見つからない"
open(sys.argv[2], "w").write(t)
PY
  cat > $d/bcondConfig.yaml <<'EOF'
inlet:  {physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {Pt: 200000.0, Tt: 300.0, k: 1.0, omega: 1000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 100000.0, Pt: 100000.0, Tt: 300.0}}
wall:   {physID: 3, kind: wall_isothermal,  outputHDFflg: 0, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 300.0}}
axis:   {physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }
EOF
  [ $T = T3_step_multicorner ] && echo 'wall2:  {physID: 6, kind: wall,             outputHDFflg: 0, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}' >> $d/bcondConfig.yaml
  ( cd $d && $CONV $T.msh m.h5 > conv.log 2>&1 ); rc=$?
  if [ $rc -ne 0 ] && grep -q "Write Input HDF5" $d/conv.log && grep -q "GPUassert: invalid argument" $d/conv.log; then rc=0; fi   # 終了時の既知の GPUassert
  echo "== $T 変換 rc=$rc" >> $W/hoop_m1.log
  python3 $HERE/hoop_verify_conv.py check $d/m.h5 --rw >> $W/hoop_m1.log 2>&1; echo "  check rc=$?" >> $W/hoop_m1.log
done
python3 $HERE/hoop_m1_an.py $W/hoop_m1.log | tee $W/hoop_m1_verdict.txt
