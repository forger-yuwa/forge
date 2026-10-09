#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.8 E2 (2026-10-09): 冷却平板 B (300 K、y₁ 3 µm) の収束場 (run_0005 の res_48000、単精度) を、
# FP64 の変換器で変換し直した同じ格子 (mesh/fp_y1_3um.msh) に restart_field.py で載せ、新バイナリ (FORGE_BIN) で cfl 2・2000 step・200 ごと:
#   run_0070_vj_point (ラインなし) / run_0071_vj_linedir_lvc0 (ライン + 方向別 + キー 5 + 値 0) / run_0072_vj_linedir_lvc2 (同 + 値 2)。forge は run_case.sh 経由。
# 置くもの: e2_src/ に run_0005 の solverConfig.yaml・bcondConfig.yaml・probe.yaml・res_48000.h5、mesh/fp_y1_3um.msh。AWS の case dir で実行する。
cd "$(dirname "$0")"
CONV=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
TOOLS=$HOME/forge-wallfit/solver_density_cuda/tools
export FORGE_CUDA_BLOCKSIZE=128 LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
unset FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER
prep() {  # $1 run, $2 deltaT に足すキー (空ならなし)
  local r=$1 keys=$2
  [ -e $r ] && { echo "$r が既にある — 止める"; return 1; }
  mkdir $r && cp mesh/fp_y1_3um.msh e2_src/bcondConfig.yaml e2_src/probe.yaml $r/
  python3 - "$r" "$keys" <<'PY'
import re, sys
r, keys = sys.argv[1], sys.argv[2]
s = open("e2_src/solverConfig.yaml").read()
s = re.sub(r"nStepOuter: \d+", "nStepOuter: 2000", s)
s = re.sub(r"outStepInterval: \d+", "outStepInterval: 200", s)
if keys:
    s = s.replace("deltaT: {", "deltaT: {" + keys + ", ", 1)
open(f"{r}/solverConfig.yaml", "w").write(s)
PY
  (cd $r && $CONV fp_y1_3um.msh mesh.h5 > convert.log 2>&1) || { echo "$r 変換に失敗"; return 1; }
  python3 $TOOLS/restart_field.py e2_src/res_48000.h5 $r/mesh.h5 --dst-run $r > $r/restart_field.log 2>&1
  tail -1 $r/restart_field.log
  rm -f $r/fp_y1_3um.msh
}
prep run_0070_vj_point "" && bash $RC $PWD/run_0070_vj_point > run_0070_vj_point/run_case_stdout.log 2>&1 &
prep run_0071_vj_linedir_lvc0 "lineImplicit: 1, lineDtDirectional: 1, implicitThermalJacobian: 5" && bash $RC $PWD/run_0071_vj_linedir_lvc0 > run_0071_vj_linedir_lvc0/run_case_stdout.log 2>&1 &
prep run_0072_vj_linedir_lvc2 "lineImplicit: 1, lineDtDirectional: 1, implicitThermalJacobian: 5, lineViscCoupling: 2" && bash $RC $PWD/run_0072_vj_linedir_lvc2 > run_0072_vj_linedir_lvc2/run_case_stdout.log 2>&1 &
wait
for r in run_0070_vj_point run_0071_vj_linedir_lvc0 run_0072_vj_linedir_lvc2; do echo "$r $(grep -o 'forge exit=[0-9]*' $r/run_case_stdout.log) $(grep -m1 'detectNaN.*Non' $r/forge_run.log | cut -c1-90)" >> e2_plate.log; done
touch e2_plate.done
