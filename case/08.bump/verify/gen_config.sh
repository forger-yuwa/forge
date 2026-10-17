#!/usr/bin/env bash
# 08.bump 検証用 config 生成。低M/高M × 陽解法/陰解法 の solverConfig.yaml を新しい run_* に書き出す。
# mesh(.h5)・bcond は既存の正準 run ディレクトリから複製する。
# usage: gen_config.sh <runDir> <loM|hiM> <exp|imp> <cfl> <nStepOuter>
set -euo pipefail
BUMP="/home/sano/work/forge/case/08.bump"
runDir="$1"; CASE="$2"; METHOD="$3"; CFL="$4"; NOUT="$5"

if [[ "$CASE" == "loM" ]]; then
  MESH="bump.h5"; VISC="0.0025"; THERM="0.0241"; INIT="bump"; SRC="$BUMP/run_slau"
else
  MESH="bump_4pct.h5"; VISC="0.0"; THERM="0.0"; INIT="bump_M1.65"; SRC="$BUMP/run_slau_mach1.65_4pctHt"
fi

if [[ "$METHOD" == "exp" ]]; then
  # 定常 (unsteady=0) の実効 CFL は cfl_pseudo (setDT_d.cu)。陽解法も cfl_pseudo で駆動。
  TI=3; NINNER=1
  CFLBLOCK="    cfl: ${CFL}
    cfl_pseudo: ${CFL}"
else
  TI=11; NINNER=20
  CFLBLOCK="    cfl: 1.0
    cfl_pseudo: ${CFL}
    implicitRelax: 1.0
    blockDPLUR: 1"
fi

mkdir -p "$runDir"
cp "$SRC/$MESH" "$runDir/"
cp "$SRC/bcondConfig.yaml" "$runDir/"
cp "$SRC/probe.yaml" "$runDir/" 2>/dev/null || true

cat > "$runDir/solverConfig.yaml" <<EOF
mesh:
  meshFormat: "hdf5"
  meshFileName: "${MESH}"
  valueFileName: "${MESH}"
gpu: 1
solver: "SLAU"
physProp:
  isCompressible: 1
  thermalMethod: 0
  viscMethod: 0
  ro : 1.4
  visc : ${VISC}
  thermCond: ${THERM}
  cp : 1005.0
  gamma : 1.4
time:
  unsteady: 0
  dualTime: 0
  last:
    control: 0
    nStepOuter: ${NOUT}
  deltaT:
    control: 1
    dt: 0.00001
${CFLBLOCK}
    dt_min: 0.00000001
    dt_max: 1.0
  outStepInterval: 2000
  outStepStart: 0
  timeIntegration: ${TI}
  nStepInner: ${NINNER}
space:
  convMethod: 2
  limiter: 2
turbulence:
  model: "none"
initial: "${INIT}"
EOF
