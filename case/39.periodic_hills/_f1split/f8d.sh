#!/usr/bin/env bash
# #8d (plan boundary-node-periodic-gradient-fix §5.1 #8d、codex result-2 M4): 1266aba1 に m3 (F1 初期値 1 を allocVariables へ、
# buildScalarDescs の初回充填を削除) だけを当てたバイナリで、run_0038 と同じ restart 入力 ic_S3_ext.h5 から 1 step (sstSigmaBlend 1 既定)。
#   m3_b1 : 1266aba1 + m3only.patch、sstSigmaBlend 1   (判定対象)
#   old_b1_rerun : 旧 1266aba1、sstSigmaBlend 1 (同じセッションでの対照。_f1split/old_b1 の再現確認、判定外)
set -uo pipefail
C=~/forge-pgrad-new/case/39.periodic_hills
T=$C/_f1split
declare -A BIN=([m3_b1]=$HOME/forge-pgrad-m3/solver_density_cuda/build/forge [old_b1_rerun]=$HOME/forge-pgrad-old/solver_density_cuda/build/forge)
for v in m3_b1 old_b1_rerun; do
  R=$T/$v
  [ -e $R/res_1.h5 ] && { echo "exists $R"; continue; }
  mkdir -p $R
  cp $C/run_0038_r1_gradfix_old_ext/{hill_r1_80x50x30.h5,ic_S3_ext.h5,bcondConfig.yaml} $R/
  printf "outStepInterval: 1\noutStepStart: 0\npoints:\nsurfaces:\n" > $R/probe.yaml
  python3 - "$C/run_0038_r1_gradfix_old_ext/solverConfig.yaml" "$R/solverConfig.yaml" 1 <<'PY'
import sys, yaml
src, dst, blend = sys.argv[1], sys.argv[2], int(sys.argv[3])
cfg = yaml.safe_load(open(src))
cfg["time"]["last"]["nStepOuter"] = 1
cfg["time"]["outStepInterval"] = 1
cfg["time"]["outStepStart"] = 0
cfg["turbulence"]["sstSigmaBlend"] = blend
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False)
PY
  FORGE_BIN=${BIN[$v]} FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128 ~/forge-pgrad-new/solver_density_cuda/tools/run_case.sh $R > $R/run_case_stdout.log 2>&1
  echo "$v rc=$? sigmaBlend=$(grep -o 'sigmaBlend=[0-9]' $R/forge_run.log | head -1) sha=$(grep forge_sha256 $R/RUN_PROVENANCE.txt)"
done
# 入力の確認 (_f1split/input_check.txt と同じ): ic_S3_ext.h5 の保存量が run_0036 res_800000 とビット一致
python3 - <<'PY'
import h5py, numpy as np
a = h5py.File("/home/ubuntu/forge-pgrad-new/case/39.periodic_hills/run_0036_r1_gradfix_old/res_800000.h5", "r")["VALUE"]
b = h5py.File("/home/ubuntu/forge-pgrad-new/case/39.periodic_hills/_f1split/m3_b1/ic_S3_ext.h5", "r")["VALUE"]
print("input:", ", ".join(k + (" bit-identical" if (np.asarray(a[k]).view(np.uint32) == np.asarray(b[k]).view(np.uint32)).all() else " DIFFERENT")
                          for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega")), "| sstF1 in ic:", "sstF1" in b)
PY
echo F8D_DONE
