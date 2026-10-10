#!/bin/bash
# plan tooling-nozzle-core-grid §5.1 #13 (b): Gref の準備と検査 → 新しい変換器 (fgeom7-fp64) で変換し直して rSurfVect と閉性を確かめる (実行中に編集しない)
set -u
cd "$(dirname "$0")/.."
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
export FORGE_CONVERTER=$(pwd)/conv_tolerant.sh
export FORGE_CUDA_BLOCKSIZE=128
unset FORGE_ALLOW_UNVERIFIED_SPECIES
G=_band_ab/core_grid
python3 core_grid_mesh.py prep Gref > $G/prep_Gref.log 2>&1; echo "rc=$? prep Gref $(date -Is)"
python3 core_grid_mesh.py prep-ic Gref > $G/prep_ic_Gref.log 2>&1; echo "rc=$? prep-ic Gref $(date -Is)"
V7=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge; echo "V7 sha256 $(sha256sum $V7 | cut -c1-16)"
d=$G/prep_Gref/new64; rm -rf $d; mkdir -p $d; cp $G/prep_Gref/*.yaml $d/
( cd $d && $V7 ../nozzle.msh nozzle.h5 > conv.log 2>&1 ); rc=$?
if [ $rc -ne 0 ] && grep -q "Write Input HDF5" $d/conv.log && grep -q "GPUassert: invalid argument" $d/conv.log; then rc=0; fi
echo "rc=$rc new64 Gref $(date -Is)"
python3 hoop_verify_conv.py check $d/nozzle.h5 --rw > $G/prep_Gref_new64_check.log 2>&1; echo "rc=$? check"
python3 hoop_verify_conv.py diff $G/prep_Gref/nozzle.h5 $d/nozzle.h5 > $G/prep_Gref_new64_diff.log 2>&1; echo "rc=$? diff"
python3 - $G/prep_Gref/nozzle.h5 $d/nozzle.h5 <<'PY'
import h5py, numpy as np, sys
a = h5py.File(sys.argv[1], "r"); b = h5py.File(sys.argv[2], "r")
print("COORD 一致:", np.array_equal(np.array(a["MESH/COORD"]), np.array(b["MESH/COORD"])), b["MESH/COORD"].dtype,
      "| wall_dist 一致:", np.array_equal(np.array(a["VALUE/wall_dist"]), np.array(b["VALUE/wall_dist"])),
      "| rSurfVect:", "PLANES/rSurfVect" in b, b["PLANES/rSurfVect"].shape if "PLANES/rSurfVect" in b else None)
PY
sha256sum $d/nozzle.h5 | cut -c1-16
df -h ~ | tail -1
echo DONE
