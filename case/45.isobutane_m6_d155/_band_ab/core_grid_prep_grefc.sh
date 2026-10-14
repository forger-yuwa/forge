#!/bin/bash
# plan tooling-nozzle-core-grid §4.36: Grefc の準備と検査 → 新しい変換器 (fgeom7-fp64) で変換し直して rSurfVect と閉性を確かめる (実行中に編集しない)
set -u
cd "$(dirname "$0")/.."
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
export FORGE_CONVERTER=$(pwd)/conv_tolerant.sh
export FORGE_CUDA_BLOCKSIZE=128
unset FORGE_ALLOW_UNVERIFIED_SPECIES
G=_band_ab/core_grid
python3 core_grid_mesh.py prep Grefc > $G/prep_Grefc.log 2>&1; echo "rc=$? prep Grefc $(date -Is)"
python3 core_grid_mesh.py prep-ic Grefc > $G/prep_ic_Grefc.log 2>&1; echo "rc=$? prep-ic Grefc $(date -Is)"
V7=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge; echo "V7 sha256 $(sha256sum $V7 | cut -c1-16)"
d=$G/prep_Grefc/new64; rm -rf $d; mkdir -p $d; cp $G/prep_Grefc/*.yaml $d/
( cd $d && $V7 ../nozzle.msh nozzle.h5 > conv.log 2>&1 ); rc=$?
if [ $rc -ne 0 ] && grep -q "Write Input HDF5" $d/conv.log && grep -q "GPUassert: invalid argument" $d/conv.log; then rc=0; fi
echo "rc=$rc new64 Grefc $(date -Is)"
python3 hoop_verify_conv.py check $d/nozzle.h5 --rw > $G/prep_Grefc_new64_check.log 2>&1; echo "rc=$? check"
python3 hoop_verify_conv.py diff $G/prep_Grefc/nozzle.h5 $d/nozzle.h5 > $G/prep_Grefc_new64_diff.log 2>&1; echo "rc=$? diff"
python3 - $G/prep_Grefc/nozzle.h5 $d/nozzle.h5 <<'PY'
import h5py, numpy as np, sys
a = h5py.File(sys.argv[1], "r"); b = h5py.File(sys.argv[2], "r")
print("COORD 一致:", np.array_equal(np.array(a["MESH/COORD"]), np.array(b["MESH/COORD"])), b["MESH/COORD"].dtype,
      "| wall_dist 一致:", np.array_equal(np.array(a["VALUE/wall_dist"]), np.array(b["VALUE/wall_dist"])),
      "| rSurfVect:", "PLANES/rSurfVect" in b, b["PLANES/rSurfVect"].shape if "PLANES/rSurfVect" in b else None)
PY
sha256sum $d/nozzle.h5 | cut -c1-16
df -h ~ | tail -1
# 区間の外の格子点が Gref (prep_Gref/new64) とビット単位で同じか (§4.36 の投入前の条件)
python3 - $G/prep_Gref/new64/nozzle.h5 $G/prep_Grefc/new64/nozzle.h5 $G/yaml-grefc.json <<'PY'
import h5py, numpy as np, sys, json
a = np.asarray(h5py.File(sys.argv[1], "r")["MESH/COORD"][:]).reshape(4719, 160, 3)
info = json.load(open(sys.argv[3]))["Grefc"]; ni = info["ni"]; nref, nnew = info["seg_cols_ref"], info["seg_cols_new"]
b = np.asarray(h5py.File(sys.argv[2], "r")["MESH/COORD"][:]).reshape(ni, 160, 3)
# 区間の始めの列 ia は、b と a が一致し続ける最後の列
same = [np.array_equal(a[i], b[i]) for i in range(min(len(a), len(b)))]
ia = same.index(False) - 1 if False in same else len(same)
tail_same = np.array_equal(a[ia + nref + 1:], b[ia + nnew + 1:])
print(f"区間の前の列 0〜{ia} が Gref とビット一致: {all(same[:ia+1])}、区間の後の列 {len(a)-ia-nref-1} 本が一致: {tail_same}、x の区間 {b[ia,-1,0]:.6g}〜{b[ia+nnew,-1,0]:.6g} m")
PY
echo DONE
