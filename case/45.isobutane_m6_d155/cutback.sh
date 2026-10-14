#!/bin/bash
# plan time_integration-implicit-thermal-jacobian §6.0 の切り戻し (run_0262) と point 側の標本の追加 (run_0263)。
# 同じ新バイナリ (thermjac_cap_fp64、35e498b1) で point cfl 4 (キー 0・ライン 0、run_0217 と同じ設定)、40000 step・5000 ごと。forge は cold_cfl.py run 経由。
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target.txt) COLD_ALT_BINARY=thermjac_cap_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE
python3 cold_cfl.py prep run_0217_ns_coldmesh_tw300_cfl4_ext3 run_0262_ns_coldmesh_tw300_cutback_point --steps 40000 --cfl 4 --out 5000 \
    --field-from run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2 > /dev/null && rm -f run_0262_ns_coldmesh_tw300_cutback_point/nozzle.msh
python3 cold_cfl.py prep run_0217_ns_coldmesh_tw300_cfl4_ext3 run_0263_ns_coldmesh_tw300_cfl4_ext4 --steps 40000 --cfl 4 --out 5000 > /dev/null \
    && rm -f run_0263_ns_coldmesh_tw300_cfl4_ext4/nozzle.msh
for r in run_0262_ns_coldmesh_tw300_cutback_point run_0263_ns_coldmesh_tw300_cfl4_ext4; do
  (python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo "$r rc=$?" >> cutback.log) &
done
wait
touch cutback.done
