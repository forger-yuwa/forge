#!/bin/bash
# plan tooling-rerun-conditions §6 (ii″)・(iv″) (事前登録 2026-10-06)
cd "$(dirname "$0")"
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
RC=../../solver_density_cuda/tools/rerun_conditions.py; QS=../../solver_density_cuda/tools/check_quasisteady.py
NS='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="full"); print("forge exit", rc); sys.exit(rc)'
EU='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged; from pathlib import Path; rc=run_staged(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)'
FW='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_forge; from pathlib import Path; rc=run_forge(Path(sys.argv[1])); print("forge exit", rc); sys.exit(rc)'
# (iv″) の Euler 参照
python3 $RC run_0086_euler_wallfit_pincal_r1_ext6k run_0134_rerun_euler_tt1500 --Tt 1500 --Y H2O=0.10 --keep-Ps --steps 6000 --out-interval 500 --cfl 2.0 | tail -1
python3 -c "$EU" run_0134_rerun_euler_tt1500 || echo "FAILED run_0134"
python3 rerun_series.py run_0134_rerun_euler_tt1500 > run_0134_rerun_euler_tt1500/series.log 2>&1 || true
# (iv″) ブロック延長 (cfl 5・6000・500 ごと、最大 4)
prev=run_0128_rerun_fullpath_ext
for k in 1 2 3 4; do
  R=run_013$((4+k))_rerun_fullpath_blk$k
  mkdir $R
  for f in nozzle.h5 nozzle.xmf bcondConfig.yaml solverConfig.yaml species_meta.yaml probe.yaml prepare_info.json MESH_QUALITY.txt wall_design.csv wall_physical.csv delta_r_initial.csv delta_r_initial.json target_axis_M.csv RERUN_CONDITIONS.json; do [ -e $prev/$f ] && cp $prev/$f $R/; done
  cp $prev/resolved_species_*.yaml $R/
  L=$(ls $prev/res_[0-9]*.h5 | sort -V | tail -1); python3 ../../solver_density_cuda/tools/restart_field.py $L $R/nozzle.h5 | tail -1
  python3 -c "$FW" $R > $R.log 2>&1 || { echo "FAILED $R"; break; }
  python3 rerun_series.py $R run_0134_rerun_euler_tt1500 > $R/series.log 2>&1
  grep -v -i warn $R/series.log | tail -6
  n=$(grep -cE "^  (delta_E|exitM_A|mdot) .* STEADY" $R/series.log); echo "blk$k STEADY count $n"
  prev=$R
  [ "$n" -ge 3 ] && break
done
echo IV2DONE
# (ii″) B3
python3 $RC run_0094_ns_c2pin_pass2_ext6k run_0133_rerun_pt08_noscale_full_cfl1 --Pt 4.4e6 --Ps 1789.6 --scale-ic none --steps 60000 --out-interval 5000 --cfl 1.0 | tail -1
python3 -c "$NS" run_0133_rerun_pt08_noscale_full_cfl1 > run_0133_rerun_pt08_noscale_full_cfl1.stdout.log 2>&1 || echo "FAILED run_0133"
python3 rerun_series.py run_0133_rerun_pt08_noscale_full_cfl1 run_0120_rerun_euler_pt08 > run_0133_rerun_pt08_noscale_full_cfl1/series.log 2>&1 || true
grep -v -i warn run_0133_rerun_pt08_noscale_full_cfl1/series.log | tail -6
echo FOLLOWUPDONE
