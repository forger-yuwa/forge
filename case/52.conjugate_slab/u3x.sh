#!/bin/bash
# plan time_integration-line-viscous-jacobian §6 U3 (組み直し、2026-10-09): 側面を周期にして上壁を x 方向に 300 m/s で動かす Couette。
# ライン + 方向別、lineViscCoupling 0/2 × cfl 5/50、5000 step・25 ごと。run_0016〜0019。forge は run_case.sh 経由。
cd "$(dirname "$0")"
export FORGE_BIN=$(cat forge_target_linevisc.txt) FORGE_CUDA_BLOCKSIZE=128
CONV=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
n=16
for cfl in 5 50; do for lvc in 0 2; do
  r=$(printf "run_%04d_u3x_lvc%d_cfl%d" $n $lvc $cfl); n=$((n+1))
  python3 u1_thermjac.py prep $r --key 0 --cfl $cfl --steps 5000 --out 25 --conv $CONV --line 1 --dir 1 --lvc $lvc --couette 300 > $r.prep.log 2>&1 \
    && bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1
  echo "$r $(grep -o 'forge exit=[0-9]*' $r/run_case_stdout.log) $(grep -m1 'periodic' $r/forge_run.log | cut -c1-80)" >> u3x.log
done; done
python3 u1_thermjac.py eval run_0016_u3x_lvc0_cfl5 run_0017_u3x_lvc2_cfl5 run_0018_u3x_lvc0_cfl50 run_0019_u3x_lvc2_cfl50 >> u3x.log 2>&1
touch u3x.done
