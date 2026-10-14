#!/bin/bash
# plan time_integration-line-viscous-jacobian §6 U3 (Poiseuille、2026-10-09 に組み直し): 側面を周期、両壁 300 K で止めたまま、体積力 f = 8μU/H² (中央 U = 300 m/s) で押す。
# ライン + 方向別、lineViscCoupling 0/2 × cfl 5/50、5000 step・25 ごと。run_0020〜0023。forge は run_case.sh 経由。
cd "$(dirname "$0")"
export FORGE_BIN=$(cat forge_target_linevisc.txt) FORGE_CUDA_BLOCKSIZE=128
CONV=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
n=20
for cfl in 5 50; do for lvc in 0 2; do
  r=$(printf "run_%04d_u3p_lvc%d_cfl%d" $n $lvc $cfl); n=$((n+1))
  python3 u1_thermjac.py prep $r --key 0 --cfl $cfl --steps 5000 --out 25 --conv $CONV --line 1 --dir 1 --lvc $lvc --poiseuille 300 > $r.prep.log 2>&1 \
    && bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1
  echo "$r $(grep -o 'forge exit=[0-9]*' $r/run_case_stdout.log) $(grep -m1 "'bodyForce'" $r/forge_run.log)" >> u3p.log
done; done
python3 u1_thermjac.py eval run_0020_u3p_lvc0_cfl5 run_0021_u3p_lvc2_cfl5 run_0022_u3p_lvc0_cfl50 run_0023_u3p_lvc2_cfl50 >> u3p.log 2>&1
touch u3p.done
