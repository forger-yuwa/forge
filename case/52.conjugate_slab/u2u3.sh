#!/bin/bash
# plan time_integration-line-viscous-jacobian §6 U2 (純伝導) と U3 (z 方向の Couette)。ライン + 方向別、lineViscCoupling 0/2 × cfl 5/50、5000 step・25 ごと。
# バイナリは forge_target_linevisc.txt (5ab83056 + typedef double)。forge は run_case.sh 経由。AWS の case dir で実行する。
cd "$(dirname "$0")"
export FORGE_BIN=$(cat forge_target_linevisc.txt) FORGE_CUDA_BLOCKSIZE=128
CONV=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
n=8
for kind in u2 u3; do for cfl in 5 50; do for lvc in 0 2; do
  r=$(printf "run_%04d_%s_lvc%d_cfl%d" $n $kind $lvc $cfl); n=$((n+1))
  extra="--line 1 --dir 1 --lvc $lvc"; [ $kind = u3 ] && extra="$extra --couette 300"
  python3 u1_thermjac.py prep $r --key 0 --cfl $cfl --steps 5000 --out 25 --conv $CONV $extra > /dev/null && bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1
  echo "$r $(grep -o 'forge exit=[0-9]*' $r/run_case_stdout.log) $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> u2u3.log
done; done; done
python3 u1_thermjac.py eval run_0008_u2_lvc0_cfl5 run_0009_u2_lvc2_cfl5 run_0010_u2_lvc0_cfl50 run_0011_u2_lvc2_cfl50 \
  run_0012_u3_lvc0_cfl5 run_0013_u3_lvc2_cfl5 run_0014_u3_lvc0_cfl50 run_0015_u3_lvc2_cfl50 >> u2u3.log 2>&1
touch u2u3.done
