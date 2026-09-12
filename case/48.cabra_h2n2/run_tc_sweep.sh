#!/usr/bin/env bash
# Cabra T_c 応答曲線スイープ (plan chemistry-cmc-tci §5.1 #6): 1045 K の収束場 (+ Q(η)) から coflow 温度だけ変えて順次回す。
#   usage: bash run_tc_sweep.sh <restart_res.h5 (run dir 相対)> <restart cmc_Q.bin (run dir 相対)> <first run number> [nstep] [Tc list]
#   dual-time (dt は環境変数 DT, 既定 1e-5 s; 定常反復は柱モードで ξ フラックスが保存されないため使わない)
#   例:    bash run_tc_sweep.sh ../run_0099_cmc_c5_fp32_cont/res_16000.h5 ../run_0099_cmc_c5_fp32_cont/cmc_Q.bin 100 12000 "1015 1030 1060 1075"
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
RES=$1; QBIN=$2; N0=$3; NSTEP=${4:-12000}; TCS=${5:-"1015 1030 1060 1075"}
PY=/home/sano/work/forge/.venv-chem/bin/python
n=$N0
for tc in $TCS; do
  d=$(printf "run_%04d_cmc_tc%d" "$n" "$tc")
  echo "== $d (T_c=$tc K)"
  rm -rf "$HERE/$d"
  (cd "$HERE" && $PY setup_cabra_case.py "$d" --chem 1 --mixfrac 1 --sdm 0 --cmc 1 --couple 7 --cmcchem 1 --cmcdt 0.5 --jac 2 --ji 5 --tci 0 \
      --cfl 0.5 --conv 1 --relax 0.5 --iccol 1 --nstep "$NSTEP" --out 500 --cmcfp32 1 --dualtime "${DT:-1e-5}" --tcof "$tc" --restart "$RES" --cmcq "$QBIN" | tail -1)
  cp "$ROOT/solver_density_cuda/tools/mechanisms/h2co_li2004_cantera.yaml" "$HERE/$d/mech.yaml"
  (cd "$ROOT" && bash solver_density_cuda/tools/run_case.sh "case/48.cabra_h2n2/$d" > "$HERE/$d.launch.log" 2>&1) || true   # run_case.sh はリポジトリルート相対
  $PY "$HERE/analyze_cmc.py" "$HERE/$d" | tail -8
  n=$((n+1))
done
