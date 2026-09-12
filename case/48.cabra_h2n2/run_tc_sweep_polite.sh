#!/usr/bin/env bash
# GPU が他の計算に使われていない (compute apps 無し) ときだけ次の T_c run を投入する版。30 分毎に確認。
set -u
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
RES=$1; QBIN=$2; N0=$3; NSTEP=$4; TCS=$5; PY=/home/sano/work/forge/.venv-chem/bin/python
n=$N0
for tc in $TCS; do
  while [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" -gt 0 ]; do
    echo "$(date "+%F %T") GPU busy (other job) -> wait 30 min"; sleep 1800
  done
  d=$(printf "run_%04d_cmc_tc%d" "$n" "$tc"); echo "$(date "+%F %T") == $d (T_c=$tc K)"
  rm -rf "$HERE/$d"
  (cd "$HERE" && $PY setup_cabra_case.py "$d" --chem 1 --mixfrac 1 --sdm 0 --cmc 1 --couple 7 --cmcchem 1 --cmcdt 0.5 --jac 2 --ji 5 --tci 0 \
      --cfl 0.5 --conv 1 --relax 0.5 --iccol 1 --nstep "$NSTEP" --out 500 --cmcfp32 1 --dualtime "${DT:-1e-5}" --tcof "$tc" --restart "$RES" --cmcq "$QBIN" | tail -1)
  cp "$HERE/mech.yaml" "$HERE/$d/mech.yaml"
  (cd "$ROOT" && bash solver_density_cuda/tools/run_case.sh "case/48.cabra_h2n2/$d" > "$HERE/$d.launch.log" 2>&1) || true
  $PY "$HERE/analyze_cmc.py" "$HERE/$d" 2>/dev/null | tail -3
  n=$((n+1))
done
echo "$(date "+%F %T") sweep done"
