#!/bin/bash
# plan tooling-rerun-conditions §6 検証 run (事前登録 2026-10-06)。粗格子 run_0094 系列でツールのゲート。
#  (i)  run_0119_rerun_ctrl         : 無変更、cfl 5・12000 step・1000 ごと、stages none
#  (ii) run_0120_rerun_euler_pt08    : Euler 対参照 (run_0086 → Pt 4.4e6・Ps 1789.6・scale-ic pt、cfl 2・6000 step・500 ごと [準定常に 10 枚以上要るため]、stages none)
#       run_0121_rerun_pt08_scale    : NS 腕 A (scale-ic pt)、cfl 5・12000・1000 ごと、stages none (実験として明示)
#       run_0122_rerun_pt08_noscale  : NS 腕 B (scale-ic none)、同上
#  (iv) run_0123_rerun_fullpath      : --Tt 1500 --Y H2O=0.10 --keep-Ps、run_staged_ns(full) (soft 3000・mid 3000・本段 cfl 5・6000 step・500 ごと)
set -e
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
RC=../../solver_density_cuda/tools/rerun_conditions.py
NS="import sys; sys.path.insert(0,'../../design'); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages=sys.argv[2]); print('forge exit', rc); sys.exit(rc)"
EU="import sys; sys.path.insert(0,'../../design'); from forge_design.evaluate.runner_axismach import run_staged; from pathlib import Path; rc=run_staged(Path(sys.argv[1]), stages='none'); print('forge exit', rc); sys.exit(rc)"
REF=run_0094_ns_c2pin_pass2_ext6k; EREF=run_0086_euler_wallfit_pincal_r1_ext6k
[ -d run_0119_rerun_ctrl ] || python3 $RC $REF run_0119_rerun_ctrl --steps 12000 --out-interval 1000
python3 $RC $EREF run_0120_rerun_euler_pt08 --Pt 4.4e6 --Ps 1789.6 --scale-ic pt --steps 6000 --out-interval 500 --cfl 2.0
python3 $RC $REF run_0121_rerun_pt08_scale --Pt 4.4e6 --Ps 1789.6 --scale-ic pt --steps 12000 --out-interval 1000
python3 $RC $REF run_0122_rerun_pt08_noscale --Pt 4.4e6 --Ps 1789.6 --scale-ic none --steps 12000 --out-interval 1000
python3 $RC $REF run_0123_rerun_fullpath --Tt 1500 --Y H2O=0.10 --keep-Ps --steps 6000 --out-interval 500
for R in run_0119_rerun_ctrl run_0121_rerun_pt08_scale run_0122_rerun_pt08_noscale; do python3 -c "$NS" $R none || echo "FAILED $R"; done
python3 -c "$EU" run_0120_rerun_euler_pt08 || echo "FAILED run_0120"
python3 -c "$NS" run_0123_rerun_fullpath full || echo "FAILED run_0123"
python3 rerun_series.py run_0119_rerun_ctrl $EREF > run_0119_rerun_ctrl/series.log 2>&1
python3 rerun_series.py run_0120_rerun_euler_pt08 > run_0120_rerun_euler_pt08/series.log 2>&1
python3 rerun_series.py run_0121_rerun_pt08_scale run_0120_rerun_euler_pt08 > run_0121_rerun_pt08_scale/series.log 2>&1
python3 rerun_series.py run_0122_rerun_pt08_noscale run_0120_rerun_euler_pt08 > run_0122_rerun_pt08_noscale/series.log 2>&1
python3 rerun_series.py run_0123_rerun_fullpath > run_0123_rerun_fullpath/series.log 2>&1
python3 rerun_series.py $EREF > /tmp/ref0086_series.log 2>&1 || true
echo ALLDONE
