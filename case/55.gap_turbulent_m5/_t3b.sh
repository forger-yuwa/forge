#!/bin/bash
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
python3 tools/make_case.py --run run_0007_T1_cav_gapref --mesh t1_m5c --series K1 --plate-thermal adiabatic \
  --ic-delta 4e-4 --soft-cfl 0.2 --mid-cfl 0.5 --soft-steps 3000 --mid-steps 3000 \
  --ramp 0.3,0.6,1.0 --ramp-steps 2500 --cfl 1.2 --main-steps 30000 --out-int 5000
long=run_0007_T1_cav_gapref_main; mkdir -p $long
cp run_0007_T1_cav_gapref/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
echo "run_0007_T1_cav_gapref" > $long/CONTINUED_FROM
python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 120000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 10000',s)
p.write_text(s)
PY
../../solver_density_cuda/tools/run_case.sh $long >/dev/null 2>&1 || true
python3 tools/y1plus.py $long
python3 tools/cavity_eval.py $long --series --qfp-from run_0005_T0_iso500_main | tail -8
echo "=== gap refine done ==="
