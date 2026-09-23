#!/bin/bash
set -e
cd "$(dirname "$0")"
# 0.211 の本段継続 (内部緩和)
long=run_0008_T1_wd0211_long
mkdir -p $long
cp run_0008_T1_wd0211/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
echo "run_0008_T1_wd0211 (main 段の収束場)" > $long/CONTINUED_FROM
python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial ../../solver_density_cuda/tools/run_case.sh $long >/dev/null 2>&1 || true
echo "=== 0.211 long done ==="
# 0.383
python3 tools/make_case.py --run run_0009_T1_wd0383 --mesh t1_wd0383_v2 --series C1 \
  --soft-cfl 0.15 --mid-cfl 0.4 --soft-steps 3000 --mid-steps 3000 \
  --ramp "0.3,0.6,1.2" --ramp-steps 1500 --cfl 2.0 --main-steps 20000 --out-int 5000
long9=run_0009_T1_wd0383_long
mkdir -p $long9
cp run_0009_T1_wd0383/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long9/
echo "run_0009_T1_wd0383 (main 段の収束場)" > $long9/CONTINUED_FROM
python3 - $long9 <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial ../../solver_density_cuda/tools/run_case.sh $long9 >/dev/null 2>&1 || true
echo "=== 0.383 long done ==="
