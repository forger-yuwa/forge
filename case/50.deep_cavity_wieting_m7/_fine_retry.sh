#!/bin/bash
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
RUN="../../solver_density_cuda/tools/run_case.sh"
while pgrep -f "forge$" >/dev/null || pgrep -f "_precond_ab.sh" >/dev/null; do sleep 30; done
echo "=== 格子細分 再投入 (soft 0.1 / mid 0.3) ==="
rm -rf run_0014_T1_wd0063_fine
python3 tools/make_case.py --run run_0014_T1_wd0063_fine --mesh t1_wd0063_fine --series A1 \
  --soft-cfl 0.1 --mid-cfl 0.3 --soft-steps 4000 --mid-steps 4000 \
  --ramp "0.2,0.5,1.0" --ramp-steps 2000 --cfl 2.0 --main-steps 30000 --out-int 5000
long=run_0014_T1_wd0063_fine_long
mkdir -p $long
cp run_0014_T1_wd0063_fine/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
echo "run_0014_T1_wd0063_fine" > $long/CONTINUED_FROM
python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
$RUN $long >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $long --series | tail -4
echo "=== 格子細分 done ==="
