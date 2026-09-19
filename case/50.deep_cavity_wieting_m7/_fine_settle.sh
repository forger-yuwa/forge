#!/bin/bash
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
r=run_0016_T1_wd0063_fine_settle
rm -rf $r; mkdir -p $r
cp run_0014_T1_wd0063_fine_long/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $r/
last=$(ls run_0014_T1_wd0063_fine_long/res_[0-9]*.h5 | sort -t_ -k2 -n | tail -1)
echo "run_0014_T1_wd0063_fine_long/$(basename $last) (DRIFTING のため 150k step 継続)" > $r/CONTINUED_FROM
python3 - $r <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 150000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 10000',s)
p.write_text(s)
PY
python3 ../../solver_density_cuda/tools/interp_field.py "$last" $r/mesh.h5
../../solver_density_cuda/tools/run_case.sh $r >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $r --series | tail -4
python3 ../../solver_density_cuda/tools/check_quasisteady.py --series-csv $r/cavity_series.csv --series-cols Qc_per_span_W_m | tail -3
echo "=== fine settle done ==="
