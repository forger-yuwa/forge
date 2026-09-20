#!/bin/bash
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
src=run_0007_T1_cav_gapref_main
r=run_0008_T1_cav_urans
rm -rf $r; mkdir -p $r
cp $src/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $r/
last=$(ls $src/res_[0-9]*.h5 | sort -t_ -k2 -n | tail -1)
echo "$src/$(basename $last) → URANS (dual-time BDF2)" > $r/CONTINUED_FROM
python3 - $r <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=s.replace('  unsteady: 0','  unsteady: 1').replace('  dualTime: 0','  dualTime: 1')
s=re.sub(r'control: 1, dt: [^,]+,','control: 0, dt: 1.0e-07,',s)
s=re.sub(r'nStepOuter: \d+','nStepOuter: 5000',s)
s=re.sub(r'outStepInterval: \d+','outStepInterval: 50',s)
s=s.replace('  nStepInner: 4','  nStepInner: 4\n  bdfOrder: 2\n  nSubIterDualTime: 20')
p.write_text(s)
PY
python3 ../../solver_density_cuda/tools/interp_field.py "$last" $r/mesh.h5
../../solver_density_cuda/tools/run_case.sh $r >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $r --series --qfp-from run_0005_T0_iso500_main | tail -10
echo "=== T3 URANS done ==="
