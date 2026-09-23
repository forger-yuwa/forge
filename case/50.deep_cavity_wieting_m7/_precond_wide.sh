#!/bin/bash
# 検証が成立した幅 (w/d=0.383) でも低マッハ前処理に感度があるかを見る
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
RUN="../../solver_density_cuda/tools/run_case.sh"
r=run_0015_T1_wd0383_precond2
rm -rf $r; mkdir -p $r
cp run_0009_T1_wd0383_long/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $r/
echo "run_0009_T1_wd0383_long/res_60000.h5 (interp_field)。lowMachPrecond 0 → 2" > $r/CONTINUED_FROM
python3 - $r <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=s.replace('blockDPLUR: 1','blockDPLUR: 1, lowMachPrecond: 2')
s=re.sub(r'cfl: [\d.]+, cfl_pseudo: [\d.]+','cfl: 0.3, cfl_pseudo: 0.3',s)
s=re.sub(r'nStepOuter: \d+','nStepOuter: 3000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 3000',s)
p.write_text(s)
PY
last=$(ls run_0009_T1_wd0383_long/res_[0-9]*.h5 | sort -t_ -k2 -n | tail -1)
python3 ../../solver_density_cuda/tools/interp_field.py "$last" $r/mesh.h5
$RUN $r >/dev/null 2>&1 || true
if ls $r/res_nan_*.h5 >/dev/null 2>&1; then echo "cfl 0.3 でも発散"; exit 0; fi
python3 - $r <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'cfl: [\d.]+, cfl_pseudo: [\d.]+','cfl: 1.0, cfl_pseudo: 1.0',s)
s=re.sub(r'nStepOuter: \d+','nStepOuter: 30000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
python3 ../../solver_density_cuda/tools/interp_field.py $r/res_3000.h5 $r/mesh.h5
rm -f $r/res_[0-9]*.h5
$RUN $r >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $r --series | tail -4
echo "=== 0.383 precond A/B done ==="
