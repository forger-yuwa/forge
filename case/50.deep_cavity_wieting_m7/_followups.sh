#!/bin/bash
# 残りの数値側の切り分け: (a) 低マッハ前処理 A/B, (b) w/d=0.524, (c) すきま格子細分
set -e
cd "$(dirname "$0")"
RUN="../../solver_density_cuda/tools/run_case.sh"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial

echo "=== (a) 低マッハ前処理 A/B (w/d=0.063, run_0006 の場から) ==="
rm -rf run_0010_T1_wd0063_precond2
mkdir -p run_0010_T1_wd0063_precond2
cp run_0006_T1_wd0063_long/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} run_0010_T1_wd0063_precond2/
echo "run_0006_T1_wd0063_long (lowMachPrecond 0 → 2 の A/B)" > run_0010_T1_wd0063_precond2/CONTINUED_FROM
python3 - run_0010_T1_wd0063_precond2 <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=s.replace('blockDPLUR: 1','blockDPLUR: 1, lowMachPrecond: 2')
s=re.sub(r'nStepOuter: \d+','nStepOuter: 20000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
$RUN run_0010_T1_wd0063_precond2 >/dev/null 2>&1 || true
python3 tools/cavity_eval.py run_0010_T1_wd0063_precond2 --series | tail -4
echo "=== (a) done ==="

echo "=== (b) w/d = 0.524 ==="
python3 gen_mesh.py --case T1 --wd 0.524 --tag t1_wd0524_v2 --nx-gap 101 >/dev/null 2>&1
python3 tools/make_case.py --run run_0011_T1_wd0524 --mesh t1_wd0524_v2 --series D1 \
  --soft-cfl 0.15 --mid-cfl 0.4 --soft-steps 3000 --mid-steps 3000 \
  --ramp "0.3,0.6,1.2" --ramp-steps 1500 --cfl 2.0 --main-steps 20000 --out-int 5000
long=run_0011_T1_wd0524_long
mkdir -p $long
cp run_0011_T1_wd0524/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
echo "run_0011_T1_wd0524" > $long/CONTINUED_FROM
python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
$RUN $long >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $long --series | tail -4
echo "=== (b) done ==="

echo "=== (c) すきま格子細分 (w/d=0.063, 開口 60 → 120 セル) ==="
python3 gen_mesh.py --case T1 --wd 0.063 --tag t1_wd0063_fine --nx-gap 121 --ny-cav 301 >/dev/null 2>&1
python3 tools/make_case.py --run run_0012_T1_wd0063_fine --mesh t1_wd0063_fine --series A1 \
  --soft-cfl 0.3 --mid-cfl 0.7 --soft-steps 2000 --mid-steps 2000 \
  --ramp "0.5,1,2" --ramp-steps 1500 --cfl 2.0 --main-steps 30000 --out-int 5000
long2=run_0012_T1_wd0063_fine_long
mkdir -p $long2
cp run_0012_T1_wd0063_fine/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long2/
echo "run_0012_T1_wd0063_fine" > $long2/CONTINUED_FROM
python3 - $long2 <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
$RUN $long2 >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $long2 --series | tail -4
echo "=== (c) done ==="
