#!/bin/bash
# 低マッハ前処理 A/B のやり直し: run_0006 の 80k 収束場から interp_field で場を渡し、CFL を上げずに投入
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
RUN="../../solver_density_cuda/tools/run_case.sh"
while pgrep -f "forge$" >/dev/null; do sleep 30; done      # 先行 run の終了待ち
rm -rf run_0013_T1_wd0063_precond2b
mkdir -p run_0013_T1_wd0063_precond2b
cp run_0006_T1_wd0063_long/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} run_0013_T1_wd0063_precond2b/
echo "run_0006_T1_wd0063_long/res_80000.h5 (interp_field で場を渡す)。lowMachPrecond 0 → 2 の A/B" > run_0013_T1_wd0063_precond2b/CONTINUED_FROM
python3 - run_0013_T1_wd0063_precond2b <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=s.replace('blockDPLUR: 1','blockDPLUR: 1, lowMachPrecond: 2')
s=re.sub(r'cfl: [\d.]+, cfl_pseudo: [\d.]+','cfl: 0.3, cfl_pseudo: 0.3',s)
s=re.sub(r'nStepOuter: \d+','nStepOuter: 3000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 3000',s)
p.write_text(s)
PY
python3 ../../solver_density_cuda/tools/interp_field.py run_0006_T1_wd0063_long/res_80000.h5 run_0013_T1_wd0063_precond2b/mesh.h5
$RUN run_0013_T1_wd0063_precond2b >/dev/null 2>&1 || true
if ls run_0013_T1_wd0063_precond2b/res_nan_*.h5 >/dev/null 2>&1; then
  echo "precond2 は cfl 0.3 でも発散 → A/B 不成立 (前処理は本ケースで使えない)"
  exit 0
fi
# 通ったら本番: 30000 step
python3 - run_0013_T1_wd0063_precond2b <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'cfl: [\d.]+, cfl_pseudo: [\d.]+','cfl: 1.0, cfl_pseudo: 1.0',s)
s=re.sub(r'nStepOuter: \d+','nStepOuter: 30000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
python3 ../../solver_density_cuda/tools/interp_field.py run_0013_T1_wd0063_precond2b/res_3000.h5 run_0013_T1_wd0063_precond2b/mesh.h5
rm -f run_0013_T1_wd0063_precond2b/res_[0-9]*.h5
$RUN run_0013_T1_wd0063_precond2b >/dev/null 2>&1 || true
python3 tools/cavity_eval.py run_0013_T1_wd0063_precond2b --series | tail -4
echo "=== precond A/B done ==="
