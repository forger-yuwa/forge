#!/bin/bash
# T3 チェーン: ① 断熱平板 (T_aw) ② 等温 500K 平板 (分母 q_fp) ③ キャビティ
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
COMMON="--series K1 --ic-delta 4e-4 --soft-cfl 0.2 --mid-cfl 0.5 --soft-steps 3000 --mid-steps 3000 --ramp 0.3,0.6,1.0 --ramp-steps 2500 --cfl 1.2 --out-int 5000"
echo "=== ① 断熱平板 (T_aw と δ99) ==="
python3 tools/make_case.py --run run_0004_T0_ad --mesh t0_m5b --plate-thermal adiabatic --main-steps 25000 $COMMON
echo "=== ② 等温 500 K 平板 (分母 q_fp) ==="
python3 tools/make_case.py --run run_0005_T0_iso500 --mesh t0_m5b --plate-thermal isothermal --main-steps 25000 $COMMON
for r in run_0004_T0_ad run_0005_T0_iso500; do
  long=${r}_main; mkdir -p $long
  cp $r/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
  echo "$r" > $long/CONTINUED_FROM
  python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 15000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
  ../../solver_density_cuda/tools/run_case.sh $long >/dev/null 2>&1 || true
  python3 tools/plate_eval.py $long --x-ref 0.34875 | tail -4
done
echo "=== ③ キャビティ (平板断熱 + すきま 500 K) ==="
python3 tools/make_case.py --run run_0006_T1_cav --mesh t1_m5b --plate-thermal adiabatic --main-steps 30000 $COMMON
long=run_0006_T1_cav_main; mkdir -p $long
cp run_0006_T1_cav/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} $long/
echo "run_0006_T1_cav" > $long/CONTINUED_FROM
python3 - $long <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 120000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 10000',s)
p.write_text(s)
PY
../../solver_density_cuda/tools/run_case.sh $long >/dev/null 2>&1 || true
python3 tools/cavity_eval.py $long --series --qfp-from run_0005_T0_iso500_main | tail -12
echo "=== T3 チェーン done ==="
