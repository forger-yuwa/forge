#!/bin/bash
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
RD=run_0006_xcheck_forge
rm -f $RD/res_*.h5 $RD/res_*.xmf $RD/*.log $RD/residual_history*
python3 - <<'PY'
import h5py, numpy as np, json
st=json.load(open('xcheck_state.json'))
with h5py.File('run_0006_xcheck_forge/mesh.h5','r+') as f:
    n=f['/VALUE/ro'].shape[0]; wd=f['/VALUE/wall_dist'][:]
    ro=np.full(n,st['ro_inf']); u=st['U_inf']*np.tanh(np.maximum(wd,0)/2e-4)
    T=np.full(n,st['T_inf']); wall=wd<=0
    u[wall]=0.0; T[wall]=st['Tw']; ro[wall]=st['p_inf']/(st['R']*st['Tw'])
    cv=st['cp']-st['R']
    f['/VALUE/ro'][:]=ro.astype(np.float32)
    f['/VALUE/roUx'][:]=(ro*u).astype(np.float32)
    f['/VALUE/roUy'][:]=np.zeros(n,np.float32); f['/VALUE/roUz'][:]=np.zeros(n,np.float32)
    f['/VALUE/roe'][:]=(ro*(cv*T+0.5*u**2)).astype(np.float32)
    k=1.5*(0.005*st['U_inf'])**2; om=st['ro_inf']*k/(10*st['mu'])
    for nm,v in (('roK',k),('roOmega',om)):
        key=f'/VALUE/{nm}'
        if key in f: f[key][:]=(ro*v).astype(np.float32)
        else: f.create_dataset(key,data=(ro*v).astype(np.float32))
PY
# 段階起動: lam -> soft -> mid -> ramp -> main
run_stage () {  # tag conv lim cfl inner nsteps turb
  python3 - "$@" <<'PY'
import sys, re, pathlib
tag,conv,lim,cfl,inner,nsteps,turb = sys.argv[1:8]
p=pathlib.Path('run_0006_xcheck_forge/solverConfig.yaml'); s=p.read_text()
s=re.sub(r'nStepOuter: \d+', f'nStepOuter: {nsteps}', s)
s=re.sub(r'outStepInterval: \d+', f'outStepInterval: {nsteps}', s)
s=re.sub(r'cfl: [0-9.]+, cfl_pseudo: [0-9.]+', f'cfl: {cfl}, cfl_pseudo: {cfl}', s)
s=re.sub(r'nStepInner: \d+', f'nStepInner: {inner}', s)
s=re.sub(r'convMethod: \d+, limiter: \d+', f'convMethod: {conv}, limiter: {lim}', s)
s=re.sub(r'model: "\w+"', f'model: "{turb}"', s)
p.write_text(s)
PY
  ../../solver_density_cuda/tools/run_case.sh $RD >/dev/null 2>&1 || true
  for f in forge_run.log residual_history.csv CONVERGENCE_VERDICT.txt; do
    [ -f "$RD/$f" ] && mv "$RD/$f" "$RD/${f%.*}_$1.${f##*.}"
  done
  [ -f "$RD/res_$6.h5" ] && python3 ../../solver_density_cuda/tools/interp_field.py "$RD/res_$6.h5" "$RD/mesh.h5" >/dev/null
  python3 - "$6" <<'PY'
import sys, h5py, numpy as np, json
st=json.load(open('xcheck_state.json'))
k=1.5*(0.005*st['U_inf'])**2; om=st['ro_inf']*k/(10*st['mu'])
with h5py.File('run_0006_xcheck_forge/mesh.h5','r+') as f:
    ro=f['/VALUE/ro'][:].astype(float)
    for nm,v in (('roK',k),('roOmega',om)):
        key=f'/VALUE/{nm}'
        if key not in f: f.create_dataset(key,data=(ro*v).astype(np.float32))
        elif float(np.max(np.abs(f[key][:])))==0.0: f[key][:]=(ro*v).astype(np.float32)
PY
  echo "  stage $1 done"
}
run_stage lam   0 0 0.2 10 3000 none
run_stage soft  0 0 0.2 10 3000 sst
run_stage mid   0 0 0.5 10 3000 sst
run_stage ramp0 1 2 0.3  6 2500 sst
run_stage ramp1 1 2 0.6  6 2500 sst
run_stage ramp2 1 2 1.0  6 2500 sst
run_stage main  1 2 1.5  4 120000 sst
echo "=== chain done ==="
