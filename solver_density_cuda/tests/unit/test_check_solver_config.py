#!/usr/bin/env python3
"""check_solver_config.py の試験: 残差では気づけない設定ミス (dual-time の緩和・不活性な S3・作動しない FCT・存在しないキー) を FAIL/WARN で止めること。"""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('csc', os.path.join(here, '..', '..', 'tools', 'check_solver_config.py'))
csc = importlib.util.module_from_spec(spec); spec.loader.exec_module(csc)

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)

def cfg(**kw):
    dT = {'dt': 8e-6, 'cfl_pseudo': kw.pop('cflp', 12.0), 'passiveScalarScheme': kw.pop('scheme', 1),
          'passiveFct': kw.pop('fct', 1), 'speciesFaceReconstruction': kw.pop('sfr', 2)}
    for k in ('implicitRelax', 'speciesImplicitCoupling'):
        if k in kw: dT[k] = kw.pop(k)
    y = {'solver': kw.pop('solver', 'SLAU'), 'mesh': {}, 'space': {'convMethod': kw.pop('conv', 1), 'limiter': 1},
         'physProp': {'species': ['A', 'B']}, 'condensation': {'condensation': kw.pop('cond', 1), 'nCondSpecies': 1},
         'time': {'unsteady': kw.pop('unsteady', 1), 'dualTime': kw.pop('dualTime', 1), 'timeIntegration': kw.pop('ti', 11),
                  'nSubIterDualTime': kw.pop('nsub', 20), 'last': {'nStepOuter': 100}, 'deltaT': dT}}
    for k, v in kw.items(): y[k] = v
    return y

def keys(res): return {k for k, _ in res}

f, w = csc.check(cfg())
check(not f and not w, f'recommended dual-time settings must be clean (fails {keys(f)}, warns {keys(w)})')

f, w = csc.check(cfg(implicitRelax=0.7))
check('time.deltaT.implicitRelax' in keys(f), 'implicitRelax < 1 in dual-time must FAIL')

f, w = csc.check(cfg(implicitRelax=0.7, unsteady=0, dualTime=0, ti=4))
check('time.deltaT.implicitRelax' not in keys(f), 'implicitRelax 0.7 in a steady run must be accepted')

f, w = csc.check(cfg(cflp=2.0))
check('time.deltaT.cfl_pseudo' in keys(w), 'small cfl_pseudo in dual-time must WARN')

f, w = csc.check(cfg(nsub=5))
check('time.nSubIterDualTime' in keys(w), 'too few sub-iterations must WARN')

f, w = csc.check(cfg(conv=0))
check('space.convMethod' in keys(f), 'speciesFaceReconstruction 2 with convMethod 0 must FAIL (S3 inactive)')

f, w = csc.check(cfg(scheme=0))
check('time.deltaT.passiveScalarScheme' in keys(f), 'condensation dual-time with scheme 0 must FAIL (no BDF term on the moments)')

f, w = csc.check(cfg(solver='KEEP'))
check('time.deltaT.passiveFct' in keys(w), 'passiveFct 1 with a non-SLAU solver must WARN (never activates)')

f, w = csc.check(cfg(sfr=0, fct=1))
check('time.deltaT.passiveFct' in keys(w), 'passiveFct 1 with speciesFaceReconstruction 0 must WARN')

f, w = csc.check(cfg(speciesImplicitCoupling=0))
check('time.deltaT.speciesImplicitCoupling' in keys(w), 'S3 with coupling 0 must WARN')

y = cfg(); y['turbulence'] = {'model': 'SST', 'kInf': 1.0, 'omegaInf': 1000.0}
f, w = csc.check(y)
check('turbulence.kInf' in keys(w) and 'turbulence.omegaInf' in keys(w), 'keys the solver never reads must WARN (kInf/omegaInf are not real keys)')

y = cfg(); y['turbulence'] = {'model': 'SST', 'kInit': 1.0, 'omegaInit': 1000.0}
f, w = csc.check(y)
check('turbulence.kInit' not in keys(w), 'the real initial-value keys must not be flagged')

# 誤配置の検出 (codex plan-2 M5 の実例)。末端名だけの照合では全部素通りしていた。
for path, val in [('lowMachPrecond', 2), ('physProp.speciesFaceReconstruction', 2), ('space.keepDissType', 1)]:
    y = cfg(); node = y
    parts = path.split('.')
    for seg in parts[:-1]: node = node.setdefault(seg, {})
    node[parts[-1]] = val
    f, w = csc.check(y)
    check(path in keys(w), f'misplaced key {path} must WARN (silently ignored where it is written)')

# 起動時に拒否されるキーは WARN でなく FAIL
y = cfg(); y['mesh']['gradLSQDegenThresh'] = 1e-6
f, w = csc.check(y)
check('mesh.gradLSQDegenThresh' in keys(f), 'a key the solver rejects at startup must FAIL')

# 読まれないが必須として全 run が書いているキーは、誤検出しない
y = cfg(); y['time']['last']['time'] = 1.0
f, w = csc.check(y)
check('time.last.time' in keys(w), 'time.last.time is not read by the solver and must WARN')

# 正しい位置のキーは素通りする (偽陽性の回帰)
y = cfg(); y['mesh']['primPack'] = 1; y['time']['deltaT']['blockDPLURDqPack'] = 1
f, w = csc.check(y)
check(not keys(f) and 'mesh.primPack' not in keys(w), 'the restored opt-in performance switches must be accepted')

y = cfg(); y['mesh']['bndFirstOrder'] = 1
f, w = csc.check(y)
check('mesh.bndFirstOrder' in keys(f), 'bndFirstOrder must FAIL (banned)')

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
