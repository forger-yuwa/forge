#!/usr/bin/env python3
"""tools/twophase_state.py の試験: 二相拡散の実効状態の判定表 (plan condensation-two-phase-default §4-2, G5) と、
check_solver_config.py・stage_manifest.py への反映。負例 (dual-time 凝縮 NS で指定 ON → FAIL 等) を主に見る。

既定 (DEFAULT) が 0 の今と、S1-c で 1 にした後の両方を、モジュールの DEFAULT を差し替えて試す。"""
import importlib.util, os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
tools = os.path.join(here, '..', '..', 'tools')


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(tools, name + '.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


tps = _load('twophase_state')
csc = _load('check_solver_config')
sm = _load('stage_manifest')

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)


def cfg(tp=None, **kw):
    """包絡内 (active) の TP carrier 凝縮 NS 定常 config。kw で 1 点ずつ外す。tp=None は condTwoPhaseDiffusion 省略。"""
    cd = {'condensation': kw.pop('cond', 1), 'nCondSpecies': kw.pop('ncond', 1), 'condModel': 1,
          'condensationSpecies': kw.pop('gas', 'H2O')}
    for k in ('condEquilibrium', 'condLimiterMode', 'condVaporMassFraction'):
        if k in kw: cd[k] = kw.pop(k)
    if tp is not None: cd['condTwoPhaseDiffusion'] = tp
    dT = {'cfl_pseudo': 2.0, 'speciesFaceReconstruction': 2, 'speciesImplicitCoupling': kw.pop('coupling', 1)}
    if 'scheme' in kw: dT['passiveScalarScheme'] = kw.pop('scheme')
    y = {'solver': 'SLAU', 'mesh': {'discretization': 'node'}, 'space': {'convMethod': 1},
         'physProp': {'thermalMethod': kw.pop('thermal', 2), 'viscMethod': kw.pop('visc', 1),
                      'species': kw.pop('species', ['N2', 'O2', 'H2O'])},
         'condensation': cd,
         'time': {'unsteady': kw.pop('unsteady', 0), 'dualTime': kw.pop('dualTime', 0), 'timeIntegration': kw.pop('ti', 11),
                  'last': {'nStepOuter': 100}, 'deltaT': dT}}
    if 'axisym' in kw: y['mesh']['isAxisymmetric'] = kw.pop('axisym')
    assert not kw, kw
    return y


PERIODIC = {'left': {'physID': 1, 'kind': 'periodic'}, 'right': {'physID': 2, 'kind': 'periodic'}}
WALLS = {'wall': {'physID': 1, 'kind': 'wall'}}

# ---- 既定値の一致 (ソルバの kCondTwoPhaseDiffusionDefault) ----
hpp = open(os.path.join(here, '..', '..', 'input', 'solverConfig.hpp')).read()
m = re.search(r'kCondTwoPhaseDiffusionDefault\s*=\s*(\d+)', hpp)
check(m is not None and int(m.group(1)) == tps.DEFAULT,
      f'twophase_state.DEFAULT ({tps.DEFAULT}) must equal solverConfig::kCondTwoPhaseDiffusionDefault ({m.group(1) if m else "?"})')

# ---- 判定表 (状態は指定値に依らない) ----
cases = [
    ('active', {}, None),
    ('inactive-a', {'cond': 0}, None),
    ('inactive-a', {'visc': 0}, None),
    ('inactive-a', {'visc': 0, 'unsteady': 1, 'dualTime': 1}, None),       # dual-time + Euler は (a)
    ('inactive-b', {'thermal': 1, 'species': [], 'gas': 'N2', 'condVaporMassFraction': 0.7671}, None),   # CPG carrier
    ('inactive-b', {'species': ['H2O'], 'gas': None}, None),                # pure
    ('inactive-b', {'condEquilibrium': 2}, None),
    ('unsupported-c', {'unsteady': 1, 'dualTime': 1}, None),
    ('unsupported-c', {'unsteady': 1, 'dualTime': 0, 'ti': 4}, None),       # RK
    ('unsupported-c', {'ti': 4}, None),                                     # 定常陽解法
    ('unsupported-c', {'coupling': 2}, None),
    ('unsupported-c', {'scheme': 0}, None),
    ('unsupported-c', {'condEquilibrium': 1}, None),
    ('unsupported-c', {'condLimiterMode': 0}, None),
    ('unsupported-c', {'ncond': 2}, None),
    ('unsupported-c', {'axisym': 1}, None),
    ('unsupported-c', {}, PERIODIC),
    ('active', {}, WALLS),
]
for want, kw, bc in cases:
    for tp in (None, 0, 1):
        y = cfg(tp, **dict(kw))
        if kw.get('gas', 'x') is None: y['condensation'].pop('condensationSpecies')
        r = tps.classify(y, bc)
        check(r['state'] == want, f'{kw} bc={bool(bc)} tp={tp}: state {r["state"]} (want {want})')
        check(r['effective'] == (1 if (want == 'active' and (tp == 1 or (tp is None and tps.DEFAULT == 1))) else 0),
              f'{kw} tp={tp}: effective {r["effective"]}')

r = tps.classify(cfg(None)); check(r['requested'] == 'omitted' and r['request_on'] == (tps.DEFAULT == 1), f'omitted: {r}')
r = tps.classify(cfg(0));    check(r['requested'] == '0' and not r['request_on'], f'explicit 0: {r}')
r = tps.classify(cfg(1));    check(r['requested'] == '1' and r['request_on'] and r['effective'] == 1, f'explicit 1: {r}')
r = tps.classify({});        check(r['state'] == 'inactive-a' and r['effective'] == 0, f'empty config: {r}')

# ---- check_solver_config (G5): 既定 0 の今と既定 1 の後 ----
K = 'condensation.condTwoPhaseDiffusion'
def keys(res): return {k for k, _ in res}

for default in (0, 1):
    tps_saved = csc._twophase_state().DEFAULT
    csc._twophase_state().DEFAULT = default
    try:
        dual_ns = dict(unsteady=1, dualTime=1)
        # dual-time + 凝縮 NS + 明示 1 → FAIL
        f, w = csc.check(cfg(1, **dual_ns))
        check(K in keys(f), f'[default {default}] dual-time condensing NS with explicit 1 must FAIL')
        # 省略は既定 ON のときだけ FAIL
        f, w = csc.check(cfg(None, **dual_ns))
        check((K in keys(f)) == (default == 1), f'[default {default}] dual-time condensing NS with the key omitted: FAIL only after the default flip')
        # dual-time + Euler は通す (明示 1 でも)
        for tp in (None, 1):
            f, w = csc.check(cfg(tp, visc=0, **dual_ns))
            check(K not in keys(f) and K not in keys(w), f'[default {default}] dual-time Euler (tp={tp}) must pass')
        # 明示 0 は既定 ON の後だけ WARN、FAIL にはしない
        f, w = csc.check(cfg(0))
        check(K not in keys(f) and ((K in keys(w)) == (default == 1)), f'[default {default}] explicit 0 in the envelope: WARN only after the default flip')
        # 明示 0 は未対応構成でも FAIL にしない (旧作用素の逃げ道)
        f, w = csc.check(cfg(0, **dual_ns))
        check(K not in keys(f), f'[default {default}] explicit 0 must open the legacy path even for dual-time')
        # 軸対称・周期 + 明示 1 → FAIL
        f, w = csc.check(cfg(1, axisym=1))
        check(K in keys(f), f'[default {default}] axisymmetric with explicit 1 must FAIL')
        f, w = csc.check(cfg(1), PERIODIC)
        check(K in keys(f), f'[default {default}] periodic with explicit 1 must FAIL')
        # 構造的に適用できない (b) + 明示 1 → WARN
        f, w = csc.check(cfg(1, condEquilibrium=2))
        check(K not in keys(f) and K in keys(w), f'[default {default}] condEquilibrium 2 with explicit 1 must WARN (inactive-b)')
        # 包絡内 + 明示 1 → 何も出ない
        f, w = csc.check(cfg(1))
        check(K not in keys(f) and K not in keys(w), f'[default {default}] explicit 1 in the envelope must be clean')
    finally:
        csc._twophase_state().DEFAULT = tps_saved

# ---- stage_manifest: OFF と ON は別区間 ----
import yaml
def stage(y, bc=None):
    return sm.stage_key(yaml.safe_dump(y), yaml.safe_dump(bc) if bc else '')
k_off, k_on = stage(cfg(0)), stage(cfg(1))
check(k_off.get('twophase_diffusion_effective') == '0' and k_on.get('twophase_diffusion_effective') == '1',
      f'manifest key must carry the effective value (off {k_off.get("twophase_diffusion_effective")}, on {k_on.get("twophase_diffusion_effective")})')
check(k_off != k_on, 'OFF and ON stages must be different segments')
man = {'stages': [{'tag': 'off', 'key': k_off}, {'tag': 'on', 'key': k_on}]}
check(len(sm.segments(man)) == 2, 'segments() must split OFF -> ON')
# 不活性な指定の違い (Euler で 0 と 1) は同じ区間
check(stage(cfg(0, visc=0)) == stage(cfg(1, visc=0)), 'inactive (Euler) 0 vs 1 must stay in one segment')
# 凝縮でない run の署名は従来どおり (キーを足さない)
y = cfg(None, cond=0); y.pop('condensation')
check('twophase_diffusion_effective' not in stage(y), 'non-condensing runs must not get the new key')

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
