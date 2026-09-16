#!/usr/bin/env python3
"""check_passive_budget.py の失敗系試験 (codex plan-7 M1): NaN, 合計超過, FCT 記録欠落, 総量 0 の補正 が FAIL になること。"""
import importlib.util, io, os, sys, math
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cpb', os.path.join(here, '..', '..', 'tools', 'check_passive_budget.py'))
cpb = importlib.util.module_from_spec(spec); spec.loader.exec_module(cpb)

def run(last, fct, clamp, expect):
    return cpb.evaluate(last, fct, clamp, 1e-6, False, expect, 0, 0, out=lambda *_: None)

base = dict(step=100, floor_abs=0.0, total=1.0, floor_rel=0.0, lim_abs=0.0, lim_signed=0.0, lim_rel=0.0)
fctok = dict(dropped=0, dropped_rel=0, pin=0, pin_rel=0, base=0, base_rel=0, bnd_signed=0, bnd_dropped=0, src=0, rem_signed=0, rem_abs=0, rem_rel=0, increment=0, max_relres=1e-7, sweeps=3, max_rh=1e-7)
fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)
check(run({'roXi': base}, {'roXi': fctok}, {}, True) is True, 'clean case must PASS')
check(run({'roXi': dict(base, floor_rel=float('nan'))}, {'roXi': fctok}, {}, True) is False, 'NaN must FAIL')
check(run({'roXi': dict(base, floor_rel=6e-7, lim_rel=6e-7)}, {'roXi': fctok}, {}, True) is False, 'sum 1.2e-6 must FAIL')
check(run({'roXi': base}, {}, {}, True) is False, 'missing FCT record must FAIL when expected')
check(run({'roXi': base}, {}, {}, False) is True, 'missing FCT record is fine when not expected')
check(run({'roXi': dict(base, total=0.0, floor_rel=1.0)}, {'roXi': fctok}, {}, True) is False, 'zero total with correction (rel=1) must FAIL')
check(run({'roXi': base}, {'roXi': dict(fctok, rem_rel=2e-6)}, {}, True) is False, 'history remainder above tol must FAIL')
check(run({'rog_0': base}, {'rog_0': fctok}, {0: dict(g=0, Q0=0, Q1=2e-6, Q2=0)}, True) is False, 'clamp budget above tol must FAIL')
print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
