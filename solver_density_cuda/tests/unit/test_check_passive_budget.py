#!/usr/bin/env python3
"""check_passive_budget.py の失敗系試験 (codex plan-7 M1, plan-8 M1/M2): **実際の log 書式**を parse → evaluate まで通し、
NaN, 合計超過, FCT 記録欠落, 総量 0 の補正, 不完全な記録, 閉合不成立, 総量と増分の不一致, 低次残差超過 が FAIL になること。"""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cpb', os.path.join(here, '..', '..', 'tools', 'check_passive_budget.py'))
cpb = importlib.util.module_from_spec(spec); spec.loader.exec_module(cpb)

FLOOR = ("[passive] step {step} floorCorr {nm:<8s} per-step(avg 10): lo 0.000e+00 hi 0.000e+00 abs 0.000e+00 | cumulative: lo 0.000000e+00 hi 0.000000e+00 abs {fabs} "
         "| total {tot} rel(abs/total) {frel} | limCorr per-step 0.000e+00 cumulative abs 0.000000e+00 signed 0.000000e+00 rel {lrel} cells 0 thetaMin(interval) 1.0000 initialTotal {init}")
FCT = ("[passive]   fctCorr {nm:<8s} cumulative: dropped antidiffusion 1.000000e-03 (rel 1.000000e-03) faces 12 prelimited 0.000000e+00 pinCorr 0.000000e+00 (rel 0.000000e+00) "
       "baseViol 0.000000e+00 (rel {base}) bndFluxSigned {bnd} bndDropped 0.000000e+00 upperViol 0.000000e+00 (rel 0.000000e+00) | budget: srcHist {src} remSigned {rem} remAbs 0.000000e+00 (rel {remrel}) increment {inc} "
       "| qL rel-residual interval-max 1.00e-07 run-max {relres} (sweeps last 3) HO residual rel interval-max 1.00e-07 run-max 1.00e-07")
CLAMP = "[passive]   clampBudget species 0 cumulative (signed/abs, rel to total): g 0.000000e+00/0.000000e+00 (0.000000e+00) Q0 0.000000e+00/0.000000e+00 (0.000000e+00) Q1 0.000000e+00/0.000000e+00 ({q1}) Q2 0.000000e+00/0.000000e+00 (0.000000e+00)"

def log(nm='roXi', step=100, tot='1.000000e+00', init='1.000000e+00', frel='0.000000e+00', lrel='0.000000e+00', fabs='0.000000e+00', fct=True,
        base='0.000000e+00', bnd='0.000000e+00', src='0.000000e+00', rem='0.000000e+00', remrel='0.000000e+00', inc='0.000000e+00', relres='1.00e-07', clamp=None, active=True):
    lines = []
    if active: lines.append('[passiveFct] active: post-step conservative FCT for 1 passive scalars (prelimit 0, sweeps 100, tol 1.0e-06)')
    lines.append(FLOOR.format(step=step, nm=nm, tot=tot, init=init, frel=frel, lrel=lrel, fabs=fabs))
    if fct: lines.append(FCT.format(nm=nm, base=base, bnd=bnd, src=src, rem=rem, remrel=remrel, inc=inc, relres=relres))
    if clamp is not None: lines.append(CLAMP.format(q1=clamp))
    return lines

def run(lines, final=100, expect='auto'):
    last, fct, clamp, nproj, ndeg, act, ls = cpb.parse_lines(lines)
    e = (expect == 'yes') or (expect == 'auto' and act)
    return cpb.evaluate(last, fct, clamp, 1e-6, 1e-4, False, e, final, out=lambda *_: None)

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)
check(run(log()) is True, 'clean case (real format, %-8s name padding) must PASS')
check(run(log(frel='nan')) is False, 'NaN must FAIL')
check(run(log(frel='6.000000e-07', lrel='6.000000e-07')) is False, 'sum 1.2e-6 must FAIL')
check(run(log(fct=False)) is False, 'missing FCT record must FAIL when [passiveFct] active')
check(run(log(fct=False, active=False)) is True, 'missing FCT record is fine when FCT never activated')
check(run(log(tot='0.000000e+00', init='0.000000e+00', frel='1.000000e+00')) is False, 'zero total with correction (rel=1) must FAIL')
check(run(log(step=91), final=100) is False, 'record at step 91 of a 100-step run must FAIL (incomplete)')
check(run(log(inc='1.000000e+00')) is False, 'increment 1 with zero boundary/source/remainder must FAIL (closure and total-vs-increment)')
check(run(log(bnd='nan')) is False, 'NaN boundary flux must FAIL')
check(run(log(relres='1.00e+00')) is False, 'low-order residual 1 must FAIL')
check(run(log(base='6.000000e-07', nm='roQ1_0', clamp='6.000000e-07')) is False, 'base 6e-7 + clamp 6e-7 on the same component must FAIL')
check(run(log(remrel='2.000000e-06')) is False, 'history remainder above tol must FAIL')
check(run(log(nm='roQ1_0', clamp='2.000000e-06')) is False, 'clamp budget above tol must FAIL')
# 閉合は成立するが総量が増分と合わない
check(run(log(tot='1.000010e+00', init='1.000000e+00', inc='1.000000e-05', src='1.000000e-05')) is True, 'consistent source-driven increase must PASS')
check(run(log(tot='1.100000e+00', init='1.000000e+00', inc='1.000000e-05', src='1.000000e-05')) is False, 'total change not matching the increment must FAIL')
print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
