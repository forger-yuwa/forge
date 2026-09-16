#!/usr/bin/env python3
"""受動種 (トレーサ・凝縮モーメント) の補正収支ゲート (plans/active/species-passive-scalar-unification.md §4.7 v5 / §6-2 / §6-6; codex plan-7 M1)。

forge_run.log の最後の `[passive]` 行群 (monitor 区間ごとに出る全期間積算) から、受動種ごとに
  floorCorr, limCorr, FCT の基点逸脱 (baseViol), ピン交換 (pinCorr), 履歴の非物理局所残り (remAbs), 実現可能性クランプの成分別 |Δ|
の**総量比の合計**が閾値 (既定 1e-6) 以下であることを判定する。FAIL 条件: 合計 > tol、いずれかが NaN/Inf、
FCT 記録が期待されるのに無い (`--expect-fct`; 既定は log の `[passiveFct] active` 行から自動判定)、最終 monitor 行が欠ける、総量 0 で補正が非ゼロ (log 側が rel=1 を出す)。
使い方: check_passive_budget.py RUN_DIR [--tol 1e-6] [--allow-lim] [--expect-fct {auto,yes,no}]
"""
import argparse, math, os, re, sys


def finite(*xs):
    return all(math.isfinite(x) for x in xs)


def parse(log):
    last, fct, clamp = {}, {}, {}
    nproj = ndeg = 0; realiz = None; cfg_fct = None; dual = None; last_step = None
    re_floor = re.compile(r'\[passive\] step (\d+) floorCorr (\S+) .*cumulative: lo (\S+) hi (\S+) abs (\S+) \| total (\S+) rel\(abs/total\) (\S+) \| limCorr per-step \S+ cumulative abs (\S+) signed (\S+) rel (\S+)')
    re_fct = re.compile(r'\[passive\]\s+fctCorr (\S+) cumulative: dropped antidiffusion (\S+) \(rel (\S+)\) faces (\S+) prelimited (\S+) pinCorr (\S+) \(rel (\S+)\) baseViol (\S+) \(rel (\S+)\) bndFluxSigned (\S+) bndDropped (\S+) \| budget: srcHist (\S+) remSigned (\S+) remAbs (\S+) \(rel (\S+)\) increment (\S+) \| interval max: qL rel-residual (\S+) \(sweeps last (\d+)\) HO residual rel (\S+)')
    re_clamp = re.compile(r'\[passive\]\s+clampBudget species (\d+) cumulative .*: g (\S+)/(\S+) \((\S+)\) Q0 (\S+)/(\S+) \((\S+)\) Q1 (\S+)/(\S+) \((\S+)\) Q2 (\S+)/(\S+) \((\S+)\)')
    re_realiz = re.compile(r'\[passive\] step (\d+) moment realizability corrections since last log: nearest-point (\d+), degenerate->monodisperse (\d+)')
    f = lambda x: float(x)
    with open(log, errors='replace') as fh:
        for line in fh:
            if line.startswith('[passiveFct] active'): cfg_fct = 1; dual = True   # FCT が実際に作動した印 (wrapper が最初の補正で出す)
            m = re_floor.search(line)
            if m:
                last[m.group(2)] = dict(step=int(m.group(1)), floor_abs=f(m.group(5)), total=f(m.group(6)), floor_rel=f(m.group(7)),
                                        lim_abs=f(m.group(8)), lim_signed=f(m.group(9)), lim_rel=f(m.group(10)))
                last_step = int(m.group(1)); continue
            m = re_fct.search(line)
            if m:
                fct[m.group(1)] = dict(dropped=f(m.group(2)), dropped_rel=f(m.group(3)), pin=f(m.group(6)), pin_rel=f(m.group(7)), base=f(m.group(8)), base_rel=f(m.group(9)),
                                       bnd_signed=f(m.group(10)), bnd_dropped=f(m.group(11)), src=f(m.group(12)), rem_signed=f(m.group(13)), rem_abs=f(m.group(14)),
                                       rem_rel=f(m.group(15)), increment=f(m.group(16)), max_relres=f(m.group(17)), sweeps=int(m.group(18)), max_rh=f(m.group(19)))
                continue
            m = re_clamp.search(line)
            if m:
                clamp[int(m.group(1))] = dict(g=f(m.group(4)), Q0=f(m.group(7)), Q1=f(m.group(10)), Q2=f(m.group(13))); continue
            m = re_realiz.search(line)
            if m: nproj += int(m.group(2)); ndeg += int(m.group(3)); realiz = int(m.group(1))
    return last, fct, clamp, nproj, ndeg, realiz, cfg_fct, dual, last_step


def evaluate(last, fct, clamp, tol, allow_lim, expect_fct, nproj, ndeg, out=print):
    ok = True
    if not last:
        out('no [passive] budget lines'); return False
    for nm, v in last.items():
        vals = [v['floor_rel'], v['lim_rel'], v['total']]
        fl = [] if finite(*vals) else ['NONFINITE']
        total = v['floor_rel'] + (0.0 if allow_lim else v['lim_rel'])
        fdesc = ''
        if expect_fct:
            fe = fct.get(nm)
            if fe is None: fl.append('NO_FCT_RECORD')
            else:
                fv = [fe['base_rel'], fe['pin_rel'], fe['rem_rel'], fe['max_relres'], fe['max_rh']]
                if not finite(*fv): fl.append('NONFINITE')
                total += fe['base_rel'] + fe['pin_rel'] + fe['rem_rel']
                fdesc = (f" | fct: dropped rel {fe['dropped_rel']:.2e} base rel {fe['base_rel']:.2e} pin rel {fe['pin_rel']:.2e} remainder rel {fe['rem_rel']:.2e}"
                         f" boundary flux(signed) {fe['bnd_signed']:.3e} qL max rel-res {fe['max_relres']:.2e} HO res {fe['max_rh']:.2e}")
        if math.isfinite(total) and total > tol: fl.append('SUM>tol')
        st = 'FAIL(' + ','.join(fl) + ')' if fl else 'ok'
        ok = ok and not fl
        out(f"  {nm:8s}: total {v['total']:.6e} floor rel {v['floor_rel']:.2e} lim rel {v['lim_rel']:.2e}{fdesc} | sum {total:.2e} -> {st}")
    for sp, c in clamp.items():
        vals = [c[k] for k in ('g', 'Q0', 'Q1', 'Q2')]
        fl = [] if finite(*vals) else ['NONFINITE']
        ssum = sum(vals)
        if math.isfinite(ssum) and ssum > tol: fl.append('SUM>tol')
        ok = ok and not fl
        out(f"  clamp species {sp}: |Δ| rel g {c['g']:.2e} Q0 {c['Q0']:.2e} Q1 {c['Q1']:.2e} Q2 {c['Q2']:.2e} | sum {ssum:.2e} -> {'FAIL(' + ','.join(fl) + ')' if fl else 'ok'}")
    out(f'  realizability corrections: nearest-point {nproj}, degenerate->monodisperse {ndeg}')
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--tol', type=float, default=1.0e-6)
    ap.add_argument('--allow-lim', action='store_true', help='limCorr (定常の起動緩和) を合否に含めない')
    ap.add_argument('--expect-fct', choices=['auto', 'yes', 'no'], default='auto')
    a = ap.parse_args()
    log = os.path.join(a.run_dir, 'forge_run.log')
    if not os.path.exists(log):
        print(f'[{a.run_dir}] NO forge_run.log'); sys.exit(2)
    last, fct, clamp, nproj, ndeg, realiz, cfg_fct, dual, last_step = parse(log)
    if not last:
        print(f'[{a.run_dir}] no [passive] budget lines (passiveScalarScheme 1 の run のみ対象)'); sys.exit(2)
    expect = (a.expect_fct == 'yes') or (a.expect_fct == 'auto' and cfg_fct == 1 and bool(dual))
    print(f'passive budget gate for {a.run_dir} (tol {a.tol:g}, last monitor step {last_step}, FCT record expected: {expect})')
    ok = evaluate(last, fct, clamp, a.tol, a.allow_lim, expect, nproj, ndeg)
    print('VERDICT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
