#!/usr/bin/env python3
"""受動種 (トレーサ・凝縮モーメント) の補正収支ゲート (plans/active/species-passive-scalar-unification.md §4.7 v5 / §6-2 / §6-6)。

forge_run.log の最後の `[passive]` 行群から、受動種ごとの全期間積算
  floorCorr (abs/total), limCorr (rel), fctCorr の baseViol (rel), pinCorr, 実現可能性クランプの成分別 |Δ| (rel)
を読み、閾値 (既定 1e-6) を超える成分があれば FAIL。非定常 (dual-time) 計算の「有界化の穴」を事後に拒否する判定 (codex plan-5/6)。
使い方: check_passive_budget.py RUN_DIR [--tol 1e-6] [--allow-lim] (定常の起動緩和 limCorr を合否に含めない)
"""
import argparse, os, re, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--tol', type=float, default=1.0e-6)
    ap.add_argument('--allow-lim', action='store_true', help='limCorr (定常の起動緩和) を合否に含めない')
    a = ap.parse_args()
    log = os.path.join(a.run_dir, 'forge_run.log')
    if not os.path.exists(log):
        print(f'[{a.run_dir}] NO forge_run.log'); sys.exit(2)
    last = {}   # name -> dict of latest values
    fct = {}
    clamp = {}
    realiz = None
    re_floor = re.compile(r'\[passive\] step (\d+) floorCorr (\S+) .*cumulative: lo (\S+) hi (\S+) abs (\S+) \| total (\S+) rel\(abs/total\) (\S+) \| limCorr per-step \S+ cumulative abs (\S+) signed (\S+) rel (\S+)')
    re_fct = re.compile(r'\[passive\]\s+fctCorr (\S+) cumulative: dropped antidiffusion (\S+) \(rel (\S+)\) faces (\S+) prelimited (\S+) pinCorr (\S+) bndExchange (\S+) baseViol (\S+) \(rel (\S+)\)')
    re_clamp = re.compile(r'\[passive\]\s+clampBudget species (\d+) cumulative .*: g (\S+)/(\S+) \((\S+)\) Q0 (\S+)/(\S+) \((\S+)\) Q1 (\S+)/(\S+) \((\S+)\) Q2 (\S+)/(\S+) \((\S+)\)')
    re_realiz = re.compile(r'\[passive\] step (\d+) moment realizability corrections since last log: nearest-point (\d+), degenerate->monodisperse (\d+)')
    nproj = ndeg = 0
    with open(log, errors='replace') as f:
        for line in f:
            m = re_floor.search(line)
            if m:
                last[m.group(2)] = dict(step=int(m.group(1)), floor_abs=float(m.group(5)), total=float(m.group(6)), floor_rel=float(m.group(7)),
                                        lim_abs=float(m.group(8)), lim_signed=float(m.group(9)), lim_rel=float(m.group(10)))
                continue
            m = re_fct.search(line)
            if m:
                fct[m.group(1)] = dict(dropped=float(m.group(2)), dropped_rel=float(m.group(3)), pin=float(m.group(6)), bnd=float(m.group(7)), base=float(m.group(8)), base_rel=float(m.group(9)))
                continue
            m = re_clamp.search(line)
            if m:
                clamp[int(m.group(1))] = dict(g=float(m.group(4)), Q0=float(m.group(7)), Q1=float(m.group(10)), Q2=float(m.group(13)))
                continue
            m = re_realiz.search(line)
            if m:
                nproj += int(m.group(2)); ndeg += int(m.group(3)); realiz = int(m.group(1))
    if not last:
        print(f'[{a.run_dir}] no [passive] budget lines (passiveScalarScheme 1 の dual-time/定常 run のみ対象)'); sys.exit(2)
    ok = True
    print(f'passive budget gate for {a.run_dir} (tol {a.tol:g}, last step {max(v["step"] for v in last.values())})')
    for nm, v in last.items():
        flags = []
        if v['floor_rel'] > a.tol: flags.append('FLOOR')
        if not a.allow_lim and v['lim_rel'] > a.tol: flags.append('LIM')
        f = fct.get(nm)
        fs = ''
        if f:
            if f['base_rel'] > a.tol: flags.append('BASE')
            pin_rel = f['pin'] / v['total'] if v['total'] > 0 else 0.0
            if pin_rel > a.tol: flags.append('PIN')
            fs = f" | fct: dropped rel {f['dropped_rel']:.2e} base-violation rel {f['base_rel']:.2e} pin rel {pin_rel:.2e} boundary exchange {f['bnd']:.2e}"
        st = 'FAIL(' + ','.join(flags) + ')' if flags else 'ok'
        ok = ok and not flags
        print(f"  {nm:8s}: total {v['total']:.6e} floor rel {v['floor_rel']:.2e} lim rel {v['lim_rel']:.2e}{fs}  -> {st}")
    for sp, c in clamp.items():
        flags = [k for k in ('g', 'Q0', 'Q1', 'Q2') if c[k] > a.tol]
        ok = ok and not flags
        print(f"  clamp species {sp}: |Δ| rel g {c['g']:.2e} Q0 {c['Q0']:.2e} Q1 {c['Q1']:.2e} Q2 {c['Q2']:.2e}  -> {'FAIL(' + ','.join(flags) + ')' if flags else 'ok'}")
    if realiz is not None:
        print(f'  realizability corrections: nearest-point {nproj}, degenerate->monodisperse {ndeg}')
    print('VERDICT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
