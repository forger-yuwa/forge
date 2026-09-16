#!/usr/bin/env python3
"""受動種 (トレーサ・凝縮モーメント) の補正収支ゲート (plans/active/species-passive-scalar-unification.md §4.7 v6 / §6-2 / §6-6; codex plan-7 M1, plan-8 M1/M2)。

forge_run.log の最後の `[passive]` 行群 (monitor 区間ごと + 終了時の全期間積算) から、受動種ごとに
  floorCorr + limCorr + FCT の基点逸脱 (baseViol) + ピン交換 (pinCorr) + 上限条件の逸脱 (upperViol) + 履歴の非物理局所残り (remAbs)
  + 実現可能性クランプの成分別 |Δ| (同じ成分に合算)
の**総量比の合計**が閾値 (既定 1e-6) 以下であることを判定する。さらに
  - 全生値が有限 (solver 側の非有限フラグも)、FCT 記録・モーメントのクランプ記録が同じ step に無い → FAIL
  - 最後の収支記録の step が run の最終 step (residual_history.csv の最終 outer step + 1) に一致しない → FAIL (不完全な記録)
  - 収支の閉合: 総増分 = −境界流束 + ソース履歴 + 残り (log の生値) の残差が総量比 tol を超える → FAIL
  - 独立照合: (全後処理後の最終総量 − 計算開始前の総量; root のみ) と 総増分 の差が総量比 tol を超える → FAIL
  - 低次陰解の受入 (全期間の最大相対線形残差) と HO 残差の全期間最大 → tol_lin (既定 1e-4) を超えたら FAIL
  - 総量 0 で補正が非ゼロ (log 側が rel=1) → FAIL
使い方: check_passive_budget.py RUN_DIR [--tol 1e-6] [--tol-lin 1e-4] [--mode auto|fct|conservative|unsteady|steady] [--no-field]
全 [passive] 行の全数値トークンを読込時点で有限性検査し (per-step 部分も)、一度でも非有限・解析不能なら FAIL を保持する。終了場は res_<nStepOuter>.h5 に固定、CSV は必須。
必須成分・FCT 作動条件・終了 step は solverConfig.yaml から確定する (passive_gate_common)。確定場の有界性・実現可能性 (check_passive_field) も併せて判定する。
"""
import argparse, csv, math, os, re, sys


def finite(*xs):
    return all(isinstance(x, (int, float)) and math.isfinite(x) for x in xs)


RE_FLOOR = re.compile(r'\[passive\] step (\d+) floorCorr (\S+)\s+.*cumulative: lo (\S+) hi (\S+) abs (\S+) \| total (\S+) rel\(abs/total\) (\S+) \| limCorr per-step \S+ cumulative abs (\S+) signed (\S+) rel (\S+) cells (\S+) thetaMin\(interval\) (\S+) initialTotal (\S+)')
RE_FCT = re.compile(r'\[passive\]\s+fctCorr (\S+)\s+cumulative: dropped antidiffusion (\S+) \(rel (\S+)\) faces (\S+) prelimited (\S+) pinCorr (\S+) \(rel (\S+)\) baseViol (\S+) \(rel (\S+)\) bndFluxSigned (\S+) bndDropped (\S+) upperViol (\S+) \(rel (\S+)\) \| budget: srcHist (\S+) remSigned (\S+) remAbs (\S+) \(rel (\S+)\) increment (\S+) \| qL rel-residual interval-max (\S+) run-max (\S+) \(sweeps last (\d+)\) HO residual rel interval-max (\S+) run-max (\S+)(?: nonfinite (\d))?')
RE_CLAMP = re.compile(r'\[passive\]\s+clampBudget species (\d+) cumulative .*: g (\S+)/(\S+) \((\S+)\) Q0 (\S+)/(\S+) \((\S+)\) Q1 (\S+)/(\S+) \((\S+)\) Q2 (\S+)/(\S+) \((\S+)\)')
RE_REALIZ = re.compile(r'\[passive\] step (\d+) moment realizability corrections since last log: nearest-point (\d+), degenerate->monodisperse (\d+)')


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float('nan')


NUM_RE = re.compile(r'(?<![A-Za-z_])[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?![A-Za-z_])')
RE_INIT = re.compile(r'\[passive\] initial total (\S+)\s+(\S+)')


def parse_lines(lines):
    """[passive] 行を全部読む。戻り (last, fct, clamp, nproj, ndeg, fct_active, last_step, problems)。
    全 [passive] 行の全数値トークンの有限性 (per-step 部分も) を読込時点で検査し、一度でも非有限・解析不能があれば problems に残す (最後まで失敗)。"""
    last, fct, clamp = {}, {}, {}
    nproj = ndeg = 0; fct_active = False; last_step = None; problems = []
    for ln, line in enumerate(lines, 1):
        if line.startswith('[passiveFct] active'):
            fct_active = True; continue
        if not line.startswith('[passive]'):
            continue
        for tok in re.findall(r'\b(?:nan|inf|-nan|-inf|NaN|Inf|-Inf)\b', line):
            problems.append(f'line {ln}: non-finite token {tok!r}')
        for tok in NUM_RE.findall(line):
            try:
                if not math.isfinite(float(tok)): problems.append(f'line {ln}: non-finite number {tok}')
            except ValueError:
                problems.append(f'line {ln}: unparsable number {tok!r}')
        m = RE_FLOOR.search(line)
        if m:
            last[m.group(2)] = dict(step=int(m.group(1)), floor_lo=fnum(m.group(3)), floor_hi=fnum(m.group(4)), floor_abs=fnum(m.group(5)), total=fnum(m.group(6)), floor_rel=fnum(m.group(7)),
                                    lim_abs=fnum(m.group(8)), lim_signed=fnum(m.group(9)), lim_rel=fnum(m.group(10)), cells=fnum(m.group(11)), thetamin=fnum(m.group(12)), initial=fnum(m.group(13)))
            last_step = int(m.group(1)); continue
        m = RE_FCT.search(line)
        if m:
            g = m.groups()
            fct[g[0]] = dict(dropped=fnum(g[1]), dropped_rel=fnum(g[2]), faces=fnum(g[3]), prelim=fnum(g[4]), pin=fnum(g[5]), pin_rel=fnum(g[6]),
                             base=fnum(g[7]), base_rel=fnum(g[8]), bnd_signed=fnum(g[9]), bnd_dropped=fnum(g[10]), upper=fnum(g[11]), upper_rel=fnum(g[12]),
                             src=fnum(g[13]), rem_signed=fnum(g[14]), rem_abs=fnum(g[15]), rem_rel=fnum(g[16]), increment=fnum(g[17]),
                             relres_int=fnum(g[18]), relres_run=fnum(g[19]), sweeps=int(g[20]), rh_int=fnum(g[21]), rh_run=fnum(g[22]), nonfinite=int(g[23]) if g[23] is not None else -1, step=last_step)
            continue
        m = RE_CLAMP.search(line)
        if m:
            g = m.groups()
            clamp[int(g[0])] = dict(g_abs=fnum(g[2]), g=fnum(g[3]), Q0_abs=fnum(g[5]), Q0=fnum(g[6]), Q1_abs=fnum(g[8]), Q1=fnum(g[9]), Q2_abs=fnum(g[11]), Q2=fnum(g[12]),
                                    g_signed=fnum(g[1]), Q0_signed=fnum(g[4]), Q1_signed=fnum(g[7]), Q2_signed=fnum(g[10]), step=last_step)
            continue
        m = RE_REALIZ.search(line)
        if m:
            nproj += int(m.group(2)); ndeg += int(m.group(3)); continue
        if RE_INIT.search(line) or 'fctDensity' in line or 'per-step' in line:
            continue
        problems.append(f'line {ln}: unrecognised [passive] line: {line.strip()[:80]}')
    return last, fct, clamp, nproj, ndeg, fct_active, last_step, problems


def clamp_component(nm):
    m = re.match(r'ro(g|Q0|Q1|Q2)_(\d+)$', nm)
    return (int(m.group(2)), m.group(1)) if m else None


def evaluate(last, fct, clamp, tol, tol_lin, mode, required, expect_fct, fct_active, final_step, out=print):
    """mode: 'fct' (FCT 作動 run: 閉合・独立照合・残差), 'conservative' (非 FCT の非定常保存試験: 総量の変化 ≤ tol), 'steady' (定常: lim は許容し floor だけ)。
    required: config から確定した必須成分名の集合。欠落・解析失敗は FAIL。"""
    ok = True
    if not last:
        out('no [passive] budget lines'); return False
    if expect_fct and not fct_active:
        out('  FAIL: FCT is configured (scheme 1, passiveFct 1, dual-time, SFR>=2, SLAU) but the log has no [passiveFct] active line'); ok = False
    missing = sorted(set(required) - set(last))
    if missing:
        out(f"  FAIL: budget records missing for required components {missing}"); ok = False
    for nm in sorted(set(last) - set(required)):
        out(f"  FAIL: unexpected component {nm} (not in the config-derived set)"); ok = False
    if final_step is None:
        out('  FAIL: final step unknown (config nStepOuter missing)'); ok = False
    for nm in sorted(required):
        v = last.get(nm)
        if v is None: continue
        fl = []
        if not finite(*v.values()): fl.append('NONFINITE')
        if final_step is not None and v['step'] != final_step: fl.append(f"INCOMPLETE(last record step {v['step']} != final {final_step})")
        if not (v['initial'] >= 0.0): fl.append('NO_INITIAL_TOTAL')
        total = v['floor_rel'] + (0.0 if mode == 'steady' else v['lim_rel'])
        cc = clamp_component(nm)
        if cc is not None:
            cl = clamp.get(cc[0])
            if cl is None: fl.append('NO_CLAMP_RECORD')
            else:
                if not finite(*cl.values()): fl.append('NONFINITE')
                if cl['step'] != v['step']: fl.append('CLAMP_RECORD_STEP_MISMATCH')
                total += cl[cc[1]]
        fdesc = ''
        fe = fct.get(nm)
        scale = max(abs(v['total']), abs(v['initial']), 1e-300)
        if mode == 'fct':
            if fe is None:
                fl.append('NO_FCT_RECORD')
            else:
                if not finite(*fe.values()): fl.append('NONFINITE')
                if fe.get('nonfinite', -1) < 0: fl.append('NONFINITE_FLAG_MISSING')
                elif fe['nonfinite'] != 0: fl.append('SOLVER_NONFINITE')
                if fe['step'] != v['step']: fl.append('FCT_RECORD_STEP_MISMATCH')
                total += fe['base_rel'] + fe['pin_rel'] + fe['rem_rel'] + fe['upper_rel']
                if not (fe['relres_run'] <= tol_lin): fl.append(f"LOWORDER_RESIDUAL({fe['relres_run']:.1e})")
                if not (fe['rh_run'] <= tol_lin): fl.append(f"HO_RESIDUAL({fe['rh_run']:.1e})")
                closure = fe['increment'] + fe['bnd_signed'] - fe['src'] - fe['rem_signed']
                if not (abs(closure) <= tol*scale): fl.append(f"CLOSURE({closure/scale:.1e})")
                indep = (v['total'] - v['initial']) - fe['increment']
                if not (abs(indep) <= tol*scale): fl.append(f"TOTAL_VS_INCREMENT({indep/scale:.1e})")
                fdesc = (f" | fct: dropped rel {fe['dropped_rel']:.2e} base {fe['base_rel']:.2e} pin {fe['pin_rel']:.2e} upper {fe['upper_rel']:.2e} remainder {fe['rem_rel']:.2e}"
                         f" boundary flux {fe['bnd_signed']:.3e} closure {closure/scale:.1e} total-vs-increment {indep/scale:.1e} qL res(run max) {fe['relres_run']:.1e} HO res(run max) {fe['rh_run']:.1e}")
        elif mode == 'conservative':
            drift = (v['total'] - v['initial'])/scale
            fdesc = f" | conservation: (final - initial)/scale {drift:.2e}"
            if not (abs(drift) <= tol): fl.append(f"NOT_CONSERVED({drift:.1e})")
        elif mode == 'unsteady':
            fdesc = f" | (final - initial)/scale {(v['total'] - v['initial'])/scale:.2e} (境界流束・ソース込みの収支記録が無いので保存は判定不能; floor/lim と場だけ判定)"
        if not (math.isfinite(total) and total <= tol): fl.append(f'SUM>tol({total:.1e})')
        st = 'FAIL(' + ','.join(fl) + ')' if fl else 'ok'
        ok = ok and not fl
        out(f"  {nm:8s}: total {v['total']:.9e} (initial {v['initial']:.9e}) floor {v['floor_rel']:.2e} lim {v['lim_rel']:.2e}{fdesc} | sum {total:.2e} -> {st}")
    return ok


def final_step_of(run_dir):
    p = os.path.join(run_dir, 'residual_history.csv')
    if not os.path.exists(p):
        return None
    last = None
    with open(p) as f:
        for row in csv.DictReader(f):
            if row.get('phase', 'outer_begin').startswith('outer'):
                try: last = int(row['step'])
                except (KeyError, ValueError): pass
    return (last + 1) if last is not None else None


def closed_and_sourceless(run_dir, cfg):
    """全 bcond が periodic かつ凝縮ソースなし (conservative モードの自動判定)。"""
    import yaml
    p = os.path.join(run_dir, 'bcondConfig.yaml')
    if not os.path.exists(p) or cfg['nCond'] > 0:
        return False
    with open(p) as fh:
        b = yaml.safe_load(fh) or {}
    kinds = [str((v or {}).get('kind', '')) for v in b.values() if isinstance(v, dict)]
    return bool(kinds) and all(k == 'periodic' for k in kinds)


def main():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from passive_gate_common import load_config, check_field, final_res
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--tol', type=float, default=1.0e-6)
    ap.add_argument('--tol-lin', type=float, default=1.0e-4, help='低次陰解・HO 残差の全期間最大相対値の許容')
    ap.add_argument('--mode', choices=['auto', 'fct', 'conservative', 'unsteady', 'steady'], default='auto',
                    help='auto: FCT 設定なら fct; 非 FCT の dual-time は閉境界・無ソースなら conservative (総量不変) さもなくば unsteady (保存は判定不能); 定常は steady')
    ap.add_argument('--no-field', action='store_true', help='確定場の有界性・実現可能性検査を省く')
    a = ap.parse_args()
    log = os.path.join(a.run_dir, 'forge_run.log')
    if not os.path.exists(log):
        print(f'[{a.run_dir}] NO forge_run.log'); sys.exit(2)
    try:
        cfg = load_config(a.run_dir)
    except Exception as e:
        print(f'[{a.run_dir}] cannot read solverConfig.yaml: {e}'); sys.exit(2)
    if cfg['passiveScalarScheme'] != 1 or not cfg['passives']:
        print(f'[{a.run_dir}] no passive scalars on the species path (scheme {cfg["passiveScalarScheme"]}, passives {cfg["passives"]}): not a gate target'); sys.exit(2)
    with open(log, errors='replace') as f:
        last, fct, clamp, nproj, ndeg, fct_active, last_step, problems = parse_lines(f)
    mode = a.mode
    if mode == 'auto':
        if cfg['fct_configured']: mode = 'fct'
        elif cfg['unsteady'] == 1 and cfg['dualTime'] == 1: mode = 'conservative' if closed_and_sourceless(a.run_dir, cfg) else 'unsteady'
        else: mode = 'steady'
    final_cfg = int(cfg['nStepOuter']) if cfg['nStepOuter'] is not None else None
    final_csv = final_step_of(a.run_dir)
    print(f'passive budget gate for {a.run_dir} (mode {mode}, tol {a.tol:g}, tol_lin {a.tol_lin:g}, required {cfg["passives"]}, last record step {last_step}, '
          f'final step config {final_cfg} csv {final_csv}, FCT configured {cfg["fct_configured"]} active {fct_active})')
    ok = True
    for pr in problems[:20]:
        print('  FAIL(log):', pr)
    if problems: ok = False
    if final_csv is None:
        print('  FAIL: residual_history.csv missing or without outer rows'); ok = False
    elif final_cfg is not None and final_cfg != final_csv:
        print(f'  FAIL: run did not complete (config nStepOuter {final_cfg} vs csv last step+1 {final_csv})'); ok = False
    if final_res(a.run_dir, cfg) is None:
        print(f'  FAIL: final field res_{final_cfg}.h5 missing'); ok = False
    ok = evaluate(last, fct, clamp, a.tol, a.tol_lin, mode, cfg['passives'], cfg['fct_configured'], fct_active, final_cfg) and ok
    print(f'  realizability corrections: nearest-point {nproj}, degenerate->monodisperse {ndeg}')
    if not a.no_field:
        fok, probs = check_field(a.run_dir, cfg)
        for pr in probs: print('  FAIL(field):', pr)
        ok = ok and fok
    print('VERDICT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
