#!/usr/bin/env python3
"""受動種ゲート (check_passive_budget.py / check_passive_field.py / case/44 analyze_moment_order.py) の end-to-end 試験 (codex plan-11):
合成 run (config + bcond + mesh/IC ダミー + log + csv + res_<n>.h5) を一時ディレクトリに作り、各ツールの main() を終了コードまで通す。
正常系 PASS と、反例 (最終場欠落・古い場のみ・CSV 欠落・列欠落・0→1→0 残差・inner 番号不足・許可外の config 差・非有限/非数値トークン・過去の nonfinite フラグ・実現不能な場・入力欠落・入口プロファイル差・outer_end 不一致・dt 非依存の非ゼロ量) の FAIL、
部分評価 (--fields) の PARTIAL と非 FCT 開放 run の INDETERMINATE (exit 3)、restart (開始時刻 0.1) と bdfOrder 省略の正常系 (codex plan-12)。"""
import os, subprocess, sys, tempfile, textwrap
import numpy as np
import h5py, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, '..', '..', 'tools')
ORDER = os.path.join(HERE, '..', '..', '..', 'case', '44.vitiated_air_wt', 'analyze_moment_order.py')
fails = 0


def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)


def cfg_dict(dt, nstep, nsub, fct=1, sfr=2, solver='SLAU', bdf=2, tracer=True, ncond=0, extra=None):
    c = {'mesh': {'meshFormat': 'hdf5', 'discretization': 'node', 'meshFileName': 'mesh.h5', 'valueFileName': 'ic.h5'}, 'gpu': 1, 'solver': solver,
         'physProp': {'isCompressible': 1, 'tracer': 'exhaust' if tracer else 'none', 'species': ['N2', 'O2']},
         'time': {'unsteady': 1, 'dualTime': 1, 'last': {'control': 0, 'nStepOuter': nstep}, 'deltaT': {'control': 0, 'dt': dt, 'cfl': 2.0, 'cfl_pseudo': 2.0, 'monitorInterval': 10,
                  'speciesFaceReconstruction': sfr, 'passiveScalarScheme': 1, 'passiveFct': fct}, 'outStepInterval': nstep, 'timeIntegration': 11, 'bdfOrder': bdf, 'nSubIterDualTime': nsub},
         'space': {'convMethod': 1, 'limiter': 1}, 'turbulence': {'model': 'none'}}
    if ncond: c['condensation'] = {'condensation': 1, 'nCondSpecies': ncond, 'condModel': 1}
    if extra:
        for k, v in extra.items():
            d = c
            ks = k.split('.')
            for kk in ks[:-1]: d = d.setdefault(kk, {})
            d[ks[-1]] = v
    return c


def make_run(root, name, dt, nstep, nsub, *, fct=1, sfr=2, solver='SLAU', bdf=2, tracer=True, ncond=0, extra=None, res_step=None, csv_ok=True,
             drop_cols=(), zero_bump=False, inner_short=False, nan_token=False, field=None, fct_log=True, bad_field=False, pert=0.0,
             restart_t0=None, broken=False, old_nonf=False, no_mesh=False, bcond='periodic', inlet_profile=None, outer_end_bad=False, omit_bdf=False, dust_field=False, tinyg_field=False, underflow_field=False):
    d = os.path.join(root, name); os.makedirs(d, exist_ok=True)
    c = cfg_dict(dt, nstep, nsub, fct, sfr, solver, bdf, tracer, ncond, extra)
    if omit_bdf: del c['time']['bdfOrder']   # 省略 = solver 既定 2 (実効設定で同一視されること)
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as f: yaml.safe_dump(c, f)
    with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as f:
        if bcond == 'periodic': f.write('a: {physID: 1, kind: periodic, ints: {type: 0, partnerBCID: 2}}\nb: {physID: 2, kind: periodic, ints: {type: 0, partnerBCID: 1}}\n')
        else: f.write('a: {physID: 1, kind: inlet_Pressure, ints: {inletProfile: 1}, floats: {Tt: 300.0, Pt: 1.0e5}}\nb: {physID: 2, kind: outlet_statPress, floats: {Ps: 1.0e5}}\n')
    if inlet_profile is not None:
        with open(os.path.join(d, 'inlet_profile_1.csv'), 'w') as f: f.write(inlet_profile)
    if not no_mesh:
        with open(os.path.join(d, 'mesh.h5'), 'wb') as f: f.write(b'MESHDUMMY')
    if restart_t0 is None:
        with open(os.path.join(d, 'ic.h5'), 'wb') as f: f.write(b'ICDUMMY')
    else:
        with h5py.File(os.path.join(d, 'ic.h5'), 'w') as h:
            ck = h.create_group('CHECKPOINT'); ck.attrs['totalTime'] = restart_t0; ck.attrs['dt'] = float(np.float32(dt)); ck.attrs['nHistoryValid'] = 2; ck.attrs['layout'] = 'x'
    dt_eff = float(np.float32(dt))
    # residual history
    cols = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe', 'rms_roY0', 'rms_roY1'] + (['rms_roXi'] if tracer else []) + [f'rms_{k}_{s}' for s in range(ncond) for k in ('rog', 'roQ2', 'roQ1', 'roQ0')]
    cols = [cc for cc in cols if cc not in drop_cols]
    if csv_ok:
        with open(os.path.join(d, 'residual_history.csv'), 'w') as f:
            f.write('step,inner,phase,' + ','.join(cols) + '\n')
            for st in range(nstep):
                base = [1e-3]*len(cols)
                f.write(f'{st},-1,outer_begin,' + ','.join(f'{v:.6e}' for v in base) + '\n')
                f.write(f'{st},0,inner_begin,' + ','.join(f'{v:.6e}' for v in base) + '\n')
                nin = (nsub - 1) if not inner_short else (nsub - 3)
                for i in range(1, nin + 1):
                    vals = [b*10**(-3.0*i/(nsub-1)) for b in base]
                    if zero_bump:
                        vals[0] = 0.0 if i == 1 else (1.0 if i == 2 else 0.0)
                    f.write(f'{st},{i},inner_iter,' + ','.join(f'{v:.6e}' for v in vals) + '\n')
                lastv = [b*10**(-3.0*nin/(nsub-1)) for b in base]
                if zero_bump: lastv[0] = 0.0
                if outer_end_bad: lastv[0] *= 2.0   # outer_end ≠ 最終 inner_iter
                f.write(f'{st},{nsub},outer_end,' + ','.join(f'{v:.6e}' for v in lastv) + '\n')
    # log
    lines = []
    if restart_t0 is not None:
        lines.append(f'[dual-time] history restored from ic.h5 (/CHECKPOINT: 10 levels, totalTime={restart_t0}, dt_file={dt_eff}, nHistoryValid=2)')
    if fct and sfr >= 2 and solver in ('SLAU', 'SLAU2') and fct_log:
        lines.append('[passiveFct] active: post-step conservative FCT for 1 passive scalars (prelimit 0, sweeps 100, tol 1.0e-06)')
    names = (['roXi'] if tracer else []) + [f'{k}_{s}' for s in range(ncond) for k in ('rog', 'roQ2', 'roQ1', 'roQ0')]
    for nm in names: lines.append(f'[passive] initial total {nm:<8s} 1.000000000000e+00 (root-only, before the first physical step)')
    pslo = 'nan' if nan_token else ('BROKEN' if broken else '0.000e+00')
    if old_nonf:
        for nm in names:
            lines.append(f'[passive] step 10 floorCorr {nm:<8s} per-step(avg 10): lo 0.000e+00 hi 0.000e+00 abs 0.000e+00 | cumulative: lo 0.000000000e+00 hi 0.000000000e+00 abs 0.000000000e+00 | total 1.000000000000e+00 rel(abs/total) 0.000000e+00 | limCorr per-step 0.000e+00 cumulative abs 0.000000000e+00 signed 0.000000000e+00 rel 0.000000e+00 cells 0 thetaMin(interval) 1.0000 initialTotal 1.000000000000e+00')
            lines.append(f'[passive]   fctCorr {nm:<8s} cumulative: dropped antidiffusion 1.000000000e-03 (rel 1.000000e-03) faces 12 prelimited 0.000000000e+00 pinCorr 0.000000000e+00 (rel 0.000000e+00) baseViol 0.000000000e+00 (rel 0.000000e+00) bndFluxSigned 0.000000000000e+00 bndDropped 0.000000000e+00 upperViol 0.000000000e+00 (rel 0.000000e+00) | budget: srcHist 0.000000000000e+00 remSigned 0.000000000000e+00 remAbs 0.000000000e+00 (rel 0.000000e+00) increment 0.000000000000e+00 | qL rel-residual interval-max 1.00e-07 run-max 1.00e-07 (sweeps last 3) HO residual rel interval-max 1.00e-07 run-max 1.00e-07 nonfinite 1')
    for nm in names:
        lines.append(f'[passive] step {nstep} floorCorr {nm:<8s} per-step(avg 10): lo {pslo} hi 0.000e+00 abs 0.000e+00 | cumulative: lo 0.000000000e+00 hi 0.000000000e+00 abs 0.000000000e+00 | total 1.000000000000e+00 rel(abs/total) 0.000000e+00 | limCorr per-step 0.000e+00 cumulative abs 0.000000000e+00 signed 0.000000000e+00 rel 0.000000e+00 cells 0 thetaMin(interval) 1.0000 initialTotal 1.000000000000e+00')
        if fct and sfr >= 2 and solver in ('SLAU', 'SLAU2') and fct_log:
            lines.append(f'[passive]   fctCorr {nm:<8s} cumulative: dropped antidiffusion 1.000000000e-03 (rel 1.000000e-03) faces 12 prelimited 0.000000000e+00 pinCorr 0.000000000e+00 (rel 0.000000e+00) baseViol 0.000000000e+00 (rel 0.000000e+00) bndFluxSigned 0.000000000000e+00 bndDropped 0.000000000e+00 upperViol 0.000000000e+00 (rel 0.000000e+00) | budget: srcHist 0.000000000000e+00 remSigned 0.000000000000e+00 remAbs 0.000000000e+00 (rel 0.000000e+00) increment 0.000000000000e+00 | qL rel-residual interval-max 1.00e-07 run-max 1.00e-07 (sweeps last 3) HO residual rel interval-max 1.00e-07 run-max 1.00e-07 nonfinite 0')
    if ncond:
        lines.append(f'[passive] step {nstep} moment realizability corrections since last log: nearest-point 0, degenerate->monodisperse 0')
        for s in range(ncond):
            lines.append(f'[passive]   clampBudget species {s} cumulative (signed/abs, rel to total): g 0.000000e+00/0.000000e+00 (0.000000e+00) Q0 0.000000e+00/0.000000e+00 (0.000000e+00) Q1 0.000000e+00/0.000000e+00 (0.000000e+00) Q2 0.000000e+00/0.000000e+00 (0.000000e+00)')
    with open(os.path.join(d, 'forge_run.log'), 'w') as f: f.write('\n'.join(lines) + '\n')
    # field
    n = 50
    rs = res_step if res_step is not None else nstep
    with h5py.File(os.path.join(d, f'res_{rs}.h5'), 'w') as h:
        V = h.create_group('VALUE'); ck = h.create_group('CHECKPOINT')
        ck.attrs['totalTime'] = (restart_t0 or 0.0) + dt_eff*nstep; ck.attrs['dt'] = dt_eff; ck.attrs['nHistoryValid'] = 2; ck.attrs['layout'] = 'x'
        fld = field if field is not None else (np.linspace(0.1, 0.9, n) + pert)
        V['ro'] = 1.0 + 0.01*fld; V['T'] = np.full(n, 250.0); V['P'] = np.full(n, 1e5)
        # 全保存量 (次数ゲートの必須評価量): 流れ・化学種は同じ摂動場、roUz は恒等 0 (2D)
        V['roUx'] = 10.0 + fld; V['roUy'] = fld; V['roUz'] = np.zeros(n); V['roe'] = 2e5 + fld; V['roY0'] = 0.7 + fld*0.1; V['roY1'] = 0.3 - fld*0.1
        if tracer: V['roXi'] = fld
        for s in range(ncond):
            g = np.full(n, 1e-3); Q0 = np.full(n, 1e14); rl = 1000.0 - 0.12*(277.0 - 250.0)
            r = np.cbrt(g/((4/3)*np.pi*rl)/Q0)
            if tinyg_field or underflow_field:
                # codex result-4 M1: 倍精度でも判定すること / 半径が underflow するセルを「判定不能」で見逃さないこと
                import math as _m
                T0 = 300.0; rl = max(1000.0 - 0.12*(277.0 - T0), 920.0)
                gg = 1e-60 if tinyg_field else 1e-320
                q3 = gg/((4/3)*_m.pi*rl); Q0 = np.ones(n); r = (q3/1.0)**(1/3) if q3 > 0 else 0.0
                x, y = (2.0, 0.1)
                Q1 = np.full(n, Q0[0]*r*x) if r > 0 else np.full(n, 1e-100)
                Q2 = np.full(n, Q0[0]*r*r*y) if r > 0 else np.full(n, 1e-200)
                for nm, v in ((f'rog_{s}', np.full(n, gg)), (f'roQ0_{s}', Q0), (f'roQ1_{s}', Q1), (f'roQ2_{s}', Q2)):
                    V[nm] = v; V[nm.replace('ro', '', 1) if nm.startswith('rog') else nm[2:]] = v
                continue
            if dust_field:   # codex result-3 M5: 微小な rog に桁違いの Q1/Q2 (絶対閾値で除外してはならない)
                g = np.full(n, 1e-31); Q0 = np.ones(n); Q1 = np.full(n, 1e8)
                V[f'rog_{s}'] = g; V[f'roQ0_{s}'] = Q0; V[f'roQ1_{s}'] = Q1; V[f'roQ2_{s}'] = np.full(n, 1e8)
                V[f'g_{s}'] = g; V[f'Q0_{s}'] = Q0; V[f'Q1_{s}'] = Q1; V[f'Q2_{s}'] = np.full(n, 1e8)
                continue
            if bad_field: Q1 = np.zeros(n)   # 特異不整合 (Q1=0, Q3>0)
            else: Q1 = Q0*r*0.9
            V[f'rog_{s}'] = g; V[f'roQ0_{s}'] = Q0; V[f'roQ1_{s}'] = Q1; V[f'roQ2_{s}'] = Q0*r*r*0.9   # x=0.9,y=0.9: 内部
            V[f'g_{s}'] = g; V[f'Q0_{s}'] = Q0; V[f'Q1_{s}'] = Q1; V[f'Q2_{s}'] = Q0*r*r*0.9
    return d


def run_tool(args):
    p = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


with tempfile.TemporaryDirectory() as td:
    B = os.path.join(TOOLS, 'check_passive_budget.py'); F = os.path.join(TOOLS, 'check_passive_field.py')
    good = make_run(td, 'good', 8e-6, 20, 40)
    rc, out = run_tool([B, good]); check(rc == 0, f'budget: clean synthetic run must PASS\n{out}')
    rc, out = run_tool([F, good]); check(rc == 0, f'field: clean synthetic run must PASS\n{out}')
    rc, out = run_tool([B, make_run(td, 'oldres', 8e-6, 20, 40, res_step=0)]); check(rc != 0, 'budget: only res_0 (no final field) must FAIL')
    rc, out = run_tool([F, make_run(td, 'oldres2', 8e-6, 20, 40, res_step=0)]); check(rc != 0, 'field: only res_0 must FAIL')
    rc, out = run_tool([B, make_run(td, 'nocsv', 8e-6, 20, 40, csv_ok=False)]); check(rc != 0, 'budget: missing csv must FAIL')
    rc, out = run_tool([B, make_run(td, 'nantok', 8e-6, 20, 40, nan_token=True)]); check(rc != 0, 'budget: nan per-step token must FAIL')
    rc, out = run_tool([B, make_run(td, 'nofctlog', 8e-6, 20, 40, fct_log=False)]); check(rc != 0, 'budget: FCT configured but no active line / record must FAIL')
    rc, out = run_tool([B, make_run(td, 'badfield', 8e-6, 20, 40, ncond=1, bad_field=True)]); check(rc != 0, 'budget: singular moment state must FAIL (field gate)')
    rc, out = run_tool([F, make_run(td, 'dustfield', 8e-6, 20, 40, ncond=1, dust_field=True)]); check(rc != 0, 'field: tiny rog with huge Q1/Q2 must FAIL (no absolute exclusion)')
    rc, out = run_tool([F, make_run(td, 'tinyg', 8e-6, 20, 40, ncond=1, tinyg_field=True), '--double']); check(rc != 0, 'field: rog 1e-60 with x=2 (violating) must FAIL in double as well')
    rc, out = run_tool([F, make_run(td, 'uflow', 8e-6, 20, 40, ncond=1, underflow_field=True), '--double']); check(rc != 0, 'field: cells whose radius underflows must FAIL (realizability cannot be established)')
    rc, out = run_tool([B, make_run(td, 'xi2', 8e-6, 20, 40, field=np.full(50, 1.5))]); check(rc != 0, 'budget: roXi/ro = 1.5 must FAIL (field gate)')
    rc, out = run_tool([B, make_run(td, 'keep', 8e-6, 20, 40, solver='KEEP')]); check(rc == 0, f'budget: KEEP (FCT not configured, closed periodic, no source) -> conservative mode PASS\n{out}')
    # order gate: normal series
    L = [make_run(td, f'lvl{i}', dt, 20*2**i, 40, pert=p) for i, (dt, p) in enumerate(((1.6e-5, 4e-4), (8e-6, 1e-4), (4e-6, 2.5e-5)))]
    N = make_run(td, 'nsub', 8e-6, 40, 80, pert=1e-4 + 1e-7)
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', N]); check(rc == 0, f'order: 2nd-order synthetic series (orders 2, ratio 1e-3) must PASS with default thresholds\n{out}')
    rc, out = run_tool([ORDER, '--fields', 'ro,roXi', '--levels', *L, '--nsub', N]); check(rc == 3 and 'PARTIAL' in out, f'order: --fields not covering the required set must be PARTIAL (exit 3)\n{out}')
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', make_run(td, 'nsub_nobdf', 8e-6, 40, 80, pert=1e-4 + 1e-7, omit_bdf=True)]); check(rc == 0, f'order: omitted bdfOrder (solver default 2) must equal explicit 2 (effective config)\n{out}')
    rc, out = run_tool([B, make_run(td, 'rs', 8e-6, 40, 40, pert=1e-4, restart_t0=0.1)]); check(rc == 0 and 'restart (restored)' in out, f'budget: restart run (t0 0.1 from checkpoint + nStepOuter*dt_eff) must PASS\n{out}')
    rc, out = run_tool([B, make_run(td, 'rs0', 8e-6, 40, 40, pert=1e-4, restart_t0=0.1, res_step=None, extra=None)]); check(rc == 0, 'budget: restart run must PASS (second instance)')
    rc, out = run_tool([B, make_run(td, 'broken', 8e-6, 20, 40, broken=True)]); check(rc != 0, 'budget: BROKEN token in the per-step part must FAIL')
    rc, out = run_tool([B, make_run(td, 'oldnonf', 8e-6, 20, 40, old_nonf=True)]); check(rc != 0, 'budget: nonfinite 1 at an earlier record must FAIL (sticky)')
    rc, out = run_tool([B, make_run(td, 'unst', 8e-6, 20, 40, solver='KEEP', bcond='inlet', inlet_profile='y,Tt\n0,300\n1,300\n')]); check(rc == 3 and 'INDETERMINATE' in out, f'budget: non-FCT open run must be INDETERMINATE (exit 3), not PASS\n{out}')
    Lm = [make_run(td, f'nm{i}', dt, 20*2**i, 40, pert=1e-4, no_mesh=True) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Nm = make_run(td, 'nmn', 8e-6, 40, 80, pert=1e-4, no_mesh=True)
    rc, out = run_tool([ORDER, '--levels', *Lm, '--nsub', Nm, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0 and 'missing' in out, f'order: mesh missing in all runs must FAIL (not silently equal)\n{out}')
    Lp = [make_run(td, f'ip{i}', dt, 20*2**i, 40, pert=p, bcond='inlet', inlet_profile='y,Tt\n0,300\n1,300\n') for i, (dt, p) in enumerate(((1.6e-5, 4e-4), (8e-6, 1e-4), (4e-6, 2.5e-5)))]
    Np = make_run(td, 'ipn', 8e-6, 40, 80, pert=1e-4 + 1e-7, bcond='inlet', inlet_profile='y,Tt\n0,300\n1,310\n')
    rc, out = run_tool([ORDER, '--levels', *Lp, '--nsub', Np]); check(rc != 0 and 'inlet_profile_1.csv' in out, f'order: differing inlet profile CSV must FAIL\n{out}')
    Le = [make_run(td, f'oe{i}', dt, 20*2**i, 40, pert=1e-4, outer_end_bad=True) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Ne = make_run(td, 'oen', 8e-6, 40, 80, pert=1e-4, outer_end_bad=True)
    rc, out = run_tool([ORDER, '--levels', *Le, '--nsub', Ne, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0 and 'outer_end' in out, f'order: outer_end != last inner_iter must FAIL\n{out}')
    # 3 水準で同じ場 (次数は nan) → 次数 FAIL するはずだが、まず config/csv/field の合格経路を "--order-lo -inf" 相当で見る: 代わりに次数だけ緩める
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', N, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9', '--expect-fct']); check(rc == 0, f'order: consistent synthetic series must PASS with relaxed thresholds\n{out}')
    Lc = [make_run(td, f'cst{i}', dt, 20*2**i, 40, pert=1e-4) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Ncst = make_run(td, 'cstn', 8e-6, 40, 80, pert=1e-4)
    rc, out = run_tool([ORDER, '--levels', *Lc, '--nsub', Ncst]); check(rc != 0 and 'not demonstrated' in out, f'order: dt-independent nonzero fields must FAIL (order not demonstrated); identically-zero roUz is tolerated\n{out}')
    Lo = [make_run(td, f'o1{i}', dt, 20*2**i, 40, pert=p) for i, (dt, p) in enumerate(((1.6e-5, 2e-4), (8e-6, 1e-4), (4e-6, 5e-5)))]; No = make_run(td, 'o1n', 8e-6, 40, 80, pert=1e-4 + 1e-7)
    rc, out = run_tool([ORDER, '--levels', *Lo, '--nsub', No]); check(rc != 0, 'order: first-order synthetic series must FAIL the BDF2 gate')
    Ld = [make_run(td, f'dc{i}', dt, 20*2**i, 40, drop_cols=('rms_roY1',)) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Nd = make_run(td, 'dcn', 8e-6, 40, 80, drop_cols=('rms_roY1',))
    rc, out = run_tool([ORDER, '--levels', *Ld, '--nsub', Nd, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0, 'order: missing required residual column must FAIL')
    Lz = [make_run(td, f'zb{i}', dt, 20*2**i, 40, zero_bump=True) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Nz = make_run(td, 'zbn', 8e-6, 40, 80, zero_bump=True)
    rc, out = run_tool([ORDER, '--levels', *Lz, '--nsub', Nz, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0, 'order: 0 -> 1 -> 0 residual must FAIL')
    Li = [make_run(td, f'is{i}', dt, 20*2**i, 40, inner_short=True) for i, dt in enumerate((1.6e-5, 8e-6, 4e-6))]; Ni = make_run(td, 'isn', 8e-6, 40, 80, inner_short=True)
    rc, out = run_tool([ORDER, '--levels', *Li, '--nsub', Ni, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0, 'order: inner_iter numbering short of nSub-1 must FAIL')
    Nc = make_run(td, 'nsubcm', 8e-6, 40, 80, extra={'space.convMethod': 0})
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', Nc, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0, 'order: nsub run with convMethod changed must FAIL')
    Nb = make_run(td, 'nsubbc', 8e-6, 40, 80)
    with open(os.path.join(Nb, 'bcondConfig.yaml'), 'a') as f: f.write('# changed\n')
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', Nb, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9']); check(rc != 0, 'order: nsub run with different bcondConfig must FAIL')
    Nk = make_run(td, 'nsubkeep', 8e-6, 40, 80, solver='KEEP')
    rc, out = run_tool([ORDER, '--levels', *L, '--nsub', Nk, '--order-lo=-1e9', '--order-hi=1e9', '--subiter-ratio=1e9', '--expect-fct']); check(rc != 0, 'order: --expect-fct with a KEEP run must FAIL')
print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
