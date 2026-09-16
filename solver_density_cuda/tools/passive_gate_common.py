#!/usr/bin/env python3
"""受動種ゲート共通部 (plan species-passive-scalar-unification §4.7 v9 / §6-2 / §6-6; codex plan-10 M1–M4)。
  - config から必須成分 (受動種名・化学種数・凝縮種数)・FCT 作動条件・終了 step・刻み・名目終了時刻を確定する
  - 確定場 (最終 res_<int>.h5) の有界性 (0 ≤ roXi/ro ≤ 1)・モーメントの非負・実現可能性 (solver と同じ ρ_l(T), 無次元 (x,y), 退化条件) を保存量から判定する
check_passive_budget.py / check_passive_field.py / case/44 analyze_moment_order.py が import する。"""
import glob, math, os, re
import numpy as np


def load_config(run_dir):
    import yaml
    p = os.path.join(run_dir, 'solverConfig.yaml')
    with open(p) as fh:
        y = yaml.safe_load(fh) or {}
    t = y.get('time', {}) or {}; dT = t.get('deltaT', {}) or {}; cd = y.get('condensation', {}) or {}; pp = y.get('physProp', {}) or {}
    last = t.get('last', {}) or {}
    species = pp.get('species', []) or []
    c = dict(
        solver=y.get('solver'), unsteady=int(t.get('unsteady', 0) or 0), dualTime=int(t.get('dualTime', 0) or 0),
        timeIntegration=int(t.get('timeIntegration', 0) or 0), bdfOrder=int(t.get('bdfOrder', 1) or 1), nsub=t.get('nSubIterDualTime'),
        nStepOuter=last.get('nStepOuter'), dt=dT.get('dt'), passiveScalarScheme=int(dT.get('passiveScalarScheme', 1) or 0),
        passiveFct=int(dT.get('passiveFct', 1) if dT.get('passiveFct', 1) is not None else 1), sfr=int(dT.get('speciesFaceReconstruction', 0) or 0),
        tracer=(pp.get('tracer', 'none') or 'none') != 'none', nCond=int(cd.get('nCondSpecies', 0) or 0) if int(cd.get('condensation', 0) or 0) else 0,
        condModel=cd.get('condModel'), nSpecies=len(species), condensation=int(cd.get('condensation', 0) or 0),
    )
    c['fct_configured'] = (c['passiveScalarScheme'] == 1 and c['passiveFct'] == 1 and c['timeIntegration'] == 11 and c['unsteady'] == 1
                           and c['dualTime'] == 1 and c['sfr'] >= 2 and c['solver'] in ('SLAU', 'SLAU2'))
    c['passives'] = (['roXi'] if c['tracer'] else []) + [f'{k}_{s}' for s in range(c['nCond']) for k in ('rog', 'roQ2', 'roQ1', 'roQ0')]
    c['nominal_time'] = (float(c['dt']) * int(c['nStepOuter'])) if (c['dt'] is not None and c['nStepOuter'] is not None) else None
    return c


def last_res(run_dir):
    fs = [f for f in glob.glob(os.path.join(run_dir, 'res_*.h5')) if re.search(r'res_(\d+)\.h5$', f)]
    if not fs:
        return None
    return max(fs, key=lambda f: int(re.search(r'res_(\d+)\.h5$', f).group(1)))


def rho_l(T, model):
    """solver と同じ液密度: H2O = h2o_rho_cond (1000 − 0.12(277 − T), 床 920), N2 = Nowak 式 (n2_rho_cond)。"""
    T = np.asarray(T, dtype=np.float64)
    if model == 'h2o':
        return np.maximum(1000.0 - 0.12*(277.0 - T), 920.0)
    Tc, rhoc = 126.192, 313.3
    tau = np.maximum(1.0 - np.minimum(T, Tc)/Tc, 0.0)
    lnr = 1.48654237*tau**0.3294 - 0.280476066*tau**(4/6) + 0.0894143085*tau**(16/6) - 0.119879866*tau**(35/6)
    return rhoc*np.exp(lnr)


def check_field(run_dir, cfg=None, out=print):
    """最終 res の確定保存量から: 成分の存在・有限性、ρ>0・T 有限、0 ≤ roXi/ro ≤ 1 (1e-6)、モーメント ≥0、実現可能性 (solver の条件)。戻り (ok, problems)。"""
    import h5py
    cfg = cfg or load_config(run_dir)
    f = last_res(run_dir)
    probs = []
    if f is None:
        return False, [f'{run_dir}: no res_<int>.h5']
    with h5py.File(f, 'r') as h:
        V = h['VALUE']
        need = ['ro', 'T'] + cfg['passives']
        for k in need:
            if k not in V: probs.append(f'{os.path.basename(f)}: field {k} missing')
        if probs:
            return False, probs
        ro = np.asarray(V['ro'], dtype=np.float64); T = np.asarray(V['T'], dtype=np.float64)
        if not np.isfinite(ro).all() or (ro <= 0).any(): probs.append('density non-finite or non-positive')
        if not np.isfinite(T).all(): probs.append('temperature non-finite')
        if cfg['tracer']:
            q = np.asarray(V['roXi'], dtype=np.float64)
            if not np.isfinite(q).all(): probs.append('roXi non-finite')
            else:
                xi = q/ro; nlo = int((xi < -1e-6).sum()); nhi = int((xi > 1 + 1e-6).sum())
                out(f'  field {os.path.basename(f)}: roXi/ro in [{xi.min():.3e}, {xi.max():.3e}] (violations <0: {nlo}, >1: {nhi})')
                if nlo or nhi: probs.append(f'tracer bounds violated ({nlo} below 0, {nhi} above 1)')
        model = 'h2o' if cfg.get('condModel') == 1 else 'n2'
        for s in range(cfg['nCond']):
            g, Q2, Q1, Q0 = (np.asarray(V[f'{k}_{s}'], dtype=np.float64) for k in ('rog', 'roQ2', 'roQ1', 'roQ0'))
            if not all(np.isfinite(a).all() for a in (g, Q0, Q1, Q2)): probs.append(f'species {s}: moments non-finite'); continue
            neg = int((g < 0).sum() + (Q0 < 0).sum() + (Q1 < 0).sum() + (Q2 < 0).sum())
            wet = (Q0 > 0) & (g > 0)
            rl = rho_l(T, model); Q3 = g/((4.0/3.0)*math.pi*rl)
            with np.errstate(all='ignore'):
                rr = np.where(wet, np.cbrt(np.where(wet, Q3/np.where(wet, Q0, 1.0), 0.0)), 0.0)
                okr = wet & np.isfinite(rr) & (rr > 0) & (rr < 1e300)
                x = np.where(okr, Q1/(Q0*rr), 0.0); y = np.where(okr, Q2/(Q0*rr*rr), 0.0)
                eps = 1e-6
                sing = okr & ((x <= 1e-30) | (y <= 1e-30))                       # solver の退化条件 (無次元)
                viol = okr & ~sing & ((x > 1 + eps) | (y < x*x*(1 - eps)) | (y*y > x*(1 + eps)))
            nv, ns = int(viol.sum()), int(sing.sum())
            out(f'  field {os.path.basename(f)} species {s} ({model}, conserved): wet {int(wet.sum())}, inequality violations {nv}, singular (x or y <= 1e-30) {ns}, negative {neg}')
            if nv or ns or neg: probs.append(f'species {s}: realizability viol {nv} singular {ns} negative {neg}')
    return (not probs), probs
