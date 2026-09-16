#!/usr/bin/env python3
"""受動種ゲート共通部 (plan species-passive-scalar-unification §4.7 v9 / §6-2 / §6-6; codex plan-10 M1–M4, plan-11 M1–M6, plan-12 M1/M3/M6/m1)。
  - solverConfig.yaml を solver と同じ既定値で正規化した実効設定にし、必須成分 (受動種・化学種・流れ・SST 列)・FCT 作動条件・終了 step・
    実効刻み (float32 に丸めた dt) と名目終了時刻 nStepOuter×dt_eff を確定する
  - 次数試験の run 同士は「dt / nStepOuter / nSubIterDualTime / outStepInterval / monitorInterval 以外の config が同一、bcond・メッシュ・IC が同一 (md5)」を要求する
  - 終了場は res_<nStepOuter>.h5 に固定 (最大番号のファイルを黙って使わない)
  - 確定場 (保存量) の有界性 (0 ≤ roXi/ro ≤ 1)・モーメント非負・実現可能性 (solver と同じ ρ_l(T), 無次元 (x,y), 退化条件) を判定する
check_passive_budget.py / check_passive_field.py / case/44 analyze_moment_order.py が import する。"""
import hashlib, math, os
import numpy as np


def _f32(x):
    return float(np.float32(x))


def load_config(run_dir):
    import yaml
    p = os.path.join(run_dir, 'solverConfig.yaml')
    with open(p) as fh:
        y = yaml.safe_load(fh) or {}
    t = y.get('time', {}) or {}; dT = t.get('deltaT', {}) or {}; cd = y.get('condensation', {}) or {}; pp = y.get('physProp', {}) or {}
    last = t.get('last', {}) or {}; sp = y.get('space', {}) or {}; tb = y.get('turbulence', {}) or {}
    species = pp.get('species', []) or []

    def req(d, k, sec):
        if k not in d or d[k] is None:
            raise ValueError(f'required key {sec}.{k} missing in solverConfig.yaml')
        return d[k]

    def opt_int(d, k, default):
        v = d.get(k, default)
        if isinstance(v, bool) or not isinstance(v, (int, float)) or int(v) != v:
            raise ValueError(f'key {k} must be an integer (got {v!r})')
        return int(v)

    c = dict(
        solver=y.get('solver'), unsteady=opt_int(t, 'unsteady', None) if 'unsteady' in t else None, dualTime=opt_int(t, 'dualTime', None) if 'dualTime' in t else None,
        timeIntegration=opt_int(t, 'timeIntegration', None) if 'timeIntegration' in t else None,
        bdfOrder=opt_int(t, 'bdfOrder', 2), nsub=opt_int(t, 'nSubIterDualTime', 20), nStepOuter=opt_int(last, 'nStepOuter', None) if 'nStepOuter' in last else None,
        outStepInterval=opt_int(t, 'outStepInterval', None) if 'outStepInterval' in t else None, monitorInterval=opt_int(dT, 'monitorInterval', 1),
        dt=dT.get('dt'), passiveScalarScheme=opt_int(dT, 'passiveScalarScheme', 1), passiveFct=opt_int(dT, 'passiveFct', 1), sfr=opt_int(dT, 'speciesFaceReconstruction', 0),
        convMethod=sp.get('convMethod'), limiter=sp.get('limiter'), turbulence=(tb.get('model') or 'none'),
        tracer=(pp.get('tracer', 'none') or 'none') != 'none', condensation=opt_int(cd, 'condensation', 0),
        condModel=cd.get('condModel'), nSpecies=len(species), speciesDBFile=pp.get('speciesDBFile'),
    )
    c['nCond'] = opt_int(cd, 'nCondSpecies', 0) if c['condensation'] else 0
    c['fct_configured'] = (c['passiveScalarScheme'] == 1 and c['passiveFct'] == 1 and c['timeIntegration'] == 11 and c['unsteady'] == 1
                           and c['dualTime'] == 1 and c['sfr'] >= 2 and c['solver'] in ('SLAU', 'SLAU2'))
    c['passives'] = (['roXi'] if c['tracer'] else []) + [f'{k}_{s}' for s in range(c['nCond']) for k in ('rog', 'roQ2', 'roQ1', 'roQ0')]
    # 必須残差列 (流れ + SST + 化学種 + 受動種)
    cols = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe']
    if str(c['turbulence']).lower() not in ('none', 'laminar', ''):
        cols += ['rms_roK', 'rms_roOmega']
    if c['nSpecies'] > 1:
        cols += [f'rms_roY{s}' for s in range(c['nSpecies'])]
    cols += ['rms_' + p for p in c['passives']]
    c['required_cols'] = cols
    c['dt_eff'] = _f32(float(c['dt'])) if c['dt'] is not None else None   # solver は dt を float32 で持つ
    c['nominal_time'] = (c['dt_eff'] * int(c['nStepOuter'])) if (c['dt_eff'] is not None and c['nStepOuter'] is not None) else None
    c['raw'] = y
    # solver と同じ既定値を埋めた実効設定 (比較用; 省略と明示既定を同一視する — codex plan-12 m1)
    import copy
    eff = copy.deepcopy(y)
    et = eff.setdefault('time', {}) or eff['time']; edT = et.setdefault('deltaT', {}) or et['deltaT']
    et.setdefault('bdfOrder', 2); et.setdefault('nSubIterDualTime', 20)
    edT.setdefault('monitorInterval', 1); edT.setdefault('passiveScalarScheme', 1); edT.setdefault('passiveFct', 1); edT.setdefault('speciesFaceReconstruction', 0)
    ep = eff.setdefault('physProp', {}) or eff['physProp']; ep.setdefault('tracer', 'none')
    ec = eff.setdefault('condensation', {}) or eff['condensation']; ec.setdefault('condensation', 0)
    if ec['condensation']: ec.setdefault('nCondSpecies', 0)
    c['effective'] = eff
    return c


def _md5(path):
    if not os.path.exists(path):
        return None
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


ALLOWED_DIFF = {('time', 'deltaT', 'dt'), ('time', 'last', 'nStepOuter'), ('time', 'nSubIterDualTime'), ('time', 'outStepInterval'),
                ('time', 'outStepStart'), ('time', 'deltaT', 'monitorInterval'), ('mesh', 'valueFileName')}


def _flatten(d, prefix=()):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flatten(v, prefix + (str(k),)))
    else:
        out[prefix] = d
    return out


def _inlet_profiles(run_dir):
    """bcondConfig.yaml で inletProfile: 1 の入口が読む inlet_profile_<physID>.csv の一覧。"""
    import yaml
    p = os.path.join(run_dir, 'bcondConfig.yaml')
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        b = yaml.safe_load(fh) or {}
    out = []
    for v in b.values():
        if isinstance(v, dict) and int((v.get('ints') or {}).get('inletProfile', 0) or 0) == 1 and v.get('physID') is not None:
            out.append(f"inlet_profile_{int(v['physID'])}.csv")
    return sorted(out)


def required_inputs(run_dir, cfg=None):
    """solver が読む入力ファイル (run dir 相対): bcondConfig, メッシュ, IC, 入口プロファイル CSV, 外部物性 DB。戻り {label: path}。"""
    cfg = cfg or load_config(run_dir)
    m = cfg['raw'].get('mesh', {}) or {}
    files = {'bcondConfig.yaml': 'bcondConfig.yaml'}
    for key in ('meshFileName', 'valueFileName'):
        files[f'mesh.{key}'] = m.get(key)
    for f in _inlet_profiles(run_dir):
        files[f] = f
    if cfg.get('speciesDBFile'):
        files['physProp.speciesDBFile'] = cfg['speciesDBFile']
    return files


def config_diff(run_a, run_b):
    """許可された差分以外の実効 solverConfig の違い、および solver が読む入力 (bcondConfig / メッシュ / IC / 入口プロファイル CSV / 物性 DB) の
    欠落と内容 (md5) の違いを列挙する (codex plan-11 M2, plan-12 M3/m1)。"""
    ca, cb = load_config(run_a), load_config(run_b)
    fa, fb = _flatten(ca['effective']), _flatten(cb['effective'])
    diffs = []
    for k in sorted(set(fa) | set(fb)):
        if k in ALLOWED_DIFF:
            continue
        if fa.get(k) != fb.get(k):
            diffs.append(f"solverConfig {'.'.join(k)}: {fa.get(k)!r} vs {fb.get(k)!r}")
    ia, ib = required_inputs(run_a, ca), required_inputs(run_b, cb)
    for label in sorted(set(ia) | set(ib)):
        pa, pb = ia.get(label), ib.get(label)
        if pa is None or pb is None:
            diffs.append(f'{label}: not referenced in both runs ({pa} vs {pb})'); continue
        ha, hb = _md5(os.path.join(run_a, pa)), _md5(os.path.join(run_b, pb))
        if ha is None: diffs.append(f'{label} missing in {run_a} ({pa})')
        if hb is None: diffs.append(f'{label} missing in {run_b} ({pb})')
        if ha is not None and hb is not None and ha != hb:
            diffs.append(f'{label} content differs ({pa} vs {pb})')
    return diffs


def restart_start_time(run_dir, cfg, out=None):
    """restart 継続 run の開始時刻 (codex plan-12 M6)。solver は checkpoint 復元時に totalTime を引き継ぐ:
      '[dual-time] history restored'                              → t0 = IC の /CHECKPOINT/totalTime (BDF2 継続)
      '[dual-time] passive FCT flux-form history missing ...'    → t0 = 同上 (時刻は引き継ぐが全系 BDF1 再開)
      '[dual-time] history NOT restored' / 行なし                 → t0 = 0
    戻り (t0, kind) kind ∈ {'restored', 'bdf1-restart', 'fresh'}。"""
    lg = os.path.join(run_dir, 'forge_run.log')
    kind = 'fresh'
    if os.path.exists(lg):
        for l in open(lg, errors='replace'):
            if l.startswith('[dual-time] history restored'): kind = 'restored'; break
            if l.startswith('[dual-time] passive FCT flux-form history missing'): kind = 'bdf1-restart'; break
            if l.startswith('[dual-time] history NOT restored'): kind = 'fresh'; break
    if kind == 'fresh':
        return 0.0, kind
    ic = (cfg['raw'].get('mesh', {}) or {}).get('valueFileName')
    p = os.path.join(run_dir, ic) if ic else None
    if p and os.path.exists(p):
        import h5py
        with h5py.File(p, 'r') as h:
            if 'CHECKPOINT' in h and 'totalTime' in h['CHECKPOINT'].attrs:
                t0 = float(h['CHECKPOINT'].attrs['totalTime'])
                if out: out(f'  restart ({kind}) from {ic}: start time {t0:.12e}')
                return t0, kind
    raise ValueError(f'{run_dir}: log says history {kind} but the IC {ic} has no /CHECKPOINT/totalTime')


def float_floor_allowance(nstep, kind='signed'):
    """float32 の状態更新の丸め床の許容 (nstep に比例; plan §6-2):
      signed: 総量の保存 (final − initial) と符号付き残り H_rem の積算: max(1e-6, nstep × 2e-8)
              (case/09 seam: run_0156 160 step +1.4e-6 = 0.9e-8/step, run_0157 no-FCT −1.2e-6, run_0159 venkat 100 step +1.9e-6 = 1.9e-8/step; 系統的な正バイアス)
      abs   : |H_rem| の絶対値積算 (毎 step のセルごとの丸め ~ulp/2 を絶対値で足す): max(1e-6, nstep × 6e-8)  (run_0156 7.2e-6/160 step = 4.5e-8/step, run_0159 4.75e-6/100 step)
    閉合 (増分 = −境界 + ソース + 残り) と独立照合は丸め床でなく tol (1e-6) のまま。"""
    return max(1.0e-6, (2.0e-8 if kind == 'signed' else 6.0e-8)*float(nstep))


def use_double_dt(cfg):
    """倍精度ビルド (flow_float=double) の run: 実効刻みを float32 に丸めない (診断用の検証 run; 生産は float)。"""
    cfg['dt_eff'] = float(cfg['dt']) if cfg['dt'] is not None else None
    cfg['nominal_time'] = (cfg['dt_eff']*int(cfg['nStepOuter'])) if (cfg['dt_eff'] is not None and cfg['nStepOuter'] is not None) else None
    return cfg


def final_res(run_dir, cfg=None):
    """終了 step の場 res_<nStepOuter>.h5 (無ければ None)。"""
    cfg = cfg or load_config(run_dir)
    if cfg['nStepOuter'] is None:
        return None
    p = os.path.join(run_dir, f'res_{int(cfg["nStepOuter"])}.h5')
    return p if os.path.exists(p) else None


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
    """終了 step の res から: 成分の存在・有限性、ρ>0・T 有限、0 ≤ roXi/ro ≤ 1 (1e-6)、モーメント ≥0、実現可能性 (solver の条件)、checkpoint の時刻・刻み。戻り (ok, problems)。"""
    import h5py
    cfg = cfg or load_config(run_dir)
    f = final_res(run_dir, cfg)
    probs = []
    if f is None:
        return False, [f'{run_dir}: final field res_{cfg.get("nStepOuter")}.h5 missing (run incomplete or outStepInterval does not hit the last step)']
    with h5py.File(f, 'r') as h:
        V = h['VALUE']
        ck = dict(h['CHECKPOINT'].attrs) if 'CHECKPOINT' in h else {}
        if cfg['unsteady'] == 1 and cfg['dualTime'] == 1:
            tt = ck.get('totalTime'); dck = ck.get('dt')
            try: t0, kind = restart_start_time(run_dir, cfg, out)
            except ValueError as e: return False, [str(e)]
            tnom = t0 + cfg['nominal_time']
            if tt is None or not math.isfinite(float(tt)) or abs(float(tt) - tnom) > 1e-9*abs(tnom):
                probs.append(f'checkpoint totalTime {tt} != start {t0:.12e} ({kind}) + nStepOuter*dt_eff {tnom:.12e}')
            if dck is None or not math.isfinite(float(dck)) or abs(float(dck) - cfg['dt_eff']) > 1e-12*cfg['dt_eff']:
                probs.append(f'checkpoint dt {dck} != effective (float32) config dt {cfg["dt_eff"]:.12e}')
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
            # float32 の正規化数の下限 (1.2e-38) 付近の「塵」は相対 1e-6 の判定ができない (書き戻しで桁が無い) → ρg > 1e-30 の湿りセルだけ判定 (solver の射影対象と同じ意味の状態)
            wet = (Q0 > 0) & (g > 1.0e-30)
            rl = rho_l(T, model); Q3 = g/((4.0/3.0)*math.pi*rl)
            with np.errstate(all='ignore'):
                rr = np.where(wet, np.cbrt(np.where(wet, Q3/np.where(wet, Q0, 1.0), 0.0)), 0.0)
                okr = wet & np.isfinite(rr) & (rr > 0) & (rr < 1e300)
                x = np.where(okr, Q1/(Q0*rr), 0.0); y = np.where(okr, Q2/(Q0*rr*rr), 0.0)
                eps = 1e-6
                sing = okr & ((x <= 1e-30) | (y <= 1e-30))
                viol = okr & ~sing & ((x > 1 + eps) | (y < x*x*(1 - eps)) | (y*y > x*(1 + eps)))
            nv, ns = int(viol.sum()), int(sing.sum())
            out(f'  field {os.path.basename(f)} species {s} ({model}, conserved, rho*g > 1e-30): wet {int(wet.sum())}, inequality violations {nv}, singular (x or y <= 1e-30) {ns}, negative {neg}')
            if nv or ns or neg: probs.append(f'species {s}: realizability viol {nv} singular {ns} negative {neg}')
    return (not probs), probs


def residual_history_check(run_dir, cfg, min_dec, out=print):
    """residual_history.csv (codex plan-11 M1, plan-12 M1): 期待する全物理 step 0..nStepOuter−1 に outer_begin / inner_begin (反復 0) / inner_iter 1..nSub−1 / outer_end、
    重複・範囲外番号なし、outer_end の残差 = 最終 inner_iter の残差 (solver は最終 sub-iter を outer_end に書く)、必須列の存在、全行の全数値が有限、
    反復 0 が 0 の列は step 内の全反復が 0 のときだけ受理、各 step の低下 (反復 0 / 最終反復 nSub−1) の最小 ≥ min_dec を全列で。
    戻り (mins dict, problems)。"""
    import csv
    p = os.path.join(run_dir, 'residual_history.csv')
    probs = []
    if not os.path.exists(p):
        return {}, [f'{run_dir}: residual_history.csv missing']
    nstep = int(cfg['nStepOuter']); nsub = int(cfg['nsub'])
    with open(p) as f:
        r = csv.DictReader(f)
        header = r.fieldnames or []
        cols = [c for c in header if c.startswith('rms_') and not c.startswith('rms_dq_')]
        for c in cfg['required_cols']:
            if c not in cols: probs.append(f'{run_dir}: required residual column {c} missing')
        outer_b, outer_e, inner = set(), {}, {}
        dup = 0
        for row in r:
            try: st = int(row['step']); inn = int(row['inner'])
            except (KeyError, ValueError, TypeError): probs.append(f'{run_dir}: malformed step/inner in csv'); break
            ph = row.get('phase', '')
            vals = {}
            for c in cols:
                try: v = float(row[c])
                except (KeyError, ValueError, TypeError): v = float('nan')
                if not math.isfinite(v): probs.append(f'{run_dir}: non-finite residual {c} at step {st} ({ph} {inn})')
                vals[c] = v
            if ph == 'outer_begin': outer_b.add(st)
            elif ph == 'outer_end':
                if st in outer_e: dup += 1
                outer_e[st] = vals
            elif ph in ('inner_begin', 'inner_iter'):
                if ph == 'inner_begin' and inn != 0: probs.append(f'{run_dir}: inner_begin with iteration {inn} != 0 at step {st}')
                if ph == 'inner_iter' and not (1 <= inn < nsub): probs.append(f'{run_dir}: inner_iter number {inn} out of range 1..{nsub-1} at step {st}')
                d = inner.setdefault(st, {})
                if inn in d: dup += 1
                d[inn] = vals
            else:
                probs.append(f'{run_dir}: unknown phase {ph!r} at step {st}')
    if dup: probs.append(f'{run_dir}: {dup} duplicated iteration rows')
    expected = set(range(nstep))
    if sorted(expected - outer_b): probs.append(f'{run_dir}: steps without outer_begin: {sorted(expected - outer_b)[:5]}')
    if sorted(expected - set(outer_e)): probs.append(f'{run_dir}: steps without outer_end: {sorted(expected - set(outer_e))[:5]}')
    if sorted(expected - set(inner)): probs.append(f'{run_dir}: steps without inner rows: {sorted(expected - set(inner))[:5]}')
    if sorted(set(inner) - expected): probs.append(f'{run_dir}: inner rows for steps beyond nStepOuter: {sorted(set(inner) - expected)[:5]}')
    for st in sorted(inner):
        idx = sorted(inner[st])
        if idx != list(range(nsub)):
            probs.append(f'{run_dir}: step {st} iteration numbering {idx[:3]}..{idx[-1:]} != 0..{nsub-1} (inner_begin + inner_iter 1..{nsub-1})'); break
        if st in outer_e and (nsub - 1) in inner[st]:
            lastv = inner[st][nsub - 1]
            for c in cols:
                a, b = lastv.get(c, float('nan')), outer_e[st].get(c, float('nan'))
                if math.isfinite(a) and math.isfinite(b) and abs(a - b) > 1e-9*max(abs(a), abs(b), 1e-300):
                    probs.append(f'{run_dir}: outer_end residual {c} at step {st} ({b:.6e}) != last inner_iter ({a:.6e})'); break
    mins = {}
    for c in cols:
        decs = []
        for st in sorted(inner):
            rows = [inner[st][i] for i in sorted(inner[st])]
            seq = [rw.get(c, float('nan')) for rw in rows]
            if any(not math.isfinite(v) for v in seq): break
            if seq[0] == 0.0:
                if any(v != 0.0 for v in seq): probs.append(f'{run_dir}: {c} is 0 at iteration 0 but nonzero later at step {st}'); break
                continue
            decs.append(math.log10(seq[0]/max(seq[-1], 1e-300)))
        if decs:
            mins[c] = (min(decs), float(np.median(decs)))
            if min(decs) < min_dec: probs.append(f'{run_dir}: sub-iter drop {c} min {min(decs):.2f} dec < {min_dec}')
    seen = set(); probs = [x for x in probs if not (x in seen or seen.add(x))]
    return mins, probs
