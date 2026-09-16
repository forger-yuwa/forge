#!/usr/bin/env python3
"""dual-time の 3 水準時間次数ゲート (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20; codex plan-8 M5, plan-9 M3/M4)。

各 run の最終 res_<整数>.h5 (境界出力 res_wall_* 等は除く) を読み、非重み付き L2 差で
  e(dt) = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt)),  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)
を量ごとに出し、次を全部検査して PASS/FAIL (exit 0/1) を返す:
  - /CHECKPOINT の totalTime が有限で全 run 一致 (相対 1e-6)、dt の刻み比が 2、solverConfig.yaml の bdfOrder が --bdf、
    nSub run の nSubIterDualTime が中間水準の 2 倍、passiveFct/speciesFaceReconstruction が全 run で同じ
  - 必要成分が揃い全て有限; 次数: BDF2 [1.7, 2.3] / BDF1 [0.7, 1.3]; sub-iter 比 ≤ --subiter-ratio (0.1) (**--nsub 必須**)
  - residual_history.csv: 全物理 step に inner_iter 行があり、必要全残差列の初回/最終 sub-iter が有限、初回 0 の列が後で非ゼロなら拒否、
    各 step の低下 (初回/最終) の最小値 ≥ --subiter-decades (2.0)
  - 実現可能性 (全凝縮種, 確定保存量, solver と同じ ρ_l(T) [H2O: 1000 − 0.12(277 − T), 床 920; N2: Nowak 式] と無次元 (x, y) = (Q1/(Q0 r), Q2/(Q0 r²)),
    r = (Q3/Q0)^{1/3}, ρQ3 = ρg/(4/3 π ρ_l), 許容 x ≤ 1(1+ε), y ≥ x²(1−ε), y² ≤ x(1+ε), ε=1e-6; 特異不整合 (Q0>0, g>0 で Q1 または Q2 が 0), 負値, 非有限も違反)
使い方: analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 --nsub RUN_dt_nsubx2 [--bdf 2] [--fields ...] [--cond-model h2o|n2]
"""
import argparse, csv, glob, math, os, re, sys
import numpy as np
import h5py
import yaml


def last_res(run):
    fs = [f for f in glob.glob(os.path.join(run, 'res_*.h5')) if re.search(r'res_(\d+)\.h5$', f)]
    if not fs:
        sys.exit(f'no res_<int>.h5 in {run}')
    return max(fs, key=lambda f: int(re.search(r'res_(\d+)\.h5$', f).group(1)))


def load(run, fields):
    f = last_res(run)
    with h5py.File(f, 'r') as h:
        V = h['VALUE']
        d = {}
        for k in fields:
            if k not in V:
                sys.exit(f'{f}: field {k} missing')
            d[k] = np.asarray(V[k], dtype=np.float64)
        ck = dict(h['CHECKPOINT'].attrs) if 'CHECKPOINT' in h else {}
        cond = {}
        for s in range(8):
            if f'g_{s}' in V and f'Q0_{s}' in V:
                cond[s] = {k: np.asarray(V[f'{k}_{s}'], dtype=np.float64) for k in ('g', 'Q0', 'Q1', 'Q2')}
        T = np.asarray(V['T'], dtype=np.float64) if 'T' in V else None
    return f, d, ck, cond, T


def cfg_of(run):
    p = os.path.join(run, 'solverConfig.yaml')
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        y = yaml.safe_load(fh)
    t = y.get('time', {}); dT = t.get('deltaT', {}); cd = y.get('condensation', {})
    return dict(bdfOrder=t.get('bdfOrder'), nsub=t.get('nSubIterDualTime'), dt=dT.get('dt'), fct=dT.get('passiveFct', 1), sfr=dT.get('speciesFaceReconstruction', 0),
                condModel=cd.get('condModel'))


def rho_l(T, model):
    if model == 'h2o':
        return np.maximum(1000.0 - 0.12*(277.0 - T), 920.0)
    Tc, rhoc = 126.192, 313.3
    tau = np.maximum(1.0 - np.minimum(T, Tc)/Tc, 0.0)
    lnr = 1.48654237*tau**0.3294 - 0.280476066*tau**(4/6) + 0.0894143085*tau**(16/6) - 0.119879866*tau**(35/6)
    return rhoc*np.exp(lnr)


def subiter_check(run, required_cols, min_dec):
    """全物理 step に inner_iter 行、必要列の初回/最終が有限、初回 0 → 後で非ゼロは拒否、低下の最小。戻り (mins dict, problems list)。"""
    p = os.path.join(run, 'residual_history.csv')
    probs = []
    if not os.path.exists(p):
        return {}, [f'{run}: no residual_history.csv']
    outer, first, lastv = set(), {}, {}
    with open(p) as f:
        r = csv.DictReader(f)
        cols = [c for c in (r.fieldnames or []) if c.startswith('rms_') and not c.startswith('rms_dq_')]
        for row in r:
            st = int(row['step'])
            if row.get('phase', '').startswith('outer'):
                outer.add(st); continue
            if row.get('phase') != 'inner_iter':
                continue
            vals = {}
            for c in cols:
                try: vals[c] = float(row[c])
                except (KeyError, ValueError, TypeError): vals[c] = float('nan')
            if st not in first: first[st] = vals
            lastv[st] = vals
    missing = sorted(outer - set(first))
    if missing: probs.append(f'{run}: {len(missing)} physical steps without inner_iter rows (e.g. {missing[:3]})')
    for c in required_cols:
        if c not in cols: probs.append(f'{run}: residual column {c} missing')
    mins = {}
    for c in cols:
        decs = []
        for st in sorted(first):
            a, b = first[st].get(c, float('nan')), lastv[st].get(c, float('nan'))
            if not (math.isfinite(a) and math.isfinite(b)):
                probs.append(f'{run}: non-finite sub-iter residual {c} at step {st}'); break
            if a == 0.0:
                if b != 0.0: probs.append(f'{run}: {c} starts at 0 and becomes nonzero at step {st}'); break
                continue
            decs.append(math.log10(a / max(b, 1e-300)))
        if decs:
            mins[c] = (min(decs), float(np.median(decs)))
            if c in required_cols and min(decs) < min_dec: probs.append(f'{run}: sub-iter drop {c} min {min(decs):.2f} dec < {min_dec}')
        elif c in required_cols:
            probs.append(f'{run}: no usable sub-iter data for {c}')
    return mins, probs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--levels', nargs=3, required=True, help='dt 2dt, dt, dt/2 の run (同一物理時刻)')
    ap.add_argument('--nsub', required=True, help='dt 水準で nSub を倍にした run (必須)')
    ap.add_argument('--fields', default='g_0,Q0_0,Q1_0,Q2_0,T,ro,P')
    ap.add_argument('--bdf', type=int, default=2)
    ap.add_argument('--order-lo', type=float, default=None); ap.add_argument('--order-hi', type=float, default=None)
    ap.add_argument('--subiter-ratio', type=float, default=0.1)
    ap.add_argument('--subiter-decades', type=float, default=2.0)
    ap.add_argument('--subiter-cols', default='rms_ro,rms_roUx,rms_roUy,rms_roe,rms_rog_0,rms_roQ2_0,rms_roQ1_0,rms_roQ0_0')
    ap.add_argument('--cond-model', choices=['auto', 'h2o', 'n2'], default='auto')
    a = ap.parse_args()
    lo = a.order_lo if a.order_lo is not None else (1.7 if a.bdf == 2 else 0.7)
    hi = a.order_hi if a.order_hi is not None else (2.3 if a.bdf == 2 else 1.3)
    fields = a.fields.split(','); req_cols = a.subiter_cols.split(',')
    runs = list(a.levels) + [a.nsub]
    data = [load(r, fields) for r in runs]
    cfgs = [cfg_of(r) for r in runs]
    bad = []
    for r, (f, _, ck, _, _), c in zip(runs, data, cfgs):
        print(f'  {r}: {os.path.basename(f)} totalTime={ck.get("totalTime")} dt={ck.get("dt")} cfg bdf={c.get("bdfOrder")} nSub={c.get("nsub")} fct={c.get("fct")} sfr={c.get("sfr")}')
    # 時刻・刻み比・設定
    tts = [d[2].get('totalTime') for d in data]; dts = [d[2].get('dt') for d in data]
    if any(t is None or not math.isfinite(float(t)) for t in tts): bad.append(f'totalTime missing/non-finite: {tts}')
    else:
        t0 = float(tts[0])
        if any(abs(float(t) - t0) > 1e-6*abs(t0) for t in tts): bad.append(f'totalTime differs across runs: {tts}')
    if any(d is None for d in dts): bad.append('missing /CHECKPOINT dt')
    else:
        r1, r2 = dts[0]/dts[1], dts[1]/dts[2]
        if abs(r1 - 2.0) > 1e-3 or abs(r2 - 2.0) > 1e-3: bad.append(f'dt ratios not 2: {r1:.4f}, {r2:.4f}')
        if abs(dts[3]/dts[1] - 1.0) > 1e-6: bad.append('nsub run dt differs from the middle level')
    for r, c in zip(runs, cfgs):
        if c.get('bdfOrder') != a.bdf: bad.append(f'{r}: bdfOrder {c.get("bdfOrder")} != --bdf {a.bdf}')
    if cfgs[1].get('nsub') and cfgs[3].get('nsub') and cfgs[3]['nsub'] != 2*cfgs[1]['nsub']: bad.append(f'nsub run nSubIterDualTime {cfgs[3]["nsub"]} != 2 x {cfgs[1]["nsub"]}')
    if len({(c.get('fct'), c.get('sfr')) for c in cfgs}) != 1: bad.append('passiveFct/speciesFaceReconstruction differ across runs')
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s}  max|q|   verdict")
    for k in fields:
        q0, q1, q2, q3 = data[0][1][k], data[1][1][k], data[2][1][k], data[3][1][k]
        if not all(np.isfinite(q).all() for q in (q0, q1, q2, q3)):
            print(f'{k:8s} non-finite'); bad.append(f'{k}: non-finite'); continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2); es = np.linalg.norm(q3 - q1)
        order = math.log2(e0/e1) if e1 > 0 and e0 > 0 else float('nan')
        ratio = es/e1 if e1 > 0 else float('inf')
        v = []
        if not (lo <= order <= hi): v.append(f'order outside [{lo},{hi}]')
        if not (ratio <= a.subiter_ratio): v.append(f'sub-iter ratio {ratio:.3f} > {a.subiter_ratio}')
        print(f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f} {es:11.4e} {ratio:7.3f}  {np.abs(q1).max():.3e}  ' + ('ok' if not v else 'FAIL: ' + '; '.join(v)))
        bad.extend(f'{k}: {x}' for x in v)
    for r in runs:
        mins, probs = subiter_check(r, req_cols, a.subiter_decades)
        if mins:
            worst = min(mins.items(), key=lambda kv: kv[1][0])
            print(f'  sub-iter drop {os.path.basename(r)}: min {worst[1][0]:.2f} dec ({worst[0]}); medians ' + ', '.join(f'{c} {v[1]:.2f}' for c, v in mins.items()))
        bad.extend(probs)
    # 実現可能性 (全種, 確定保存量, solver と同じ ρ_l(T))
    for r, (f, _, _, cond, T), c in zip(runs, data, cfgs):
        model = a.cond_model if a.cond_model != 'auto' else ('h2o' if c.get('condModel') == 1 else 'n2')
        for s, m in cond.items():
            g, Q0, Q1, Q2 = m['g'], m['Q0'], m['Q1'], m['Q2']
            nonfin = int((~np.isfinite(g)).sum() + (~np.isfinite(Q0)).sum() + (~np.isfinite(Q1)).sum() + (~np.isfinite(Q2)).sum())
            neg = int((g < 0).sum() + (Q0 < 0).sum() + (Q1 < 0).sum() + (Q2 < 0).sum())
            wet = (Q0 > 0) & (g > 0)
            sing = int(((Q1 <= 0) | (Q2 <= 0))[wet].sum())
            rl = rho_l(T, model) if T is not None else np.full_like(g, 1000.0)
            Q3 = g/((4.0/3.0)*math.pi*rl)
            with np.errstate(all='ignore'):
                rr = np.cbrt(np.where(wet, Q3/np.where(wet, Q0, 1.0), 0.0))
                ok_r = wet & (rr > 0)
                x = np.where(ok_r, Q1/(Q0*rr), 0.0); y = np.where(ok_r, Q2/(Q0*rr*rr), 0.0)
                eps = 1e-6
                viol = ok_r & ((x > 1 + eps) | (y < x*x*(1 - eps)) | (y*y > x*(1 + eps)))
            nv = int(viol.sum())
            print(f'  realizability {os.path.basename(r)} species {s} ({model}): wet {int(wet.sum())}, inequality violations {nv}, singular (Q1 or Q2 = 0) {sing}, negative {neg}, non-finite {nonfin}')
            if nv or sing or neg or nonfin: bad.append(f'{r} species {s}: realizability viol {nv} singular {sing} negative {neg} non-finite {nonfin}')
    if bad:
        print('VERDICT: FAIL'); [print('  -', b) for b in bad]; sys.exit(1)
    print('VERDICT: PASS'); sys.exit(0)


if __name__ == '__main__':
    main()
