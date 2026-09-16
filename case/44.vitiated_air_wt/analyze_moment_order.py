#!/usr/bin/env python3
"""dual-time の 3 水準時間次数ゲート (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20; codex plan-8 M5)。

各 run の最終 res_<整数>.h5 (境界出力 res_wall_* 等は除く; 同一物理時刻を /CHECKPOINT/totalTime で照合) を読み、非重み付き L2 差で
  e(dt) = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt)),  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)
を量ごとに出し、次を検査して PASS/FAIL (exit 0/1) を返す:
  - 刻み比が 2 (/CHECKPOINT/dt), 全 run の totalTime が一致 (相対 1e-6), 必要成分が揃い全て有限
  - 次数: BDF2 は [--order-lo, --order-hi] (既定 1.7–2.3), BDF1 は 0.7–1.3 (--bdf 1)
  - sub-iter 比 ≤ --subiter-ratio (既定 0.1; --nsub 指定時)
  - 各 run の各物理 step で sub-iter 残差低下 (inner_iter 行の 初回/最終 rms) の最小値 ≥ --subiter-decades (既定 2.0) を全列で
  - 実現可能性: 最終場で Q1² ≤ Q0 Q2 (1+1e-6), Q2² ≤ Q1 Q3 (1+1e-6) (Q3 = g/((4/3)π ρ_l), --rho-l 既定 1000 kg/m³ [H2O]; N2 は 810), 違反ノード 0 (Q0>0, g>0 のみ)
使い方:
  analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 [--nsub RUN_dt_nsub2x] [--bdf 2] [--fields g_0,Q0_0,Q1_0,Q2_0,T,ro,P]
"""
import argparse, csv, glob, math, os, re, sys
import numpy as np
import h5py


def last_res(run):
    fs = [f for f in glob.glob(os.path.join(run, 'res_*.h5')) if re.search(r'res_(\d+)\.h5$', f)]
    if not fs:
        sys.exit(f'no res_<int>.h5 in {run}')
    return max(fs, key=lambda f: int(re.search(r'res_(\d+)\.h5$', f).group(1)))


def load(run, fields):
    f = last_res(run)
    with h5py.File(f, 'r') as h:
        d = {}
        for k in fields:
            if k not in h['VALUE']:
                sys.exit(f'{f}: field {k} missing')
            d[k] = np.asarray(h['VALUE'][k], dtype=np.float64)
        ck = dict(h['CHECKPOINT'].attrs) if 'CHECKPOINT' in h else {}
        g = np.asarray(h['VALUE']['g_0'], dtype=np.float64) if 'g_0' in h['VALUE'] else None
        Qs = [np.asarray(h['VALUE'][k], dtype=np.float64) for k in ('Q0_0', 'Q1_0', 'Q2_0')] if 'Q0_0' in h['VALUE'] else None
    return f, d, ck, g, Qs


def subiter_decades(run):
    """residual_history.csv の各物理 step について inner_iter 行の 初回/最終 rms の低下桁 (列ごとの最小値)。"""
    p = os.path.join(run, 'residual_history.csv')
    if not os.path.exists(p):
        return None
    first, lastv, cols = {}, {}, None
    with open(p) as f:
        r = csv.DictReader(f)
        cols = [c for c in r.fieldnames if c.startswith('rms_') and not c.startswith('rms_dq_')]
        for row in r:
            if row.get('phase') != 'inner_iter':
                continue
            st = int(row['step'])
            vals = {c: float(row[c]) for c in cols if row.get(c) not in (None, '')}
            if st not in first: first[st] = vals
            lastv[st] = vals
    mins = {}
    for c in cols:
        decs = []
        for st in first:
            a, b = first[st].get(c), lastv[st].get(c)
            if a is None or b is None or a <= 0.0:
                continue
            decs.append(math.log10(a / max(b, 1e-300)))
        if decs:
            mins[c] = (min(decs), float(np.median(decs)))
    return mins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--levels', nargs=3, required=True, help='dt 2dt, dt, dt/2 の run (同一物理時刻)')
    ap.add_argument('--nsub', default=None, help='dt 水準で nSub を倍にした run')
    ap.add_argument('--fields', default='g_0,Q0_0,Q1_0,Q2_0,T,ro,P')
    ap.add_argument('--bdf', type=int, default=2)
    ap.add_argument('--order-lo', type=float, default=None); ap.add_argument('--order-hi', type=float, default=None)
    ap.add_argument('--subiter-ratio', type=float, default=0.1)
    ap.add_argument('--subiter-decades', type=float, default=2.0)
    ap.add_argument('--rho-l', type=float, default=1000.0)
    a = ap.parse_args()
    lo = a.order_lo if a.order_lo is not None else (1.7 if a.bdf == 2 else 0.7)
    hi = a.order_hi if a.order_hi is not None else (2.3 if a.bdf == 2 else 1.3)
    fields = a.fields.split(',')
    runs = list(a.levels) + ([a.nsub] if a.nsub else [])
    data = [load(r, fields) for r in runs]
    bad = []
    for r, (f, _, ck, _, _) in zip(runs, data):
        print(f'  {r}: {os.path.basename(f)} totalTime={ck.get("totalTime")} dt={ck.get("dt")}')
    # 時刻・刻み比
    tts = [d[2].get('totalTime') for d in data]; dts = [d[2].get('dt') for d in data]
    if any(t is None for t in tts):
        bad.append('missing /CHECKPOINT totalTime in some run')
    else:
        t0 = tts[0]
        if any(abs(t - t0) > 1e-6*abs(t0) for t in tts): bad.append(f'totalTime differs across runs: {tts}')
    if any(d is None for d in dts[:3]):
        bad.append('missing /CHECKPOINT dt')
    else:
        r1, r2 = dts[0]/dts[1], dts[1]/dts[2]
        if abs(r1 - 2.0) > 1e-3 or abs(r2 - 2.0) > 1e-3: bad.append(f'dt ratios not 2: {r1:.4f}, {r2:.4f}')
        if a.nsub and abs(dts[3]/dts[1] - 1.0) > 1e-6: bad.append('nsub run dt differs from the middle level')
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s}  max|q|   verdict")
    for k in fields:
        q0, q1, q2 = data[0][1][k], data[1][1][k], data[2][1][k]
        if not (np.isfinite(q0).all() and np.isfinite(q1).all() and np.isfinite(q2).all()):
            print(f'{k:8s} non-finite'); bad.append(f'{k}: non-finite'); continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2)
        order = math.log2(e0/e1) if e1 > 0 and e0 > 0 else float('nan')
        line = f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f}'
        v = []
        if not (lo <= order <= hi): v.append(f'order outside [{lo},{hi}]')
        if a.nsub:
            es = np.linalg.norm(data[3][1][k] - q1)
            ratio = es/e1 if e1 > 0 else float('inf')
            line += f' {es:11.4e} {ratio:7.3f}'
            if not (ratio <= a.subiter_ratio): v.append(f'sub-iter ratio {ratio:.3f} > {a.subiter_ratio}')
        else:
            line += ' '*20
        line += f'  {np.abs(q1).max():.3e}  ' + ('ok' if not v else 'FAIL: ' + '; '.join(v))
        print(line)
        bad.extend(f'{k}: {x}' for x in v)
    # sub-iter 低下 (各 run)
    for r in runs:
        m = subiter_decades(r)
        if m is None:
            bad.append(f'{r}: no residual_history.csv'); continue
        worst = min(m.items(), key=lambda kv: kv[1][0]) if m else None
        if worst:
            print(f'  sub-iter drop {os.path.basename(r)}: min over steps/columns {worst[1][0]:.2f} dec ({worst[0]}), medians ' +
                  ', '.join(f'{c} {v[1]:.2f}' for c, v in m.items()))
            if worst[1][0] < a.subiter_decades: bad.append(f'{r}: sub-iter drop {worst[1][0]:.2f} dec ({worst[0]}) < {a.subiter_decades}')
    # 実現可能性 (最終場)
    for r, (f, _, _, g, Qs) in zip(runs, data):
        if g is None or Qs is None:
            continue
        Q0, Q1, Q2 = Qs; m = (Q0 > 0) & (g > 0)
        Q3 = g/((4.0/3.0)*math.pi*a.rho_l)
        v1 = int(np.sum(Q1[m]**2 > Q0[m]*Q2[m]*(1 + 1e-6) + 1e-300)); v2 = int(np.sum(Q2[m]**2 > Q1[m]*Q3[m]*(1 + 1e-6) + 1e-300))
        print(f'  realizability {os.path.basename(r)}: wet nodes {int(m.sum())}, Q1^2<=Q0Q2 violations {v1}, Q2^2<=Q1Q3 violations {v2}')
        if v1 or v2: bad.append(f'{r}: realizability violations {v1}/{v2}')
    if bad:
        print('VERDICT: FAIL'); [print('  -', b) for b in bad]; sys.exit(1)
    print('VERDICT: PASS'); sys.exit(0)


if __name__ == '__main__':
    main()
