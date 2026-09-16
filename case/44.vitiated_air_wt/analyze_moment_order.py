#!/usr/bin/env python3
"""dual-time の 3 水準時間次数ゲート (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20; codex plan-8 M5, plan-9 M3/M4, plan-10 M2–M4/m1)。

各 run の最終 res_<整数>.h5 を読み、非重み付き L2 差で
  e(dt) = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt)),  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)
を量ごとに出し、次を全部検査して PASS/FAIL (exit 0/1) を返す:
  - config (solverConfig.yaml): 必須キー (dt, nStepOuter, nSubIterDualTime, bdfOrder) の存在・型・範囲、bdfOrder == --bdf、nSub run は nSub が 2 倍で他は同一、
    solver/unsteady/dualTime/timeIntegration/passiveScalarScheme/speciesFaceReconstruction/passiveFct が 4 run で同一、刻み比 2 (config)
  - 名目終了時刻 dt×nStepOuter (double) が 4 run で一致 (相対 1e-12); checkpoint の totalTime はそれと 1e-9 相対で一致、checkpoint dt == config dt
  - --expect-fct: config で FCT が有効 (scheme 1, passiveFct 1, dual-time, SFR≥2, SLAU) かつ各 run の forge_run.log に `[passiveFct] active`
  - 必要成分が揃い全て有限; 次数: BDF2 [1.7, 2.3] / BDF1 [0.7, 1.3]; sub-iter 比 ≤ --subiter-ratio (0.1)
  - residual_history.csv: 期待する全物理 step (0..nStepOuter−1) に outer 行と inner_iter 行、全輸送列 (rms_* から rms_dq_* を除く全列) の全 inner 行が有限、
    初回 0 の列が後で非ゼロなら拒否 (全反復ゼロは可)、各 step の低下 (初回/最終) の最小値 ≥ --subiter-decades (2.0) を**全列**で
  - 確定場 (保存量) の有界性・非負・実現可能性 (solver と同じ ρ_l(T)・無次元・退化条件; passive_gate_common.check_field) を全 run で
使い方: analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 --nsub RUN_dt_nsubx2 [--bdf 2] [--expect-fct] [--fields ...]
"""
import argparse, csv, math, os, sys
import numpy as np
import h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'solver_density_cuda', 'tools'))
from passive_gate_common import load_config, last_res, check_field  # noqa: E402


def load(run, fields):
    f = last_res(run)
    if f is None:
        sys.exit(f'no res_<int>.h5 in {run}')
    with h5py.File(f, 'r') as h:
        V = h['VALUE']
        d = {}
        for k in fields:
            if k not in V:
                sys.exit(f'{f}: field {k} missing')
            d[k] = np.asarray(V[k], dtype=np.float64)
        ck = dict(h['CHECKPOINT'].attrs) if 'CHECKPOINT' in h else {}
    return f, d, ck


def subiter_check(run, nstep, min_dec):
    """全物理 step に outer 行と inner_iter 行、全輸送列の全 inner 行が有限、初回 0 → 後で非ゼロは拒否、低下の最小 (全列)。"""
    p = os.path.join(run, 'residual_history.csv')
    probs = []
    if not os.path.exists(p):
        return {}, [f'{run}: no residual_history.csv']
    outer, first, lastv, cols = set(), {}, {}, None
    with open(p) as f:
        r = csv.DictReader(f)
        cols = [c for c in (r.fieldnames or []) if c.startswith('rms_') and not c.startswith('rms_dq_')]
        for row in r:
            try: st = int(row['step'])
            except (KeyError, ValueError): probs.append(f'{run}: malformed step in csv'); break
            ph = row.get('phase', '')
            if ph.startswith('outer'):
                outer.add(st); continue
            if ph != 'inner_iter':
                continue
            vals = {}
            for c in cols:
                try: v = float(row[c])
                except (KeyError, ValueError, TypeError): v = float('nan')
                if not math.isfinite(v): probs.append(f'{run}: non-finite sub-iter residual {c} at step {st}')
                vals[c] = v
            if st not in first: first[st] = vals
            lastv[st] = vals
    expected = set(range(nstep)) if nstep else set()
    mo = sorted(expected - outer); mi = sorted(expected - set(first))
    if mo: probs.append(f'{run}: {len(mo)} physical steps without outer rows (e.g. {mo[:3]})')
    if mi: probs.append(f'{run}: {len(mi)} physical steps without inner_iter rows (e.g. {mi[:3]})')
    if not cols: probs.append(f'{run}: no rms_ columns')
    mins = {}
    for c in cols:
        decs = []
        for st in sorted(first):
            a, b = first[st].get(c, float('nan')), lastv[st].get(c, float('nan'))
            if not (math.isfinite(a) and math.isfinite(b)):
                break
            if a == 0.0:
                if b != 0.0: probs.append(f'{run}: {c} starts at 0 and becomes nonzero at step {st}'); break
                continue
            decs.append(math.log10(a / max(b, 1e-300)))
        if decs:
            mins[c] = (min(decs), float(np.median(decs)))
            if min(decs) < min_dec: probs.append(f'{run}: sub-iter drop {c} min {min(decs):.2f} dec < {min_dec}')
    # 重複した非有限報告をまとめる
    seen = set(); probs = [x for x in probs if not (x in seen or seen.add(x))]
    return mins, probs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--levels', nargs=3, required=True, help='dt 2dt, dt, dt/2 の run (同一物理時刻)')
    ap.add_argument('--nsub', required=True, help='dt 水準で nSub を倍にした run (必須)')
    ap.add_argument('--fields', default='g_0,Q0_0,Q1_0,Q2_0,T,ro,P')
    ap.add_argument('--bdf', type=int, default=2)
    ap.add_argument('--expect-fct', action='store_true', help='FCT 有効試験: config で FCT が有効で log に [passiveFct] active があること')
    ap.add_argument('--order-lo', type=float, default=None); ap.add_argument('--order-hi', type=float, default=None)
    ap.add_argument('--subiter-ratio', type=float, default=0.1)
    ap.add_argument('--subiter-decades', type=float, default=2.0)
    a = ap.parse_args()
    lo = a.order_lo if a.order_lo is not None else (1.7 if a.bdf == 2 else 0.7)
    hi = a.order_hi if a.order_hi is not None else (2.3 if a.bdf == 2 else 1.3)
    fields = a.fields.split(',')
    runs = list(a.levels) + [a.nsub]
    bad = []
    cfgs = []
    for r in runs:
        try: cfgs.append(load_config(r))
        except Exception as e: sys.exit(f'{r}: cannot read solverConfig.yaml ({e})')
    # config の存在・型・範囲
    for r, c in zip(runs, cfgs):
        for k in ('dt', 'nStepOuter', 'nsub', 'bdfOrder'):
            v = c.get(k)
            if v is None or not isinstance(v, (int, float)) or not math.isfinite(float(v)) or float(v) <= 0: bad.append(f'{r}: config {k} missing/invalid ({v})')
        if c.get('bdfOrder') != a.bdf: bad.append(f'{r}: bdfOrder {c.get("bdfOrder")} != --bdf {a.bdf}')
        if not (c['unsteady'] == 1 and c['dualTime'] == 1 and c['timeIntegration'] == 11): bad.append(f'{r}: not a dual-time run')
        if a.expect_fct:
            if not c['fct_configured']: bad.append(f'{r}: --expect-fct but FCT not configured (scheme {c["passiveScalarScheme"]}, passiveFct {c["passiveFct"]}, sfr {c["sfr"]}, solver {c["solver"]})')
            lg = os.path.join(r, 'forge_run.log')
            act = os.path.exists(lg) and any(l.startswith('[passiveFct] active') for l in open(lg, errors='replace'))
            if not act: bad.append(f'{r}: --expect-fct but no [passiveFct] active line in forge_run.log')
    if bad:
        print('VERDICT: FAIL'); [print('  -', b) for b in bad]; sys.exit(1)
    keys = ('solver', 'unsteady', 'dualTime', 'timeIntegration', 'passiveScalarScheme', 'sfr', 'passiveFct', 'condModel', 'nCond', 'tracer')
    for k in keys:
        if len({c.get(k) for c in cfgs}) != 1: bad.append(f'config {k} differs across runs: {[c.get(k) for c in cfgs]}')
    dts = [float(c['dt']) for c in cfgs]; nst = [int(c['nStepOuter']) for c in cfgs]; nsb = [int(c['nsub']) for c in cfgs]
    if abs(dts[0]/dts[1] - 2.0) > 1e-9 or abs(dts[1]/dts[2] - 2.0) > 1e-9: bad.append(f'config dt ratios not 2: {dts}')
    if dts[3] != dts[1] or nst[3] != nst[1]: bad.append('nsub run must have the same dt and nStepOuter as the middle level')
    if nsb[3] != 2*nsb[1] or len({nsb[0], nsb[1], nsb[2]}) != 1: bad.append(f'nSubIterDualTime must be equal on the 3 levels and doubled on the nsub run: {nsb}')
    tnom = [d*n for d, n in zip(dts, nst)]
    if any(abs(t - tnom[0]) > 1e-12*abs(tnom[0]) for t in tnom): bad.append(f'nominal end times differ: {tnom}')
    data = [load(r, fields) for r in runs]
    for r, (f, _, ck), c, t in zip(runs, data, cfgs, tnom):
        tt = ck.get('totalTime'); dck = ck.get('dt')
        print(f'  {r}: {os.path.basename(f)} totalTime={tt} dt={dck} nominal {t:.12e} cfg bdf={c["bdfOrder"]} nSub={c["nsub"]} fct={c["passiveFct"]} sfr={c["sfr"]}')
        if tt is None or not math.isfinite(float(tt)): bad.append(f'{r}: checkpoint totalTime missing/non-finite')
        elif abs(float(tt) - t) > 1e-9*abs(t): bad.append(f'{r}: checkpoint totalTime {float(tt):.12e} != nominal {t:.12e}')
        if dck is None or not math.isfinite(float(dck)) or abs(float(dck) - float(c['dt'])) > 1e-12*float(c['dt']): bad.append(f'{r}: checkpoint dt {dck} != config dt {c["dt"]}')
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s}  max|q|   verdict")
    for k in fields:
        q0, q1, q2, q3 = data[0][1][k], data[1][1][k], data[2][1][k], data[3][1][k]
        if not all(np.isfinite(q).all() for q in (q0, q1, q2, q3)) or not all(q.shape == q0.shape for q in (q1, q2, q3)):
            print(f'{k:8s} non-finite or shape mismatch'); bad.append(f'{k}: non-finite or shape mismatch'); continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2); es = np.linalg.norm(q3 - q1)
        order = math.log2(e0/e1) if e1 > 0 and e0 > 0 else float('nan')
        ratio = es/e1 if e1 > 0 else float('inf')
        v = []
        if not (lo <= order <= hi): v.append(f'order outside [{lo},{hi}]')
        if not (ratio <= a.subiter_ratio): v.append(f'sub-iter ratio {ratio:.3f} > {a.subiter_ratio}')
        print(f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f} {es:11.4e} {ratio:7.3f}  {np.abs(q1).max():.3e}  ' + ('ok' if not v else 'FAIL: ' + '; '.join(v)))
        bad.extend(f'{k}: {x}' for x in v)
    for r, c in zip(runs, cfgs):
        mins, probs = subiter_check(r, int(c['nStepOuter']), a.subiter_decades)
        if mins:
            worst = min(mins.items(), key=lambda kv: kv[1][0])
            print(f'  sub-iter drop {os.path.basename(r)}: min {worst[1][0]:.2f} dec ({worst[0]}); medians ' + ', '.join(f'{cc} {vv[1]:.2f}' for cc, vv in mins.items()))
        bad.extend(probs)
    for r, c in zip(runs, cfgs):
        ok, probs = check_field(r, c)
        bad.extend(f'{r}: {x}' for x in probs)
    if bad:
        print('VERDICT: FAIL'); [print('  -', b) for b in bad]; sys.exit(1)
    print('VERDICT: PASS'); sys.exit(0)


if __name__ == '__main__':
    main()
