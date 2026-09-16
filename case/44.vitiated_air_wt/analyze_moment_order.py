#!/usr/bin/env python3
"""dual-time の 3 水準時間次数ゲート (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20; codex plan-8 M5, plan-9 M3/M4, plan-10 M2–M4/m1, plan-11 M1/M2/M4/M5)。

各 run の終了 step の場 res_<nStepOuter>.h5 を読み、非重み付き L2 差で
  e(dt) = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt)),  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)
を量ごとに出し、次を全部検査して PASS/FAIL (exit 0/1) を返す (実体は solver_density_cuda/tools/passive_gate_common.py):
  - config: solver と同じ既定値で正規化した実効設定。3 水準 + nSub run は dt / nStepOuter / nSubIterDualTime / outStepInterval / monitorInterval / valueFileName
    以外の solverConfig が同一で、bcondConfig・メッシュ・IC (md5) も同一。bdfOrder == --bdf。刻み比 2 (float32 実効 dt)、nSub run は nSub が 2 倍で dt/nStepOuter は同一
  - 名目終了時刻 nStepOuter × dt_eff (dt_eff = float32(dt): solver は dt を float32 で持つ) が一致; checkpoint の totalTime/dt がそれと一致 (1e-9 / 1e-12)
  - --expect-fct: config で FCT が有効かつ各 run の forge_run.log に `[passiveFct] active`
  - 必要成分が揃い全て有限; 次数: BDF2 [1.7, 2.3] / BDF1 [0.7, 1.3]; sub-iter 比 ≤ --subiter-ratio (0.1)
  - residual_history.csv: 全物理 step に outer_begin/outer_end と inner_iter 1..nSub−1、config 由来の必須列 (流れ・SST・化学種・受動種) の存在、全行の全数値が有限、
    初回 0 の列は step 内の全 inner 行が 0 のときだけ受理、各 step の低下 (初回/最終 inner) の最小値 ≥ --subiter-decades (2.0) を**全列**で
  - 確定場ゲート (0 ≤ roXi/ro ≤ 1、モーメント非負、solver と同じ実現可能性) を全 run で
使い方: analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 --nsub RUN_dt_nsubx2 [--bdf 2] [--expect-fct] [--fields ...]
"""
import argparse, math, os, sys
import numpy as np
import h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'solver_density_cuda', 'tools'))
from passive_gate_common import load_config, final_res, check_field, config_diff, residual_history_check  # noqa: E402


def load(run, cfg, fields):
    f = final_res(run, cfg)
    if f is None:
        sys.exit(f'{run}: final field res_{cfg.get("nStepOuter")}.h5 missing')
    with h5py.File(f, 'r') as h:
        V = h['VALUE']
        d = {}
        for k in fields:
            if k not in V:
                sys.exit(f'{f}: field {k} missing')
            d[k] = np.asarray(V[k], dtype=np.float64)
        ck = dict(h['CHECKPOINT'].attrs) if 'CHECKPOINT' in h else {}
    return f, d, ck


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
        except Exception as e: print('VERDICT: FAIL'); print(f'  - {r}: config error: {e}'); sys.exit(1)
    for r, c in zip(runs, cfgs):
        for k in ('dt', 'nStepOuter', 'nsub', 'bdfOrder', 'unsteady', 'dualTime', 'timeIntegration', 'outStepInterval'):
            v = c.get(k)
            if v is None or not isinstance(v, (int, float)) or not math.isfinite(float(v)) or (k not in ('unsteady', 'dualTime') and float(v) <= 0):
                bad.append(f'{r}: config {k} missing/invalid ({v})')
        if c.get('bdfOrder') != a.bdf: bad.append(f'{r}: bdfOrder {c.get("bdfOrder")} != --bdf {a.bdf}')
        if not (c['unsteady'] == 1 and c['dualTime'] == 1 and c['timeIntegration'] == 11): bad.append(f'{r}: not a dual-time run')
        if a.expect_fct:
            if not c['fct_configured']: bad.append(f'{r}: --expect-fct but FCT not configured (scheme {c["passiveScalarScheme"]}, passiveFct {c["passiveFct"]}, sfr {c["sfr"]}, solver {c["solver"]})')
            lg = os.path.join(r, 'forge_run.log')
            act = os.path.exists(lg) and any(l.startswith('[passiveFct] active') for l in open(lg, errors='replace'))
            if not act: bad.append(f'{r}: --expect-fct but no [passiveFct] active line in forge_run.log')
    if bad:
        print('VERDICT: FAIL'); [print('  -', b) for b in bad]; sys.exit(1)
    # run 同士の同一性 (許可された差分以外)
    for r in runs[1:]:
        for d in config_diff(runs[0], r): bad.append(f'{r} vs {runs[0]}: {d}')
    dts = [c['dt_eff'] for c in cfgs]; nst = [int(c['nStepOuter']) for c in cfgs]; nsb = [int(c['nsub']) for c in cfgs]
    if abs(dts[0]/dts[1] - 2.0) > 1e-12 or abs(dts[1]/dts[2] - 2.0) > 1e-12: bad.append(f'effective dt ratios not 2: {dts}')
    if dts[3] != dts[1] or nst[3] != nst[1]: bad.append('nsub run must have the same dt and nStepOuter as the middle level')
    if nsb[3] != 2*nsb[1] or len({nsb[0], nsb[1], nsb[2]}) != 1: bad.append(f'nSubIterDualTime must be equal on the 3 levels and doubled on the nsub run: {nsb}')
    tnom = [c['nominal_time'] for c in cfgs]
    if any(abs(t - tnom[0]) > 1e-12*abs(tnom[0]) for t in tnom): bad.append(f'nominal end times differ: {tnom}')
    data = [load(r, c, fields) for r, c in zip(runs, cfgs)]
    for r, (f, _, ck), c, t in zip(runs, data, cfgs, tnom):
        print(f'  {r}: {os.path.basename(f)} totalTime={ck.get("totalTime")} dt={ck.get("dt")} nominal {t:.12e} cfg bdf={c["bdfOrder"]} nSub={c["nsub"]} fct={c["passiveFct"]} sfr={c["sfr"]}')
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s}  max|q|   verdict")
    for k in fields:
        q0, q1, q2, q3 = data[0][1][k], data[1][1][k], data[2][1][k], data[3][1][k]
        if not all(np.isfinite(q).all() for q in (q0, q1, q2, q3)) or not all(q.shape == q0.shape for q in (q1, q2, q3)):
            print(f'{k:8s} non-finite or shape mismatch'); bad.append(f'{k}: non-finite or shape mismatch'); continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2); es = np.linalg.norm(q3 - q1)
        if e0 == 0.0 and e1 == 0.0 and es == 0.0:
            print(f'{k:8s} {e0:11.4e} {e1:11.4e} {"exact":>7s} {es:11.4e} {"-":>7s}  {np.abs(q1).max():.3e}  ok (no dt dependence)'); continue
        order = math.log2(e0/e1) if e1 > 0 and e0 > 0 else float('nan')
        ratio = es/e1 if e1 > 0 else float('inf')
        v = []
        if not (lo <= order <= hi): v.append(f'order outside [{lo},{hi}]')
        if not (ratio <= a.subiter_ratio): v.append(f'sub-iter ratio {ratio:.3f} > {a.subiter_ratio}')
        print(f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f} {es:11.4e} {ratio:7.3f}  {np.abs(q1).max():.3e}  ' + ('ok' if not v else 'FAIL: ' + '; '.join(v)))
        bad.extend(f'{k}: {x}' for x in v)
    for r, c in zip(runs, cfgs):
        mins, probs = residual_history_check(r, c, a.subiter_decades)
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
