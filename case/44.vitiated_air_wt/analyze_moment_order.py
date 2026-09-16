#!/usr/bin/env python3
"""dual-time の 3 水準時間次数ゲート (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20; codex plan-8 M5, plan-9 M3/M4, plan-10 M2–M4/m1, plan-11 M1/M2/M4/M5)。

各 run の終了 step の場 res_<nStepOuter>.h5 を読み、非重み付き L2 差で
  e(dt) = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt)),  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)
を量ごとに出し、次を全部検査して PASS/FAIL (exit 0/1) を返す (実体は solver_density_cuda/tools/passive_gate_common.py):
  - config: solver と同じ既定値で正規化した実効設定。3 水準 + nSub run は dt / nStepOuter / nSubIterDualTime / outStepInterval / monitorInterval / valueFileName
    以外の solverConfig が同一で、bcondConfig・メッシュ・IC (md5) も同一。bdfOrder == --bdf。刻み比 2 (float32 実効 dt)、nSub run は nSub が 2 倍で dt/nStepOuter は同一
  - 名目終了時刻 nStepOuter × dt_eff (dt_eff = float32(dt): solver は dt を float32 で持つ) が一致; checkpoint の totalTime/dt がそれと一致 (1e-9 / 1e-12)
  - --expect-fct: config で FCT が有効かつ各 run の forge_run.log に `[passiveFct] active`
  - 評価量は config の全保存量 (ro,roUx,roUy,roUz,roe,roY*,受動種) が既定・必須; --fields で必須集合を覆わなければ PARTIAL (exit 3)
  - 必要成分が揃い全て有限; 次数: BDF2 [1.7, 2.3] / BDF1 [0.7, 1.3] (--expect-fct では下限 --fct-order-lo [既定 1.3]: 非線形リミッタは極値近傍で L2 次数を下げる); sub-iter 比 ≤ --subiter-ratio (0.1); dt 非依存の非ゼロ量は「次数未実証」で FAIL (恒等 0 だけ情報なしで通す)
  - residual_history.csv: 全物理 step に outer_begin/inner_begin(0)/inner_iter 1..nSub−1/outer_end (outer_end = 最終 inner)、config 由来の必須列 (流れ・SST・化学種・受動種) の存在、全行の全数値が有限、
    初回 0 の列は step 内の全 inner 行が 0 のときだけ受理、各 step の低下 (反復 0/最終反復) の最小値 ≥ --subiter-decades (2.0) を**全列**で
  - 確定場ゲート (0 ≤ roXi/ro ≤ 1、モーメント非負、solver と同じ実現可能性) を全 run で
使い方: analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 --nsub RUN_dt_nsubx2 [--bdf 2] [--expect-fct] [--fields ...]
終了コード: 0 PASS / 1 FAIL / 3 PARTIAL (--fields が必須集合を覆わない、または --noise 比較でノイズ律速の量がある)
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
    ap.add_argument('--noise', default=None, help='dt 水準で nSub を +1 した run (同一バイナリのノイズ床; float32 では凝縮モーメントが ulp の数百倍のノイズを持つ [plan §9 2026-09-17])。'
                    '与えると nSub 倍増差は max(--subiter-ratio × e(dt), 2 × ノイズ) で判定し、e(dt) < 5 × ノイズ の量は次数を判定せず NOISE-LIMITED (PARTIAL 扱い) にする')
    ap.add_argument('--fields', default=None, help='評価量 (省略時は config から全保存量: ro,roUx,roUy,roUz,roe,roY*,受動種)。必須集合を覆わない指定は PARTIAL (exit 3) で正式 PASS にしない')
    ap.add_argument('--bdf', type=int, default=2)
    ap.add_argument('--double', action='store_true', help='倍精度ビルド (flow_float=double) の run: 実効刻みを float32 に丸めない (診断用; 生産は float)')
    ap.add_argument('--expect-fct', action='store_true', help='FCT 有効試験: config で FCT が有効で log に [passiveFct] active があること')
    ap.add_argument('--order-lo', type=float, default=None); ap.add_argument('--order-hi', type=float, default=None)
    ap.add_argument('--fct-order-lo', type=float, default=1.3,
                    help='--expect-fct のときの下限 (既定 1.3)。FCT は非線形リミッタなので極値の近くで L2 の観測次数が下がる '
                         '(case/44 `run_0417`-`0420` vs 対照 `run_0425`-`0428`: 同一設定で passiveFct 0 なら全モーメント 1.99-2.01)。'
                         'スキーム自体の 2 次は **passiveFct 0 の対照系列**を既定閾値で通して示すこと')
    ap.add_argument('--subiter-ratio', type=float, default=0.1)
    ap.add_argument('--subiter-decades', type=float, default=2.0)
    a = ap.parse_args()
    lo_def = a.order_lo if a.order_lo is not None else (1.7 if a.bdf == 2 else 0.7)
    hi = a.order_hi if a.order_hi is not None else (2.3 if a.bdf == 2 else 1.3)
    runs = list(a.levels) + [a.nsub] + ([a.noise] if a.noise else [])
    bad = []
    cfgs = []
    for r in runs:
        try: cfgs.append(load_config(r))
        except Exception as e: print('VERDICT: FAIL'); print(f'  - {r}: config error: {e}'); sys.exit(1)
    # 必須評価量 = config の全保存量 (流れ + 化学種 + 受動種; codex plan-12 M4)。--fields はその部分集合/追加を許すが、必須集合を覆わなければ PARTIAL
    c0 = cfgs[0]
    required_fields = ['ro', 'roUx', 'roUy', 'roUz', 'roe'] + ([f'roY{s}' for s in range(c0['nSpecies'])] if c0['nSpecies'] > 1 else []) + list(c0['passives'])
    # FCT 作動試験の下限緩和は **BDF2 の受動種だけ** (codex result-4 M2: 流れ・化学種や BDF1 まで緩めない)
    relaxed = set(c0['passives']) if (a.expect_fct and a.bdf == 2 and a.order_lo is None) else set()
    def bounds(k):
        return (a.fct_order_lo if k in relaxed else lo_def), hi
    if relaxed:
        print(f'note: FCT 作動試験: 受動種 {sorted(relaxed)} の下限のみ {a.fct_order_lo} (非線形リミッタによる極値近傍の低下を許す; 流れ・化学種は [{lo_def}, {hi}])。'
              f'スキームの 2 次は passiveFct 0 の対照系列を既定閾値で通して示すこと')
    fields = a.fields.split(',') if a.fields else list(required_fields)
    partial = sorted(set(required_fields) - set(fields))
    if a.double:
        for c in cfgs:
            c['dt_eff'] = float(c['dt']); c['nominal_time'] = c['dt_eff']*int(c['nStepOuter'])
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
    if a.noise and (dts[4] != dts[1] or nst[4] != nst[1] or nsb[4] != nsb[1] + 1): bad.append(f'--noise run must have the same dt/nStepOuter as the middle level and nSub+1 (got dt {dts[4]} nStepOuter {nst[4]} nSub {nsb[4]})')
    tnom = [c['nominal_time'] for c in cfgs]
    if any(abs(t - tnom[0]) > 1e-12*abs(tnom[0]) for t in tnom): bad.append(f'nominal end times differ: {tnom}')
    data = [load(r, c, fields) for r, c in zip(runs, cfgs)]
    for r, (f, _, ck), c, t in zip(runs, data, cfgs, tnom):
        print(f'  {r}: {os.path.basename(f)} totalTime={ck.get("totalTime")} dt={ck.get("dt")} nominal {t:.12e} cfg bdf={c["bdfOrder"]} nSub={c["nsub"]} fct={c["passiveFct"]} sfr={c["sfr"]}')
    vol = None
    with h5py.File(data[1][0], 'r') as h:
        if 'volume' in h['VALUE']: vol = np.asarray(h['VALUE']['volume'], dtype=np.float64)
    noise_limited = []
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s} {'noise':>9s}  max|q|   verdict")
    for k in fields:
        q0, q1, q2, q3 = data[0][1][k], data[1][1][k], data[2][1][k], data[3][1][k]
        qn = data[4][1][k] if a.noise else None
        if not all(np.isfinite(q).all() for q in (q0, q1, q2, q3)) or not all(q.shape == q0.shape for q in (q1, q2, q3)) or (qn is not None and (not np.isfinite(qn).all() or qn.shape != q0.shape)):
            print(f'{k:8s} non-finite or shape mismatch'); bad.append(f'{k}: non-finite or shape mismatch'); continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2); es = np.linalg.norm(q3 - q1)
        en = np.linalg.norm(qn - q1) if qn is not None else 0.0
        if e0 == 0.0 and e1 == 0.0 and es == 0.0:
            # dt 非依存 (全水準で同一) は時間次数を実証しない: 全水準で恒等的に 0 の量 (2D の roUz 等) だけ情報なしとして通す (codex plan-12 M4)
            if all(not np.any(q) for q in (q0, q1, q2, q3)):
                print(f'{k:8s} {e0:11.4e} {e1:11.4e} {"zero":>7s} {es:11.4e} {"-":>7s}  {np.abs(q1).max():.3e}  ok (identically zero on all levels: no order information)'); continue
            print(f'{k:8s} {e0:11.4e} {e1:11.4e} {"exact":>7s} {es:11.4e} {"-":>7s}  {np.abs(q1).max():.3e}  FAIL: nonzero field identical on all levels (order not demonstrated)')
            bad.append(f'{k}: dt-independent nonzero field, order not demonstrated'); continue
        order = math.log2(e0/e1) if e1 > 0 and e0 > 0 else float('nan')
        ratio = es/e1 if e1 > 0 else float('inf')
        v = []
        if qn is not None and e1 < 5.0*en:
            # 同一バイナリの nSub+1 差 (ノイズ床) に対して e(dt) が 5 倍未満: 次数も nSub 比も判定不能 (float32 のモーメント; 倍精度で再判定する)
            noise_limited.append(k)
            print(f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f} {es:11.4e} {ratio:7.3f} {en:9.2e}  {np.abs(q1).max():.3e}  NOISE-LIMITED (e(dt) < 5 x noise floor): not judged'); 
        else:
            klo, khi = bounds(k)
            if not (klo <= order <= khi): v.append(f'order outside [{klo},{khi}]')
            if not (ratio <= a.subiter_ratio or (qn is not None and es <= 2.0*en)): v.append(f'sub-iter ratio {ratio:.3f} > {a.subiter_ratio}' + (f' and nSub diff > 2 x noise {en:.2e}' if qn is not None else ''))
            print(f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f} {es:11.4e} {ratio:7.3f} {en:9.2e}  {np.abs(q1).max():.3e}  ' + ('ok' if not v else 'FAIL: ' + '; '.join(v)))
        bad.extend(f'{k}: {x}' for x in v)
        # 診断 (判定には使わない): 体積重み総量の次数、e(dt)² の上位 1 % セルへの集中率とその補集合の L2 次数 (前線 [核生成 onset] の非平滑性の指標)
        d0, d1 = q0 - q1, q1 - q2
        if vol is not None and vol.shape == q1.shape:
            t0, t1, t2 = float((q0*vol).sum()), float((q1*vol).sum()), float((q2*vol).sum())
            otot = math.log2(abs(t0 - t1)/abs(t1 - t2)) if (t1 != t2 and t0 != t1) else float('nan')
        else:
            otot = float('nan')
        n1 = max(1, len(d1)//100); idx = np.argsort(np.abs(d1))[::-1]; m = np.ones(len(d1), bool); m[idx[:n1]] = False
        share = float((d1[~m]**2).sum()/max((d1**2).sum(), 1e-300)); oc = math.log2(np.linalg.norm(d0[m])/np.linalg.norm(d1[m])) if np.linalg.norm(d1[m]) > 0 and np.linalg.norm(d0[m]) > 0 else float('nan')
        ol1 = math.log2(np.abs(d0).sum()/np.abs(d1).sum()) if np.abs(d1).sum() > 0 and np.abs(d0).sum() > 0 else float('nan')
        print(f'{"":8s} diag: L1 order {ol1:.2f}, volume-weighted total order {otot:.2f}, top-1% cells hold {100*share:.0f}% of e(dt)^2, L2 order on the other 99% {oc:.2f}')
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
    if partial or noise_limited:
        print(f'VERDICT: PARTIAL (' + ('; '.join(x for x in [f'fields {fields} do not cover the config-derived required set; missing {partial}' if partial else '', f'noise-limited (not judged): {noise_limited}' if noise_limited else ''] if x)) + ')'); sys.exit(3)
    print('VERDICT: PASS'); sys.exit(0)


if __name__ == '__main__':
    main()
