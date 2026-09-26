#!/usr/bin/env python3
"""起点の残差床に対する継続 run の収束判定 (plateau 許容 + 床比)。plan gradient-scalar-lsq-unification §6 S2「共通規則」・§5.1 #2f。

用途: 収束済み (ただし残差がプラトーで止まっている) 起点から `restart_field.py` で継続した A/B の双子が、
作用素の切り替え後に「起点と同じ残差床へ戻ったか」を判定する。`check_convergence.py --from-floor` は
**全期間のピーク**も床の倍率以内を要求し、かつ参照 run が通常判定 PASS でないと REFUSED になるので、
切り替え直後の跳ねがあり起点が plateau の S2 には使えない (codex plan M4)。本ツールはその代替:

  (i)  `check_convergence.analyze` で NaN/Inf が無く、RISING 列が無い (plateau = STALLED は可。低下桁数は問わない)。
  (ii) 各残差列の末尾平均 (末尾 --tail の |値| 平均) が、起点 run の末尾平均の --factor 倍以内 (既定 1.5)。
       切り替え直後のピークは判定に使わない (記録のみ)。**再進入 step** (以後ずっと factor 倍以内に留まる最初の step)
       を列ごとに記録する。列の対応は `--from-floor` と同じ規則 (起点に無い列・起点で all-zero の列が非ゼロなら FAIL)。
  比較量の定常性 (iii) は `check_quasisteady.py` で別に判定する (本ツールは残差だけ)。

使い方:
  python3 tools/check_floor_ratio.py --start REF_RUN RUN [RUN ...] [--factor 1.5] [--tail 0.2]
  python3 tools/check_floor_ratio.py --selftest
VERDICT: PASS / NOT BACK ON FLOOR / RISING / DIVERGED / 判定不能 (判定不能は合格ではない)。
"""
import argparse
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_convergence as CC  # noqa: E402


def csv_path(rd):
    return rd if rd.endswith('.csv') else os.path.join(rd, 'residual_history.csv')


def judge(path, floor, zero_cols, factor, tail_frac, ref_path):
    """戻り値 (verdict, laststep, report)。report は列 -> (msg, ok)。"""
    res = CC.analyze(path, 0.0, tail_frac)   # min_drop 0: 低下桁数は問わず trend と NaN だけ使う
    if res is None:
        return '判定不能 (空の残差ファイル)', '?', {}
    laststep, rep_cc, _, any_nan, _, _ = res
    rows, cols = CC.load_series(path)
    steps = [r.get('step', str(i)) for i, r in enumerate(rows)]
    report = {}
    ok_in = True
    if CC.transition_active(ref_path) != CC.transition_active(path):
        report['(起点)'] = ('起点 run と遷移モデルの有無が違う (別の方程式系)  <-- 判定不能', False); ok_in = False
    for k in ('(入力)',):
        if k in rep_cc and not rep_cc[k][1]:
            report[k] = rep_cc[k]; ok_in = False
    for c in floor:
        if c not in cols:
            report[c] = ('起点に非ゼロの床があるのに対象に列が無い  <-- COLUMN MISMATCH', False); ok_in = False
    rising = False
    back = True
    for c, ser in cols.items():
        if not any(v != 0.0 for v in ser):
            if c in floor:
                # 起点で活動していた列が対象で全ゼロ = 方程式が消えた / 出力されていない (codex result-1 M1 の追加例)
                report[c] = ('起点で活動していた列が対象で all-zero  <-- COLUMN MISMATCH', False); ok_in = False
            else:
                report[c] = ('all-zero (both inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False); continue
        if c not in floor:
            why = '起点で all-zero の列が対象で非ゼロ' if c in zero_cols else '起点に無い列'
            report[c] = (f'{why}  <-- COLUMN MISMATCH', False); ok_in = False; continue
        cc_msg = rep_cc.get(c, ("", True))[0]
        if "判定不能" in cc_msg:
            # 通常判定が列単位で判定不能 (末尾窓の代表値が 0 など) なら床比でも合格にしない (codex result-1 M1)
            report[c] = (f"{cc_msg.strip()}  (check_convergence の列判定)", False); ok_in = False; continue
        a = ser[int(len(ser) * (1 - tail_frac)):]
        r_tail = (sum(abs(x) for x in a) / len(a)) / floor[c]
        peak_i = max(range(len(ser)), key=lambda i: abs(ser[i]))
        r_peak = abs(ser[peak_i]) / floor[c]
        lim = factor * floor[c]
        last_out = max((i for i, x in enumerate(ser) if abs(x) > lim), default=-1)
        reentry = 'never left' if last_out < 0 else (steps[last_out + 1] if last_out + 1 < len(ser) else 'not re-entered')
        trend = rep_cc.get(c, ('', True))[0]
        is_rising = 'RISING' in trend
        rising |= is_rising
        col_ok = (r_tail <= factor) and not is_rising
        back &= r_tail <= factor
        report[c] = (f"tail/floor={r_tail:6.3f} peak/floor={r_peak:8.2f} @step {steps[peak_i]} re-entry={reentry}"
                     f"{'  <-- RISING' if is_rising else ''}{'' if r_tail <= factor else '  <-- NOT BACK ON FLOOR'}", col_ok)
    if any_nan:
        v = 'DIVERGED (NaN/Inf)'
    elif not ok_in:
        v = '判定不能 (列・方程式系の不一致、または列単位の判定不能)'
    elif rising:
        v = 'RISING'
    elif not back:
        v = f'NOT BACK ON FLOOR (末尾平均 > {factor}× 起点床)'
    else:
        v = f'PASS (plateau 許容、末尾平均 ≤ {factor}× 起点床、RISING なし)'
    return v, laststep, report


def ref_problems(ref_path, floor):
    """起点系列の検査 (codex result-1 M1): 必須列・NaN/Inf・床が正の有限値か。問題の文字列リスト (空なら可)。"""
    out = []
    _, cols = CC.load_series(ref_path)
    need = CC.REQUIRED_COLS + (CC.TRANSITION_COLS if CC.transition_active(ref_path) else ())
    miss = [c for c in need if c not in cols]
    if miss:
        out.append("必須列が無い " + ", ".join(miss))
    for c, ser in cols.items():
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            out.append(f"{c} に NaN/Inf")
    for c, f in floor.items():
        if not (math.isfinite(f) and f > 0.0):
            out.append(f"{c} の床 {f!r} が正の有限値でない")
    return out


def run(start, targets, factor, tail_frac):
    ref_path = csv_path(start)
    if not os.path.exists(ref_path):
        print(f"[{start}] NO residual_history.csv (起点)  <-- 判定不能"); return False
    floor, zero_cols = CC.reference_floor(ref_path, tail_frac)
    bad = ref_problems(ref_path, floor)
    if bad:
        print(f"[{start}] 起点の残差系列が床の基準にならない: {'; '.join(bad)}  <-- 判定不能")
        return False
    print(f"起点床: {start} (末尾 {tail_frac:.0%} の |値| 平均、判定は通常基準を問わない)")
    for c in sorted(floor):
        print(f"  {c:14s}: {floor[c]:.3e}")
    if zero_cols:
        print(f"  inactive: {sorted(zero_cols)}")
    all_ok = True
    for rd in targets:
        p = csv_path(rd)
        if not os.path.exists(p):
            print(f"\n[{rd}] NO residual_history.csv  <-- 判定不能"); all_ok = False; continue
        v, laststep, report = judge(p, floor, zero_cols, factor, tail_frac, ref_path)
        print(f"\n=== {rd}  [last step {laststep}]  -> {v} ===")
        for c, (msg, _) in report.items():
            print(f"  {c:14s}: {msg}")
        all_ok &= v.startswith('PASS')
    print(f"\nOVERALL: {'ALL PASS' if all_ok else 'CHECK FAILURES ABOVE'}")
    return all_ok


def _write(path, series, extra_cols=()):
    cols = list(CC.REQUIRED_COLS) + list(extra_cols)
    with open(path, 'w') as f:
        f.write('step,phase,' + ','.join(cols) + '\n')
        for i, v in enumerate(series):
            f.write(f"{i},outer_end," + ','.join(str(v) for _ in cols) + '\n')


def selftest():
    """合成系列で 5 通りの判定を確かめる (判定ロジックを触ったら回す)。"""
    import random
    random.seed(0)
    d = tempfile.mkdtemp(prefix='floor_ratio_')
    base = [1e-6 * (1 + 0.3 * random.random()) for _ in range(400)]
    cases = {
        'plateau_same': (base, 'PASS'),
        'kick_then_back': ([1e-3] * 20 + base[20:], 'PASS'),
        'higher_floor': ([2e-6 * (1 + 0.3 * random.random()) for _ in range(400)], 'NOT BACK'),
        'rising': ([1e-6 * (1 + i / 40.0) for i in range(400)], None),   # RISING か NOT BACK のどちらでも不合格
        'nan': (base[:399] + [float('nan')], 'DIVERGED'),
    }
    ref = os.path.join(d, 'ref.csv'); _write(ref, base)
    floor, zc = CC.reference_floor(ref, 0.2)
    ok = True
    for name, (ser, want) in cases.items():
        p = os.path.join(d, name + '.csv'); _write(p, ser)
        v, _, _ = judge(p, floor, zc, 1.5, 0.2, ref)
        good = (v.startswith(want) if want else not v.startswith('PASS'))
        ok &= good
        print(f"{name:16s} -> {v}  [{'ok' if good else 'NG'}]")
    # 列の不一致: 起点に無い列が対象で非ゼロ
    p = os.path.join(d, 'extra.csv'); _write(p, base, extra_cols=('rms_roXi',))
    v, _, _ = judge(p, floor, zc, 1.5, 0.2, ref)
    good = v.startswith('判定不能'); ok &= good
    print(f"{'extra_column':16s} -> {v}  [{'ok' if good else 'NG'}]")
    # codex result-1 M1 の 2 例: 末尾が 0 の系列 / 参照が Inf
    p = os.path.join(d, 'tailzero.csv'); _write(p, [1e-6] * 320 + [0.0] * 80)
    v, _, _ = judge(p, floor, zc, 1.5, 0.2, ref)
    good = not v.startswith('PASS'); ok &= good
    print(f"{'tail_zero':16s} -> {v}  [{'ok' if good else 'NG'}]")
    refinf = os.path.join(d, 'refinf.csv'); _write(refinf, [float('inf')] * 400)
    fi, _ = CC.reference_floor(refinf, 0.2)
    good = bool(ref_problems(refinf, fi)) and not run(refinf, [os.path.join(d, 'plateau_same.csv')], 1.5, 0.2)
    ok &= good
    print(f"{'ref_inf':16s} -> {'判定不能 (拒否)' if good else 'PASS してしまう'}  [{'ok' if good else 'NG'}]")
    p = os.path.join(d, 'allzero.csv'); _write(p, [0.0] * 400)
    v, _, _ = judge(p, floor, zc, 1.5, 0.2, ref)
    good = not v.startswith('PASS'); ok &= good
    print(f"{'target_allzero':16s} -> {v}  [{'ok' if good else 'NG'}]")
    print('SELFTEST', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='*')
    ap.add_argument('--start', help='起点 run (restart 元の収束場を出した run、または その residual_history.csv)')
    ap.add_argument('--factor', type=float, default=1.5)
    ap.add_argument('--tail', type=float, default=0.2)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not a.start or not a.runs:
        ap.error('--start と対象 run が要る')
    sys.exit(0 if run(a.start, a.runs, a.factor, a.tail) else 1)
