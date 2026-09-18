#!/usr/bin/env python3
"""check_field_regress.py の試験 (plan config-key-pruning §6.2')。

重点は **NaN を含む候補が PASS にならないこと** (codex plan-3 M2: `max(0.0, NaN)` が `0.0` になるため、
比較の中では異常が消える)。欠落・形状不一致も数値判定の前に落ちること。
"""
import importlib.util, os, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cfr', os.path.join(here, '..', '..', 'tools', 'check_field_regress.py'))
cfr = importlib.util.module_from_spec(spec); spec.loader.exec_module(cfr)

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)

Q = ['ro', 'roUx']
def fld(**kw):
    out = {q: np.array([1.0, 2.0, 3.0]) for q in Q}
    out.update(kw); return out

# 健全な 3 反復 + 候補は素通りする
d = {r: fld() for r in ['r1', 'r2', 'r3', 'cand']}
check(cfr.data_check(d, ['r1', 'r2', 'r3', 'cand'], Q) == [], 'clean data must produce no problems')

# 候補に NaN
d['cand'] = fld(ro=np.array([1.0, np.nan, 3.0]))
p = cfr.data_check(d, ['r1', 'r2', 'r3', 'cand'], Q)
check(len(p) == 1 and 'cand' in p[0] and '非有限' in p[0], f'NaN in the candidate must be reported, got {p}')

# 反復側に Inf (ノイズ床を壊す)
d = {r: fld() for r in ['r1', 'r2', 'r3', 'cand']}
d['r2'] = fld(roUx=np.array([1.0, np.inf, 3.0]))
p = cfr.data_check(d, ['r1', 'r2', 'r3', 'cand'], Q)
check(len(p) == 1 and 'r2' in p[0], f'Inf in a repeat must be reported, got {p}')

# 必須量の欠落
d = {r: fld() for r in ['r1', 'r2', 'cand']}
del d['cand']['roUx']
p = cfr.data_check(d, ['r1', 'r2', 'cand'], Q)
check(len(p) == 1 and 'roUx' in p[0], f'a missing required quantity must be reported, got {p}')

# 形状不一致
d = {r: fld() for r in ['r1', 'r2', 'cand']}
d['cand'] = fld(ro=np.array([1.0, 2.0]))
p = cfr.data_check(d, ['r1', 'r2', 'cand'], Q)
check(any('形状' in x for x in p), f'a shape mismatch must be reported, got {p}')

# ノルムの向き: 差がゼロなら 0、正規化は基準側
a = np.array([1.0, 2.0, 3.0])
check(cfr.norms(a, a) == (0.0, 0.0), 'identical fields must give zero norms')
n2, ni = cfr.norms(a, a + 1.0)
check(abs(ni - 1.0/3.0) < 1e-12, f'relative L-inf must normalize by the reference maximum, got {ni}')

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
