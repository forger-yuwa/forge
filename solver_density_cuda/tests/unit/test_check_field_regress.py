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

# --- main() を通した反例 (codex result M1): 任意量の NaN・欠落が PASS になっていた ---
import subprocess, tempfile, os
try:
    import h5py
    HAVE_H5 = True
except Exception:
    HAVE_H5 = False

if HAVE_H5:
    tool = os.path.join(here, '..', '..', 'tools', 'check_field_regress.py')

    def make_run(d, step, extra=None, drop=()):
        os.makedirs(d, exist_ok=True)
        with h5py.File(os.path.join(d, f'res_{step}.h5'), 'w') as f:
            g = f.create_group('VALUE')
            vals = {q: np.array([1.0, 2.0, 3.0]) for q in ['ro', 'roUx', 'roUy', 'roUz', 'roe', 'P', 'T']}
            vals['roK'] = np.array([0.1, 0.2, 0.3])
            if extra: vals.update(extra)
            for q in drop: vals.pop(q, None)
            for q, v in vals.items():
                g.create_dataset(q, data=v)

    def run_tool(*runs_and_cands):
        reps, cand = runs_and_cands[:-1], runs_and_cands[-1]
        r = subprocess.run([sys.executable, tool, '--repeat', *reps, '--candidate', cand],
                           capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    with tempfile.TemporaryDirectory() as tmp:
        reps = [os.path.join(tmp, f'r{i}') for i in (1, 2, 3)]
        for d in reps: make_run(d, 10)
        # (a) 候補の任意量 roK が NaN
        c = os.path.join(tmp, 'c_nan'); make_run(c, 10, extra={'roK': np.array([0.1, np.nan, 0.3])})
        rc, out = run_tool(*reps, c)
        check(rc == 2 and 'roK' in out, f'NaN in an optional quantity must stop the run (rc={rc})')
        check('VERDICT: PASS' not in out, 'a NaN must never be reported as PASS')
        # (b) 候補から roK が欠落
        c2 = os.path.join(tmp, 'c_missing'); make_run(c2, 10, drop=('roK',))
        rc, out = run_tool(*reps, c2)
        check(rc == 2 and 'roK' in out, f'a quantity missing from one side must stop the run (rc={rc})')
        # (c) 健全なら通る
        c3 = os.path.join(tmp, 'c_ok'); make_run(c3, 10)
        rc, out = run_tool(*reps, c3)
        check(rc == 0 and 'VERDICT: PASS' in out, f'clean data must still pass (rc={rc})')
else:
    print('  (h5py が無いので main() 経由の反例は省略)')

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
