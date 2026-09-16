#!/usr/bin/env python3
"""dual-time の 3 水準時間次数 (plan species-passive-scalar-unification §6-6, §5.1 #12/#18/#20)。

各 run の最終 res_*.h5 (同一物理時刻) を読み、非重み付き L2 差で
  e(dt)   = ||q(dt) − q(dt/2)||,  観測次数 p = log2(e(2dt)/e(dt))
  sub-iter 比 = ||q(nSub×2) − q(nSub)|| / e(最小水準)   (≤0.1 を要求)
を量ごとに出す。使い方:
  analyze_moment_order.py --levels RUN_2dt RUN_dt RUN_dt/2 [--nsub RUN_dt_nsub2x] [--fields g_0,Q0_0,T,ro]
"""
import argparse, glob, os, re, sys
import numpy as np
import h5py


def last_res(run):
    fs = glob.glob(os.path.join(run, 'res_*.h5'))
    if not fs:
        sys.exit(f'no res_*.h5 in {run}')
    return max(fs, key=lambda f: int(re.search(r'res_(\d+)\.h5', f).group(1)))


def load(run, fields):
    f = last_res(run)
    with h5py.File(f, 'r') as h:
        d = {k: np.asarray(h['VALUE'][k], dtype=np.float64) for k in fields}
    return f, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--levels', nargs=3, required=True, help='dt 2dt, dt, dt/2 の run (同一物理時刻)')
    ap.add_argument('--nsub', default=None, help='dt 水準で nSub を倍にした run')
    ap.add_argument('--fields', default='g_0,Q0_0,Q1_0,Q2_0,T,ro,P')
    a = ap.parse_args()
    fields = a.fields.split(',')
    files, data = zip(*[load(r, fields) for r in a.levels])
    for r, f in zip(a.levels, files):
        print(f'  {r}: {os.path.basename(f)}')
    nsub = load(a.nsub, fields) if a.nsub else None
    if nsub:
        print(f'  nSub x2: {a.nsub}: {os.path.basename(nsub[0])}')
    bad = False
    print(f"{'field':8s} {'e(2dt)':>11s} {'e(dt)':>11s} {'order':>7s} {'nSub diff':>11s} {'ratio':>7s}  max|q|")
    for k in fields:
        q0, q1, q2 = data[0][k], data[1][k], data[2][k]
        if np.isnan(q0).any() or np.isnan(q1).any() or np.isnan(q2).any():
            print(f'{k:8s} NaN present'); bad = True; continue
        e0 = np.linalg.norm(q0 - q1); e1 = np.linalg.norm(q1 - q2)
        order = np.log2(e0 / e1) if e1 > 0 else float('inf')
        line = f'{k:8s} {e0:11.4e} {e1:11.4e} {order:7.3f}'
        if nsub:
            es = np.linalg.norm(nsub[1][k] - q1)
            ratio = es / e1 if e1 > 0 else float('inf')
            line += f' {es:11.4e} {ratio:7.3f}'
            if ratio > 0.1:
                line += '  <-- sub-iter error not resolved (>0.1)'; bad = True
        else:
            line += ' ' * 20
        line += f'  {np.abs(q1).max():.3e}'
        print(line)
    print('VERDICT:', 'FAIL' if bad else 'see orders above')


if __name__ == '__main__':
    main()
