#!/usr/bin/env python3
"""V2a 音響反射の判定 (plan boundary-node-farfield-characteristic §6 V2a)。
評価点 (probe 0、x = 0.8 m) の圧力時系列で、
  入射振幅 A = max_t |P_long − P∞|、反射 = max_t |P_short − P_long| / A  (全時間窓 [0, t_end]。長領域の右端からの反射は窓の外)
  python3 eval_v2a.py reflect SHORT LONG [--limit 0.05 | --min 0.9]    反射の判定 (上限 / 対照の下限)
  python3 eval_v2a.py same RUN_A RUN_B --ref LONG [--limit 0.01]      時間精度: 2 run の差の最大 / A ≤ limit
時刻は TotalTime で照合し、刻みの違う run 同士は線形補間で揃える。
"""
import sys
import numpy as np

P_INF = 101325.0


def load(run, i=0):
    d = np.genfromtxt(f"{run}/point_probe_{i}.out", delimiter=",", names=True)
    return np.asarray(d["TotalTime"], float), np.asarray(d["P"], float)


def on(t, tp, p):
    return np.interp(t, tp, p)


mode = sys.argv[1]
lim = float(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else None
mn = float(sys.argv[sys.argv.index("--min") + 1]) if "--min" in sys.argv else None
if mode == "reflect":
    s, l = sys.argv[2], sys.argv[3]
    ts, ps = load(s); tl, pl = load(l)
    tmax = min(ts[-1], tl[-1])
    t = ts[ts <= tmax]
    A = np.max(np.abs(pl[tl <= tmax] - P_INF))
    d = np.abs(ps[ts <= tmax] - on(t, tl, pl))
    r = d.max() / A
    print(f"{s} vs {l}: 入射振幅 {A:.4g} Pa、反射 max|ΔP| {d.max():.4g} Pa (t = {t[d.argmax()]:.4e} s) → 反射率 {r:.4%}")
    if lim is not None:
        print(f"VERDICT: {'PASS' if r <= lim else 'FAIL'} (反射率 {r:.4%} ≤ {lim:.2%})")
    if mn is not None:
        print(f"VERDICT: {'PASS' if r >= mn else 'FAIL'} (対照: 反射率 {r:.4%} ≥ {mn:.0%})")
elif mode == "same":
    a, b = sys.argv[2], sys.argv[3]
    ref = sys.argv[sys.argv.index("--ref") + 1]
    ta, pa = load(a); tb, pb = load(b); tl, pl = load(ref)
    tmax = min(ta[-1], tb[-1])
    A = np.max(np.abs(pl[tl <= tmax] - P_INF))
    t = ta[ta <= tmax]
    d = np.abs(pa[ta <= tmax] - on(t, tb, pb))
    r = d.max() / A
    print(f"{a} vs {b}: max|ΔP| {d.max():.4g} Pa (t = {t[d.argmax()]:.4e} s) / 入射振幅 {A:.4g} = {r:.4%}")
    if lim is not None:
        print(f"VERDICT: {'PASS' if r <= lim else 'FAIL'} ({r:.4%} ≤ {lim:.2%})")
