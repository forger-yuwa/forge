"""#4pjg 補: 射影集合の節点で、実現可能域の各境界の余裕 H を段ごとにどちらへ動かすかを出す。

境界 (Q3 ∝ g、ρ は斉次 0 次なので保存量のまま使える):
  H_up  = ln Q1 + ln g − 2 ln Q2      (y² ≤ x  ⇔ Q2² ≤ Q1·Q3; 0 半径への集積側)
  H_low = ln Q0 + ln Q2 − 2 ln Q1     (x² ≤ y  ⇔ 単分散側)
  H_x1  = 2 ln Q0 + ln g − 3 ln Q1    (x ≤ 1)
段 s の寄与 dH_s/dt = Σ_c a_c R_c,s/(V q_c) (a_c は上の係数、R は記録点の差)。H < 0 側へ押す段が射影と戦っている段。
使い方: python3 boundary_push.py <ON tp_operator.h5> <proj.npy>
"""
import sys
import h5py
import numpy as np

f = h5py.File(sys.argv[1], "r")
mask = np.load(sys.argv[2])
labels = [l.decode() if isinstance(l, bytes) else l for l in f["res/labels"][:]]
snap = f["res/snap"][:].astype(np.float64)
comps = [c.decode() if isinstance(c, bytes) else c for c in f.attrs["components"]]
st = f["node/state"][:].astype(np.float64)
V = f["node/volume"][:].astype(np.float64)
ci = {c: comps.index(c) for c in ["g", "Q2", "Q1", "Q0"]}
L = {l: labels.index(l) for l in labels}
stages = {"adv": ("cm_zero", "cm_adv"), "diff": ("cm_adv", "tp_diff"), "src": ("tp_diff", "cond_src"),
          "after_src": ("cond_src", "final")}
idx = np.asarray(mask, dtype=np.int64)
q = {c: st[ci[c], idx] for c in ci}
pos = np.all([q[c] > 0 for c in ci], axis=0)
print(f"射影集合 {idx.size} 節点のうち全成分正 {int(pos.sum())} (他は除外)")
idx = idx[pos]
q = {c: st[ci[c], idx] for c in ci}
# どの境界に張り付いているか (ρ_l は 1000 kg/m³ 近似; 判別用)。x = Q1/(Q0 r), y = Q2/(Q0 r²), r³ = Q3/Q0
rho_l = 1000.0
r = np.cbrt(q["g"] / (4.0 / 3.0 * np.pi * rho_l) / q["Q0"])
x = q["Q1"] / (q["Q0"] * r); y = q["Q2"] / (q["Q0"] * r * r)
print(f"tight (1e-3 以内): y=√x {np.mean(np.abs(y*y/x-1)<1e-3):.3f}, y=x² {np.mean(np.abs(x*x/y-1)<1e-3):.3f}, x=1 {np.mean(np.abs(x-1)<1e-3):.3f}")
print(f"  median x {np.median(x):.4f}, y {np.median(y):.4f}, y²/x {np.median(y*y/x):.6f}")
Q3 = q["g"]
Hs = {"H_up": {"Q1": 1, "g": 1, "Q2": -2}, "H_low": {"Q0": 1, "Q2": 1, "Q1": -2}, "H_x1": {"Q0": 2, "g": 1, "Q1": -3}}
for h, co in Hs.items():
    print(f"\n{h}: coefficients {co}")
    print(f"  {'stage':10s} {'median dH/dt':>14s} {'frac<0':>8s} {'p10':>12s} {'p90':>12s}")
    tot = np.zeros(idx.size)
    for s, (a, b) in stages.items():
        R = {c: snap[L[b], ci[c], idx] - snap[L[a], ci[c], idx] for c in ci}
        d = sum(a_ * R[c] / (V[idx] * q[c]) for c, a_ in co.items())
        tot += d
        print(f"  {s:10s} {np.median(d):+14.4e} {np.mean(d < 0):8.3f} {np.percentile(d,10):+12.4e} {np.percentile(d,90):+12.4e}")
    print(f"  {'total':10s} {np.median(tot):+14.4e} {np.mean(tot < 0):8.3f} {np.percentile(tot,10):+12.4e} {np.percentile(tot,90):+12.4e}")
br = f["source/slots"][list(map(lambda s: s.decode() if isinstance(s, bytes) else s, f["source/slot_names"][:])).index("branch"), idx]
print("\nbranch on set:", dict(zip(*np.unique(br, return_counts=True))))
