#!/usr/bin/env python3
"""V2c 判定 (plan boundary-node-farfield-characteristic §6 V2c): 下面 (評価線 = ランプ 0.2–0.7 m + 水平部 –1.2 m、z = 5 mm の節点列) の圧力で
  B (farfield) 合格: max|p_B − p_C| ≤ 0.02 Δp_shock、 対照: max|p_A − p_C| ≥ 0.1 Δp_shock
Δp_shock = 斜め衝撃波理論 (M 2.5、θ 10°、γ 1.4) の p2 − p1。各 run の最終 res_*.h5 を使う。
  python3 eval_v2c.py RUN_A RUN_B RUN_C
"""
import glob, math, os, re, sys
import h5py, numpy as np

XC, TH, DZ = 0.2, math.radians(10.0), 0.005
XE_OF = {"v1": 0.7, "v2": 1.8}


def oblique(M, th, g=1.4):
    lo, hi = math.asin(1.0 / M) + 1e-9, math.radians(64.0)   # 弱い解の枝
    f = lambda b: math.tan(th) - 2.0 / math.tan(b) * (M * M * math.sin(b) ** 2 - 1.0) / (M * M * (g + math.cos(2 * b)) + 2.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        (lo, hi) = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    b = 0.5 * (lo + hi); Mn = M * math.sin(b)
    return b, 1.0 + 2.0 * g / (g + 1.0) * (Mn * Mn - 1.0)


def ramp_p(run):
    fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(run + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
    st, last = fs[-1]
    with h5py.File(glob.glob(run + "/ramp.h5")[0]) as m:
        xyz = np.array(m["MESH/COORD"]).reshape(-1, 3)
    with h5py.File(last) as f:
        P = f["VALUE/P"][:].astype(float)
        nonfin = sum(int(np.sum(~np.isfinite(f["VALUE"][k][:]))) for k in f["VALUE"].keys())
    g = open(run + "/GEOM.txt").read().strip() if os.path.exists(run + "/GEOM.txt") else "v1"
    ybx = (np.clip(xyz[:, 0], XC, XE_OF[g]) - XC) * math.tan(TH)
    sel = (np.abs(xyz[:, 1] - ybx) < 1e-6) & (np.abs(xyz[:, 2] - DZ) < 1e-6) & (xyz[:, 0] >= XC - 1e-9)
    x = xyz[sel, 0]; o = np.argsort(x)
    return st, np.round(x[o], 6), P[sel][o], nonfin


b, pr = oblique(2.5, TH)
P1 = 101325.0
dps = (pr - 1.0) * P1
ra, rb, rc = sys.argv[1:4]
(sa, xa, pa, na), (sb, xb, pb, nb), (sc, xc, pc, nc) = ramp_p(ra), ramp_p(rb), ramp_p(rc)
assert np.array_equal(xa, xc) and np.array_equal(xb, xc), "評価線の節点が一致しない"
dB, dA = np.abs(pb - pc), np.abs(pa - pc)
print(f"斜め衝撃波理論: β {math.degrees(b):.2f}°、p2/p1 {pr:.4f}、Δp_shock {dps:.1f} Pa。評価線 {len(xc)} 節点 (x {xc[0]:.3f}–{xc[-1]:.3f} m)")
print(f"  ランプ面の圧力 (C、角の直後 x 0.3 m): {np.interp(0.3, xc, pc):.1f} Pa (理論 p2 {pr * P1:.1f})")
print(f"  B (farfield) vs C: max|Δp| {dB.max():.1f} Pa = {dB.max() / dps:.4f} Δp (x {xc[dB.argmax()]:.3f} m)  step {sb} / {sc}、非有限 {nb}/{nc}")
print(f"  A (slip)     vs C: max|Δp| {dA.max():.1f} Pa = {dA.max() / dps:.4f} Δp (x {xc[dA.argmax()]:.3f} m)  step {sa}、非有限 {na}")
okB, okA = dB.max() <= 0.02 * dps and nb == 0 and nc == 0, dA.max() >= 0.1 * dps
print(f"VERDICT: {'PASS' if okB and okA else 'FAIL'} (B ≤ 0.02 Δp: {'OK' if okB else 'NG'}、対照 A ≥ 0.1 Δp: {'OK' if okA else 'NG'})")
np.savetxt(sys.argv[4] if len(sys.argv) > 4 else "v2c_ramp_pressure.csv", np.c_[xc, pa, pb, pc], delimiter=",", header="x,p_A,p_B,p_C", comments="")
