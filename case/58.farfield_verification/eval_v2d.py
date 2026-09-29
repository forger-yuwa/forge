#!/usr/bin/env python3
"""V2d の判定 (plan boundary-node-farfield-characteristic §6 V2d、判定 (b) = 短領域 vs 長領域 = 境界が加える誤差)。
  contact SHORT LONG : V2d-1。評価点 (右端から 5 セル内側) で
       max_t |P_short − P∞|/P∞ ≤ 1e-3、 max_t |T_short − T_long| ≤ 0.02 × (600 − 220) K、 max_t |Y_short − Y_long| ≤ 0.002
       (Y は 10 step ごとの res_*.h5 から評価点の節点値を読む。T・P はプローブ)
  nopulse SHORT      : V2d-2 パルスなし対照。評価点 (x 0.8 m) で max_t |P − P∞| ≤ 0.01 × 1e-3 P∞
  acoustic SHORT LONG: V2d-2。反射 = max_t |P_short − P_long| / 入射振幅 ≤ 0.05 (入射振幅 = max_t |P_long − P∞|)
"""
import glob, re, sys
import h5py, numpy as np

P0 = 2851.0


def probe(run, i=0):
    d = np.genfromtxt(f"{run}/point_probe_{i}.out", delimiter=",", names=True)
    return np.asarray(d["TotalTime"], float), np.asarray(d["P"], float), np.asarray(d["T"], float)


def y_series(run, xeval):
    fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(run + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
    with h5py.File(run + "/chan.h5") as m:
        xyz = np.array(m["MESH/COORD"]).reshape(-1, 3)
    node = int(np.argmin((xyz[:, 0] - xeval) ** 2 + (xyz[:, 1] - 0.005) ** 2 + (xyz[:, 2] - 0.005) ** 2))
    t, Y = [], []
    for st, p in fs:
        with h5py.File(p) as f:
            V = f["VALUE"]
            if "roY0" not in V:
                return None
            Y.append(float(V["roY0"][node] / V["ro"][node]))
            t.append(st)
    return np.array(t), np.array(Y)


mode = sys.argv[1]
if mode == "contact":
    s, l = sys.argv[2], sys.argv[3]
    ts, Ps, Ts = probe(s); tl, Pl, Tl = probe(l)
    tm = min(ts[-1], tl[-1]); k = ts <= tm
    dP = np.max(np.abs(Ps[k] - P0)) / P0
    dT = np.max(np.abs(Ts[k] - np.interp(ts[k], tl, Tl)))
    Tamp = np.max(Tl) - 220.0
    out = [f"max|P−P∞|/P∞ {dP:.3e} (≤ 1e-3)", f"max|ΔT| {dT:.3f} K (≤ {0.02 * 380:.1f} K、塊の到達振幅 {Tamp:.1f} K)"]
    ok = dP <= 1e-3 and dT <= 0.02 * 380.0
    ys, yl = y_series(s, 0.975), y_series(l, 0.975)
    if ys is not None:
        n = min(len(ys[0]), len(yl[0]))
        assert np.array_equal(ys[0][:n], yl[0][:n])
        dY = np.max(np.abs(ys[1][:n] - yl[1][:n]))
        out.append(f"max|ΔY| {dY:.2e} (≤ 0.002、塊の到達 Y {np.max(yl[1]):.4f})")
        ok = ok and dY <= 0.002
    print(f"{s} vs {l}: " + "、".join(out))
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
elif mode == "nopulse":
    s = sys.argv[2]
    ts, Ps, _ = probe(s)
    d = np.max(np.abs(Ps - P0))
    print(f"{s}: パルスなし max|P − P∞| {d:.4e} Pa (≤ {0.01 * 1e-3 * P0:.4e})")
    print(f"VERDICT: {'PASS' if d <= 0.01 * 1e-3 * P0 else 'FAIL'}")
elif mode == "acoustic":
    s, l = sys.argv[2], sys.argv[3]
    ts, Ps, _ = probe(s); tl, Pl, _ = probe(l)
    tm = min(ts[-1], tl[-1]); k = ts <= tm
    A = np.max(np.abs(Pl[tl <= tm] - P0))
    d = np.abs(Ps[k] - np.interp(ts[k], tl, Pl))
    print(f"{s} vs {l}: 入射振幅 {A:.4g} Pa、max|ΔP| {d.max():.4g} Pa (t {ts[k][d.argmax()]:.4e} s) → 反射率 {d.max() / A:.4%}")
    print(f"VERDICT: {'PASS' if d.max() / A <= 0.05 else 'FAIL'} (≤ 5 %)")
