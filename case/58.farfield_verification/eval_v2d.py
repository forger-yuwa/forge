#!/usr/bin/env python3
"""V2d の判定 (plan boundary-node-farfield-characteristic §6 V2d、判定 (b) = 短領域 vs 長領域 = 境界が加える誤差)。
  contact SHORT LONG : V2d-1。評価点 (右端から 5 セル内側) で
       max_t |P_short − P∞|/P∞ ≤ 1e-3、 max_t |T_short − T_long| ≤ 0.02 × (600 − 220) K、 max_t |Y_short − Y_long| ≤ 0.002
       (Y は 10 step ごとの res_*.h5 から評価点の節点値を読む。T・P はプローブ)
  nopulse SHORT      : V2d-2 パルスなし対照。評価点 (x 0.8 m) で max_t |P − P∞| ≤ 0.01 × 1e-3 P∞
  acoustic SHORT LONG: V2d-2 (旧)。反射 = max_t |P_short − P_long| / 入射振幅 ≤ 0.05 (入射振幅 = max_t |P_long − P∞|)
  acoustic_win SHORT LONG: V2d-2 (再設計、plan §5.1 #3b)。入射窓 = t_inc ± 3σ/(c+u) で入射振幅 max|P_long − P∞|、
       反射到達窓 = t_ref ± 3σ/(c−u) で max|P_short − P_long|。t_inc = (0.8 − 0.4)/(c+u)、t_ref = (1 − 0.4)/(c+u) + (1 − 0.8)/(c−u)、
       σ = 0.04/2.355、c・u は内部の値 (IC_FROM.txt から)。反射 ≤ 0.05
"""
import re as _re
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
elif mode == "acoustic_win":
    s, l = sys.argv[2], sys.argv[3]
    info = open(f"{s}/IC_FROM.txt").read()
    U = float(_re.search(r"U ([-0-9.eE+]+) \(dir", info).group(1)); ci = float(_re.search(r"内部 T [0-9.]+ Y [0-9.]+ c ([0-9.eE+]+)", info).group(1))
    sig = 0.04 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    t_inc = 0.4 / (ci + U); t_ref = 0.6 / (ci + U) + 0.2 / (ci - U)
    w_inc, w_ref = 3 * sig / (ci + U), 3 * sig / (ci - U)
    ts, Ps, _ = probe(s); tl, Pl, _ = probe(l)
    ki = (tl >= t_inc - w_inc) & (tl <= t_inc + w_inc)
    kr = (ts >= t_ref - w_ref) & (ts <= t_ref + w_ref)
    if ts[-1] < t_ref + w_ref or tl[-1] < t_ref + w_ref:
        raise SystemExit(f"run が反射到達窓の終わり ({t_ref + w_ref:.4e} s) まで回っていない")
    A = np.max(np.abs(Pl[ki] - P0))
    d = np.abs(Ps[kr] - np.interp(ts[kr], tl, Pl))
    print(f"{s} vs {l}: c_i {ci:.5g} u {U:.5g}、入射窓 [{t_inc - w_inc:.3e}, {t_inc + w_inc:.3e}] 振幅 {A:.4g} Pa、"
          f"反射到達窓 [{t_ref - w_ref:.3e}, {t_ref + w_ref:.3e}] max|ΔP| {d.max():.4g} Pa → 反射率 {d.max() / A:.4%}")
    print(f"VERDICT: {'PASS' if d.max() / A <= 0.05 else 'FAIL'} (≤ 5 %)")
elif mode == "same":
    # 時間精度: 2 run の評価点 P の差の最大 / 入射振幅 (長領域の入射窓の振幅、P∞ = 2851 Pa) ≤ 0.01
    a, b, l = sys.argv[2], sys.argv[3], sys.argv[4]
    info = open(f"{l}/IC_FROM.txt").read()
    U = float(_re.search(r"U ([-0-9.eE+]+) \(dir", info).group(1)); ci = float(_re.search(r"内部 T [0-9.]+ Y [0-9.]+ c ([0-9.eE+]+)", info).group(1))
    sig = 0.04 / (2.0 * np.sqrt(2.0 * np.log(2.0))); t_inc = 0.4 / (ci + U); w_inc = 3 * sig / (ci + U)
    # 固定評価区間 [0, 反射到達窓の末尾] (codex diagnose 2026-10-03 M4): 3 run とも区間を完全に覆い、時刻が単調増加、
    # 全値が有限、入射振幅が有限かつ正であることを必須にする (満たさなければ判定不能で終了コード 2)。以前は短い方の終了時刻までしか比べなかった
    t_ref = 0.6 / (ci + U) + 0.2 / (ci - U); w_ref = 3 * sig / (ci - U); t_end = t_ref + w_ref
    tl, Pl, _ = probe(l); ta, Pa, _ = probe(a); tb, Pb, _ = probe(b)
    bad = []
    for nm, t, P in ((a, ta, Pa), (b, tb, Pb), (l, tl, Pl)):
        if len(t) < 2 or t[-1] < t_end:
            bad.append(f"{nm} が評価区間の末尾 {t_end:.4e} s まで無い (最終 {t[-1] if len(t) else float('nan'):.4e})")
        if np.any(np.diff(t) <= 0):
            bad.append(f"{nm} の時刻が単調増加でない")
        if not (np.all(np.isfinite(t)) and np.all(np.isfinite(P))):
            bad.append(f"{nm} に非有限値")
    ki = (tl >= t_inc - w_inc) & (tl <= t_inc + w_inc)
    A = float(np.max(np.abs(Pl[ki] - P0))) if ki.any() else float("nan")
    if not (np.isfinite(A) and A > 0):
        bad.append(f"入射振幅が有限・正でない ({A})")
    if bad:
        print("判定不能:\n  " + "\n  ".join(bad)); print("VERDICT: UNDECIDABLE"); sys.exit(2)
    k = ta <= t_end
    d = np.abs(Pa[k] - np.interp(ta[k], tb, Pb))
    print(f"{a} vs {b}: 評価区間 [0, {t_end:.4e}] s ({int(k.sum())} 点)、max|ΔP| {d.max():.4g} Pa (t {ta[k][d.argmax()]:.4e}) / 入射振幅 {A:.4g} Pa = {d.max() / A:.3%}")
    rel = d.max() / A
    if len(sys.argv) > 5 and sys.argv[5] == "--nsub":
        # nSub 感度 (codex diagnose 2026-10-03 ①): D_N ≤ 0.002 → 反復数依存の説明は弱い、> 0.01 → 「nSub 20 で十分」を棄却、中間は保留
        print(f"D_N = {rel:.4%} → " + ("≤ 0.2 %: 20→40 の感度は小さい" if rel <= 0.002 else ("> 1 %: nSub 20 で十分を棄却" if rel > 0.01 else "判定保留 (0.2–1 %)")))
    else:
        print(f"VERDICT: {'PASS' if rel <= 0.01 else 'FAIL'} (≤ 1 %)")
