#!/usr/bin/env python3
"""R2 (plan boundary-node-periodic-gradient-fix §6 R2) の保存量・KE・エントロピー履歴。

    python3 r2_conservation.py --keep OLD_KEEP NEW_KEEP --slau OLD_SLAU NEW_SLAU [--png r2_ke_entropy.png]

積分の定義は `plot_ke_entropy_history.py` / `analyze_conservation.py` と同じ
(K=Σ½ρ|u|²V, S=Σρ c_v ln(P/ρ^γ) V, E=Σ roe V, M=Σ ρ V)。これに全運動量 P_i=Σ ρu_i V を足し、
node の周期重複 DOF を**一意 DOF** に揃える:

- `plot_ke_entropy_history.py` の `cv_weight` (境界 CV を多重度 2^(境界方向数) で割る) は、res の `volume` が
  重複 DOF ごとに**合併体積**を持つ場合の補正。部分体積を持つ場合は重み 1 が正しい。
  どちらかを決め打ちせず、**Σ V w が箱の体積 (2π)^3 と相対 1e-6 で一致する方**を採り、どちらも一致しなければ止める。

合格 (§6 R2、測る前に固定):
- KEEP (非粘性、既存閾値): |K/K0−1| ≲ 1 %、|ΔS/S0| ≲ 1e-4、全運動量 ≲ 1e-6 (ρ0U0V 正規化)。
- SLAU (Re=1600、新 run で判定): |M−M0|/M0 ≤ 1e-6、|P_i|/(ρ0U0V) ≤ 1e-6、|E−E0|/E0 ≤ 1e-5。
- KE・S 履歴の旧新差は記録のみ。
全 snapshot の最大値で判定する (最終値も併記)。
"""
import argparse
import csv
import glob
import os
import re

import h5py
import numpy as np

GAMMA, CP = 1.4, 0.4
CV = CP / GAMMA
DT = 0.007
L = 2.0 * np.pi
VBOX = L ** 3
RHO0, U0 = 1.0, 0.4


def step_of(p):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(p))
    return int(m.group(1)) if m else -1


def cv_weight(coord, n):
    if coord.size != 3 * n:
        return np.ones(n)
    c = coord.reshape(-1, 3)
    m = np.ones(n)
    for k in range(3):
        onb = (np.abs(c[:, k]) < 1e-5) | (np.abs(c[:, k] - L) < 1e-5)
        m = m * np.where(onb, 2.0, 1.0)
    return 1.0 / m


def weights(f):
    V = f["VALUE/volume"][:].astype(np.float64)
    w = cv_weight(f["MESH/COORD"][:].astype(np.float64), V.size)
    for name, ww in (("unit", np.ones_like(V)), ("cv_weight", w)):
        if abs(np.sum(V * ww) / VBOX - 1.0) < 1e-6:
            return V * ww, name
    raise SystemExit("ΣV が箱の体積と合わない: raw %.9g / weighted %.9g / box %.9g"
                     % (np.sum(V), np.sum(V * w), VBOX))


def totals(p):
    with h5py.File(p, "r") as f:
        Vw, how = weights(f)
        g = {k: f["VALUE/" + k][:].astype(np.float64) for k in ("ro", "Ux", "Uy", "Uz", "P", "roe")}
    ro = g["ro"]
    K = np.sum(0.5 * ro * (g["Ux"] ** 2 + g["Uy"] ** 2 + g["Uz"] ** 2) * Vw)
    S = np.sum(ro * CV * np.log(g["P"] / ro ** GAMMA) * Vw)
    E = np.sum(g["roe"] * Vw)
    M = np.sum(ro * Vw)
    P = [np.sum(ro * g[c] * Vw) for c in ("Ux", "Uy", "Uz")]
    return K, S, E, M, P, how


def series(run):
    fs = sorted([p for p in glob.glob(os.path.join(run, "res_*.h5")) if step_of(p) >= 0], key=step_of)
    rows, how = [], set()
    for p in fs:
        K, S, E, M, P, h = totals(p)
        how.add(h)
        rows.append([step_of(p), step_of(p) * DT, K, S, E, M] + P)
    a = np.array(rows)
    K0, S0, E0, M0 = a[0, 2], a[0, 3], a[0, 4], a[0, 5]
    out = np.column_stack([a[:, 0], a[:, 1], a[:, 2] / K0, (a[:, 3] - S0) / abs(S0), (a[:, 4] - E0) / E0,
                           (a[:, 5] - M0) / M0, a[:, 6:9] / (RHO0 * U0 * VBOX)])
    with open(os.path.join(run, "r2_conservation.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "time", "K_K0", "dS_S0", "dE_E0", "dM_M0", "Px_n", "Py_n", "Pz_n"])
        w.writerows(out.tolist())
    return out, ",".join(sorted(how))


def judge(label, val_max, val_fin, lim):
    ok = val_max <= lim
    print("    %-28s max %.3e  final %.3e  (<= %.0e)  %s" % (label, val_max, val_fin, lim, "PASS" if ok else "FAIL"))
    return ok


def report(run, kind):
    a, how = series(run)
    print("\n## %s  [%s, %d snapshots, t_end=%.3f, 体積重み=%s]" % (run, kind, len(a), a[-1, 1], how))
    pm = np.max(np.abs(a[:, 6:9]), axis=1)
    ok = True
    if kind == "keep":
        ok &= judge("|K/K0-1|", np.max(np.abs(a[:, 2] - 1)), abs(a[-1, 2] - 1), 1e-2)
        ok &= judge("|dS/S0|", np.max(np.abs(a[:, 3])), abs(a[-1, 3]), 1e-4)
        ok &= judge("|P_i|/(rho0 U0 V)", np.max(pm), pm[-1], 1e-6)
    else:
        ok &= judge("|M-M0|/M0", np.max(np.abs(a[:, 5])), abs(a[-1, 5]), 1e-6)
        ok &= judge("|P_i|/(rho0 U0 V)", np.max(pm), pm[-1], 1e-6)
        ok &= judge("|E-E0|/E0", np.max(np.abs(a[:, 4])), abs(a[-1, 4]), 1e-5)
        print("    (記録) K/K0 final %.6f, dS/S0 final %.3e" % (a[-1, 2], a[-1, 3]))
    print("  VERDICT(%s): %s" % (os.path.basename(run.rstrip("/")), "PASS" if ok else "FAIL"))
    return a


def diff(tag, a_old, a_new):
    n = min(len(a_old), len(a_new))
    if not np.array_equal(a_old[:n, 0], a_new[:n, 0]):
        print("  %s: snapshot の step が旧新で揃っていない" % tag)
        return
    dk = np.abs(a_new[:n, 2] - a_old[:n, 2])
    ds = np.abs(a_new[:n, 3] - a_old[:n, 3])
    i = int(np.argmax(dk))
    print("  %s 旧新差 (記録のみ): max|ΔK/K0| %.3e (t=%.2f)、final %.3e; max|Δ(dS/S0)| %.3e、final %.3e"
          % (tag, dk[i], a_old[i, 1], dk[-1], np.max(ds), ds[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", nargs=2, metavar=("OLD", "NEW"))
    ap.add_argument("--slau", nargs=2, metavar=("OLD", "NEW"))
    ap.add_argument("--png", default=None)
    a = ap.parse_args()
    data = {}
    for kind, pair in (("keep", a.keep), ("slau", a.slau)):
        if not pair:
            continue
        ao, an = report(pair[0], kind), report(pair[1], kind)
        diff(kind.upper(), ao, an)
        data[kind] = (pair, ao, an)
    if a.png and data:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for kind, (pair, ao, an) in data.items():
            for lab, arr, ls in (("old 1266aba1", ao, "--"), ("new 565959c7", an, "-")):
                nm = "%s %s (%s)" % (kind.upper(), lab, os.path.basename(pair[0 if ls == "--" else 1].rstrip("/")))
                ax[0].plot(arr[:, 1], arr[:, 2], ls, lw=1.4, label=nm)
                ax[1].plot(arr[:, 1], arr[:, 3], ls, lw=1.4, label=nm)
        ax[0].set_xlabel("t"); ax[0].set_ylabel("K/K0"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=7)
        ax[1].set_xlabel("t"); ax[1].set_ylabel("(S-S0)/|S0|"); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=7)
        fig.suptitle("R2 Taylor-Green 32^3 node: old vs new gradient (t_c=2.5, end 10 t_c)")
        fig.tight_layout()
        fig.savefig(a.png, dpi=110)
        print("\nsaved", a.png)


if __name__ == "__main__":
    main()
