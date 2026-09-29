#!/usr/bin/env python3
"""T4-0b-TF (thermoFloat 1/0 A/B) の集計 (case/60)。

    python3 tools/tf_ab.py RUN_A RUN_B [--windows 3000:4000,4000:5000]

1. 残差履歴: 窓ごとの RMS 平均と B/A 比 (同じメッシュなので √Σres² 比と同じ)、窓間の中央値変化。
2. 残差場 (最終 res_*.h5): 自由流域 (y 0.03–0.79 m) と第一内部列 (y = y1) の成分別 √Σres²、B/A。
3. 固定節点 (A の第一内部列 res_roOmega 上位 5 点 + 自由流 3 点) の ρ・P・roK・roOmega の 100 step 毎時系列の振幅。
4. 第一内部列の ω 収支 omg_trans + omg_prod − omg_dest + omg_cross と res_roOmega (最終場、上位 5 点)。
"""
import argparse, glob, re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py


def hist(run):
    d = pd.read_csv(Path(run) / "residual_history.csv")
    if d.shape[1] > 2 and d.iloc[:, 2].dtype == object:
        d = d[d.iloc[:, 2] == "outer_end"]
    return d


def snaps(run):
    ss = sorted(int(re.search(r"res_(\d+)\.h5$", f).group(1)) for f in glob.glob(str(Path(run) / "res_[0-9]*.h5")))
    return [(s, Path(run) / f"res_{s}.h5") for s in ss]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--windows", default="3000:4000,4000:5000")
    a = ap.parse_args()
    wins = [tuple(int(v) for v in w.split(":")) for w in a.windows.split(",")]
    hA, hB = hist(a.A), hist(a.B)
    cols = [c for c in hA.columns if c.startswith("rms_") and not c.startswith("rms_dq") and c != "rms_roUz"]
    for lo, hi in wins:
        mA = hA[(hA.iloc[:, 0] >= lo) & (hA.iloc[:, 0] < hi)][cols].mean()
        mB = hB[(hB.iloc[:, 0] >= lo) & (hB.iloc[:, 0] < hi)][cols].mean()
        print(f"\n[残差履歴 step {lo}–{hi} の平均]\n" + pd.DataFrame({"A": mA, "B": mB, "B/A": mB / mA}).to_string(float_format=lambda v: f"{v:.3g}"))
    w0, w1 = wins[-2], wins[-1]
    for nm, h in (("A", hA), ("B", hB)):
        m0 = h[(h.iloc[:, 0] >= w0[0]) & (h.iloc[:, 0] < w0[1])][cols].median()
        m1 = h[(h.iloc[:, 0] >= w1[0]) & (h.iloc[:, 0] < w1[1])][cols].median()
        print(f"[{nm}] 窓間の中央値変化 %: " + "  ".join(f"{c[4:]} {v:+.1f}" for c, v in ((m1 / m0 - 1) * 100).items()))

    SA, SB = snaps(a.A), snaps(a.B)
    with h5py.File(SA[-1][1], "r") as h:
        X = np.asarray(h["MESH/COORD"]).reshape(-1, 3)
    x, y = X[:, 0], X[:, 1]
    y1 = np.min(y[y > 0])
    reg = {"自由流域 (y 0.03–0.79)": (y >= 0.03) & (y <= 0.79), "第一内部列 (y = y1)": np.isclose(y, y1),
           "壁近傍 (y1 < y < 1 mm)": (y > y1 * 1.01) & (y < 1e-3), "全域": np.ones_like(x, bool)}
    print(f"\n[残差場 √Σres² (最終 step A {SA[-1][0]} / B {SB[-1][0]})]  y1 = {y1:.3e}")
    fA, fB = h5py.File(SA[-1][1], "r"), h5py.File(SB[-1][1], "r")
    ks = sorted(k for k in fA["VALUE"] if k.startswith("res_") and k in fB["VALUE"])
    for k in ks:
        rA, rB = np.asarray(fA["VALUE"][k], float), np.asarray(fB["VALUE"][k], float)
        s = "  ".join(f"{n}: B/A {np.sqrt((rB[m]**2).sum())/np.sqrt((rA[m]**2).sum()):.3g}" for n, m in reg.items())
        print(f"  {k:12s} A 全域 {np.sqrt((rA**2).sum()):.3e}  {s}")
    # 固定節点
    rO = np.asarray(fA["VALUE"]["res_roOmega"], float)
    row = np.where(np.isclose(y, y1))[0]
    top = row[np.argsort(abs(rO[row]))[::-1][:5]]
    def near(px, py):
        return int(np.argmin((x - px) ** 2 + (y - py) ** 2))
    fixed = [(f"壁列 x={x[i]:.3f}", i) for i in top] + [(f"自由流 ({px},{py})", near(px, py)) for px, py in ((2.156, 0.7498), (1.0, 0.6), (2.0, 0.4))]
    print(f"\n[固定節点の時系列 (step {SA[0][0]}–{SA[-1][0]}、{len(SA)} 枚): 振幅 (max−min)/|平均|  A → B]")
    for nm, i in fixed:
        out = []
        for q in ("ro", "P", "roK", "roOmega"):
            amp = []
            for S in (SA, SB):
                v = np.array([h5py.File(f, "r")["VALUE"][q][i] for _, f in S], float)
                amp.append((v.max() - v.min()) / abs(v.mean()) if v.mean() else float("nan"))
            out.append(f"{q} {amp[0]:.2e}→{amp[1]:.2e}")
        print(f"  {nm:26s} " + "  ".join(out))
    # ω 収支
    print("\n[第一内部列の ω 収支 (最終場、A の res_roOmega 上位 5 点)]")
    for nm, f in (("A", fA), ("B", fB)):
        V = f["VALUE"]
        if "omg_trans" not in V:
            print(f"  {nm}: omg_* なし"); continue
        for i in top:
            t, p, d_, c = (float(V[k][i]) for k in ("omg_trans", "omg_prod", "omg_dest", "omg_cross"))
            print(f"  {nm} x={x[i]:.3f}: trans {t:+.4e} prod {p:+.4e} dest {d_:+.4e} cross {c:+.4e}  和 {t+p-d_+c:+.3e}  res_roOmega {float(V['res_roOmega'][i]):+.3e}")


if __name__ == "__main__":
    main()
