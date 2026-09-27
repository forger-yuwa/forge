#!/usr/bin/env python3
"""T4-0b-H (上端の高さ A/B) の集計 (case/60)。

    python3 tools/height_ab.py RUN_A RUN_B [--windows 12000:16000,16000:20000]

1. 残差履歴: 窓ごとの RMS 平均と √Σres² (= RMS × √N、N は DOF 数) の平均、B/A 比。
   RMS は節点数で割るので、節点を足した B は同じ Σres² でも RMS が小さく出る (codex 2026-09-27)。
2. 残差場 (最終 res_*.h5 の res_*): 領域別の √Σres² と最大値、B/A 比。
3. 点プローブ (point_probe_*.out): 末尾窓の P・T の振幅 (max−min) と標準偏差、B/A 比。
"""
import argparse, glob, re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py

REGIONS = {
    "旧上端付近 (x 2.6–2.8, y 0.45–0.50)": lambda x, y: (x >= 2.6) & (y >= 0.45) & (y <= 0.5 + 1e-9),
    "B の追加域 (y > 0.5)": lambda x, y: y > 0.5 + 1e-9,
    "新上端 (y ≥ 0.76)": lambda x, y: y >= 0.76,
    "出口列 (x ≥ 2.79)": lambda x, y: x >= 2.79,
    "比較域 (x 1.07–2.55)": lambda x, y: (x >= 1.07) & (x <= 2.55),
    "壁近傍 x≥2.6 (y < 1 mm)": lambda x, y: (x >= 2.6) & (y < 1e-3),
    "全域": lambda x, y: np.ones_like(x, bool),
}


def hist(run):
    d = pd.read_csv(Path(run) / "residual_history.csv")
    if d.shape[1] > 2 and d.iloc[:, 2].dtype == object:
        d = d[d.iloc[:, 2] == "outer_end"]
    return d


def ndof(run):
    with h5py.File(Path(run) / "mesh.h5", "r") as h:
        return len(h["VALUE/ro"])


def last_res(run):
    s = max(int(re.search(r"res_(\d+)\.h5$", f).group(1)) for f in glob.glob(str(Path(run) / "res_[0-9]*.h5")))
    return s, Path(run) / f"res_{s}.h5"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--windows", default="12000:16000,16000:20000")
    a = ap.parse_args()
    wins = [tuple(int(v) for v in w.split(":")) for w in a.windows.split(",")]
    hA, hB = hist(a.A), hist(a.B)
    nA, nB = ndof(a.A), ndof(a.B)
    cols = [c for c in hA.columns if c.startswith("rms_") and not c.startswith("rms_dq") and c != "rms_roUz"]
    print(f"DOF: A {nA}  B {nB}  (√N 比 B/A {np.sqrt(nB/nA):.4f})")
    for lo, hi in wins:
        mA = hA[(hA.iloc[:, 0] >= lo) & (hA.iloc[:, 0] < hi)][cols].mean()
        mB = hB[(hB.iloc[:, 0] >= lo) & (hB.iloc[:, 0] < hi)][cols].mean()
        t = pd.DataFrame({"RMS_A": mA, "RMS_B": mB, "RMS B/A": mB / mA,
                          "√Σres² B/A": (mB * np.sqrt(nB)) / (mA * np.sqrt(nA))})
        print(f"\n[残差履歴 step {lo}–{hi} の平均]\n" + t.to_string(float_format=lambda v: f"{v:.3g}"))
    # 窓間の中央値変化 (過渡判定)
    if len(wins) >= 2:
        for nm, h in (("A", hA), ("B", hB)):
            w0, w1 = wins[-2], wins[-1]
            m0 = h[(h.iloc[:, 0] >= w0[0]) & (h.iloc[:, 0] < w0[1])][cols].median()
            m1 = h[(h.iloc[:, 0] >= w1[0]) & (h.iloc[:, 0] < w1[1])][cols].median()
            ch = (m1 / m0 - 1) * 100
            print(f"\n[{nm}] 窓間の中央値変化 % ({w0}→{w1}): " + "  ".join(f"{c[4:]} {v:+.1f}" for c, v in ch.items()))
    # 残差場
    print("\n[残差場 √Σres² (最終スナップショット)]")
    out = {}
    for nm, run in (("A", a.A), ("B", a.B)):
        s, f = last_res(run)
        with h5py.File(f, "r") as h:
            X = np.asarray(h["MESH/COORD"]).reshape(-1, 3); V = h["VALUE"]
            ks = sorted(k for k in V if k.startswith("res_"))
            out[nm] = (s, {k: np.asarray(V[k], float) for k in ks}, X[:, 0], X[:, 1])
    sA, rA, xA, yA = out["A"]; sB, rB, xB, yB = out["B"]
    print(f"  step A {sA} / B {sB}")
    for k in sorted(set(rA) & set(rB)):
        print(f"  {k}")
        for reg, f in REGIONS.items():
            mA_, mB_ = f(xA, yA), f(xB, yB)
            vA = np.sqrt((rA[k][mA_] ** 2).sum()) if mA_.any() else float("nan")
            vB = np.sqrt((rB[k][mB_] ** 2).sum()) if mB_.any() else float("nan")
            print(f"    {reg:34s} A {vA:.3e}  B {vB:.3e}  B/A {vB/vA if vA else float('nan'):.3g}")
        iA, iB = np.argmax(abs(rA[k])), np.argmax(abs(rB[k]))
        print(f"    最大: A {abs(rA[k][iA]):.3e} @({xA[iA]:.3f},{yA[iA]:.4f})  B {abs(rB[k][iB]):.3e} @({xB[iB]:.3f},{yB[iB]:.4f})")
    # 点プローブ
    lo, hi = wins[-1]
    print(f"\n[点プローブ step {lo}–{hi}: 振幅 (max−min) / 平均]")
    for i in range(6):
        row = []
        for run in (a.A, a.B):
            p = Path(run) / f"point_probe_{i}.out"
            if not p.exists():
                row.append(None); continue
            d = pd.read_csv(p, sep=r"[\s,]+", engine="python", comment="#")
            d = d[(d.iloc[:, 0] >= lo) & (d.iloc[:, 0] < hi)]
            row.append(d)
        if row[0] is None or row[1] is None:
            continue
        s = []
        for c in [c for c in row[0].columns if c in ("P", "T")]:
            aA = (row[0][c].max() - row[0][c].min()) / abs(row[0][c].mean())
            aB = (row[1][c].max() - row[1][c].min()) / abs(row[1][c].mean())
            s.append(f"{c}: A {aA:.2e} B {aB:.2e} B/A {aB/aA if aA else float('nan'):.3g}")
        print(f"  probe {i}: " + "  ".join(s))


if __name__ == "__main__":
    main()
