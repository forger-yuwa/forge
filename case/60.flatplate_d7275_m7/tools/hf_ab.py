#!/usr/bin/env python3
"""T4-0b-HF (SLAU 面エンタルピー float/double A/B) の集計 (case/60)。

    python3 tools/hf_ab.py RUN_A RUN_B [--windows 3000:4000,4000:5000]

tf_ab.py との違い (codex 2026-09-27): 領域別の √Σres² を最終 1 枚でなく保存した全スナップショット
(100 step 毎) の時系列で出し、窓ごとの中央値・窓間変化・B/A を示す。点プローブは 1 step 毎の連続反復で振幅を取る。
"""
import argparse, glob, re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py


def snaps(run):
    ss = sorted(int(re.search(r"res_(\d+)\.h5$", f).group(1)) for f in glob.glob(str(Path(run) / "res_[0-9]*.h5")))
    return [(s, Path(run) / f"res_{s}.h5") for s in ss]


def hist(run):
    d = pd.read_csv(Path(run) / "residual_history.csv")
    if d.shape[1] > 2 and d.iloc[:, 2].dtype == object:
        d = d[d.iloc[:, 2] == "outer_end"]
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--windows", default="3000:4000,4000:5000")
    a = ap.parse_args()
    wins = [tuple(int(v) for v in w.split(":")) for w in a.windows.split(",")]
    SA = snaps(a.A)
    with h5py.File(SA[-1][1], "r") as h:
        X = np.asarray(h["MESH/COORD"]).reshape(-1, 3)
        keys = sorted(k for k in h["VALUE"] if k.startswith("res_"))
    x, y = X[:, 0], X[:, 1]
    y1 = np.min(y[y > 0])
    reg = {"自由流域": (y >= 0.03) & (y <= 0.79), "第一内部列": np.isclose(y, y1), "全域": np.ones_like(x, bool)}
    rows = []
    for nm, run in (("A", a.A), ("B", a.B)):
        for s, f in snaps(run):
            with h5py.File(f, "r") as h:
                for k in keys:
                    r = np.asarray(h["VALUE"][k], float)
                    for rn, m in reg.items():
                        rows.append((nm, s, k, rn, float(np.sqrt((r[m] ** 2).sum()))))
    df = pd.DataFrame(rows, columns=["arm", "step", "field", "region", "val"])
    print(f"スナップショット: A {df[df.arm=='A'].step.nunique()} 枚 / B {df[df.arm=='B'].step.nunique()} 枚")
    for rn in reg:
        print(f"\n[{rn} の √Σres²: 窓ごとの中央値、B/A、窓間変化 %]")
        out = []
        for k in keys:
            line = {"field": k}
            for lo, hi in wins:
                for arm in ("A", "B"):
                    v = df[(df.arm == arm) & (df.field == k) & (df.region == rn) & (df.step >= lo) & (df.step < hi if hi != wins[-1][1] else df.step <= hi)].val
                    line[f"{arm}{lo//1000}-{hi//1000}k"] = v.median()
            lo, hi = wins[-1]; plo, phi = wins[-2]
            line["B/A(末尾)"] = line[f"B{lo//1000}-{hi//1000}k"] / line[f"A{lo//1000}-{hi//1000}k"]
            for arm in ("A", "B"):
                line[f"Δ{arm}%"] = (line[f"{arm}{lo//1000}-{hi//1000}k"] / line[f"{arm}{plo//1000}-{phi//1000}k"] - 1) * 100
            out.append(line)
        print(pd.DataFrame(out).set_index("field").to_string(float_format=lambda v: f"{v:.3g}"))
    # 全域 RMS 履歴
    hA, hB = hist(a.A), hist(a.B)
    cols = [c for c in hA.columns if c.startswith("rms_") and not c.startswith("rms_dq") and c != "rms_roUz"]
    lo, hi = wins[-1]
    mA = hA[(hA.iloc[:, 0] >= lo) & (hA.iloc[:, 0] < hi)][cols].mean()
    mB = hB[(hB.iloc[:, 0] >= lo) & (hB.iloc[:, 0] < hi)][cols].mean()
    print(f"\n[全域 RMS 履歴 step {lo}–{hi}]\n" + pd.DataFrame({"A": mA, "B": mB, "B/A": mB / mA}).to_string(float_format=lambda v: f"{v:.3g}"))
    # 点プローブ (1 step 毎)
    print(f"\n[点プローブ step {lo}–{hi} (連続反復): 振幅 (max−min)/|平均|  A → B  (出力 6 桁、分解能 ~1e-6)]")
    for i in range(6):
        vals = []
        for run in (a.A, a.B):
            p = Path(run) / f"point_probe_{i}.out"
            d = pd.read_csv(p, sep=r"[\s,]+", engine="python")
            d = d[(d.iloc[:, 0] >= lo) & (d.iloc[:, 0] < hi)]
            vals.append({c: (d[c].max() - d[c].min()) / abs(d[c].mean()) for c in ("T", "P")})
        print(f"  probe {i}: T {vals[0]['T']:.2e}→{vals[1]['T']:.2e}  P {vals[0]['P']:.2e}→{vals[1]['P']:.2e}")


if __name__ == "__main__":
    main()
