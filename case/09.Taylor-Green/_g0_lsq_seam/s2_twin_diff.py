#!/usr/bin/env python3
"""S2 双子の場の差と共通ゲート (plan gradient-scalar-lsq-unification §6 S2)。

  python3 s2_twin_diff.py GG_RUN LSQ_RUN [--step N] [--fields ro roUy roY0 ...] [--mesh NAME]

- 同じ step の res_N.h5 (既定 = 両方にある最大の N) で、各 VALUE の最大絶対差・場の最大値に対する相対差・位置を出す。
- 共通ゲート: NaN/Inf の有無、Y 系があれば max|ΣY−1| と min Y、Xi があれば [min, max]、凝縮モーメントの最小値。
"""
import argparse
import glob
import os
import re

import h5py
import numpy as np


def steps(run):
    return sorted(int(re.search(r"res_(\d+)\.h5$", p).group(1)) for p in glob.glob(os.path.join(run, "res_*.h5"))
                  if re.search(r"/res_\d+\.h5$", p))


def load(p):
    with h5py.File(p, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def gates(tag, V, n):
    out = []
    bad = [k for k, a in V.items() if a.dtype.kind == "f" and not np.all(np.isfinite(a[:n]))]
    out.append(f"{tag}: NaN/Inf 配列 {bad if bad else 'なし'}")
    ys = sorted(k for k in V if re.fullmatch(r"Y\d+", k))
    if ys:
        s = sum(V[k][:n].astype(np.float64) for k in ys)
        out.append(f"{tag}: max|ΣY−1| = {np.abs(s - 1).max():.3e}、min Y = {min(V[k][:n].min() for k in ys):.3e}")
    if "Xi" in V:
        out.append(f"{tag}: Xi ∈ [{V['Xi'][:n].min():.3e}, {V['Xi'][:n].max():.6f}]")
    mom = sorted(k for k in V if k.startswith("roQ") or k.startswith("rog_"))
    if mom:
        out.append(f"{tag}: モーメント最小 " + ", ".join(f"{k} {V[k][:n].min():.3e}" for k in mom))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gg"); ap.add_argument("lsq")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--fields", nargs="*", default=None)
    a = ap.parse_args()
    st = a.step if a.step is not None else max(set(steps(a.gg)) & set(steps(a.lsq)))
    G = load(os.path.join(a.gg, f"res_{st}.h5")); L = load(os.path.join(a.lsq, f"res_{st}.h5"))
    n = min(len(G["ro"]), len(L["ro"]))
    print(f"# 双子の差 lsq − gg @ res_{st}\n- gg : {a.gg}\n- lsq: {a.lsq}")
    for line in gates("gg", G, n) + gates("lsq", L, n):
        print("- " + line)
    names = a.fields or sorted(k for k in G if k in L and G[k].dtype.kind == "f")
    print("| 場 | max\\|Δ\\| | max\\|Δ\\| / max\\|gg\\| | 位置 (index) |")
    print("| --- | --- | --- | --- |")
    for k in names:
        g = G[k][:n].astype(np.float64); l = L[k][:n].astype(np.float64)
        d = np.abs(l - g)
        i = int(np.nanargmax(d)) if d.size else -1
        sc = np.nanmax(np.abs(g)) if d.size else 0.0
        print(f"| {k} | {d[i]:.3e} | {d[i] / sc if sc > 0 else float('nan'):.3e} | {i} |")


if __name__ == "__main__":
    main()
