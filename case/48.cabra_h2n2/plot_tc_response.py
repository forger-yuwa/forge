#!/usr/bin/env python3
"""Cabra T_c 応答曲線: 各 run の火炎基部 H/d (条件付き T_Q>1300 K の最上流 x/d と 平均 T>1200 K の最上流 x/d) を末尾スナップショットの
平均±範囲で集計し、実験 (T_c=1045 K で H/d≈10) と並べて図にする。  python3 plot_tc_response.py out.png run_dir:Tc [run_dir:Tc ...] [--tail 3]"""
import sys, glob, re, argparse, h5py, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("runs", nargs="+"); ap.add_argument("--tail", type=int, default=3); a = ap.parse_args()
D = 0.00457; rows = []
for spec in a.runs:
    run, tc = spec.rsplit(":", 1); tc = float(tc)
    with h5py.File(f"{run}/cabra.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
    fs = sorted([f for f in glob.glob(f"{run}/res_*.h5") if re.fullmatch(r".*/res_\d+\.h5", f)], key=lambda p: int(re.findall(r"res_(\d+)", p)[0]))[-a.tail:]
    hq, hm = [], []
    for f in fs:
        with h5py.File(f) as h: V = h["VALUE"]; T = V["T"][:]; TQ = V["cmc_TQmax"][:]
        N = len(T); x = xyz[:N, 0]; y = xyz[:N, 1]; ok = (x > 0.5 * D)   # リップ壁ノード (Dirichlet で初期場が残る) を除外
        s1 = ok & (TQ > 1300); s2 = ok & (T > 1200)
        hq.append(x[s1].min() / D if s1.any() else np.nan); hm.append(x[s2].min() / D if s2.any() else np.nan)
    rows.append((tc, run, np.nanmean(hq), np.nanmin(hq), np.nanmax(hq), np.nanmean(hm), np.nanmin(hm), np.nanmax(hm), [int(re.findall(r"res_(\d+)", f)[0]) for f in fs]))
rows.sort()
print(" T_c   H/d(T_Q>1300) [min..max]   H/d(T>1200) [min..max]   steps   run")
for r in rows: print(f"{r[0]:5.0f}  {r[2]:5.1f} [{r[3]:.1f}..{r[4]:.1f}]        {r[5]:5.1f} [{r[6]:.1f}..{r[7]:.1f}]      {r[8]}  {r[1]}")
fig, ax = plt.subplots(figsize=(6, 4))
tcs = [r[0] for r in rows]
ax.errorbar(tcs, [r[2] for r in rows], yerr=[[r[2] - r[3] for r in rows], [r[4] - r[2] for r in rows]], fmt="o-", capsize=3, label="forge CMC: T_Q>1300 K")
ax.errorbar(tcs, [r[5] for r in rows], yerr=[[r[5] - r[6] for r in rows], [r[7] - r[5] for r in rows]], fmt="s--", capsize=3, label="forge CMC: mean T>1200 K")
ax.plot([1045], [10], "k*", ms=12, label="exp (Cabra 2002) H/d≈10 @1045 K")
# 文献の感度目安: Cabra 実験 (Cabra 2002, Fig.): T_c 1030 K で H/d≈15〜20, 1060 K で ≈5〜6 程度 (目視; ±30 K 調整が査読標準)
ax.set_xlabel("coflow T_c [K]"); ax.set_ylabel("lift-off height H/d"); ax.set_ylim(0, 40); ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(a.out, dpi=120); print("fig:", a.out)
