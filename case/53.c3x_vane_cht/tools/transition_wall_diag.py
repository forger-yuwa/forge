#!/usr/bin/env python3
r"""遷移モデルつき翼 run の「どこで・どうやって乱流になったか」を壁に沿って出す診断。

壁節点ごとに、壁から `--band` [m] 以内の場の節点を最近傍で割り当て、その中の $\gamma_{eff}$ **最小** (境界層の中が層流なら小さい。帯が自由流 ($\gamma=1$) を含んでも影響されない) と $\mu_t/\mu$ 最大を取る。
帯は境界層より薄く取ること (既定 0.1 mm。厚すぎると自由流の $\mu_t/\mu$ を拾う)。
壁せん断の向き (流れ方向の符号) から剥離域も出す。横軸は `compare_h.py` と同じ弧長 $s/S$ (負圧面 +、正圧面 −)。

usage: transition_wall_diag.py RUN_DIR [--band 5e-4] [--phys 5] [--csv out.csv]
"""
import argparse, glob, os, sys
from pathlib import Path
import numpy as np, h5py
from scipy.spatial import cKDTree
sys.path.insert(0, str(Path(__file__).parent))
from compare_h import arc_map

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--band", type=float, default=1e-4); ap.add_argument("--phys", type=int, default=5); ap.add_argument("--csv", default=None)
a = ap.parse_args()
fw = sorted(glob.glob(os.path.join(a.run, f"res_wall_{a.phys}_[0-9]*.h5")), key=lambda p: int(p[:-3].split("_")[-1]))[-1]
fv = sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda p: int(os.path.basename(p)[4:-3]))[-1]
with h5py.File(fw, "r") as h:
    Cw = h["MESH/COORD"][:].reshape(-1, 3)[:, :2]; tx, ty = h["VALUE/twall_x"][:], h["VALUE/twall_y"][:]
with h5py.File(fv, "r") as h:
    C = h["MESH/COORD"][:].reshape(-1, 3)[:, :2]; V = {k: h["VALUE"][k][:].astype(float) for k in ("gammaEff", "gammaTr", "vis_turb", "vis_lam", "wall_dist", "k", "Ux", "Uy") if k in h["VALUE"]}
s, ss = arc_map(Cw); x = np.where(ss, s, -s)
m = (V["wall_dist"] > 0) & (V["wall_dist"] < a.band)
_, idx = cKDTree(Cw).query(C[m])
ge = np.full(len(Cw), 2.0); mr = np.zeros(len(Cw))
np.minimum.at(ge, idx, V.get("gammaEff", np.ones_like(V["wall_dist"]))[m]); np.maximum.at(mr, idx, (V["vis_turb"] / V["vis_lam"])[m])
# 壁せん断の流れ方向成分: 前縁から後縁へ向かう接線 (負圧面・正圧面それぞれ s が増える向き) との内積
o = np.argsort(x); xs = x[o]; t = np.gradient(Cw[o], axis=0); t /= np.linalg.norm(t, axis=1)[:, None] + 1e-30
sgn = np.where(xs >= 0, 1.0, -1.0)[:, None]            # x (= ±s/S) の増加方向は負圧面で下流、正圧面で上流
tau_s = (np.stack([tx, ty], 1)[o] * t * sgn).sum(1)
ge, mr = ge[o], mr[o]
print(f"{a.run}  ({os.path.basename(fv)}, band {a.band*1e3:.2f} mm)")
for name, lo, hi in (("正圧面", -0.9, -0.05), ("負圧面", 0.05, 0.9)):
    r = (xs > lo) & (xs < hi); xr = np.abs(xs[r]); k = np.argsort(xr); xr = xr[k]; g = ge[r][k]; q = mr[r][k]; tt = tau_s[r][k]
    on = xr[np.argmax((g > 0.5) & (xr > 0.08))] if ((g > 0.5) & (xr > 0.08)).any() else float("nan"); sep = xr[tt < 0]
    print(f"  {name}: 層内の gamma_eff 最小が 0.5 を超える s/S = {on:.3f} | mu_t/mu>10 になる s/S = {xr[np.argmax((q > 10) & (xr > 0.08))] if ((q > 10) & (xr > 0.08)).any() else float('nan'):.3f} | 逆流 (壁せん断が負) の範囲 s/S = " + (f"{sep.min():.3f}–{sep.max():.3f}" if len(sep) else "なし") + f" | 層内 gamma_eff 最小の範囲 {g[xr>0.08].min():.3f}–{g[xr>0.08].max():.2f}")
if a.csv:
    np.savetxt(a.csv, np.c_[xs, ge, mr, tau_s], delimiter=",", header="s_signed,gamma_eff_min,mut_over_mu_max,tau_streamwise", comments="")
