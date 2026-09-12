#!/usr/bin/env python3
"""Arthur ノズル凝縮 run の onset 解析 (plan condensation-air §4.4/§6):
  - dry 基準 (restart_dry.h5 の保存量から p_dry = (γ−1)(ρe − ½ρ|u|²)) との中心線 p/p_dry − 1 が 1 % を超える最初の x を onset とし、
    その点の dry 状態 (p, T) と、g>1e-4 の最初の x を出す。
  - 膨張率 Ṗ = −(u/p) dp/dx [1/s] を dry 中心線で評価 (Daum & Gyarmathy の相関パラメータ)。
  - Daum & Gyarmathy 最小 onset 曲線 (daum_gyarmathy_min_onset_n2.csv) の T_onset(p) と比べた ΔT を出す。
usage: onset_analysis.py RUN_DIR [--res res_N.h5] [--dry restart_dry.h5] [--gamma 1.4] [--R 296.8]"""
import sys, os, glob, argparse, h5py, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--res", default=None); ap.add_argument("--dry", default="restart_dry.h5")
ap.add_argument("--gamma", type=float, default=1.4); ap.add_argument("--R", type=float, default=296.8); ap.add_argument("--thr", type=float, default=0.01)
a = ap.parse_args()
def latest(run):
    fs = [f for f in glob.glob(os.path.join(run, "res_*.h5")) if "nan" not in f]; return max(fs, key=lambda f: int(os.path.basename(f)[4:-3]))
res = a.res or latest(a.run); f = h5py.File(res, "r"); V = f["VALUE"]
mesh = h5py.File(os.path.join(a.run, "arthur_nozzle.h5"), "r"); c = np.array(mesh["CELLS/centCoords"]).reshape(-1, 3)[:len(V["P"])]
d = h5py.File(os.path.join(a.run, a.dry), "r")["VALUE"]
ro = np.array(d["ro"], float); u = np.array(d["roUx"], float)/ro; v = np.array(d["roUy"], float)/ro; w = np.array(d["roUz"], float)/ro; roe = np.array(d["roe"], float)
p_dry = (a.gamma - 1)*(roe - 0.5*ro*(u*u + v*v + w*w)); T_dry = p_dry/(ro*a.R)
P = np.array(V["P"], float); g = np.array(V["g_0"], float) if "g_0" in V else np.zeros_like(P)
# 中心線: source-flow メッシュはセル中心 x が列ごとに揺れるので x をビン分け (~400 列) し、各ビンで |y| 最小のセルを取る
nb = 400; edges = np.linspace(c[:, 0].min(), c[:, 0].max() + 1e-9, nb + 1); ib = np.clip(np.digitize(c[:, 0], edges) - 1, 0, nb - 1)
C = np.array([cand[np.argmin(np.abs(c[cand, 1]))] for cand in (np.where(ib == k)[0] for k in range(nb)) if len(cand)])
order = np.argsort(c[C, 0]); C = C[order]
x = c[C, 0]; pd = p_dry[C]; Td = T_dry[C]; ud = u[C]; ratio = P[C]/pd - 1.0; gc = g[C]
Pdot = -ud/pd*np.gradient(pd, x)   # [1/s] (膨張で正)
def first(mask):
    i = np.where(mask)[0]; return int(i[0]) if len(i) else None
io = first(ratio > a.thr); ig = first(gc > 1e-4)
on = np.genfromtxt(os.path.join(a.run, "..", "daum_gyarmathy_min_onset_n2.csv"), delimiter=",", comments="#", skip_header=1)
def T_dg(p): return np.interp(np.log10(p), np.log10(on[:, 0]), on[:, 1])
print(f"res={os.path.basename(res)}  dry={a.dry}  cells on centerline={len(C)}")
for lab, i in (("Δp/p_dry > %.0f %%" % (a.thr*100), io), ("g > 1e-4", ig)):
    if i is None: print(f"  onset ({lab}): not reached"); continue
    print(f"  onset ({lab}): x={x[i]*1e3:.2f} mm = {x[i]/0.0254:.2f} in,  p_dry={pd[i]:.1f} Pa, T_dry={Td[i]:.2f} K, Pdot={Pdot[i]:.3g} 1/s;  D&G min onset T({pd[i]:.0f} Pa)={T_dg(pd[i]):.1f} K -> ΔT = {Td[i]-T_dg(pd[i]):+.1f} K (正=D&G より高温=過冷却が浅い)")
print(f"  exit: x={x[-1]/0.0254:.2f} in  p_dry={pd[-1]:.1f} Pa T_dry={Td[-1]:.1f} K  g_exit={gc[-1]:.4f}  max(p/p_dry-1)={ratio.max():.3f}  Pdot range along nozzle [{Pdot[np.isfinite(Pdot)].min():.3g}, {Pdot[np.isfinite(Pdot)].max():.3g}]")
