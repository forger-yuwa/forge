#!/usr/bin/env python3
r"""残差の床が単精度の丸めかどうかを**次元の合った量**で測る (codex result-2 M4)。

`rms_ro` は密度ではなく「CV を出入りする質量流束の和」[kg/s] なので、密度で割っても無次元にならない。ここでは節点ごとに
  丸めの目安  e_i = eps32 · Σ_f ( |ρ_f u_f·S_f| + p_f |S_f| / c_f )
    第 1 項: 面の質量流束そのものの丸め。第 2 項: SLAU の質量流束は圧力差の散逸項 Δp/c̄·S を含み、Δp = p_R − p_L は
    **絶対圧どうしの引き算**なので eps·p の誤差を持つ (この case は pRef ゲージ無し)。速度が 0 に近い壁近傍では第 2 項が支配的。
を作り、level 2 出力の残差場 |res_ro_i| と比べる。流束の和を float32 で積むと、和の各項の大きさ × eps 程度の誤差は避けられない。
比 |res_ro_i| / e_i が O(1) なら、その節点の質量残差は丸めで説明できる大きさ。エネルギーも同様に Σ|ρ u·S| · H で見る。
面の値は両側の節点値の算術平均 (目安なので再構成やリミッタは入れない)。

usage: rounding_floor.py RUN_DIR      (output.level 2 の run)
"""
import glob, os, sys
import h5py, numpy as np
run = sys.argv[1]
f = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda p: int(os.path.basename(p)[4:-3]))[-1]
with h5py.File(os.path.join(run, "mesh.h5"), "r") as m:
    st = m["PLANES/STRUCT"][:]; sv = m["PLANES/surfVect"][:].reshape(-1, 3); ncell = len(m["CELLS/volume"])
with h5py.File(f, "r") as h:
    V = {k: h["VALUE"][k][:].astype(float) for k in ("ro", "Ux", "Uy", "Uz", "h0", "P", "sonic", "res_ro", "res_roe")}
c0 = np.empty(len(sv), int); c1 = np.empty(len(sv), int); p = 0
for ip in range(len(sv)):
    nn = st[p]; p += 1 + nn; nc = st[p]; c0[ip] = st[p + 1]; c1[ip] = st[p + 2] if nc > 1 else -1; p += 1 + nc
inner = (c1 >= 0) & (c1 < ncell) & (c0 < ncell)
a, b = c0[inner], c1[inner]; S = sv[inner]
rf = 0.5 * (V["ro"][a] + V["ro"][b]); uf = 0.5 * (np.c_[V["Ux"][a], V["Uy"][a], V["Uz"][a]] + np.c_[V["Ux"][b], V["Uy"][b], V["Uz"][b]])
mf = np.abs(rf * np.einsum("ij,ij->i", uf, S)) + 0.5 * (V["P"][a] / V["sonic"][a] + V["P"][b] / V["sonic"][b]) * np.linalg.norm(S, axis=1); hf = 0.5 * (V["h0"][a] + V["h0"][b])
sm = np.zeros(ncell); se = np.zeros(ncell)
np.add.at(sm, a, mf); np.add.at(sm, b, mf); np.add.at(se, a, mf * hf); np.add.at(se, b, mf * hf)
eps = 1.1920929e-7
print(f"{f}  ({inner.sum()} 内部双対面, {ncell} 節点)")
for name, res, s in (("質量 res_ro", V["res_ro"], sm), ("エネルギー res_roe", V["res_roe"], se)):
    ok = s > 0; r = np.abs(res[ok]) / (eps * s[ok])
    print(f"  {name:<18} rms|res| = {np.sqrt(np.mean(res**2)):.3e}   rms(eps·Σ(|ρu·S|+p|S|/c)) = {np.sqrt(np.mean((eps*s[ok])**2)):.3e}   比 |res|/(eps·Σ(|ρu·S|+p|S|/c)): 中央値 {np.median(r):.2f}, 90 % 点 {np.percentile(r, 90):.2f}, 99 % 点 {np.percentile(r, 99):.2f}, 最大 {r.max():.1f}")
