#!/usr/bin/env python3
"""case/45 node 双対面: 集約面 r̄·ΣS_k と区間ごとの Σ r_k S_k の閉性比較 (codex 第 1 仮説の反証試験)。

変換器 (gmshReader.hpp buildMedianDual 2D) を再現する:
  - primal セル重心 G = 頂点平均 (gmshReader.hpp:1232-1244)
  - 内部双対面: エッジ ip=(A,B) の隣接セルごとに区間 M->G、法線 rotate(-90) を A->B 向きに揃えて和、
    重心は区間中点の区間長加重 (gmshReader.hpp:1475-1531)
  - 境界半割面: 各 bcond の境界エッジを両端ノードへ 1/2 ずつ、重心 (3N+O)/4 を長さ加重で集約
    (gmshReader.hpp:1641-1699)
  - ソルバ: S_r = S * max(pcy, 1e-20)、A_planar = dualVolume (variables.cpp:751-780)
"""
import sys
import numpy as np
import h5py

# usage: python3 dual_segment_closure.py <nozzle.h5 (FP64 の変換器の格子)> <段 ② のダンプ stage2.h5 (FP64)> <出力 .npz>
MESH_H5, DUMP_H5, OUT_NPZ = sys.argv[1], sys.argv[2], sys.argv[3]
RT = 0.076807
NJ = 121

f = h5py.File(MESH_H5, 'r')
X = f['MESH/COORD'][:].reshape(-1, 3)
nN = X.shape[0]
x, y = X[:, 0], X[:, 1]

# ---- primal quads (VIZMESH/CONNE を走査) ----
cw = f['VIZMESH/CONNE'][:]
quads = []
p = 0
while p < len(cw):
    code = cw[p]
    nn = {5: 4, 4: 3}[int(code)]
    quads.append(cw[p + 1:p + 1 + nn])
    p += 1 + nn
quads = np.array(quads, dtype=np.int64)
assert quads.shape[1] == 4, 'triangles present; script assumes all quads'
nQ = len(quads)
G = X[quads].mean(axis=1)  # 頂点平均

# ---- solver planes ----
ps = f['PLANES/STRUCT'][:]
nP = f['PLANES/surfArea'].shape[0]
SV = f['PLANES/surfVect'][:].reshape(-1, 3)
PC = f['PLANES/centCoords'][:].reshape(-1, 3)
pA = np.empty(nP, np.int64); pB = np.full(nP, -1, np.int64)
p = 0
for ip in range(nP):
    nn = ps[p]; nodes = ps[p + 1:p + 1 + nn]; p += 1 + nn
    nc = ps[p]; cells = ps[p + 1:p + 1 + nc]; p += 1 + nc
    pA[ip] = cells[0]
    if nc > 1:
        pB[ip] = cells[1]
nInt = f['MESH'].attrs['nNormalPlanes']
assert np.all(pB[:nInt] >= 0) and np.all(pB[nInt:] < 0)

# ---- 内部面: 区間の再構成 ----
# 各 quad の 4 エッジ (n_k, n_{k+1})
ea = quads.reshape(-1)                                  # 区間の端 a
eb = np.roll(quads, -1, axis=1).reshape(-1)             # 区間の端 b
gq = np.repeat(np.arange(nQ), 4)
lo = np.minimum(ea, eb); hi = np.maximum(ea, eb)
key = lo * nN + hi
# plane key -> ip
pkey = np.minimum(pA[:nInt], pB[:nInt]) * nN + np.maximum(pA[:nInt], pB[:nInt])
order = np.argsort(pkey)
pos = np.searchsorted(pkey[order], key)
assert np.all(pkey[order][pos] == key), 'edge not found among planes'
seg_ip = order[pos]
A = pA[seg_ip]; B = pB[seg_ip]
Mx = 0.5 * (x[A] + x[B]); My = 0.5 * (y[A] + y[B])
sx = G[gq, 0] - Mx; sy = G[gq, 1] - My
nx = sy.copy(); ny = -sx
eABx = x[B] - x[A]; eABy = y[B] - y[A]
flip = nx * eABx + ny * eABy < 0
nx[flip] *= -1; ny[flip] *= -1
L = np.hypot(nx, ny)
rk = My + 0.5 * sy                        # 区間中点の y
# 集約
vx = np.bincount(seg_ip, nx, nInt); vy = np.bincount(seg_ip, ny, nInt)
cyw = np.bincount(seg_ip, 0.5 * sy * L, nInt); wsum = np.bincount(seg_ip, L, nInt)
MyP = 0.5 * (y[pA[:nInt]] + y[pB[:nInt]])
rbar = MyP + cyw / wsum
exact_f = np.bincount(seg_ip, rk * ny, nInt)          # Σ_k r_k S_{k,y}
nseg = np.bincount(seg_ip, minlength=nInt)

print('internal faces: segments per face', np.bincount(nseg))
print('max |vy - SV_y| / max|SV| =', np.max(np.abs(vy - SV[:nInt, 1])) / np.max(np.abs(SV[:nInt, 1])),
      ' max|vx-SV_x|=', np.max(np.abs(vx - SV[:nInt, 0])))
print('max |rbar - PC_y| =', np.max(np.abs(rbar - PC[:nInt, 1])))

rface = np.maximum(PC[:, 1], 1e-20)
cur_f = SV[:nInt, 1] * rface[:nInt]                   # 現行 r̄·ΣS (保存値)
dlt_f = cur_f - exact_f                               # 集約誤差 (面)

# ---- 境界半割面 ----
nB = int(f['MESH'].attrs['nBconds'])
bvec = np.zeros((nP - nInt, 2)); bcy = np.zeros(nP - nInt); bexact = np.zeros(nP - nInt)
bnode_of = pA[nInt:]
# 境界エッジ -> 隣接 quad (外向き判定用)
qkey_order = np.argsort(key)
key_sorted = key[qkey_order]
for b in f['BCONDS']:
    g = f['BCONDS'][b]
    ipl = g['iPlanes'][:]; ipl = ipl[ipl >= 0]
    vbn = g['vizBfaceNodes'][:].reshape(-1, 2)
    # node -> half-plane index in this bcond
    node2hp = {int(pA[i]): int(i - nInt) for i in ipl}
    for (a, bb) in vbn:
        k = min(a, bb) * nN + max(a, bb)
        j = np.searchsorted(key_sorted, k)
        q = gq[qkey_order[j]]
        ex = x[bb] - x[a]; ey = y[bb] - y[a]
        s0, s1 = ey, -ex
        mx, my = 0.5 * (x[a] + x[bb]), 0.5 * (y[a] + y[bb])
        if s0 * (mx - G[q, 0]) + s1 * (my - G[q, 1]) < 0:   # 外向き
            s0, s1 = -s0, -s1
        w = 0.5 * np.hypot(s0, s1)
        for N, O in ((a, bb), (bb, a)):
            h = node2hp[int(N)]
            bvec[h, 0] += 0.5 * s0; bvec[h, 1] += 0.5 * s1
            r = 0.25 * (3 * y[N] + y[O])
            bcy[h] += w * r
            bexact[h] += r * 0.5 * s1
bw = np.zeros(nP - nInt)
# 重み和 (再計算)
for b in f['BCONDS']:
    g = f['BCONDS'][b]
    ipl = g['iPlanes'][:]; ipl = ipl[ipl >= 0]
    node2hp = {int(pA[i]): int(i - nInt) for i in ipl}
    for (a, bb) in g['vizBfaceNodes'][:].reshape(-1, 2):
        w = 0.5 * np.hypot(x[bb] - x[a], y[bb] - y[a])
        bw[node2hp[int(a)]] += w; bw[node2hp[int(bb)]] += w
print('boundary: max|bvec - SV| =', np.max(np.abs(bvec - SV[nInt:, :2])),
      ' max|rbar_b - PC_y| =', np.max(np.abs(bcy / bw - PC[nInt:, 1])))
cur_b = SV[nInt:, 1] * rface[nInt:]
dlt_b = cur_b - bexact

# ---- ノード集計 ----
cur = np.bincount(pA[:nInt], cur_f, nN) - np.bincount(pB[:nInt], cur_f, nN) + np.bincount(bnode_of, cur_b, nN)
exa = np.bincount(pA[:nInt], exact_f, nN) - np.bincount(pB[:nInt], exact_f, nN) + np.bincount(bnode_of, bexact, nN)
agg = np.bincount(pA[:nInt], dlt_f, nN) - np.bincount(pB[:nInt], dlt_f, nN) + np.bincount(bnode_of, dlt_b, nN)
absr = np.bincount(pA[:nInt], np.abs(cur_f), nN) + np.bincount(pB[:nInt], np.abs(cur_f), nN) + np.bincount(bnode_of, np.abs(cur_b), nN)
Avol = f['CELLS/volume'][:]

d = h5py.File(DUMP_H5, 'r')
Ap = d['axisym/A_planar'][:nN]
raw = d['axisym/closure_raw64_y'][:nN]
sumd = d['axisym/closure_sum64_y'][:nN]
print('A_planar(dump) vs CELLS/volume: max rel', np.max(np.abs(Ap - Avol) / Avol))
print('current(recon) vs closure_raw64_y: max |diff|/A', np.max(np.abs(cur - raw) / Ap))
print('closure_sum64_y vs raw64: max |diff|/A', np.max(np.abs(sumd - raw) / Ap))

# 双対多角形の面積 (sub-CV quads A,M1,G,M2) を独立に再計算
poly = np.zeros(nN)
for k in range(4):
    a = quads[:, k]; nb1 = quads[:, (k + 1) % 4]; nb0 = quads[:, (k - 1) % 4]
    m1 = 0.5 * (X[a] + X[nb1]); m2 = 0.5 * (X[a] + X[nb0])
    P = np.stack([X[a], m1, G, m2], axis=1)[:, :, :2] - X[a][:, None, :2]
    cr = P[:, :, 0] * np.roll(P[:, :, 1], -1, 1) - np.roll(P[:, :, 0], -1, 1) * P[:, :, 1]
    poly += np.bincount(a, 0.5 * np.abs(cr.sum(1)), nN)
print('dual polygon area vs A_planar: max rel', np.max(np.abs(poly - Ap) / Ap))

def_cur = Ap - cur          # 観測された欠損 (現行)
def_exa = Ap - exa          # 区間ごとの厳密形
# def_cur = def_exa - agg
idx = np.arange(nN); j = idx % NJ; i = idx // NJ
mask_x = np.ones(nN, bool)
print()
print(' j |   max|def_cur|/A    rms     | max|def_exa|/A    rms     | resid=|def_cur+agg|/|def_cur| | corr(def_cur,-agg) | sign agree')
for jj in [0, 1, 2, 3, 4, 5, 6, 7, 8, 60, 118, 119, 120]:
    s = j == jj
    rc = def_cur[s] / Ap[s]; re = def_exa[s] / Ap[s]; ra = -agg[s] / Ap[s]
    res = np.linalg.norm(rc - ra) / max(np.linalg.norm(rc), 1e-300)
    c = np.corrcoef(rc, ra)[0, 1] if np.std(rc) > 0 and np.std(ra) > 0 else float('nan')
    big = np.abs(rc) > 0.1 * np.max(np.abs(rc))
    sa = np.mean(np.sign(rc[big]) == np.sign(ra[big])) if big.any() else float('nan')
    print(f'{jj:3d} | {np.max(np.abs(rc)):.3e} {np.sqrt(np.mean(rc**2)):.3e} | {np.max(np.abs(re)):.3e} {np.sqrt(np.mean(re**2)):.3e} |'
          f' {res:.3e} | {c:.6f} | {sa:.3f}  mean def_cur/A {np.mean(rc):+.3e}')

for jj in [2, 3, 4]:
    s = np.where(j == jj)[0]
    rc = def_cur[s] / Ap[s]
    top = s[np.argsort(-np.abs(rc))[:5]]
    print(f'j={jj} largest: ' + ', '.join(f'x/rt={x[n]/RT:+.3f} r/rt={y[n]/RT:.4f} def/A={def_cur[n]/Ap[n]:+.3e} (-agg/A={-agg[n]/Ap[n]:+.3e})' for n in top))
# x 分布の概略 (j=3)
s = np.where(j == 3)[0]
rc = def_cur[s] / Ap[s]
for lo_, hi_ in [(-20, -5), (-5, -2), (-2, -0.5), (-0.5, 0.5), (0.5, 2), (2, 10), (10, 100)]:
    m = (x[s] / RT >= lo_) & (x[s] / RT < hi_)
    if m.any():
        print(f'j=3 x/rt in [{lo_},{hi_}): n={m.sum()} max|def|/A={np.max(np.abs(rc[m])):.3e} mean={np.mean(rc[m]):+.3e}')
print('global: max|def_exa|/A =', np.max(np.abs(def_exa / Ap)), ' max|def_cur|/A =', np.max(np.abs(def_cur / Ap)),
      ' at j', j[np.argmax(np.abs(def_cur / Ap))])
np.savez(OUT_NPZ, def_cur=def_cur, def_exa=def_exa, agg=agg, Ap=Ap, x=x, y=y)

# ---- 追加: 集約誤差を面の種類で分ける (radial edge = 同じ i、axial edge = 同じ j) ----
radial = np.abs(pA[:nInt] - pB[:nInt]) == 1
for name, m in (('radial-edge faces (normal ~r)', radial), ('axial-edge faces (normal ~x)', ~radial)):
    a = np.bincount(pA[:nInt][m], dlt_f[m], nN) - np.bincount(pB[:nInt][m], dlt_f[m], nN)
    for jj in (2, 3, 4):
        s = j == jj
        print(f'{name:32s} j={jj}: max|agg|/A={np.max(np.abs(a[s]/Ap[s])):.3e}')
# 丸め水準の目安: Σ|r S_y| / A (j=120)
s = j == 120
print('j=120: median Σ|rS_y|/A =', np.median(absr[s] / Ap[s]), ' max =', np.max(absr[s] / Ap[s]),
      ' -> eps*ratio max =', 2.2e-16 * np.max(absr[s] / Ap[s]))
k = np.argmax(np.abs(def_exa[s] / Ap[s]))
print('j=120 worst def_exa node: Σ|rS|/A =', (absr[s] / Ap[s])[k], ' def_exa/A =', (def_exa[s] / Ap[s])[k])
m = ~radial
a_ax = np.bincount(pA[:nInt][m], dlt_f[m], nN) - np.bincount(pB[:nInt][m], dlt_f[m], nN)
a_bd = np.bincount(bnode_of, dlt_b, nN)
for jj in (2, 3, 4):
    s = np.where(j == jj)[0]
    n = s[np.argmax(np.abs(a_ax[s] / Ap[s]))]
    print(f'j={jj} axial-max node i={i[n]} x/rt={x[n]/RT:+.3f}: ax={a_ax[n]/Ap[n]:+.3e} bnd={a_bd[n]/Ap[n]:+.3e} total agg={agg[n]/Ap[n]:+.3e} def_cur={def_cur[n]/Ap[n]:+.3e}')

# ==== plan axisymmetric-freestream-hoop-gauge §4.3 (事前登録): 面のメトリックの A/B (x・y 両成分) ====
# A: 今の r̄_f·ΣS_k、B: 区間ごとの W_f = Σ_k r_k S_k。接続・向き・境界の所属・A_planar は共通。
exact_fx = np.bincount(seg_ip, rk * nx, nInt)                     # Σ_k r_k S_{k,x} (内部面)
bexact_x = np.zeros(nP - nInt)
for b in f['BCONDS']:
    g = f['BCONDS'][b]
    ipl = g['iPlanes'][:]; ipl = ipl[ipl >= 0]
    node2hp = {int(pA[i]): int(i - nInt) for i in ipl}
    for (a, bb) in g['vizBfaceNodes'][:].reshape(-1, 2):
        k = min(a, bb) * nN + max(a, bb)
        jq = np.searchsorted(key_sorted, k); q = gq[qkey_order[jq]]
        ex = x[bb] - x[a]; ey = y[bb] - y[a]; s0, s1 = ey, -ex
        mx, my = 0.5 * (x[a] + x[bb]), 0.5 * (y[a] + y[bb])
        if s0 * (mx - G[q, 0]) + s1 * (my - G[q, 1]) < 0: s0, s1 = -s0, -s1
        for N, O in ((a, bb), (bb, a)):
            bexact_x[node2hp[int(N)]] += 0.25 * (3 * y[N] + y[O]) * 0.5 * s0
def node_sum(fi, fb):
    return np.bincount(pA[:nInt], fi, nN) - np.bincount(pB[:nInt], fi, nN) + np.bincount(bnode_of, fb, nN)
def node_abs(fi, fb):
    return np.bincount(pA[:nInt], np.abs(fi), nN) + np.bincount(pB[:nInt], np.abs(fi), nN) + np.bincount(bnode_of, np.abs(fb), nN)
curx_f = SV[:nInt, 0] * rface[:nInt]; curx_b = SV[nInt:, 0] * rface[nInt:]
M = {"A (r̄·ΣS)": (curx_f, curx_b, cur_f, cur_b), "B (Σ r_k S_k)": (exact_fx, bexact_x, exact_f, bexact)}
eps = np.finfo(np.float64).eps
print("\n== §4.3 の A/B: E_x = |Σ±W_x|/A_planar、E_y = |Σ±W_y − A_planar|/A_planar")
res43 = {}
for name, (fx_, bx_, fy_, by_) in M.items():
    Sx = node_sum(fx_, bx_); Sy = node_sum(fy_, by_); Ax = node_abs(fx_, bx_); Ay = node_abs(fy_, by_)
    Ex = np.abs(Sx) / Ap; Ey = np.abs(Sy - Ap) / Ap
    m28 = (j >= 2) & (j <= 8)
    tolx = 100 * eps * (Ap + Ax); toly = 100 * eps * (Ap + Ay)
    gx = int(np.count_nonzero(np.abs(Sx) > tolx)); gy = int(np.count_nonzero(np.abs(Sy - Ap) > toly))
    print(f"  {name}: j 2〜8 の最大 E_x {Ex[m28].max():.3e}・E_y {Ey[m28].max():.3e} | 全域の最大 E_x {Ex.max():.3e} (j {j[np.argmax(Ex)]})・E_y {Ey.max():.3e} (j {j[np.argmax(Ey)]})"
          f" | 100·ε64·(A + Σ|W|) を超える CV: x {gx}・y {gy}")
    for jj in (0, 1, 2, 3, 4, 8, 60, 119, 120):
        s_ = j == jj
        print(f"     j {jj:3d}: max E_x {Ex[s_].max():.3e}  max E_y {Ey[s_].max():.3e}  一定の U の質量の残差 max|Σ±W_x|/Σ|W_x| {np.max(np.abs(Sx[s_]) / np.maximum(Ax[s_], 1e-300)):.3e}")
    res43[name] = dict(Ex28=float(Ex[m28].max()), Ey28=float(Ey[m28].max()), over_x=gx, over_y=gy)
b = res43["B (Σ r_k S_k)"]
ok = b["Ex28"] <= 1e-10 and b["Ey28"] <= 1e-10 and b["over_x"] == 0 and b["over_y"] == 0
print(f"== VERDICT §4.3: {'B が合格 (j 2〜8 の E_x・E_y ≤ 1e-10、全域で丸めの規模以内)' if ok else 'B が不合格 (集約の誤差だけという仮説を棄却、境界の所属・向き・面積を点検し直す)'}")

# ==== plan axisymmetric-freestream-hoop-gauge §4.8 (事前登録): 境界の半割面の端点を丸めた中点 M̃ にそろえる CPU の A/B ====
# A: 今の半割面 (ベクトル 0.5·(B−A)、r = (3N+O)/4 の y) = §4.3 の B。B: N → M̃ = fl(0.5·(x_a + x_b)) の区間、S = rot(M̃ − N) を外向き、r = (N_y + M̃_y)/2。
bal_x = np.zeros(nP - nInt); bal_y = np.zeros(nP - nInt); dSmax = 0.0
for b in f['BCONDS']:
    g = f['BCONDS'][b]
    ipl = g['iPlanes'][:]; ipl = ipl[ipl >= 0]
    node2hp = {int(pA[i]): int(i - nInt) for i in ipl}
    for (a, bb) in g['vizBfaceNodes'][:].reshape(-1, 2):
        k = min(a, bb) * nN + max(a, bb)
        jq = np.searchsorted(key_sorted, k); q = gq[qkey_order[jq]]
        ex = x[bb] - x[a]; ey = y[bb] - y[a]; s0, s1 = ey, -ex
        mxe, mye = 0.5 * (x[a] + x[bb]), 0.5 * (y[a] + y[bb])            # 丸めた中点 M̃ (内部の区間と同じ演算)
        outward = (s0 * (mxe - G[q, 0]) + s1 * (mye - G[q, 1])) >= 0
        for N in (a, bb):
            dx_, dy_ = mxe - x[N], mye - y[N]
            h0, h1 = dy_, -dx_
            if (h0 * s0 + h1 * s1 >= 0) != outward: h0, h1 = -h0, -h1   # 辺の外向きの法線と同じ向きに
            r_ = 0.5 * (y[N] + mye)
            hp = node2hp[int(N)]
            bal_x[hp] += r_ * h0; bal_y[hp] += r_ * h1
print("\n== §4.8 の A/B: 境界の半割面の端点 (A 今 / B 丸めた中点にそろえる)。内部の面はどちらも区間ごとの W")
ver48 = {}
for name, (bx_, by_) in (("A (今の半割面)", (bexact_x, bexact)), ("B (M̃ にそろえる)", (bal_x, bal_y))):
    Sx = node_sum(exact_fx, bx_); Sy = node_sum(exact_f, by_); Ax = node_abs(exact_fx, bx_); Ay = node_abs(exact_f, by_)
    Ex = np.abs(Sx) / Ap; Ey = np.abs(Sy - Ap) / Ap; m28 = (j >= 2) & (j <= 8)
    ox = int(np.count_nonzero(np.abs(Sx) > 100 * eps * (Ap + Ax))); oy = int(np.count_nonzero(np.abs(Sy - Ap) > 100 * eps * (Ap + Ay)))
    print(f"  {name}: 超えた CV x {ox}・y {oy} | j 2〜8 の最大 E_x {Ex[m28].max():.3e}・E_y {Ey[m28].max():.3e} | 全域の最大 E_x {Ex.max():.3e}・E_y {Ey.max():.3e}")
    ver48[name] = (ox, oy, float(Ex[m28].max()), float(Ey[m28].max()))
print(f"  (記録) 半割面の W の変化の最大 |ΔW|/|W_A|: x {np.max(np.abs(bal_x - bexact_x) / np.maximum(np.abs(bexact_x), 1e-300)):.3e}・y {np.max(np.abs(bal_y - bexact) / np.maximum(np.abs(bexact), 1e-300)):.3e}")
a_, b_ = ver48["A (今の半割面)"], ver48["B (M̃ にそろえる)"]
if (a_[0], a_[1]) != (2664, 1534): v = f"比較は無効 (A が §4.4 の FAIL を再現しない: {a_[:2]})"
elif b_[0] == 0 and b_[1] == 0 and b_[2] <= 1e-10 and b_[3] <= 1e-10: v = "支持 (端点をそろえると全 CV で今の閾値に収まる)"
else: v = "棄却 (B にも超過が残る)"
print(f"== VERDICT §4.8: {v}")
