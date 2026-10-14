#!/usr/bin/env python3
"""plan axisymmetric-freestream-hoop-gauge §4.6 の 0: §4.4 の全域の FAIL (B = Σ r_k S_k の閉性) の高精度の内訳。

入力の binary64 の節点座標 (/MESH/COORD) を厳密な有理数とみなし、`fractions.Fraction` で次の 3 つを作る。
  (a) 一貫した双対多角形: 中点 M = (A+B)/2・四角形の重心 G = 頂点平均を厳密に作り、区間 M→G (内部) と
      節点→エッジ中点 (境界の半割面、変換器と同じく 0.5·(B−A)、r = (3N+O)/4) から W_a = Σ r_k S_k と
      双対面積 area_a を厳密に作る。閉性は恒等式なので 0 になるはず (検査する)。
  (b) 今のコードの端点: 変換器と同じ binary64 の演算で端点と区間のベクトルを作る
      (M̃ = 0.5·(x_A + x_B)、G̃ = 頂点の逐次和/4、内部 s̃ = G̃ − M̃、境界 0.5·fl(B − A)。
      gmshReader.hpp:1475-1531, 1641-1699)。r_k はその区間の厳密な中点の半径、積と和は厳密。
  binary64: dual_segment_closure.py の B と同じ演算 (rk = My + 0.5·sy、W = rk·n、面ごと・CV ごとの bincount)。

符号つきの欠損 D (x: Σ±W_x、y: Σ±W_y − A_planar、A_planar = /CELLS/volume) を次に分ける (厳密に D_obs と一致する):
  D_obs = D_a + E1 (端点の生成) + E2 (binary64 の r・積・集約) + E3 (A_planar の double の面積)
    D_a = ΣW_a − area_a                       (構成の誤り。0 でなければ接続・向き・所属の問題)
    E1  = (ΣW_b − area_b) − D_a               (端点の丸めによる多角形の不整合)
    E2  = ΣW_64 − ΣW_b                        (E2a: r と積の丸め、E2b: 面・CV の和の丸め)
    E3  = area_b − A_planar                   (y のみ。area_b は丸めた端点の多角形 (N, M̃1, G̃, M̃2) の厳密な面積)
端点の機構の予測: 境界の半割面の端 N + 0.5·fl(O − N) (= 厳密な中点 M) と、内部の区間の始点 M̃ = fl(M) の
食い違い g = M̃ − M (|g_c| ≤ u·|M_c|、u = 2^-53) を閉じる区間の W を W_close とすると E1 ≈ −ΣW_close。
CV ごとの上界: B1 (端点) = u·Σ_b (|M̃_y| + |s_h|)·(|M̃_c'| + |s_h|)、B2 (演算) = (n + 4)·u·Σ_k|W_k|、
B3 (面積) = 6u·Σ_q Σ_k (|p_x,k p_y,k+1| + |p_x,k+1 p_y,k|)/2 + n_q·u·A。

usage: python3 dual_segment_closure_hp.py <nozzle.h5 (FP64 の変換器)> <stage2_fp64.h5 | -> <出力 .npz> [--nctl 300] [--seed 0]
  標本は「§4.4 で超えた CV の全部」∪「境界の CV の全部」∪「内部の CV の無作為 nctl 個」∪「壁の隣 (j 119) の無作為 100 個」。
"""
import argparse
import time
from fractions import Fraction as Fr

import h5py
import numpy as np

ap_ = argparse.ArgumentParser()
ap_.add_argument('mesh'); ap_.add_argument('dump'); ap_.add_argument('out')
ap_.add_argument('--nctl', type=int, default=300); ap_.add_argument('--seed', type=int, default=0)
args = ap_.parse_args()
NJ = 121
eps = np.finfo(np.float64).eps
u = eps / 2
t0 = time.time()

f = h5py.File(args.mesh, 'r')
X = f['MESH/COORD'][:].reshape(-1, 3)
nN = X.shape[0]
x, y = X[:, 0], X[:, 1]
cw = f['VIZMESH/CONNE'][:]
assert len(cw) % 5 == 0 and np.all(cw[::5] == 5), 'quad 以外を含む格子は対象外'
quads = cw.reshape(-1, 5)[:, 1:].astype(np.int64)
nQ = len(quads)
G = X[quads].mean(axis=1)                                   # dual_segment_closure.py と同じ
Gseq = (((X[quads[:, 0]] + X[quads[:, 1]]) + X[quads[:, 2]]) + X[quads[:, 3]]) / 4.0   # 変換器 (cc += x_n; cc/4)
print(f'G: numpy mean と変換器の逐次和のビット不一致 {np.count_nonzero(G != Gseq)} / {G.size}')

# ---- planes (STRUCT を一括で解く: 内部 [2 n0 n1 2 c0 c1]、境界 [1 n 1 c]) ----
ps = f['PLANES/STRUCT'][:]
nInt = int(f['MESH'].attrs['nNormalPlanes']); nP = int(f['MESH'].attrs['nPlanes'])
ii = ps[:6 * nInt].reshape(-1, 6).astype(np.int64)
bb_ = ps[6 * nInt:].reshape(-1, 4).astype(np.int64)
assert np.all(ii[:, 0] == 2) and np.all(ii[:, 3] == 2) and np.all(bb_[:, 0] == 1) and np.all(bb_[:, 2] == 1)
assert len(bb_) == nP - nInt
pA = np.concatenate([ii[:, 4], bb_[:, 3]]); pB = np.concatenate([ii[:, 5], np.full(nP - nInt, -1)])
SV = f['PLANES/surfVect'][:].reshape(-1, 3)

# ---- 内部面の区間 (dual_segment_closure.py と同じ演算) ----
ea = quads.reshape(-1); eb = np.roll(quads, -1, axis=1).reshape(-1)
gq = np.repeat(np.arange(nQ), 4)
key = np.minimum(ea, eb) * nN + np.maximum(ea, eb)
pkey = np.minimum(pA[:nInt], pB[:nInt]) * nN + np.maximum(pA[:nInt], pB[:nInt])
order = np.argsort(pkey)
pos = np.searchsorted(pkey[order], key)
assert np.all(pkey[order][pos] == key)
seg_ip = order[pos]
A_ = pA[seg_ip]; B_ = pB[seg_ip]
Mx = 0.5 * (x[A_] + x[B_]); My = 0.5 * (y[A_] + y[B_])
sx = G[gq, 0] - Mx; sy = G[gq, 1] - My
nx = sy.copy(); ny = -sx
eABx = x[B_] - x[A_]; eABy = y[B_] - y[A_]
flip = nx * eABx + ny * eABy < 0
nx[flip] *= -1; ny[flip] *= -1
rk = My + 0.5 * sy
vx = np.bincount(seg_ip, nx, nInt); vy = np.bincount(seg_ip, ny, nInt)
print(f'内部面: Σ 区間の n と /PLANES/surfVect のビット不一致 x {np.count_nonzero(vx != SV[:nInt, 0])}・y {np.count_nonzero(vy != SV[:nInt, 1])}'
      ' (0 なら M̃・G̃・s̃ は変換器と同じ値)')
exact_fx = np.bincount(seg_ip, rk * nx, nInt); exact_fy = np.bincount(seg_ip, rk * ny, nInt)

# ---- 境界の半割面 (dual_segment_closure.py と同じ演算・同じ和の順序) ----
bnode_of = pA[nInt:]
qkey_order = np.argsort(key); key_sorted = key[qkey_order]
bex_x = np.zeros(nP - nInt); bex_y = np.zeros(nP - nInt)
hrec = []   # (N, O, a, bb, s0, s1, q, bcond)
bnodes = {}
for bname in f['BCONDS']:
    g = f['BCONDS'][bname]
    ipl = g['iPlanes'][:]; ipl = ipl[ipl >= 0]
    node2hp = {int(pA[i]): int(i - nInt) for i in ipl}
    bnodes[bname] = set(node2hp)
    for (a, bb) in g['vizBfaceNodes'][:].reshape(-1, 2):
        k = min(a, bb) * nN + max(a, bb)
        jq = np.searchsorted(key_sorted, k); q = int(gq[qkey_order[jq]])
        ex = x[bb] - x[a]; ey = y[bb] - y[a]
        s0, s1 = ey, -ex
        mx, my = 0.5 * (x[a] + x[bb]), 0.5 * (y[a] + y[bb])
        if s0 * (mx - G[q, 0]) + s1 * (my - G[q, 1]) < 0:
            s0, s1 = -s0, -s1
        for N, O in ((a, bb), (bb, a)):
            h = node2hp[int(N)]
            r = 0.25 * (3 * y[N] + y[O])
            bex_y[h] += r * 0.5 * s1
            bex_x[h] += r * 0.5 * s0
            hrec.append((int(N), int(O), int(a), int(bb), float(s0), float(s1), q, bname))
print(f'境界の半割面の記録 {len(hrec)}')

def node_sum(fi, fb):
    return np.bincount(pA[:nInt], fi, nN) - np.bincount(pB[:nInt], fi, nN) + np.bincount(bnode_of, fb, nN)

def node_abs(fi, fb):
    return np.bincount(pA[:nInt], np.abs(fi), nN) + np.bincount(pB[:nInt], np.abs(fi), nN) + np.bincount(bnode_of, np.abs(fb), nN)

Avol = f['CELLS/volume'][:]
if args.dump != '-':
    Ap = h5py.File(args.dump, 'r')['axisym/A_planar'][:nN]
    print(f'A_planar(dump) と /CELLS/volume のビット不一致 {np.count_nonzero(Ap != Avol)}')
else:
    Ap = Avol
Sx = node_sum(exact_fx, bex_x); Sy = node_sum(exact_fy, bex_y)
Ax = node_abs(exact_fx, bex_x); Ay = node_abs(exact_fy, bex_y)
tolx = 100 * eps * (Ap + Ax); toly = 100 * eps * (Ap + Ay)
ox = np.abs(Sx) > tolx; oy = np.abs(Sy - Ap) > toly
idx = np.arange(nN); jn = idx % NJ; inn = idx // NJ; NI = nN // NJ
print(f'§4.4 の再現: 100·ε64·(A + Σ|W|) を超える CV x {ox.sum()}・y {oy.sum()} (§4.4 は 2,664 / 1,534)')

# ---- 領域の分類 (境界の所属から) ----
isw = np.zeros(nN, bool); isi = np.zeros(nN, bool); iso = np.zeros(nN, bool); isa = np.zeros(nN, bool)
kind = {k: f['BCONDS'][k].attrs['bcondKind'] for k in f['BCONDS']}
for k, s in bnodes.items():
    arr = np.fromiter(s, np.int64)
    kd = kind[k]
    (isw if kd.startswith('wall') else isi if kd.startswith('inlet') else iso if kd.startswith('outlet') else isa)[arr] = True
region = np.full(nN, 'interior', dtype=object)
region[isa] = 'axis'; region[iso] = 'outlet'; region[isi] = 'inlet'; region[isw] = 'wall'
region[isw & (isi | iso)] = 'corner'
region[isa & (isi | iso)] = 'axis-corner'
for nm, o in (('x', ox), ('y', oy)):
    rr, cc = np.unique(region[o], return_counts=True)
    print(f'  超えた CV ({nm}) の領域: ' + ', '.join(f'{a} {b}' for a, b in zip(rr, cc)))

# ---- 標本 ----
rng = np.random.default_rng(args.seed)
isb = isw | isi | iso | isa
interior_ok = np.flatnonzero(~isb & ~ox & ~oy & (jn >= 1) & (jn <= 118))
near_wall = np.flatnonzero(~isb & ~ox & ~oy & (jn == 119))
ctl = np.concatenate([rng.choice(interior_ok, min(args.nctl, len(interior_ok)), replace=False),
                      rng.choice(near_wall, min(100, len(near_wall)), replace=False)])
sample = np.unique(np.concatenate([np.flatnonzero(ox | oy), np.flatnonzero(isb), ctl]))
print(f'標本 {len(sample)} CV (超えた CV {np.count_nonzero(ox | oy)}・境界 {isb.sum()}・対照 {len(ctl)})')

# 索引: 区間 → 節点 (A 側 +、B 側 −)、節点 → 四角形、節点 → 半割面
def csr(keys):
    o = np.argsort(keys, kind='stable'); ks = keys[o]
    return o, ks
oA, kA = csr(A_); oB, kB = csr(B_)
qn = quads.reshape(-1); qi = np.repeat(np.arange(nQ), 4); qk = np.tile(np.arange(4), nQ)
oq, kq = csr(qn)
hby = {}
for t, rec in enumerate(hrec):
    hby.setdefault(rec[0], []).append(t)

F = lambda v: Fr(float(v))
Xf = {}
def XF(n):
    if n not in Xf:
        Xf[n] = (F(x[n]), F(y[n]))
    return Xf[n]
def GF(q):   # 厳密な頂点平均
    vs = [XF(int(v)) for v in quads[q]]
    return (sum(v[0] for v in vs) / 4, sum(v[1] for v in vs) / 4)
def shoelace(P):
    s = 0
    for k in range(len(P)):
        a, b = P[k], P[(k + 1) % len(P)]
        s += a[0] * b[1] - b[0] * a[1]
    return s / 2

res = {k: [] for k in ('node', 'Dobs_x', 'Dobs_y', 'Da_x', 'Da_y', 'E1_x', 'E1_y', 'E2a_x', 'E2a_y', 'E2b_x', 'E2b_y', 'E3_y',
                       'pred1_x', 'pred1_y', 'B1_x', 'B1_y', 'B2_x', 'B2_y', 'B3_y', 'tol_x', 'tol_y',
                       'areaA_rel', 'areaB_rel', 'Wby_rel', 'flipmis', 'subinexact', 'A')}
nflipmis_total = 0; nsubinexact = 0
for n in sample:
    n = int(n)
    Wa = [Fr(0), Fr(0)]; Wb = [Fr(0), Fr(0)]; W64e = [Fr(0), Fr(0)]
    absW = [0.0, 0.0]; nterm = 0; flipmis = 0; subinex = 0
    for oarr, karr, sgn in ((oA, kA, 1), (oB, kB, -1)):
        lo = np.searchsorted(karr, n, 'left'); hi = np.searchsorted(karr, n, 'right')
        for s in oarr[lo:hi]:
            a_, b_ = int(A_[s]), int(B_[s]); q = int(gq[s])
            xa, xb = XF(a_), XF(b_)
            # (a) 厳密な端点
            M = ((xa[0] + xb[0]) / 2, (xa[1] + xb[1]) / 2); Gq = GF(q)
            sax, say = Gq[0] - M[0], Gq[1] - M[1]
            nax, nay = say, -sax
            fa = nax * (xb[0] - xa[0]) + nay * (xb[1] - xa[1]) < 0
            if fa != bool(flip[s]):
                flipmis += 1
            if fa:
                nax, nay = -nax, -nay
            ra = (M[1] + Gq[1]) / 2
            Wa[0] += sgn * ra * nax; Wa[1] += sgn * ra * nay
            # (b) 今のコードの端点 (binary64 の M̃・G̃・s̃)、積と和は厳密
            sbx, sby = F(sx[s]), F(sy[s])
            if sbx != F(G[q, 0]) - F(Mx[s]) or sby != F(G[q, 1]) - F(My[s]):
                subinex += 1
            nbx, nby = (sby, -sbx) if not flip[s] else (-sby, sbx)
            rb = F(My[s]) + sby / 2
            Wb[0] += sgn * rb * nbx; Wb[1] += sgn * rb * nby
            w64x, w64y = rk[s] * nx[s], rk[s] * ny[s]
            W64e[0] += sgn * F(w64x); W64e[1] += sgn * F(w64y)
            absW[0] += abs(w64x); absW[1] += abs(w64y); nterm += 1
    pred = [Fr(0), Fr(0)]; B1 = [0.0, 0.0]
    for t in hby.get(n, []):
        N, O, a, bb, s0, s1, q, bname = hrec[t]
        xN, xO = XF(N), XF(O); xa_, xbb = XF(a), XF(bb)
        # (a)
        exa, eya = xbb[0] - xa_[0], xbb[1] - xa_[1]
        Ma = ((xa_[0] + xbb[0]) / 2, (xa_[1] + xbb[1]) / 2)
        Gq = GF(q)
        t0a, t1a = eya, -exa
        if t0a * (Ma[0] - Gq[0]) + t1a * (Ma[1] - Gq[1]) < 0:
            t0a, t1a = -t0a, -t1a
        fa = (t0a, t1a) != (eya, -exa)                       # (a) で反転したか
        ey64_, ex64_ = F(y[bb] - y[a]), F(x[bb] - x[a])
        f64 = (F(s1) != -ex64_) if ex64_ != 0 else (F(s0) != ey64_)   # binary64 で反転したか
        if fa != f64:
            flipmis += 1
        ra = (3 * xN[1] + xO[1]) / 4
        Wa[0] += ra * t0a / 2; Wa[1] += ra * t1a / 2
        # (b): 半割面は N → N + 0.5·fl(O − N)、法線 0.5·(s0, s1)
        sg = 1 if N == a else -1
        ey64 = F(y[bb] - y[a]); ex64 = F(x[bb] - x[a])
        if ex64 != xbb[0] - xa_[0] or ey64 != xbb[1] - xa_[1]:
            subinex += 1
        rb = xN[1] + sg * ey64 / 4
        hb = (F(s0) / 2, F(s1) / 2)
        Wb[0] += rb * hb[0]; Wb[1] += rb * hb[1]
        r64 = 0.25 * (3 * y[N] + y[O])
        w64x, w64y = r64 * 0.5 * s0, r64 * 0.5 * s1
        W64e[0] += F(w64x); W64e[1] += F(w64y)
        absW[0] += abs(w64x); absW[1] += abs(w64y); nterm += 1
        # 端点の機構の予測: 半割面の端 Me = N + 0.5·fl(O−N) と内部の区間の始点 M̃ = fl(0.5·(a+bb)) を閉じる区間
        Me = (xN[0] + sg * ex64 / 2, xN[1] + sg * ey64 / 2)
        Mt = (F(0.5 * (x[a] + x[bb])), F(0.5 * (y[a] + y[bb])))
        sccw = (-hb[1], hb[0])                       # 外向き法線 (n_x, n_y) の CCW の向き (−n_y, n_x)
        fwd = sccw[0] * (Me[0] - xN[0]) + sccw[1] * (Me[1] - xN[1]) > 0
        P0, P1 = (Me, Mt) if fwd else (Mt, Me)
        sc = (P1[0] - P0[0], P1[1] - P0[1]); rc = (P0[1] + P1[1]) / 2
        pred[0] -= rc * sc[1]; pred[1] -= rc * (-sc[0])
        Mf = (float(Mt[0]), float(Mt[1])); shl = 0.5 * float(np.hypot(s0, s1))
        B1[0] += u * (abs(Mf[1]) + shl) * (abs(Mf[1]) + shl)
        B1[1] += u * (abs(Mf[1]) + shl) * (abs(Mf[0]) + shl)
    # 面積: (a) 厳密な端点、(b) 丸めた端点 (変換器の面積の式が表す多角形)
    lo = np.searchsorted(kq, n, 'left'); hi = np.searchsorted(kq, n, 'right')
    areaA = Fr(0); areaB = Fr(0); B3 = 0.0
    for t in oq[lo:hi]:
        q = int(qi[t]); k = int(qk[t])
        v1 = int(quads[q, (k + 1) % 4]); v0 = int(quads[q, (k - 1) % 4])
        xN, x1, x0 = XF(n), XF(v1), XF(v0); Gq = GF(q)
        M1 = ((xN[0] + x1[0]) / 2, (xN[1] + x1[1]) / 2); M2 = ((xN[0] + x0[0]) / 2, (xN[1] + x0[1]) / 2)
        areaA += abs(shoelace([xN, M1, Gq, M2]))
        M1t = (F(0.5 * (x[n] + x[v1])), F(0.5 * (y[n] + y[v1]))); M2t = (F(0.5 * (x[n] + x[v0])), F(0.5 * (y[n] + y[v0])))
        Gt = (F(G[q, 0]), F(G[q, 1]))
        Pb = [xN, M1t, Gt, M2t]
        areaB += abs(shoelace(Pb))
        pr = [(float(p[0] - xN[0]), float(p[1] - xN[1])) for p in Pb]
        B3 += 6 * u * sum(abs(pr[k2][0] * pr[(k2 + 1) % 4][1]) + abs(pr[(k2 + 1) % 4][0] * pr[k2][1]) for k2 in range(4)) / 2
    B3 += (hi - lo) * u * Ap[n]
    Apn = F(Ap[n])
    Dobs = (F(Sx[n]), F(Sy[n]) - Apn)
    Da = (Wa[0], Wa[1] - areaA)
    E1 = (Wb[0] - Da[0], (Wb[1] - areaB) - Da[1])
    E2a = (W64e[0] - Wb[0], W64e[1] - Wb[1])
    E2b = (F(Sx[n]) - W64e[0], F(Sy[n]) - W64e[1])
    E3 = areaB - Apn
    assert Da[0] + E1[0] + E2a[0] + E2b[0] == Dobs[0]
    assert Da[1] + E1[1] + E2a[1] + E2b[1] + E3 == Dobs[1]
    An = float(Ap[n])
    for kk, v in (('node', n), ('Dobs_x', float(Dobs[0])), ('Dobs_y', float(Dobs[1])), ('Da_x', float(Da[0])), ('Da_y', float(Da[1])),
                  ('E1_x', float(E1[0])), ('E1_y', float(E1[1])), ('E2a_x', float(E2a[0])), ('E2a_y', float(E2a[1])),
                  ('E2b_x', float(E2b[0])), ('E2b_y', float(E2b[1])), ('E3_y', float(E3)),
                  ('pred1_x', float(pred[0])), ('pred1_y', float(pred[1])), ('B1_x', B1[0]), ('B1_y', B1[1]),
                  ('B2_x', (nterm + 4) * u * absW[0]), ('B2_y', (nterm + 4) * u * absW[1]), ('B3_y', B3),
                  ('tol_x', float(tolx[n])), ('tol_y', float(toly[n])),
                  ('areaA_rel', float((areaA - Apn) / Apn)), ('areaB_rel', float((areaB - Apn) / Apn)),
                  ('Wby_rel', float((Wb[1] - Apn) / Apn)), ('flipmis', flipmis), ('subinexact', subinex), ('A', An)):
        res[kk].append(v)
R = {k: np.array(v) for k, v in res.items()}
print(f'厳密演算 {time.time() - t0:.0f} s')
nd = R['node']; reg = region[nd]; fx = ox[nd]; fy = oy[nd]

# ---- 報告 ----
print('\n== 構成の検査 (参照 (a))')
print(f'  D_a ≠ 0 の CV: x {np.count_nonzero(R["Da_x"] != 0)}・y {np.count_nonzero(R["Da_y"] != 0)} / {len(nd)} (0 なら (a) は厳密に閉じる)')
print(f'  向きの判定の (a)/binary64 の食い違い: {int(R["flipmis"].sum())}・区間のベクトルが binary64 で非厳密: {int(R["subinexact"].sum())}')

def q(v):
    return f'{np.median(v):+.3f} [{np.min(v):+.3f}, {np.max(v):+.3f}]' if len(v) else '—'

for comp in ('x', 'y'):
    print(f'\n== 成分 {comp}: 符号つきの欠損 D_obs の内訳 (E/D_obs の中央値 [最小, 最大])、上界の比')
    fails = fx if comp == 'x' else fy
    groups = [(f'超えた CV・{r}', fails & (reg == r)) for r in ('wall', 'inlet', 'outlet', 'corner', 'axis', 'axis-corner', 'interior')]
    groups += [('超えていない境界の CV', ~fails & (reg != 'interior')), ('対照 (内部)', ~fails & (reg == 'interior'))]
    for nm, m in groups:
        if not m.any():
            continue
        D = R[f'Dobs_{comp}'][m]; A = R['A'][m]
        nz = D != 0
        e1 = R[f'E1_{comp}'][m]; e2 = R[f'E2a_{comp}'][m] + R[f'E2b_{comp}'][m]
        e3 = R['E3_y'][m] if comp == 'y' else np.zeros(m.sum())
        B = R[f'B1_{comp}'][m] + R[f'B2_{comp}'][m] + (R['B3_y'][m] if comp == 'y' else 0)
        pr = R[f'pred1_{comp}'][m]
        print(f'  {nm}: n {m.sum()}  |D_obs|/A 中央値 {np.median(np.abs(D) / A):.2e}・最大 {np.max(np.abs(D) / A):.2e}')
        if nz.any():
            print(f'     端点 E1/D {q(e1[nz] / D[nz])}  演算 E2/D {q(e2[nz] / D[nz])}  面積 E3/D {q(e3[nz] / D[nz])}')
        pm = pr != 0
        if pm.any():
            print(f'     E1/予測 (−ΣW_close) {q(e1[pm] / pr[pm])}')
        print(f'     |D_obs|/上界 (B1+B2+B3) 中央値 {np.median(np.abs(D) / B):.3f}・最大 {np.max(np.abs(D) / B):.3f};'
              f'  |E1|/B1 最大 {np.max(np.abs(e1) / np.maximum(R[f"B1_{comp}"][m], 1e-300)):.3f}'
              f'  |E2|/B2 最大 {np.max(np.abs(e2) / R[f"B2_{comp}"][m]):.3f}'
              + (f'  |E3|/B3 最大 {np.max(np.abs(e3) / R["B3_y"][m]):.3f}' if comp == 'y' else '')
              + f';  B1/§4.4 の閾値 中央値 {np.median(R[f"B1_{comp}"][m] / R[f"tol_{comp}"][m]):.1f}・|D|/閾値 中央値 {np.median(np.abs(D) / R[f"tol_{comp}"][m]):.1f}')
m = reg == 'wall'
print(f'\n== 面積 (壁の CV): |area_a − A_planar|/A 中央値 {np.median(np.abs(R["areaA_rel"][m])):.2e}・最大 {np.max(np.abs(R["areaA_rel"][m])):.2e};'
      f' |area_b − A_planar|/A 最大 {np.max(np.abs(R["areaB_rel"][m])):.2e}; |ΣW_b,y − A_planar|/A 中央値 {np.median(np.abs(R["Wby_rel"][m])):.2e}')
m = reg == 'interior'
print(f'   内部の CV: |area_a − A|/A 最大 {np.max(np.abs(R["areaA_rel"][m])):.2e}・|area_b − A|/A 最大 {np.max(np.abs(R["areaB_rel"][m])):.2e}'
      f'・E1 = 0 の CV x {np.count_nonzero(R["E1_x"][m] == 0)}・y {np.count_nonzero(R["E1_y"][m] == 0)} / {m.sum()}')
np.savez(args.out, region=reg.astype(str), fail_x=fx, fail_y=fy, **R)
print(f'保存 {args.out}  ({time.time() - t0:.0f} s)')
