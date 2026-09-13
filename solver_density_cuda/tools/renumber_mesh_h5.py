#!/usr/bin/env python3
"""forge 入力 h5 (変換済みメッシュ + VALUE) の節点を RCM (reverse Cuthill–McKee) で再番号付けして
gather の局所性 (L2 ヒット率) を上げる。node (median-dual) では節点 = CV なので節点順がそのまま CV 順。
双対面は読み込み時に (min,max) 節点順で組まれるので、節点順を変えるだけで面順も追従する。

  python3 renumber_mesh_h5.py IN.h5 OUT.h5 [--check]

置換対象: /MESH/COORD (節点座標), /MESH/CONNE 系 (要素→節点), /VALUE/* (節点長の配列), 境界の節点参照。
h5 のレイアウトはファイルを読んで自動判定する (未知の節点長データセットは全て並べ替える)。
順列は /MESH/RENUMBER_PERM (new→old) に保存し、res_*.h5 を旧番号へ戻す/移すのに使える。
"""
import sys, argparse, numpy as np, h5py
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("dst"); ap.add_argument("--check", action="store_true")
ap.add_argument("--list", action="store_true", help="データセット一覧だけ表示")
a = ap.parse_args()

src = h5py.File(a.src, "r")
def walk(g, pre=""):
    for k, v in g.items():
        p = pre + "/" + k
        if isinstance(v, h5py.Dataset): yield p, v
        else: yield from walk(v, p)
dsets = list(walk(src))
if a.list:
    for p, v in dsets: print(p, v.shape, v.dtype)
    sys.exit(0)
coord = src["MESH/COORD"][...]
n = coord.size // 3 if coord.ndim == 1 else coord.shape[0]
print("nodes", n)
# 節点隣接: 要素連結から (各要素の節点を全結合)。要素配列は 1 始まり/0 始まりを判定
conne = None
for p, v in dsets:
    if p.endswith("/CONNE") and v.ndim >= 1 and v.size > n:
        conne = (p, v[...]); break
if conne is None:
    raise SystemExit("CONNE dataset not found; use --list")
cp, ce = conne
ce = np.asarray(ce)
print("using", cp, ce.shape, ce.dtype, "min", ce.min(), "max", ce.max())
base = 1 if ce.min() >= 1 and ce.max() == n else 0
rows, cols = [], []
if ce.ndim == 2:
    k = ce.shape[1]
    for i in range(k):
        for j in range(i + 1, k):
            rows.append(ce[:, i] - base); cols.append(ce[:, j] - base)
else:
    raise SystemExit("CONNE is 1-D (packed); extend the reader")
r = np.concatenate(rows); c = np.concatenate(cols)
m = (r >= 0) & (c >= 0) & (r < n) & (c < n)
r, c = r[m], c[m]
A = coo_matrix((np.ones(len(r), np.int8), (r, c)), shape=(n, n)).tocsr()
A = A + A.T
perm = reverse_cuthill_mckee(A, symmetric_mode=True)     # new index i holds old node perm[i]
inv = np.empty(n, np.int64); inv[perm] = np.arange(n)     # old -> new
# 帯域 (局所性) の前後比較
def bandwidth(A, order):
    coo = A.tocoo(); return int(np.abs(order[coo.row] - order[coo.col]).max()), float(np.abs(order[coo.row] - order[coo.col]).mean())
print("bandwidth max/mean before", bandwidth(A, np.arange(n)), "after", bandwidth(A, inv))
dst = h5py.File(a.dst, "w")
def copy_attrs(s, d):
    for k, v in s.attrs.items(): d.attrs[k] = v
for p, v in dsets:
    data = v[...]
    grp = dst.require_group(p.rsplit("/", 1)[0]) if "/" in p.strip("/") else dst
    name = p.rsplit("/", 1)[1]
    if p == "/MESH/COORD":
        data = (data.reshape(-1, 3)[perm]).reshape(data.shape)
    elif p == cp or (name == "CONNE" and data.size > n and data.min() >= base and data.max() - base < n and p.startswith("/MESH")):
        data = inv[data - base] + base
    elif data.ndim == 1 and data.shape[0] == n and p.startswith("/VALUE"):
        data = data[perm]
    elif data.ndim == 1 and data.shape[0] == n and not p.startswith("/VALUE"):
        print("  node-length non-VALUE dataset permuted:", p); data = data[perm]
    elif np.issubdtype(data.dtype, np.integer) and data.size > 0 and data.min() >= base and data.max() - base < n and p.startswith("/MESH") and "BCOND" in p.upper():
        data = inv[data - base] + base
    d = grp.create_dataset(name, data=data, compression=v.compression, compression_opts=v.compression_opts)
    copy_attrs(v, d)
for k, v in src.attrs.items(): dst.attrs[k] = v
for gname in ["MESH", "VALUE"]:
    if gname in src:
        for k, v in src[gname].attrs.items(): dst[gname].attrs[k] = v
dst.create_dataset("MESH/RENUMBER_PERM", data=perm.astype(np.int64))
dst.close(); print("wrote", a.dst)
