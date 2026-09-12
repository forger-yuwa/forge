#!/usr/bin/env python3
"""旧番号の res_*.h5 (または変換 h5 の VALUE) を、RCM 再番号付けした新メッシュ h5 の VALUE へ index 移植する。
  python3 permute_res_h5.py SRC_res.h5 DST_mesh.h5
DST の /MESH/RENUMBER_PERM (new->old) を使い、DST/VALUE/k[new] = SRC/VALUE/k[perm[new]] を書く (wall_dist は移植しない)。
同一メッシュ・同一 gmsh 節点順が前提 (旧 h5 が再番号付け無しで変換されたもの)。"""
import sys, h5py, numpy as np
src = h5py.File(sys.argv[1], "r"); dst = h5py.File(sys.argv[2], "r+")
n = len(dst["VALUE/ro"])
perm = dst["MESH/RENUMBER_PERM"][...] if "MESH/RENUMBER_PERM" in dst else np.arange(n)   # 無ければ恒等 (index コピー)
assert len(perm) == n, (len(perm), n)
nk = 0; created = []
for k in src["VALUE"]:
    if k == "wall_dist" or len(src["VALUE"][k]) != n: continue
    if k not in dst["VALUE"]:
        # 保存量 (ro*: roY0/roY1 など, 変換 h5 に無い化学種) は作る。派生量は作らない。
        if not k.startswith("ro"): continue
        dst["VALUE"].create_dataset(k, data=src["VALUE"][k][...][perm]); created.append(k); nk += 1; continue
    dst["VALUE"][k][:] = src["VALUE"][k][...][perm]; nk += 1
if created: print("created", created)
print("permuted", nk, "fields into", sys.argv[2])
