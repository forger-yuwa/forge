#!/usr/bin/env python3
"""R7b 4-1 (plan tooling-nozzle-sern-chain §5.1 R7b ④): `L_sw_exact` で x の station を写しただけの格子へ、元格子の保存量を
**節点番号で**写す。座標は少し動くので restart_field (座標の完全一致) も r7a_twin_restart (座標一致) も使えない。

検査してから写す:
- 節点数・hex 接続 (MESH/CONNE) が同一
- 各節点の境界タグの組が同一 (板の内外 = 双子の取り違えが無い)
- 座標のずれは x 成分だけ、かつ |Δx| がその節点の x 方向の局所間隔の半分以下
保存量 (ro, roUx, roUy, roUz, roe, roK, roOmega, roY*) を index コピー、wall_dist は新格子の値を残す。

  python3 r7b_index_restart.py SRC_res.h5 SRC_mesh.h5 DST_input.h5
"""
import sys
import h5py
import numpy as np


def tags(f, n):
    t = [[] for _ in range(n)]
    for k in f["BCONDS"].keys():
        for i in np.unique(np.asarray(f[f"BCONDS/{k}/iCells"])):
            t[int(i)].append(int(k))
    return [tuple(sorted(x)) for x in t]


def main():
    src_res, src_mesh, dst = sys.argv[1:4]
    with h5py.File(src_mesh, "r") as a, h5py.File(dst, "r") as b:
        ca = np.asarray(a["MESH/COORD"], float).reshape(-1, 3); cb = np.asarray(b["MESH/COORD"], float).reshape(-1, 3)
        if ca.shape != cb.shape:
            raise SystemExit(f"節点数が違う: {ca.shape[0]} vs {cb.shape[0]}")
        if not np.array_equal(np.asarray(a["MESH/CONNE"]), np.asarray(b["MESH/CONNE"])):
            raise SystemExit("hex 接続 (MESH/CONNE) が違う")
        n = ca.shape[0]
        if tags(a, n) != tags(b, n):
            raise SystemExit("境界タグの組が違う節点がある (双子の取り違えの恐れ)")
    d = cb - ca
    if np.abs(d[:, 1:]).max() > 0.0:
        raise SystemExit(f"y/z 座標が動いている (最大 {np.abs(d[:, 1:]).max():.3e})")
    xs = np.unique(np.round(ca[:, 0], 12)); dx = np.diff(xs)
    k = np.clip(np.searchsorted(xs, np.round(ca[:, 0], 12)), 1, len(xs) - 1)
    local = np.minimum(dx[k - 1], dx[np.minimum(k, len(dx) - 1)])
    ratio = np.abs(d[:, 0]) / local
    moved = d[:, 0] != 0.0
    if ratio.max() > 0.5:
        raise SystemExit(f"x のずれが局所間隔の半分を超える節点がある (最大比 {ratio.max():.3f})")
    print(f"節点 {n}: 接続・境界タグ同一、動いた節点 {int(moved.sum())} (x {ca[moved,0].min():.5f}–{ca[moved,0].max():.5f} m)、"
          f"最大 |Δx| {np.abs(d[:,0]).max():.3e} m (局所間隔比 {ratio.max():.3f})")
    with h5py.File(src_res, "r") as s, h5py.File(dst, "r+") as f:
        names = [q for q in s["VALUE"] if (q.startswith("ro")) and q in f["VALUE"]]
        for q in names:
            f["VALUE/" + q][...] = np.asarray(s["VALUE/" + q]).astype(f["VALUE/" + q].dtype)
        bad = [q for q in names if not np.array_equal(np.asarray(f["VALUE/" + q]), np.asarray(s["VALUE/" + q]).astype(f["VALUE/" + q].dtype))]
    if bad:
        raise SystemExit(f"写した保存量が一致しない: {bad}")
    print(f"VERDICT: OK ({len(names)} 量を index コピー: {names}; wall_dist は新格子の値)")


if __name__ == "__main__":
    main()
