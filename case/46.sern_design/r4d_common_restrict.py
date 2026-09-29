#!/usr/bin/env python3
"""広い格子 (z_append) の場を狭い格子へ切り出す restart (farfield plan §5.1 #4a、codex diagnose 2026-09-29)。
`r4d_common_restart.py` の逆向き: DST (狭い格子) の**全節点**が SRC (広い格子) に座標で一致する前提で、保存量を直接コピーする。
座標の一致する双子節点 (カウル・側壁のスリット) は所属する境界タグの集合で見分ける。側方遠方面 side_far のタグは比較から除く
(狭い格子の side_far 面は、広い格子では内部の面なのでタグが付かない)。
DST に対応の無い節点があれば失敗させる。最後に DST の保存量が SRC の対応節点とビット一致することを検査する。

  python3 r4d_common_restrict.py SRC_res.h5 SRC_mesh.h5 DST_input.h5
"""
import sys
import h5py, numpy as np

CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]
SIDE_FAR = 10   # PHYS_SERN3D["side_far"]


def keys(f):
    c = np.asarray(f["MESH/COORD"]).reshape(-1, 3)
    tags = [[] for _ in range(c.shape[0])]
    for k in f["BCONDS"].keys():
        if int(k) == SIDE_FAR:
            continue
        for n in np.unique(np.asarray(f[f"BCONDS/{k}/iCells"])):
            tags[int(n)].append(int(k))
    return [(tuple(np.asarray(c[i], dtype=np.float64).tolist()), tuple(sorted(tags[i]))) for i in range(c.shape[0])]


def main():
    src_res, src_mesh, dst = sys.argv[1:4]
    with h5py.File(src_mesh, "r") as f:
        ks = keys(f)
    idx = {}
    for i, k in enumerate(ks):
        if k in idx:
            raise SystemExit(f"SRC に座標・タグとも同じ節点が 2 つある: {k}")
        idx[k] = i
    with h5py.File(dst, "r") as f:
        kd = keys(f)
    miss = [k for k in kd if k not in idx]
    if miss:
        raise SystemExit(f"DST の {len(miss)} 節点が SRC に無い (例 {miss[:3]})")
    m = np.array([idx[k] for k in kd])
    with h5py.File(src_res, "r") as s, h5py.File(dst, "r+") as d:
        for q in CONS:
            v = s["VALUE"][q][:]
            d["VALUE"][q][...] = v[m].astype(d["VALUE"][q].dtype)
        ok = all(np.array_equal(d["VALUE"][q][:], s["VALUE"][q][:][m]) for q in CONS)
    print(f"DST {len(kd)} 節点 ← SRC {len(ks)} 節点 (保存量 {len(CONS)} 量)")
    print("VERDICT: " + ("OK (DST が SRC の対応節点とビット一致)" if ok else "NG"))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
