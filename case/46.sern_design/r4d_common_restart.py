#!/usr/bin/env python3
"""R4d (plan tooling-nozzle-sern-3d §5.1 R4d) の B 用 restart: 共通領域は保存量をそのまま写し、追加した外側だけ最近傍にする。

`mesh3d.z_append` の格子は共通領域の座標・接続が元格子と同一だが、節点番号は並びが変わる。
`interp_field.py` は (1) 原始量から保存量を組み直すので丸めでずれ、(2) 座標の一致する双子節点 (カウル・側壁のスリット) を
同じ元節点に写す (反対側の値を拾う)。そこで節点を座標で対応づけ、座標の一致する双子だけ**所属する境界タグの集合** (side_far を除く) で見分けて、共通領域の保存量を
index コピーする。対応の無い節点 (追加した外側) は DST に既にある値 (interp_field の最近傍) を残す。
最後に共通領域が SRC とビット一致することを検査し、しなければ失敗させる。

  python3 r4d_common_restart.py SRC_res.h5 SRC_mesh.h5 DST_input.h5
  (SRC_mesh は SRC と同じ格子の入力 h5。BCONDS の取得に使う)
"""
import sys
import h5py
import numpy as np

CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]


def keys(f):
    c = np.asarray(f["MESH/COORD"]).reshape(-1, 3)
    tags = [[] for _ in range(c.shape[0])]
    for k in f["BCONDS"].keys():
        for n in np.unique(np.asarray(f[f"BCONDS/{k}/iCells"])):
            tags[int(n)].append(int(k))
    return [(tuple(np.asarray(c[i], dtype=np.float64).tolist()), tuple(sorted(tags[i]))) for i in range(c.shape[0])]


def main():
    src_res, src_mesh, dst = sys.argv[1:4]
    SIDE_FAR = 10   # PHYS_SERN3D["side_far"]: A の遠方面は B では内部になるので照合から外す
    with h5py.File(src_mesh, "r") as m:
        ks = keys(m)
    with h5py.File(dst, "r") as d:
        kd = keys(d)
    by_c = {}
    for i, (c, t) in enumerate(ks):
        by_c.setdefault(c, []).append((tuple(x for x in t if x != SIDE_FAR), i))
    src_of = np.full(len(kd), -1)
    used = set()
    for j, (c, t) in enumerate(kd):
        cand = by_c.get(c)
        if not cand:
            continue
        if len(cand) == 1:
            src_of[j] = cand[0][1]
        else:   # 座標の一致する双子: タグ (side_far を除く) で見分ける
            tt = tuple(x for x in t if x != SIDE_FAR)
            hit = [i for (tg, i) in cand if tg == tt]
            if len(hit) != 1:
                raise SystemExit(f"双子節点をタグで一意に見分けられない: 座標 {c} DST タグ {t} 候補 {cand}")
            src_of[j] = hit[0]
    common = src_of >= 0
    if len(set(src_of[common].tolist())) != int(common.sum()):
        raise SystemExit("SRC の同じ節点に DST の複数節点が対応した")
    print(f"DST {len(kd)} 節点のうち共通 {int(common.sum())} (SRC {len(ks)})、追加 {int((~common).sum())}")
    if int(common.sum()) != len(ks):
        raise SystemExit("SRC の全節点が DST に対応しない")
    with h5py.File(src_res, "r") as s, h5py.File(dst, "r+") as d:
        for n in CONS:
            if n not in s["VALUE"] or n not in d["VALUE"]:
                raise SystemExit(f"保存量 {n} が無い")
            sv = np.asarray(s["VALUE/" + n]); dv = np.asarray(d["VALUE/" + n])
            dv[common] = sv[src_of[common]].astype(dv.dtype)
            d["VALUE/" + n][...] = dv
        bad = [n for n in CONS if not np.array_equal(np.asarray(d["VALUE/" + n])[common], np.asarray(s["VALUE/" + n])[src_of[common]].astype(d["VALUE/" + n].dtype))]
    if bad:
        raise SystemExit(f"共通領域がビット一致しない: {bad}")
    print(f"VERDICT: OK (共通領域 {int(common.sum())} 節点の保存量 {len(CONS)} 量が SRC とビット一致、追加 {int((~common).sum())} 節点は DST の値のまま)")


if __name__ == "__main__":
    main()
