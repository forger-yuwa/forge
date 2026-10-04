#!/usr/bin/env python3
"""R7a 写像対照 (plan tooling-nozzle-sern-chain §5.1 R7a、codex diagnose 2026-10-05): 側壁を延ばした B 格子へ、
A 格子 (run_1017 と同じ) の保存量を**双子節点の内外を保存して**写す。

- B の節点を座標で A に対応づける。A 側で座標の一致する既存の双子 (厚さ 0 の板の両側) は、所属する境界タグの組で見分ける。
- B で新しく分かれた節点 (延ばした側壁の両側、A では 1 節点) は、両側ともその A 節点の値 (従来の初期値と同じ)。
- 対応の無い B 節点は無い想定 (あれば失敗)。`wall_dist` は B の値を残す。保存量は index コピー (丸めなし)。
- 既存の interp_field による B の初期場 (CMP_H5) と比べ、変わった節点数・位置・量別最大差を出す。

  python3 r7a_twin_restart.py SRC_res.h5 SRC_mesh.h5 DST_input.h5 [CMP_H5]
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
    return c, [(tuple(np.asarray(c[i], dtype=np.float64).tolist()), tuple(sorted(tags[i]))) for i in range(c.shape[0])]


def main():
    src_res, src_mesh, dst = sys.argv[1:4]
    cmp_h5 = sys.argv[4] if len(sys.argv) > 4 else None
    with h5py.File(src_mesh, "r") as m:
        cs, ks = keys(m)
    with h5py.File(dst, "r") as d:
        cd, kd = keys(d)
    by_c = {}
    for i, (c, t) in enumerate(ks):
        by_c.setdefault(c, []).append((t, i))
    src_of = np.full(len(kd), -1)
    n_twin = n_split = 0
    for j, (c, t) in enumerate(kd):
        cand = by_c.get(c)
        if not cand:
            continue
        if len(cand) == 1:
            src_of[j] = cand[0][1]
        else:   # A 側の既存の双子: タグの組で見分ける
            hit = [i for (tg, i) in cand if tg == t]
            if len(hit) != 1:
                raise SystemExit(f"双子節点をタグで一意に見分けられない: 座標 {c} DST タグ {t} 候補 {cand}")
            src_of[j] = hit[0]; n_twin += 1
    if (src_of < 0).any():
        raise SystemExit(f"A に対応の無い B 節点 {int((src_of < 0).sum())} 個 (想定外)")
    cnt = np.bincount(src_of, minlength=len(ks))
    n_split = int((cnt > 1).sum())
    if (cnt == 0).any():
        raise SystemExit(f"B に対応の無い A 節点 {int((cnt == 0).sum())} 個 (想定外)")
    split_x = cs[cnt > 1][:, 0]
    print(f"B {len(kd)} 節点 ← A {len(ks)} 節点: A の既存の双子をタグで見分けた B 節点 {n_twin}、B で分かれた (A の 1 節点を共有) A 節点 {n_split} 個 "
          f"(x {split_x.min():.6f}–{split_x.max():.6f} m)" if n_split else f"B {len(kd)} ← A {len(ks)}、双子 {n_twin}、分岐 0")
    with h5py.File(src_res, "r") as s, h5py.File(dst, "r+") as d:
        for n in CONS:
            d["VALUE/" + n][...] = np.asarray(s["VALUE/" + n])[src_of].astype(d["VALUE/" + n].dtype)
        bad = [n for n in CONS if not np.array_equal(np.asarray(d["VALUE/" + n]), np.asarray(s["VALUE/" + n])[src_of].astype(d["VALUE/" + n].dtype))]
    if bad:
        raise SystemExit(f"写した保存量が SRC と一致しない: {bad}")
    if cmp_h5:
        with h5py.File(cmp_h5, "r") as c0, h5py.File(dst, "r") as d:
            changed = np.zeros(len(kd), bool)
            for n in CONS:
                a = np.asarray(c0["VALUE/" + n], float); b = np.asarray(d["VALUE/" + n], float)
                diff = np.abs(a - b); changed |= diff > 0
                w = int(np.argmax(diff))
                print(f"  {n:8s}: 従来 IC との最大差 {diff.max():.4e} (節点 {w}, x {cd[w,0]:.5f} y {cd[w,1]:.5f} z {cd[w,2]:.5f})、差のある節点 {int((diff > 0).sum())}")
            xs = cd[changed]
            print(f"  従来 IC と違う節点 {int(changed.sum())} 個" + (f" (x {xs[:,0].min():.4f}–{xs[:,0].max():.4f}, z {xs[:,2].min():.4f}–{xs[:,2].max():.4f})" if changed.any() else ""))
    print("VERDICT: OK (全 B 節点に A の保存量を index コピー、双子の内外を保存)")


if __name__ == "__main__":
    main()
