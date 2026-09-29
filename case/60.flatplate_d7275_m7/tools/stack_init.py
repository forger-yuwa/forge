#!/usr/bin/env python3
"""上に層を積んだメッシュへ、元メッシュの場を節点座標の一致で「そのまま」コピーする (case/60 T4-0b-H)。

    python3 tools/stack_init.py SRC_RES.h5 SRC_MESH.h5 DST_MESH.h5 [--free-x -0.1 --free-y 0.5]

- 元メッシュの全節点は DST に同じ座標で存在する前提 (gen_mesh.py --top-layers)。一致しない節点があれば失敗する。
- 保存量 (VALUE/ro, roUx, roUy, roUz, roe, roK, roOmega, roY*) を SRC の dtype のままビット単位でコピーする
  (interp_field.py のように原始量から組み直さない)。
- 追加した節点には、SRC の入口列上端の節点 (--free-x, --free-y) の保存量をそのまま入れる (入口の自由流)。
  旧上端の擾乱を最近傍で新上端へ引き延ばさない (codex 2026-09-27)。
- wall_dist は DST のものを残す。
"""
import argparse, sys
import numpy as np
import h5py

KEYS_FIXED = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src_res"); ap.add_argument("src_mesh"); ap.add_argument("dst_mesh")
    ap.add_argument("--free-x", type=float, default=-0.1)
    ap.add_argument("--free-y", type=float, default=0.5)
    a = ap.parse_args()
    with h5py.File(a.src_mesh, "r") as h:
        ns = np.asarray(h["MESH/COORD"]).reshape(-1, 3)
    with h5py.File(a.dst_mesh, "r") as h:
        nd = np.asarray(h["MESH/COORD"]).reshape(-1, 3)
    # 節点 = DOF (node 離散) の対応を座標の完全一致で作る
    key = lambda c: [tuple(r) for r in c]
    pos = {k: i for i, k in enumerate(key(nd))}
    try:
        idx = np.array([pos[k] for k in key(ns)])
    except KeyError:
        sys.exit("SRC の節点に DST で座標が完全一致しないものがある")
    if len(np.unique(idx)) != len(idx):
        sys.exit("対応が一対一でない")
    new = np.ones(len(nd), bool); new[idx] = False
    fi = int(np.argmin((ns[:, 0] - a.free_x) ** 2 + (ns[:, 1] - a.free_y) ** 2))
    print(f"SRC {len(ns)} 節点 → DST {len(nd)} 節点 (追加 {new.sum()})、自由流の元節点 #{fi} {ns[fi]}")
    with h5py.File(a.src_res, "r") as s, h5py.File(a.dst_mesh, "r+") as d:
        keys = KEYS_FIXED + sorted(k for k in s["VALUE"] if k.startswith("roY"))
        for k in keys:
            v = np.asarray(s["VALUE"][k])
            if len(v) != len(ns):
                sys.exit(f"{k}: SRC の長さ {len(v)} ≠ 節点数 {len(ns)}")
            out = np.empty(len(nd), dtype=v.dtype)
            out[idx] = v
            out[new] = v[fi]
            if k in d["VALUE"]:
                del d["VALUE"][k]
            d["VALUE"].create_dataset(k, data=out)
        # 検査: 共通節点がビット一致
        bad = [k for k in keys if not np.array_equal(np.asarray(d["VALUE"][k])[idx], np.asarray(s["VALUE"][k]))]
        ro = np.asarray(d["VALUE"]["ro"], float)
        sy = sum(np.asarray(d["VALUE"][k], float) for k in keys if k.startswith("roY")) / ro if any(k.startswith("roY") for k in keys) else None
        print(f"コピー {len(keys)} 量 (dtype {np.asarray(s['VALUE']['ro']).dtype})、共通節点のビット不一致: {bad or 'なし'}、"
              f"NaN {sum(int(np.isnan(np.asarray(d['VALUE'][k], float)).sum()) for k in keys)}"
              + (f"、ΣY [{sy.min():.15f}, {sy.max():.15f}]" if sy is not None else ""))
        if bad:
            sys.exit(1)


if __name__ == "__main__":
    main()
