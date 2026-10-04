#!/usr/bin/env python3
r"""`solid_mesh_to_h5.py` の自己検証 (CHT Phase 2, plan §6 V4b (a) a5)。

C++ 側はこの HDF5 だけを読んで固体を組むので、**変換で問題が変わっていないこと**を
ここで閉じておく。並べ替え (RCM) の帳簿が 1 か所ずれるだけで、孔が別の場所に付いたり
界面が入れ替わったりするが、**解はもっともらしいまま**なので後段では気づけない。

試験:
  T1 帳簿      : PERM が置換であること、逆写像で元の npz の座標・三角形・孔辺に戻ること
  T2 界面      : IFACE/NODES が IFACE/EDGES の節点集合と一致し、座標が元の界面と集合として同じ
  T3 解の同値  : 元の npz から組んだ作用素と、h5 から組んだ作用素で、**同じ界面温度に対する
                 全節点温度が PERM 越しに一致する** (max|ΔT| ≤ 1e-9 K)
  T4 帯幅      : RCM 後の帯幅が h5 の属性と一致し、元より小さい

使い方:
  python3 solver_density_cuda/tools/test_solid_mesh_to_h5.py \
      --solid case/53.c3x_vane_cht/solid_published.json \
      --npz   case/53.c3x_vane_cht/mesh/solid_c3x.npz \
      --h5    case/53.c3x_vane_cht/mesh/solid_c3x.h5
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from solid_fem2d import Fem2DOperator  # noqa: E402
from solid_mesh_to_h5 import bandwidth  # noqa: E402

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = ""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def build_from_npz(spec, d):
    nodes, tris = np.asarray(d["nodes"], float), np.asarray(d["tris"], int)
    outer = np.asarray(d["outer_edges"], int)
    iface = np.array(sorted(set(outer.ravel().tolist())), int)
    hole_keys = sorted([k for k in d.files if k.startswith("hole")], key=lambda s_: int(s_[4:]))
    hs = spec["holes"]
    robin = []
    for i, hk in enumerate(hole_keys):
        hp = hs[i] if len(hs) > 1 else hs[0]
        for (n0, n1) in np.asarray(d[hk], int):
            robin.append((int(n0), int(n1), float(hp["h"]), float(hp["T_c"])))
    k_solid = (np.asarray(spec["k_table"][0], float), np.asarray(spec["k_table"][1], float)) \
        if "k_table" in spec else float(spec["k_solid"])
    return Fem2DOperator(nodes, tris, iface, [tuple(e) for e in outer], robin, k_solid)


def build_from_h5(f):
    nodes = np.asarray(f["MESH/COORD"][:], float)
    tris = np.asarray(f["MESH/TRIS"][:], int)
    iface = np.asarray(f["IFACE/NODES"][:], int)
    outer = np.asarray(f["IFACE/EDGES"][:], int)
    re_, rh, rtc = f["ROBIN/EDGES"][:], f["ROBIN/H"][:], f["ROBIN/TC"][:]
    robin = [(int(a), int(b), float(h), float(t)) for (a, b), h, t in zip(re_, rh, rtc)]
    kT, kV = np.asarray(f["SOLID/K_T"][:], float), np.asarray(f["SOLID/K_V"][:], float)
    k_solid = (kT, kV) if len(kT) > 1 else float(kV[0])
    return Fem2DOperator(nodes, tris, iface, [tuple(e) for e in outer], robin, k_solid)


def solve_all(op, T_iface):
    """$k_s(T)$ を内部温度と自己整合させてから全節点温度を返す。"""
    u = None
    for _ in range(30):
        op.assemble(T_iface)
        u_new = op.recover_interior(T_iface)
        if u is not None and np.max(np.abs(u_new - u)) < 1e-12:
            return u_new
        u = u_new
    return u


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solid", required=True)
    ap.add_argument("--npz", required=True)
    ap.add_argument("--h5", required=True)
    a = ap.parse_args()

    spec = json.loads(Path(a.solid).read_text())
    d = np.load(a.npz)
    f = h5py.File(a.h5, "r")
    perm = np.asarray(f["MESH/PERM"][:], int)
    N = len(perm)
    inv = np.empty(N, int)
    inv[perm] = np.arange(N)

    print(f"[test_solid_mesh_to_h5] {a.h5}  N={N}")

    # T1 帳簿
    check("T1a PERM は置換", sorted(perm.tolist()) == list(range(N)))
    check("T1b 座標が戻る",
          np.array_equal(np.asarray(f["MESH/COORD"][:], float)[inv],
                         np.asarray(d["nodes"], float)[:, :2]))
    check("T1c 三角形が戻る",
          np.array_equal(perm[np.asarray(f["MESH/TRIS"][:], int)], np.asarray(d["tris"], int)))
    hole_keys = sorted([k for k in d.files if k.startswith("hole")], key=lambda s_: int(s_[4:]))
    holes_orig = np.vstack([np.asarray(d[k], int) for k in hole_keys])
    check("T1d 孔の辺が戻る",
          np.array_equal(perm[np.asarray(f["ROBIN/EDGES"][:], int)], holes_orig))

    # T2 界面
    ifn = np.asarray(f["IFACE/NODES"][:], int)
    ife = np.asarray(f["IFACE/EDGES"][:], int)
    check("T2a IFACE/NODES = 辺の節点集合", set(ifn.tolist()) == set(ife.ravel().tolist()))
    xy_h5 = np.asarray(f["IFACE/COORD"][:], float)[:, :2]
    outer_orig = np.asarray(d["outer_edges"], int)
    xy_np = np.asarray(d["nodes"], float)[sorted(set(outer_orig.ravel().tolist())), :2]
    check("T2b 界面座標が集合として一致",
          np.allclose(np.sort(xy_h5, axis=0), np.sort(xy_np, axis=0), atol=1e-12))

    # T3 解の同値
    op_n = build_from_npz(spec, d)
    op_h = build_from_h5(f)
    rng = np.random.default_rng(0)
    T0 = 500.0 + 150.0 * rng.random(len(op_n.iface))          # 非一様な界面温度
    u_n = solve_all(op_n, T0)
    # h5 側の界面順は座標で対応づける (並べ替えで順序が変わるため)
    pos = {tuple(np.round(p, 12)): i for i, p in enumerate(op_n.coords[:, :2])}
    idx = np.array([pos[tuple(np.round(p, 12))] for p in op_h.coords[:, :2]], int)
    u_h = solve_all(op_h, T0[idx])
    dmax = float(np.max(np.abs(u_h[inv] - u_n)))
    check("T3 全節点温度が PERM 越しに一致", dmax <= 1e-9, f"max|ΔT| = {dmax:.3e} K")

    # T4 帯幅
    bw_h5 = bandwidth(N, np.asarray(f["MESH/TRIS"][:], int))
    bw_np = bandwidth(N, np.asarray(d["tris"], int))
    check("T4 帯幅が属性と一致し元より小さい",
          bw_h5 == int(f.attrs["bandwidth"]) and bw_h5 < bw_np, f"{bw_np} -> {bw_h5}")

    print(f"\nVERDICT: {'PASS (all)' if not FAILS else 'FAIL: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
