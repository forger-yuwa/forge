#!/usr/bin/env python3
r"""固体メッシュ (`fem2d`) を **ソルバが読む HDF5** に変換する (CHT Phase 2)。

plan boundary-conjugate-heat-transfer §4.6a / §5.1 #67 ①。
ソルバ内 `fem2d` は固体全節点系を**バンド Cholesky** で解くので、
**帯幅を小さくする並べ替え (RCM) をここで済ませて**から書く (C++ 側は並べ替えをしない)。

入力は外部ループが使っている固体 JSON (`mesh_npz` / `holes` / `k_table` or `k_solid`) で、
**外部ループと同じ素材から同じ問題を作る**ことを保証する (移植の同値試験 §6 V4b(a) の前提)。

出力 (すべて倍精度・0 始まり index):

    MESH/COORD      (N,2)  節点座標 [m]            … RCM 並べ替え済み
    MESH/TRIS       (M,3)  三角形の節点 index      … 同上
    MESH/PERM       (N,)   new[i] = old index      … 元の npz に戻すための帳簿
    IFACE/NODES     (n,)   界面 (ガス側表面) 節点 index
    IFACE/COORD     (n,3)  その座標 (x, y, 0)      … 流体の壁節点と 1 対 1 に突き合わせる
    IFACE/EDGES     (ne,2) 界面の辺
    ROBIN/EDGES     (nr,2) 孔の辺
    ROBIN/H         (nr,)  各辺の h [W/m2K]
    ROBIN/TC        (nr,)  各辺の T_c [K]
    SOLID/K_T, K_V  (nk,)  k_s(T) のテーブル (定数のときは 1 点)

属性 (`/` に付ける): `n_nodes`, `n_tris`, `n_iface`, `bandwidth`, `iface_sha1`, `source_npz`。
`iface_sha1` は **IFACE/COORD を丸めてから取った SHA1** で、ソルバ起動時に
「この固体メッシュはこの流体壁のために作られたか」を 1 行で照合するために使う。

使い方:
  python3 solver_density_cuda/tools/solid_mesh_to_h5.py \
      --solid case/53.c3x_vane_cht/solid_published.json \
      --out   case/53.c3x_vane_cht/mesh/solid_c3x.h5
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import reverse_cuthill_mckee


def adjacency(n_nodes: int, tris: np.ndarray) -> sp.csr_matrix:
    """三角形から節点隣接行列 (対称・0/1) を作る。"""
    e = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    rows = np.concatenate([e[:, 0], e[:, 1], np.arange(n_nodes)])
    cols = np.concatenate([e[:, 1], e[:, 0], np.arange(n_nodes)])
    data = np.ones(len(rows))
    return sp.csr_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes))


def bandwidth(n_nodes: int, tris: np.ndarray) -> int:
    e = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    return int(np.abs(e[:, 0].astype(np.int64) - e[:, 1].astype(np.int64)).max())


def write_solid_h5(out, nodes, tris, outer_edges, robin_e, robin_h, robin_tc,
                   kT, kV, rcm=True, source="(generated)"):
    """固体 FE の HDF5 を書く (RCM 並べ替えつき)。**変換器と試験が同じ経路を通る**ようにする
    ため関数にしてある。戻り値は (perm, bandwidth_before, bandwidth_after)。"""
    nodes = np.asarray(nodes, float)[:, :2]
    tris = np.asarray(tris, int)
    outer_edges = np.asarray(outer_edges, int)
    robin_e = np.asarray(robin_e, int).reshape(-1, 2)
    N = len(nodes)
    bw0 = bandwidth(N, tris)

    perm = np.arange(N) if not rcm else np.asarray(
        reverse_cuthill_mckee(adjacency(N, tris).tocsr(), symmetric_mode=True), int)
    inv = np.empty(N, int)
    inv[perm] = np.arange(N)

    nodes_p, tris_p = nodes[perm], inv[tris]
    outer_p = inv[outer_edges]
    robin_p = inv[robin_e] if len(robin_e) else robin_e
    bw1 = bandwidth(N, tris_p)

    iface = np.array(sorted(set(outer_p.ravel().tolist())), int)
    iface_xyz = np.zeros((len(iface), 3))
    iface_xyz[:, :2] = nodes_p[iface]
    key = np.round(iface_xyz[np.lexsort((iface_xyz[:, 1], iface_xyz[:, 0]))], 9)
    sha = hashlib.sha1(key.tobytes()).hexdigest()

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out, "w") as f:
        f.create_dataset("MESH/COORD", data=nodes_p)
        f.create_dataset("MESH/TRIS", data=tris_p.astype(np.int32))
        f.create_dataset("MESH/PERM", data=perm.astype(np.int32))
        f.create_dataset("IFACE/NODES", data=iface.astype(np.int32))
        f.create_dataset("IFACE/COORD", data=iface_xyz)
        f.create_dataset("IFACE/EDGES", data=outer_p.astype(np.int32))
        f.create_dataset("ROBIN/EDGES", data=robin_p.astype(np.int32))
        f.create_dataset("ROBIN/H", data=np.asarray(robin_h, float))
        f.create_dataset("ROBIN/TC", data=np.asarray(robin_tc, float))
        f.create_dataset("SOLID/K_T", data=np.asarray(kT, float))
        f.create_dataset("SOLID/K_V", data=np.asarray(kV, float))
        for k, v in (("n_nodes", N), ("n_tris", len(tris)), ("n_iface", len(iface)),
                     ("bandwidth", bw1)):
            f.attrs[k] = int(v)
        f.attrs["iface_sha1"] = sha
        # **中身全体のハッシュ** — 界面座標だけでは内部節点の順序・接続・物性・冷却条件を識別できない
        # (codex result 2 巡目 M3: 内部 2 節点を入れ替えても界面ハッシュと節点数は同じで、
        #  保存温度の対応が 214 K ずれた)。再開時はこちらを照合する。
        h = hashlib.sha1()
        for arr in (nodes_p, tris_p, iface, iface_xyz, outer_p, robin_p,
                    np.asarray(robin_h, float), np.asarray(robin_tc, float),
                    np.asarray(kT, float), np.asarray(kV, float)):
            h.update(np.ascontiguousarray(arr).tobytes())
        f.attrs["content_sha1"] = h.hexdigest()
        f.attrs["source_npz"] = str(source)
    return perm, bw0, bw1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solid", required=True, help="固体 JSON (mesh_npz / holes / k_table|k_solid)")
    ap.add_argument("--out", required=True, help="出力 HDF5")
    ap.add_argument("--npz", default=None,
                    help="JSON の mesh_npz を上書きする。**JSON は絶対パスを持つことがあり、別 worktree の\n"
                         "未コミットのメッシュを黙って読んでしまう**ので、worktree 内で回すときは明示する")
    ap.add_argument("--no-rcm", action="store_true", help="RCM を掛けない (帯幅の比較用)")
    a = ap.parse_args()

    spec = json.loads(Path(a.solid).read_text())
    if "mesh_npz" not in spec:
        sys.exit("[solid_mesh_to_h5] この JSON は fem2d 用ではない (mesh_npz が無い)")
    npz_path = Path(a.npz) if a.npz else (Path(a.solid).parent / spec["mesh_npz"]) if not Path(spec["mesh_npz"]).is_absolute() \
        else Path(spec["mesh_npz"])
    if not a.npz and not npz_path.exists():
        npz_path = Path(spec["mesh_npz"])          # JSON が実行ディレクトリ基準で書いている場合
    d = np.load(npz_path)

    nodes = np.asarray(d["nodes"], float)[:, :2]
    tris = np.asarray(d["tris"], int)
    outer_edges = np.asarray(d["outer_edges"], int)
    hole_keys = sorted([k for k in d.files if k.startswith("hole")], key=lambda s_: int(s_[4:]))

    hs = spec["holes"]
    if len(hs) not in (1, len(hole_keys)):
        sys.exit(f"[solid_mesh_to_h5] JSON の holes が {len(hs)} 個、メッシュの孔は {len(hole_keys)} 個")
    robin_e, robin_h, robin_tc = [], [], []
    for i, hk in enumerate(hole_keys):
        hp = hs[i] if len(hs) > 1 else hs[0]
        for (n0, n1) in np.asarray(d[hk], int):
            robin_e.append((int(n0), int(n1)))
            robin_h.append(float(hp["h"]))
            robin_tc.append(float(hp["T_c"]))
    robin_e = np.asarray(robin_e, int).reshape(-1, 2)

    if "k_table" in spec:
        kT = np.asarray(spec["k_table"]["T"] if isinstance(spec["k_table"], dict) else spec["k_table"][0], float)
        kV = np.asarray(spec["k_table"]["k"] if isinstance(spec["k_table"], dict) else spec["k_table"][1], float)
    else:
        kT, kV = np.array([300.0]), np.array([float(spec["k_solid"])])

    # **書き出しは write_solid_h5 に一本化する** (2 経路あると必ずずれる。
    #  実際 content_sha1 を足したとき片方に入っていなかった)。
    out = Path(a.out)
    perm, bw0, bw1 = write_solid_h5(out, nodes, tris, outer_edges, robin_e, robin_h, robin_tc,
                                    kT, kV, rcm=not a.no_rcm, source=str(npz_path))
    with h5py.File(out, "r+") as f:
        f.attrs["source_sha256"] = hashlib.sha256(npz_path.read_bytes()).hexdigest()
        sha = f.attrs["iface_sha1"]
        content = f.attrs["content_sha1"]
        N = int(f.attrs["n_nodes"])
        n_iface = int(f.attrs["n_iface"])
        n_robin = len(f["ROBIN/H"])

    print(f"[solid_mesh_to_h5] {out}")
    print(f"  nodes {N}  tris {len(tris)}  iface {n_iface}  robin edges {n_robin}")
    print(f"  帯幅  {bw0} -> {bw1}" + ("  (RCM 無効)" if a.no_rcm else "  (RCM)"))
    print(f"  バンド Cholesky の作業量 ~ N*b^2 = {N * bw1 * bw1:.3e} flop")
    print(f"  k_s   {kV.min():.3f} .. {kV.max():.3f} W/mK  ({len(kT)} 点)")
    print(f"  iface_sha1 {sha}")
    print(f"  content_sha1 {content}")
    print(f"  source     {npz_path}")


if __name__ == "__main__":
    main()
