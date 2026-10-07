#!/usr/bin/env python3
"""mark_zero_thickness_edges.py (厚さ 0 の板の自由端の重み w の前処理) の単体試験。GPU 不要。

    python3 solver_density_cuda/tools/test_mark_zero_thickness_edges.py

合成の小さい node 変換 h5 (MESH/COORD・PLANES/STRUCT・BCONDS の vizBface・VIZMESH・CELLS/volume) で:
  (1) 3D 共有辺: 2 タグの面が辺を共有する節点だけが E (板の内部の座標一致の別 ID は結合しない)
  (2) 座標一致の別 ID だけ (共有 ID なし) → E が空でエラー
  (3) 有限厚 (上下の面が底面タグで隔てられる) → E が空でエラー
  (4) 端なし (2 タグが離れている) → エラー
  (5) 2D: 2 タグの線が共有する端点
  (6) BFS: 内部面グラフの距離 (経路・格子・スリットの両側)
  (7) タグ誤記・physID 無し・接続なし・rings<1・同じ physID → エラー
  (8) 書込み: w=0 が S、1 がそれ以外、属性 (署名・ハッシュ・タグ・rings・数)、XDMF (Center='Node')
  (9) 格子署名: gzip で保存し直しても同じ / 節点番号を付け替えると変わる / /VALUE・/AUX を変えても同じ
  (10) 照合 (verify): 正しい h5 は通る、w を 1 値改変すると拒否、番号を付け替えた格子に古い /AUX を残すと拒否
  (11) --control-ones は w ≡ 1 で署名つき
"""
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mark_zero_thickness_edges as m  # noqa: E402

TOOL = os.path.join(HERE, "mark_zero_thickness_edges.py")
fail = 0


def check(name, ok, detail=""):
    global fail
    print(("ok   " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    fail += 0 if ok else 1


def raises(fn, *a, **k):
    try:
        fn(*a, **k)
    except m.EdgeMarkError as e:
        return str(e)
    return None


def write_h5(path, coord, edges, bfaces, kinds=None, compression=None, values=True):
    """合成の node 変換 h5。edges = 内部面 (節点対)、bfaces = {physID: [面の節点列, ...]}。"""
    n = coord.shape[0]
    kw = {"compression": compression} if compression else {}
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        g.attrs["nNodes"] = n; g.attrs["nCells"] = n
        g.attrs["nNormalPlanes"] = len(edges); g.attrs["nPlanes"] = len(edges)
        g.attrs["nBconds"] = len(bfaces); g.attrs["nBPlanes"] = 0
        g.create_dataset("COORD", data=coord.astype(np.float32).reshape(-1), **kw)
        st = np.array([[2, a, b, 2, a, b] for a, b in edges], dtype=np.int32).reshape(-1)
        f.create_dataset("PLANES/STRUCT", data=st, **kw)
        f.create_dataset("CELLS/volume", data=np.ones(n, dtype=np.float32), **kw)
        v = f.create_group("VIZMESH"); v.attrs["nVizCells"] = 1; v.attrs["vizCONNE_dim"] = 9
        v.create_dataset("CONNE", data=np.array([9, 0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int32))
        for pid, faces in bfaces.items():
            b = f.create_group(f"BCONDS/{pid}")
            b.attrs["bcondKind"] = (kinds or {}).get(pid, "wall")
            if faces is None:
                continue
            b.create_dataset("vizBfaceSizes", data=np.array([len(x) for x in faces], dtype=np.int32), **kw)
            b.create_dataset("vizBfaceNodes", data=np.array([i for x in faces for i in x], dtype=np.int32), **kw)
        if values:
            f.create_dataset("VALUE/ro", data=np.linspace(1.0, 2.0, n).astype(np.float32), **kw)


# ---------------------------------------------------------------------------
# 3D の合成格子: nx×ny×nz の格子節点 + 板 (j=1, i<=ITE) の下面の複製節点 (座標一致の別 ID)
# ---------------------------------------------------------------------------
NX, NY, NZ, JP, ITE = 6, 4, 3, 1, 2


def gid(i, j, k):
    return (i * NY + j) * NZ + k


def build_3d(share_te=True, finite=False):
    """戻り: coord, edges, bfaces (5: 上面 A, 6: 下面 B, 7: 底面 C (有限厚のみ), 8: 遠い壁)。
    share_te=True: 後縁の線 (i=ITE) は上下で同じ ID (厚さ 0 の自由端)。False: 後縁も複製 (座標一致の別 ID)。
    finite=True: 下面は j=0 の別の節点列 (有限厚)、後縁は底面タグ 7 で閉じる。"""
    coord = [(i, j, k) for i in range(NX) for j in range(NY) for k in range(NZ)]
    ndup = {}
    if not finite:
        for i in range(ITE + 1):
            if i == ITE and share_te:
                continue
            for k in range(NZ):
                ndup[(i, k)] = len(coord); coord.append((i, JP, k))   # 座標一致の別 ID (下面側)
    coord = np.array(coord, dtype=float)

    def lower(i, k):
        if finite:
            return gid(i, JP - 1, k)
        return ndup.get((i, k), gid(i, JP, k))

    edges = set()
    def add(a, b):
        edges.add((min(a, b), max(a, b)))
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                a = gid(i, j, k)
                if i + 1 < NX: add(a, gid(i + 1, j, k))
                if k + 1 < NZ: add(a, gid(i, j, k + 1))
                if j + 1 < NY:
                    # スリット: 板の範囲 (i<=ITE) では j=JP の上面節点と j=JP-1 を結ばず、下面の複製と結ぶ
                    if (not finite) and j + 1 == JP and i <= ITE and (i, k) in ndup:
                        add(a, ndup[(i, k)])
                    else:
                        add(a, gid(i, j + 1, k))
    for (i, k), d in ndup.items():   # 下面の複製の面内の辺
        if (i + 1, k) in ndup: add(d, ndup[(i + 1, k)])
        elif i + 1 <= ITE: add(d, gid(i + 1, JP, k))
        if (i, k + 1) in ndup: add(d, ndup[(i, k + 1)])
    A, B = [], []
    for i in range(ITE):
        for k in range(NZ - 1):
            A.append([gid(i, JP, k), gid(i + 1, JP, k), gid(i + 1, JP, k + 1), gid(i, JP, k + 1)])
            B.append([lower(i, k), lower(i, k + 1), lower(i + 1, k + 1), lower(i + 1, k)])
    bf = {5: A, 6: B, 8: [[gid(NX - 1, 0, 0), gid(NX - 1, 1, 0), gid(NX - 1, 1, 1), gid(NX - 1, 0, 1)]]}
    if finite:
        bf[7] = [[gid(ITE, JP, k), gid(ITE, JP, k + 1), gid(ITE, JP - 1, k + 1), gid(ITE, JP - 1, k)] for k in range(NZ - 1)]
    return coord, sorted(edges), bf


with tempfile.TemporaryDirectory() as td:
    P = lambda name: os.path.join(td, name)
    pairs = [("cowl_in", 5), ("cowl_out", 6)]

    # (1) 3D 共有辺
    c, e, bf = build_3d(share_te=True)
    write_h5(P("s3.h5"), c, e, bf)
    w, info = m.build_weight(P("s3.h5"), pairs, 1)
    te = sorted(gid(ITE, JP, k) for k in range(NZ))
    check("(1) 3D: E = 後縁の線の節点 (共有辺の両端)", info["E"].tolist() == te and info["dim"] == 3, f"E={info['E'].tolist()}")
    dups = [i for i in range(NX * NY * NZ, c.shape[0])]
    check("(1) 3D: 板の内部の座標一致の別 ID は E に入らない", not set(dups) & set(info["E"].tolist()))
    # 3D の辺の向き・面の巡回の違いは関係しない (B の巡回を逆にしても同じ)
    ka = m._edge_keys(m._faces_by_size(np.array([4]), np.array([0, 1, 2, 3])), 10)
    kb = m._edge_keys(m._faces_by_size(np.array([4]), np.array([3, 2, 1, 0])), 10)
    check("(1) 辺のキーは巡回の向きに依らない", np.array_equal(ka, kb))

    # (2) 座標一致の別 ID だけ (後縁も複製)
    c2, e2, bf2 = build_3d(share_te=False)
    write_h5(P("dup.h5"), c2, e2, bf2)
    msg = raises(m.build_weight, P("dup.h5"), pairs, 2)
    check("(2) 座標一致の別 ID だけ → E が空でエラー", msg is not None and "E が空" in msg, msg or "")

    # (3) 有限厚
    c3, e3, bf3 = build_3d(finite=True)
    write_h5(P("thick.h5"), c3, e3, bf3)
    msg = raises(m.build_weight, P("thick.h5"), pairs, 2)
    check("(3) 有限厚 → E が空でエラー", msg is not None and "E が空" in msg, msg or "")

    # (4) 端なし (遠い壁 8 と上面 5)
    msg = raises(m.build_weight, P("s3.h5"), [("cowl_in", 5), ("far", 8)], 2)
    check("(4) 端なし (共有辺なし) → エラー", msg is not None and "E が空" in msg, msg or "")

    # (5) 2D: 線の共有端点
    c5 = np.array([[x, 0.0, 0.0] for x in range(4)] + [[x, 0.0, 0.0] for x in range(3)], dtype=float)  # 4..6 は座標一致の別 ID
    e5 = [(0, 1), (1, 2), (2, 3), (4, 5), (5, 3), (0, 4), (1, 5)]
    bf5 = {5: [[0, 1], [1, 2], [2, 3]], 6: [[4, 5], [5, 3]]}
    write_h5(P("s2.h5"), c5, e5, bf5)
    w5, info5 = m.build_weight(P("s2.h5"), pairs, 1)
    check("(5) 2D: E = 共有端点だけ (座標一致の 0/4・1/5 は結合しない)", info5["E"].tolist() == [3] and info5["dim"] == 2,
          f"E={info5['E'].tolist()}")
    check("(5) 2D: rings 1 の S = {3 と隣接 2, 5}", sorted(info5["S"].tolist()) == [2, 3, 5], f"S={sorted(info5['S'].tolist())}")
    msg = raises(m.edge_nodes, (np.array([2, 4]), np.array([0, 1, 0, 1, 2, 3])), (np.array([2]), np.array([1, 2])), 4)
    check("(5) 線と多角形の混在 → エラー", msg is not None, msg or "")

    # (6) BFS
    path = np.array([[i, i + 1] for i in range(9)], dtype=np.int64)
    ip, ix = m.adjacency(path, 10)
    d = m.bfs_rings(ip, ix, np.array([0]), 3)
    check("(6) BFS 経路: 距離 0..3 と範囲外 -1", d.tolist() == [0, 1, 2, 3] + [-1] * 6, str(d.tolist()))
    d = m.bfs_rings(ip, ix, np.array([4, 5]), 1)
    check("(6) BFS 複数の種", sorted(np.nonzero(d >= 0)[0].tolist()) == [3, 4, 5, 6])
    # 3D スリットの両側に広がる (上面と下面の複製の両方を含む)
    w3, info3 = m.build_weight(P("s3.h5"), pairs, 2)
    S = set(info3["S"].tolist())
    check("(6) 3D rings 2: S は板の上側 (j=JP+1) と下面の複製の両方を含む",
          gid(ITE, JP + 1, 0) in S and any(dd in S for dd in dups), f"|S|={len(S)}")
    ip3, ix3 = m.adjacency(np.array(e, dtype=np.int64), c.shape[0])
    dist = m.bfs_rings(ip3, ix3, info3["E"], 2)
    check("(6) 3D rings 2: S = 距離 ≤ 2 の節点", S == set(np.nonzero(dist >= 0)[0].tolist()))

    # (7) 入力の拒否
    msg = raises(m.build_weight, P("s3.h5"), [("cowl_in", 5), ("x", 99)], 2)
    check("(7) physID が h5 に無い → エラー", msg is not None and "99" in msg, msg or "")
    write_h5(P("noviz.h5"), c, e, {5: bf[5], 6: None})
    msg = raises(m.build_weight, P("noviz.h5"), pairs, 2)
    check("(7) 元の境界面の接続 (vizBfaceNodes) 無し → エラー", msg is not None and "vizBfaceNodes" in msg, msg or "")
    msg = raises(m.build_weight, P("s3.h5"), pairs, 0)
    check("(7) rings 0 → エラー", msg is not None and "rings" in msg, msg or "")
    with open(P("bcondConfig.yaml"), "w") as fh:
        fh.write("cowl_in: {physID: 5, kind: wall}\ncowl_out: {physID: 6, kind: wall}\n")
    msg = raises(m.resolve_tags, None, "cowl_inn,cowl_out", P("bcondConfig.yaml"))
    check("(7) タグ誤記 → エラー", msg is not None and "cowl_inn" in msg, msg or "")
    check("(7) --tags A,B と A B の両形式", m.resolve_tags(None, "cowl_in,cowl_out", P("bcondConfig.yaml")) == pairs
          and m.resolve_tags(None, ["cowl_in", "cowl_out"], P("bcondConfig.yaml")) == pairs)
    msg = raises(m.resolve_tags, ["a=5", "b=5"], None, None)
    check("(7) 同じ physID → エラー", msg is not None, msg or "")
    msg = raises(m.resolve_tags, ["a=5"], None, None)
    check("(7) タグが 1 つ → エラー", msg is not None, msg or "")
    r = subprocess.run([sys.executable, TOOL, P("s3.h5"), "--tags", "cowl_in", "cowl_outt", "--bcond-config", P("bcondConfig.yaml")],
                       capture_output=True, text=True)
    check("(7) CLI: タグ誤記は非零終了で /AUX を書かない", r.returncode != 0 and "AUX" not in h5py.File(P("s3.h5"), "r"),
          r.stderr.strip()[-80:])
    # 非 node (VIZMESH なし) の h5
    write_h5(P("cellish.h5"), c, e, bf)
    with h5py.File(P("cellish.h5"), "a") as f:
        del f["VIZMESH"]
    msg = raises(m.build_weight, P("cellish.h5"), pairs, 2)
    check("(7) node 変換でない h5 → エラー", msg is not None and "node" in msg, msg or "")

    # (8) 書込み
    r = subprocess.run([sys.executable, TOOL, P("s3.h5"), "--tag", "cowl_in=5", "--tag", "cowl_out=6", "--rings", "2"],
                       capture_output=True, text=True)
    check("(8) CLI 書込み 成功", r.returncode == 0, (r.stdout + r.stderr).strip()[-200:])
    with h5py.File(P("s3.h5"), "r") as f:
        ds = f["AUX/w_recon_vel"]; wv = ds[()]; at = dict(ds.attrs)
    check("(8) w: S で 0、それ以外で 1、float32", wv.dtype == np.float32 and set(np.nonzero(wv == 0)[0].tolist()) == S
          and np.all((wv == 0) | (wv == 1)))
    need = {"generator", "generator_version", "mode", "tags", "tag_physids", "rings", "n_edge_nodes", "n_marked_nodes",
            "n_nodes", "mesh_signature", "mesh_signature_version", "field_sha256", "created"}
    check("(8) 属性が揃う", need <= set(at), str(sorted(need - set(at))))
    check("(8) 属性の値 (タグ・rings・数・ハッシュ)", at["tags"] == "cowl_in,cowl_out" and at["tag_physids"] == "cowl_in:5,cowl_out:6"
          and int(at["rings"]) == 2 and int(at["n_edge_nodes"]) == NZ and int(at["n_marked_nodes"]) == len(S)
          and at["field_sha256"] == m.field_hash(wv) and at["mesh_signature"] == m.mesh_signature(P("s3.h5")))
    xmf = P("s3_w_recon_vel.xmf")
    txt = open(xmf).read() if os.path.exists(xmf) else ""
    check("(8) XDMF: /AUX/w_recon_vel を Center='Node' で参照", "Center='Node'" in txt and "s3.h5:AUX/w_recon_vel" in txt
          and "s3.h5:VIZMESH/CONNE" in txt)

    # (9) 格子署名
    sig = m.mesh_signature(P("s3.h5"))
    write_h5(P("s3_gz.h5"), c, e, bf, compression="gzip")
    check("(9) gzip で保存し直しても署名は同じ", m.mesh_signature(P("s3_gz.h5")) == sig)
    write_h5(P("s3_noval.h5"), c, e, bf, values=False)
    check("(9) /VALUE が無くても (= 場を変えても) 署名は同じ", m.mesh_signature(P("s3_noval.h5")) == sig)
    check("(9) /AUX を書いた後も署名は同じ", m.mesh_signature(P("s3.h5")) == sig)
    # 節点番号の付け替え (座標・接続を整合させた別番号の同じ格子)
    perm = np.random.default_rng(1).permutation(c.shape[0])        # 新番号 k ← 旧番号 perm[k]
    inv = np.empty_like(perm); inv[perm] = np.arange(perm.size)    # 旧 → 新
    e_r = [(int(inv[a]), int(inv[b])) for a, b in e]
    bf_r = {pid: [[int(inv[i]) for i in x] for x in fs] for pid, fs in bf.items()}
    write_h5(P("s3_renum.h5"), c[perm], e_r, bf_r)
    sig_r = m.mesh_signature(P("s3_renum.h5"))
    check("(9) 節点番号を付け替えると署名が変わる", sig_r != sig)
    w_r, info_r = m.build_weight(P("s3_renum.h5"), pairs, 2)
    check("(9) 付け替えた格子でも E・S は同じ節点 (番号は写像)", sorted(inv[info3["S"]].tolist()) == sorted(info_r["S"].tolist())
          and sorted(inv[info3["E"]].tolist()) == sorted(info_r["E"].tolist()))
    check("(9) field_hash は float32 の値で決まる", m.field_hash(np.ones(3)) == m.field_hash(np.ones(3, dtype=np.float32))
          and m.field_hash(np.ones(3)) != m.field_hash(np.array([1, 1, 0.5])))

    # (10) 照合
    ok = m.verify(P("s3.h5"))
    check("(10) 正しい h5 は verify を通る", ok["n_w_lt1"] == len(S))
    with h5py.File(P("s3.h5"), "r") as f, h5py.File(P("s3_renum.h5"), "a") as g:   # 古い /AUX をそのまま残す
        f.copy(f["AUX"], g, "AUX")
    msg = raises(m.verify, P("s3_renum.h5"))
    check("(10) 番号を付け替えた格子 + 古い /AUX → 拒否", msg is not None and "mesh_signature" in msg, msg or "")
    with h5py.File(P("s3_gz.h5"), "a") as g, h5py.File(P("s3.h5"), "r") as f:
        f.copy(f["AUX"], g, "AUX")
    check("(10) gzip で保存し直した同じ格子 + /AUX は通る", m.verify(P("s3_gz.h5"))["mesh_signature"] == sig)
    with h5py.File(P("s3_gz.h5"), "a") as g:
        ds = g["AUX/w_recon_vel"]; v = ds[()]; v[0] = np.float32(0.5) if v[0] != 0.5 else np.float32(1.0); ds[...] = v
    msg = raises(m.verify, P("s3_gz.h5"))
    check("(10) w を 1 値だけ改変 → 拒否 (field_sha256)", msg is not None and "field_sha256" in msg, msg or "")

    # (11) 対照 (w ≡ 1)
    write_h5(P("ones.h5"), c, e, bf)
    r = subprocess.run([sys.executable, TOOL, P("ones.h5"), "--control-ones", "--no-xdmf"], capture_output=True, text=True)
    with h5py.File(P("ones.h5"), "r") as f:
        wv = f["AUX/w_recon_vel"][()]; mode = f["AUX/w_recon_vel"].attrs["mode"]
    check("(11) --control-ones: w ≡ 1・署名つき・verify を通る", r.returncode == 0 and np.all(wv == 1.0) and mode == "control_ones"
          and m.verify(P("ones.h5"))["n_w_lt1"] == 0)
    r = subprocess.run([sys.executable, TOOL, P("ones.h5"), "--control-ones", "--field", "other"], capture_output=True, text=True)
    check("(11) --field は w_recon_vel だけ", r.returncode != 0)

print("RESULT:", "PASS" if fail == 0 else f"FAIL ({fail})")
sys.exit(1 if fail else 0)
