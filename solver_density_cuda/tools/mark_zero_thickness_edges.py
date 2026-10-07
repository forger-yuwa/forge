#!/usr/bin/env python3
"""厚さ 0 の板の自由端の近傍に、速度再構成の重み w (節点ごと、0/1) を付ける前処理。

    python3 solver_density_cuda/tools/mark_zero_thickness_edges.py MESH.h5 --tag cowl_in=5 --tag cowl_out=6 [--rings 2]
    python3 solver_density_cuda/tools/mark_zero_thickness_edges.py MESH.h5 --tags cowl_in,cowl_out --bcond-config bcondConfig.yaml
    python3 solver_density_cuda/tools/mark_zero_thickness_edges.py MESH.h5 --control-ones     # 対照用 (w ≡ 1)
    python3 solver_density_cuda/tools/mark_zero_thickness_edges.py MESH.h5 --verify           # 書かずに照合だけ

plan: plans/active/convection-zero-thickness-edge-reconstruction.md §4 (端集合 E・S_2)。
ソルバは端を判定しない。この道具が node 変換済みの入力 h5 (convertGmshToForge の出力) に
`/AUX/w_recon_vel` (float32, 節点 = node の CV 順, ghost を含まない) を書き、ソルバは
`space.zeroThicknessEdgeVelocity: {enabled: 1}` のとき SLAU の内部面で節点 i 側の速度を
u_f = u_i + w_i (u_f,old − u_i) にする (w=1 で従来、w=0 で節点値)。

**端集合 E** (座標一致だけで別 ID の節点を結合しない): 指定した 2 タグ (境界条件名、physID を明示) の
**元の境界面の接続** `/BCONDS/<physID>/vizBfaceNodes` (+ `vizBfaceSizes`) から作る。node 変換後の境界面
(PLANES の半割面) は 1 節点なので使えない。
  - 3D (境界面が 3 節点以上): 2 タグの面が**共有する辺** (面の隣接 2 節点) の両端の節点。
  - 2D (境界面が 2 節点の線): 2 タグの線が**共有する端点**。
**S** = 内部面 (`/PLANES/STRUCT` の先頭 nNormalPlanes 面の隣接 CV 対) の接続グラフで E から距離 rings 以下の節点。
w = 0 (S)、1 (それ以外)。初版は 0/1 だけを書く。

拒否 (エラー終了): node 変換の h5 でない / タグの physID が h5 に無い / 元の境界面の接続が無い・空 /
2 タグで次元が違う / E が空 / rings < 1。

**格子との結び付け**: 属性 `mesh_signature` (版 `zte-mesh-sig-v1`) と `field_sha256` を書き、ソルバは読込時に
両方を再計算して照合する (不一致は起動時エラー)。格子署名は HDF5 ファイルのバイト列でなく、ソルバが読む配列の値から作る
(圧縮・チャンク配置では変わらず、/VALUE・/AUX は含めないので restart_field/interp_field で場を写しても変わらない)。
署名の定義はこのファイルの `mesh_signature()` と `solver_density_cuda/input/zeroThicknessEdge.cpp` の `meshSignature`
で同一 (下の docstring)。配線 (stage_key 等) は `mesh_signature(h5)` と `field_hash(w)` を import して使う。

/VALUE でなく /AUX に置く理由: restart_field.py / interp_field.py は /VALUE だけを書き換えるので、同じ格子の継続では
重みがそのまま残り、格子を変えたら新しい h5 に /AUX が無い (ソルバが起動時に止める) → 道具の回し直しを強制できる。

出力: `/AUX/w_recon_vel` と属性、ParaView 用の `<h5 の stem>_w_recon_vel.xmf` (VIZMESH の primal セル + Center='Node')。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import os
import sys

import h5py
import numpy as np

GENERATOR = "mark_zero_thickness_edges.py"
GENERATOR_VERSION = "1"
SIG_VERSION = "zte-mesh-sig-v1"
FIELD = "w_recon_vel"
AUX_PATH = "AUX/" + FIELD


class EdgeMarkError(RuntimeError):
    """入力の不備 (タグ無し・接続なし・E が空など)。CLI はこれで非零終了する。"""


# ---------------------------------------------------------------------------
# 格子署名とハッシュ (ソルバ input/zeroThicknessEdge.cpp の meshSignature と同一の定義)
# ---------------------------------------------------------------------------
def _line(name: str, dtype: str, shape: str, data: bytes) -> str:
    return f"{name} {dtype} {shape} {hashlib.sha256(data).hexdigest()}\n"


def read_internal_face_cells(f: h5py.File) -> np.ndarray:
    """内部面 ip < nNormalPlanes の (iCells[0], iCells[1]) を int64 (nNormalPlanes, 2) で返す。
    `/PLANES/STRUCT` は面ごとに [nNodes, nodes..., nCells, cells...] (mesh.cpp の readMesh と同じ読み方)。"""
    nN = int(f["MESH"].attrs["nNormalPlanes"])
    s = np.asarray(f["PLANES/STRUCT"][()], dtype=np.int64)
    # 速い経路: node 変換の内部面は [2, a, b, 2, c0, c1] の 6 語が並ぶ (先頭から帰納的に検査するので誤読しない)
    if s.size >= 6 * nN:
        s6 = s[:6 * nN].reshape(nN, 6)
        if nN == 0 or (np.all(s6[:, 0] == 2) and np.all(s6[:, 3] == 2)):
            return np.ascontiguousarray(s6[:, 4:6])
    out = np.empty((nN, 2), dtype=np.int64)
    p = 0
    for ip in range(nN):
        nn = int(s[p]); p += 1 + nn
        nc = int(s[p]); p += 1
        if nc != 2:
            raise EdgeMarkError(f"内部面 {ip} の隣接 CV が {nc} 個 (2 であるべき)")
        out[ip, 0] = s[p]; out[ip, 1] = s[p + 1]
        p += nc
    return out


def read_bface(f: h5py.File, pid: int):
    """physID の元の境界面接続 (sizes, flat nodes) を int64 で返す。無ければ (None, None)。"""
    g = f[f"BCONDS/{pid}"]
    if "vizBfaceSizes" not in g or "vizBfaceNodes" not in g:
        return None, None
    return (np.asarray(g["vizBfaceSizes"][()], dtype=np.int64).reshape(-1),
            np.asarray(g["vizBfaceNodes"][()], dtype=np.int64).reshape(-1))


def discretization_of(f: h5py.File) -> str:
    """h5 が node 変換 (median-dual) なら "node"。VIZMESH が在り CV 数 = 節点数のときだけ node とみなす。"""
    a = f["MESH"].attrs
    if "VIZMESH" in f and int(a["nCells"]) == int(a["nNodes"]):
        return "node"
    return "cell"


def mesh_signature(h5) -> str:
    """格子署名 (版 zte-mesh-sig-v1、16 進 64 桁)。h5 はパスか開いた h5py.File。

    項目 (この順)。各項目の行 = "<名前> <型> <形 (カンマ区切り)> <little-endian の値のバイト列の SHA-256>\\n":
      discretization               str  (長さ)        "node" / "cell"
      dim                          i8   (1)           vizBface の最大節点数が 2 → 2、3 以上 → 3、どの bcond にも無ければ 0
      coord                        f8   (nNodes,3)    /MESH/COORD を倍精度にしたもの (節点順)
      internal_face_cells          i8   (nNormalPlanes,2)
      bcond/<physID>/vizBfaceSizes i8   (n)           physID の昇順。無ければ長さ 0
      bcond/<physID>/vizBfaceNodes i8   (m)
    署名 = SHA-256("zte-mesh-sig-v1\\n" + 全行)。/VALUE・/AUX・ファイルのバイト列は使わない。"""
    if not isinstance(h5, h5py.File):
        with h5py.File(h5, "r") as f:
            return mesh_signature(f)
    f = h5
    disc = discretization_of(f).encode()
    text = SIG_VERSION + "\n"
    text += _line("discretization", "str", str(len(disc)), disc)
    pids = sorted(int(k) for k in f["BCONDS"].keys())
    faces = {pid: read_bface(f, pid) for pid in pids}
    max_face = 0
    for sz, _ in faces.values():
        if sz is not None and sz.size:
            max_face = max(max_face, int(sz.max()))
    dim = 0 if max_face == 0 else (2 if max_face == 2 else 3)
    text += _line("dim", "i8", "1", np.asarray([dim], dtype="<i8").tobytes())
    nNodes = int(f["MESH"].attrs["nNodes"])
    coord = np.asarray(f["MESH/COORD"][()]).reshape(-1)
    if coord.size != 3 * nNodes:
        raise EdgeMarkError(f"/MESH/COORD の長さ {coord.size} が 3 x nNodes ({3 * nNodes}) でない")
    text += _line("coord", "f8", f"{nNodes},3", coord.astype("<f8").tobytes())
    fc = read_internal_face_cells(f)
    text += _line("internal_face_cells", "i8", f"{fc.shape[0]},2", fc.astype("<i8").tobytes())
    for pid in pids:
        sz, nd = faces[pid]
        if sz is None:
            sz = np.zeros(0, dtype=np.int64); nd = np.zeros(0, dtype=np.int64)
        text += _line(f"bcond/{pid}/vizBfaceSizes", "i8", str(sz.size), sz.astype("<i8").tobytes())
        text += _line(f"bcond/{pid}/vizBfaceNodes", "i8", str(nd.size), nd.astype("<i8").tobytes())
    return hashlib.sha256(text.encode()).hexdigest()


def field_hash(w) -> str:
    """重みの SHA-256 (float32 little-endian のバイト列)。ソルバの field_sha256 照合と同じ。"""
    return hashlib.sha256(np.ascontiguousarray(np.asarray(w, dtype="<f4")).tobytes()).hexdigest()


# ---------------------------------------------------------------------------
# 端集合 E と距離 rings の集合 S
# ---------------------------------------------------------------------------
def _faces_by_size(sizes: np.ndarray, nodes: np.ndarray):
    """可変長の面接続を {節点数: (nf, s) 配列} に分ける。"""
    if sizes.size == 0:
        return {}
    if int(sizes.sum()) != nodes.size:
        raise EdgeMarkError(f"vizBfaceSizes の総和 {int(sizes.sum())} と vizBfaceNodes の長さ {nodes.size} が違う")
    off = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    out = {}
    for s in np.unique(sizes):
        idx = off[sizes == s]
        out[int(s)] = nodes[idx[:, None] + np.arange(int(s))[None, :]]
    return out


def _edge_keys(by_size: dict, n_nodes: int) -> np.ndarray:
    """多角形の辺 (隣接 2 節点、閉じる) を min*N+max の int64 キーで返す (重複除去済み)。"""
    keys = []
    for s, F in by_size.items():
        if s < 3:
            raise EdgeMarkError(f"3D の辺の抽出に {s} 節点の面が混ざっている")
        a = F; b = np.roll(F, -1, axis=1)
        lo = np.minimum(a, b).reshape(-1); hi = np.maximum(a, b).reshape(-1)
        keys.append(lo * np.int64(n_nodes) + hi)
    return np.unique(np.concatenate(keys)) if keys else np.zeros(0, dtype=np.int64)


def edge_nodes(faces_a, faces_b, n_nodes: int):
    """2 タグの元の境界面接続 (sizes, nodes) から端集合 E (節点 ID の昇順 int64) と次元 (2/3) を返す。
    3D は共有辺の両端、2D は共有端点。座標は見ない (同じ ID のときだけ共有とみなす)。"""
    (sa, na), (sb, nb) = faces_a, faces_b
    if sa.size == 0 or sb.size == 0:
        raise EdgeMarkError("指定タグの元の境界面が空")
    da = 2 if int(sa.max()) == 2 else 3
    db = 2 if int(sb.max()) == 2 else 3
    if (da == 2 and int(sa.min()) != 2) or (db == 2 and int(sb.min()) != 2):
        raise EdgeMarkError("2 節点の線と多角形が同じタグに混在している")
    if da != db:
        raise EdgeMarkError(f"2 タグで次元が違う ({da}D と {db}D)")
    if da == 2:
        E = np.intersect1d(np.unique(na), np.unique(nb))
    else:
        ka = _edge_keys(_faces_by_size(sa, na), n_nodes)
        kb = _edge_keys(_faces_by_size(sb, nb), n_nodes)
        shared = np.intersect1d(ka, kb)
        E = np.unique(np.concatenate([shared // n_nodes, shared % n_nodes])) if shared.size else np.zeros(0, dtype=np.int64)
    return E.astype(np.int64), da


def adjacency(face_cells: np.ndarray, n_nodes: int):
    """内部面の隣接 CV 対から CSR の (indptr, indices) を作る (無向)。"""
    a = face_cells[:, 0]; b = face_cells[:, 1]
    if face_cells.size and (min(a.min(), b.min()) < 0 or max(a.max(), b.max()) >= n_nodes):
        raise EdgeMarkError("内部面の隣接 CV が節点数の範囲外 (ghost を含む面が混ざっている)")
    src = np.concatenate([a, b]); dst = np.concatenate([b, a])
    order = np.argsort(src, kind="stable")
    indices = dst[order]
    indptr = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(np.bincount(src, minlength=n_nodes), out=indptr[1:])
    return indptr, indices


def bfs_rings(indptr: np.ndarray, indices: np.ndarray, seeds: np.ndarray, rings: int) -> np.ndarray:
    """seeds から内部面グラフで距離 rings 以下の節点の距離 (それ以外は -1) を返す。"""
    n = indptr.size - 1
    dist = np.full(n, -1, dtype=np.int64)
    front = np.unique(np.asarray(seeds, dtype=np.int64))
    dist[front] = 0
    for r in range(1, rings + 1):
        if front.size == 0:
            break
        st = indptr[front]; ln = indptr[front + 1] - st
        tot = int(ln.sum())
        if tot == 0:
            break
        rep = np.repeat(st - np.concatenate([[0], np.cumsum(ln)[:-1]]), ln)
        nb = indices[rep + np.arange(tot)]
        nb = np.unique(nb)
        nb = nb[dist[nb] < 0]
        dist[nb] = r
        front = nb
    return dist


# ---------------------------------------------------------------------------
# タグ → physID
# ---------------------------------------------------------------------------
def resolve_tags(tag_args, tags_csv, bcond_config):
    """--tag NAME=ID (明示) か --tags A,B (または A B) + --bcond-config (run の bcond 設定) から [(name, physID), ...]。"""
    pairs = []
    for t in tag_args or []:
        if "=" not in t:
            raise EdgeMarkError(f"--tag は NAME=physID の形 ('{t}')")
        name, pid = t.split("=", 1)
        try:
            pairs.append((name.strip(), int(pid)))
        except ValueError:
            raise EdgeMarkError(f"--tag の physID が整数でない ('{t}')")
    if tags_csv:
        if pairs:
            raise EdgeMarkError("--tag と --tags は併用しない")
        if not bcond_config:
            raise EdgeMarkError("--tags には --bcond-config (run の bcondConfig.yaml) が要る")
        import yaml
        with open(bcond_config) as fh:
            bc = yaml.safe_load(fh) or {}
        names = [s.strip() for t in (tags_csv if isinstance(tags_csv, (list, tuple)) else [tags_csv])
                 for s in str(t).split(",") if s.strip()]
        for name in names:
            if name not in bc or not isinstance(bc[name], dict) or "physID" not in bc[name]:
                raise EdgeMarkError(f"タグ '{name}' が {bcond_config} に無い (境界条件名の誤記?)")
            pairs.append((name, int(bc[name]["physID"])))
    if len(pairs) != 2:
        raise EdgeMarkError(f"タグはちょうど 2 つ (板の上下の壁面) 指定する (指定 {len(pairs)})")
    if pairs[0][1] == pairs[1][1]:
        raise EdgeMarkError("2 タグの physID が同じ")
    return pairs


# ---------------------------------------------------------------------------
# 作成・書き込み・照合
# ---------------------------------------------------------------------------
def build_weight(h5path: str, pairs, rings: int):
    """(w, info) を返す (書かない)。"""
    if rings < 1:
        raise EdgeMarkError("rings は 1 以上の整数")
    with h5py.File(h5path, "r") as f:
        if discretization_of(f) != "node":
            raise EdgeMarkError(f"{h5path} は node 変換の h5 でない (/VIZMESH が無いか CV 数 ≠ 節点数)。"
                                " discretization: node で convertGmshToForge を回した h5 に使う")
        nN = int(f["MESH"].attrs["nNodes"])
        faces = []
        kinds = []
        for name, pid in pairs:
            if f"BCONDS/{pid}" not in f:
                raise EdgeMarkError(f"タグ '{name}' の physID {pid} が {h5path} の /BCONDS に無い")
            sz, nd = read_bface(f, pid)
            if sz is None:
                raise EdgeMarkError(f"タグ '{name}' (physID {pid}) に元の境界面の接続 (vizBfaceNodes) が無い (古い変換器の h5?)")
            if sz.size == 0:
                raise EdgeMarkError(f"タグ '{name}' (physID {pid}) の元の境界面が 0 枚")
            if nd.size and (nd.min() < 0 or nd.max() >= nN):
                raise EdgeMarkError(f"タグ '{name}' の境界面接続の節点 ID が範囲外")
            faces.append((sz, nd))
            kinds.append(str(f[f"BCONDS/{pid}"].attrs.get("bcondKind", "?")))
        E, dim = edge_nodes(faces[0], faces[1], nN)
        if E.size == 0:
            raise EdgeMarkError(f"タグ {pairs[0][0]}・{pairs[1][0]} が共有する{'辺' if dim == 3 else '端点'}が無い (E が空)。"
                                " 有限厚の板・別 ID の座標一致 (共有していない) なら処置の対象外")
        fc = read_internal_face_cells(f)
        indptr, indices = adjacency(fc, nN)
        dist = bfs_rings(indptr, indices, E, rings)
        S = np.nonzero(dist >= 0)[0]
        coord = np.asarray(f["MESH/COORD"][()], dtype=np.float64).reshape(-1, 3)
        vol = np.asarray(f["CELLS/volume"][()], dtype=np.float64) if "CELLS/volume" in f else None
        sig = mesh_signature(f)
    w = np.ones(nN, dtype=np.float32)
    w[S] = 0.0
    inS = np.zeros(nN, dtype=bool); inS[S] = True
    n_faces = int(np.count_nonzero(inS[fc[:, 0]] | inS[fc[:, 1]]))
    info = {"n_nodes": nN, "dim": dim, "E": E, "S": S, "dist": dist, "kinds": kinds,
            "E_bbox": (coord[E].min(0), coord[E].max(0)), "S_bbox": (coord[S].min(0), coord[S].max(0)),
            "S_volume": float(vol[S].sum()) if vol is not None else float("nan"), "n_faces": n_faces,
            "mesh_signature": sig, "ring_counts": [int(np.count_nonzero(dist == r)) for r in range(rings + 1)]}
    return w, info


def write_weight(h5path: str, w: np.ndarray, attrs: dict, xdmf: bool = True) -> str:
    """/AUX/w_recon_vel を (上書きで) 書き、ParaView 用 XDMF を出す。XDMF のパスを返す。"""
    with h5py.File(h5path, "a") as f:
        g = f.require_group("AUX")
        if FIELD in g:
            del g[FIELD]
        ds = g.create_dataset(FIELD, data=np.asarray(w, dtype="<f4"))
        for k, v in attrs.items():
            ds.attrs[k] = v
        nviz = int(f["VIZMESH"].attrs["nVizCells"]); vdim = int(f["VIZMESH/CONNE"].shape[0])
        nN = int(f["MESH"].attrs["nNodes"])
        prec = int(f["MESH/COORD"].dtype.itemsize)
    if not xdmf:
        return ""
    base = os.path.splitext(os.path.basename(h5path))[0]
    xp = os.path.join(os.path.dirname(os.path.abspath(h5path)), f"{base}_{FIELD}.xmf")
    h5n = os.path.basename(h5path)
    with open(xp, "w") as o:
        o.write("<?xml version='1.0' ?>\n<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n<Xdmf>\n  <Domain>\n")
        o.write(f"    <Grid Name='{FIELD}'>\n")
        o.write(f"      <Topology Type='Mixed' NumberOfElements='{nviz}'>\n")
        o.write(f"        <DataItem Format='HDF' DataType='Int' Dimensions='{vdim}'>\n          {h5n}:VIZMESH/CONNE\n        </DataItem>\n")
        o.write("      </Topology>\n      <Geometry Type='XYZ'>\n")
        o.write(f"        <DataItem Format='HDF' DataType='Float' Precision='{prec}' Dimensions='{3 * nN}'>\n          {h5n}:MESH/COORD\n        </DataItem>\n")
        o.write("      </Geometry>\n")
        o.write(f"      <Attribute Name='{FIELD}' Center='Node'>\n")
        o.write(f"        <DataItem Format='HDF' DataType='Float' Precision='4' Dimensions='{nN}'>\n          {h5n}:{AUX_PATH}\n        </DataItem>\n")
        o.write("      </Attribute>\n    </Grid>\n  </Domain>\n</Xdmf>\n")
    return xp


def verify(h5path: str) -> dict:
    """/AUX/w_recon_vel の値と属性を、ソルバと同じ規則で照合する。違反は EdgeMarkError。"""
    with h5py.File(h5path, "r") as f:
        if AUX_PATH not in f:
            raise EdgeMarkError(f"{h5path} に /{AUX_PATH} が無い")
        ds = f[AUX_PATH]
        if ds.dtype != np.dtype("<f4") or ds.ndim != 1:
            raise EdgeMarkError(f"/{AUX_PATH} は float32 の 1 次元配列でない ({ds.dtype}, {ds.shape})")
        w = ds[()]
        at = dict(ds.attrs)
        nC = int(f["MESH"].attrs["nCells"])
        if w.size != nC:
            raise EdgeMarkError(f"/{AUX_PATH} の長さ {w.size} が CV 数 {nC} と違う")
        if not np.all(np.isfinite(w)) or w.min() < 0.0 or w.max() > 1.0:
            raise EdgeMarkError(f"/{AUX_PATH} に非有限または [0,1] の外の値がある")
        sig = mesh_signature(f)
    def _s(v):
        return v.decode() if isinstance(v, bytes) else str(v)
    if _s(at.get("field_sha256", "")) != field_hash(w):
        raise EdgeMarkError("field_sha256 が値と一致しない (値が書き換えられている)")
    if _s(at.get("mesh_signature_version", "")) != SIG_VERSION:
        raise EdgeMarkError(f"mesh_signature_version が {SIG_VERSION} でない")
    if _s(at.get("mesh_signature", "")) != sig:
        raise EdgeMarkError("mesh_signature が格子と一致しない (別の格子・節点番号の違う格子の重み)")
    return {"mesh_signature": sig, "field_sha256": field_hash(w), "n_w_lt1": int(np.count_nonzero(w < 1.0))}


def _fmt3(v):
    return "(" + ", ".join(f"{x:.6g}" for x in v) + ")"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("h5", help="node 変換済みの入力 h5 (convertGmshToForge の出力)。/AUX/w_recon_vel を書く")
    ap.add_argument("--tag", action="append", help="NAME=physID (2 回。板の上下の壁面、例 --tag cowl_in=5 --tag cowl_out=6)")
    ap.add_argument("--tags", nargs="+", help="NAME,NAME または NAME NAME (physID は --bcond-config から引く)")
    ap.add_argument("--bcond-config", help="run の bcondConfig.yaml (--tags のとき。既定: h5 と同じ場所の bcondConfig.yaml)")
    ap.add_argument("--field", default=FIELD, help=f"書くデータセット名 (初版は {FIELD} 固定。ソルバはこの名前だけを読む)")
    ap.add_argument("--rings", type=int, default=2, help="E から内部面グラフで何リングまでを w=0 にするか (既定 2)")
    ap.add_argument("--control-ones", action="store_true", help="対照用: タグなしで w ≡ 1 を書く (処置は恒等)")
    ap.add_argument("--verify", action="store_true", help="書かずに既存の /AUX/w_recon_vel を照合する")
    ap.add_argument("--dry-run", action="store_true", help="作って要約を出すだけで書かない")
    ap.add_argument("--no-xdmf", action="store_true", help="XDMF を出さない")
    a = ap.parse_args(argv)
    try:
        if a.field != FIELD:
            raise EdgeMarkError(f"--field は {FIELD} だけ (ソルバは /AUX/{FIELD} だけを読む)")
        if a.tags and not a.bcond_config:
            cand = os.path.join(os.path.dirname(os.path.abspath(a.h5)), "bcondConfig.yaml")
            if not os.path.exists(cand):
                raise EdgeMarkError(f"--tags には --bcond-config が要る ({cand} も無い)")
            a.bcond_config = cand
        if a.verify:
            r = verify(a.h5)
            print(f"[mark_zero_thickness_edges] OK: {a.h5} /{AUX_PATH} mesh_signature {r['mesh_signature'][:16]} "
                  f"field_sha256 {r['field_sha256'][:16]} w<1 {r['n_w_lt1']}")
            return 0
        now = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if a.control_ones:
            if a.tag or a.tags:
                raise EdgeMarkError("--control-ones とタグは併用しない")
            with h5py.File(a.h5, "r") as f:
                if discretization_of(f) != "node":
                    raise EdgeMarkError(f"{a.h5} は node 変換の h5 でない")
                nN = int(f["MESH"].attrs["nNodes"])
                sig = mesh_signature(f)
            w = np.ones(nN, dtype=np.float32)
            attrs = {"generator": GENERATOR, "generator_version": GENERATOR_VERSION, "mode": "control_ones",
                     "tags": "", "tag_physids": "", "rings": 0, "n_edge_nodes": 0, "n_marked_nodes": 0, "n_nodes": nN,
                     "mesh_signature": sig, "mesh_signature_version": SIG_VERSION, "field_sha256": field_hash(w), "created": now}
            print(f"[mark_zero_thickness_edges] 対照 (w ≡ 1): 節点 {nN}、mesh_signature {sig[:16]}")
        else:
            pairs = resolve_tags(a.tag, a.tags, a.bcond_config)
            w, info = build_weight(a.h5, pairs, a.rings)
            attrs = {"generator": GENERATOR, "generator_version": GENERATOR_VERSION, "mode": "edges",
                     "tags": ",".join(n for n, _ in pairs), "tag_physids": ",".join(f"{n}:{p}" for n, p in pairs),
                     "rings": int(a.rings), "n_edge_nodes": int(info["E"].size), "n_marked_nodes": int(info["S"].size),
                     "n_nodes": int(info["n_nodes"]), "edge_dim": int(info["dim"]),
                     "mesh_signature": info["mesh_signature"], "mesh_signature_version": SIG_VERSION,
                     "field_sha256": field_hash(w), "created": now}
            print(f"[mark_zero_thickness_edges] {info['dim']}D、タグ {attrs['tag_physids']} (kind {','.join(info['kinds'])})")
            print(f"  E (共有{'辺の節点' if info['dim'] == 3 else '端点'}) {info['E'].size} 節点: 範囲 {_fmt3(info['E_bbox'][0])} – {_fmt3(info['E_bbox'][1])}")
            print(f"  S (rings {a.rings}) {info['S'].size} 節点 (距離別 {info['ring_counts']}): 範囲 {_fmt3(info['S_bbox'][0])} – {_fmt3(info['S_bbox'][1])}")
            print(f"  対象の内部面 {info['n_faces']}、S の体積の合計 {info['S_volume']:.6g}、mesh_signature {info['mesh_signature'][:16]}")
            if info["E"].size <= 20:
                print(f"  E の節点 ID: {info['E'].tolist()}")
        if a.dry_run:
            print("[mark_zero_thickness_edges] --dry-run: 書かない")
            return 0
        xp = write_weight(a.h5, w, attrs, xdmf=not a.no_xdmf)
        verify(a.h5)
        print(f"[mark_zero_thickness_edges] 書いた: {a.h5}:/{AUX_PATH} (field_sha256 {attrs['field_sha256'][:16]})"
              + (f"、ParaView: {xp}" if xp else ""))
        return 0
    except EdgeMarkError as e:
        print(f"[mark_zero_thickness_edges] ERROR: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
