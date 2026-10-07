#!/usr/bin/env python3
"""壁が終わる所 (壁終端の線) の格子の幾何を、流れを回さずに点検する (読むだけ、合否の判定はしない)。

    python3 solver_density_cuda/tools/check_wall_end_geometry.py MESH.h5 --bcond-config bcondConfig.yaml \\
        [--pair cowl_in,cowl_out] [--pair sidewall_in,sidewall_out] [--layers 3] [--top 20] [--sort kink]
    python3 solver_density_cuda/tools/check_wall_end_geometry.py MESH.h5 --wall cowl_in=5 --wall cowl_out=6 --wall ramp=4

plan: plans/active/convection-zero-thickness-edge-reconstruction.md §4.2「他の壁終端の幾何点検」・§5.1 #2d。
SERN g3 のカウル後縁で、壁に沿った格子線が後縁の 1 セルで約 36° 折れていた (TE_GRID_GEOM_M10.txt)。
同じ型の折れが他の「壁が終わる所」にもないかを、幾何と接続だけで全 span について集計する。

入力 (node 変換済みの入力 h5 = convertGmshToForge の出力、discretization: node):
  /MESH/COORD                      節点座標 (節点 = node の CV)
  /PLANES/STRUCT の先頭 nNormalPlanes 面   内部面 (双対面) の隣接節点対 = 格子の辺 = forge の LSQ 近傍 (内部隣接のみ)
  /VIZMESH/CONNE                   primal セル (3D は hex、2D は quad。gmsh の節点順)。端の周りの流体角と格子線のたどりに使う
  /BCONDS/<physID>/vizBfaceSizes・vizBfaceNodes   元の境界面 (3D は面、2D は線) の接続。属性 bcondKind
  node 変換後の境界面 (PLANES の半割面) は 1 節点なので使わない (gmshReader.hpp の replacePrimalWithDual が退避した
  vizBfaceNodes を使う)。
壁: --bcond-config の kind が wall / wall_isothermal (ソルバの isWallKind と同じ) のタグ。--walls A,B,... で絞る。
  --wall NAME=physID (複数回) で明示したときは、そのタグだけが壁。どちらも無ければ h5 の bcondKind で選び、名前は "pid<N>"。

定義:
  端の要素   3D は壁面の辺、2D は壁の線の端点。各要素について
             F = その要素を含む元の境界面 (全タグ)、
             φ = 要素の周りの流体角 = 要素を含む primal セルの二面角 (2D は内角) の和 (内部の辺なら 360°、平らな壁なら 180°)。
  分類       knife  (φ ≥ --knife-deg、既定 300°): 厚さ 0 の板の自由端 (カウル後縁・露出したカウル側端・側壁後縁)
             convex (--convex-deg < φ < knife、既定 225°): 凸の角 (ランプ後縁 ∩ 機体ベース、カウル下面 ∩ 側壁外面、機体の角)。
                    45° 未満の壁の折れ (ランプの膨張の角など、φ 180〜225°) は壁の曲がりであって終端ではないので入れない
             open   (--open-deg ≤ φ ≤ convex かつ F に壁でない面がある、既定 160°): 壁が壁でない境界 (出口・遠方) へほぼ同じ面で続く所
             それ以外 (凹の角、対称面・入口・出口との交線、壁どうしの平らなつなぎ目) は終端としない。
             **タグだけで「隣が壁でない辺」とはしない**: 側壁後縁は sidewall_in/out の両面が壁、ランプ後縁は vehicle_base (壁) に
             接する。逆に壁と対称面・入口の交線は隣が壁でないが終端ではない。流体角なら両方を正しく分けられる。
  --pair A,B A と B の面が**同じ節点 ID で**共有する辺 (2D は端点) の数と、そのうち knife の数を出す (plan §4 の端集合 E と
             同じ定義。座標一致の別 ID は共有とみなさない)。0 本ならエラー (タグの誤記や有限厚と区別するため)。
  終端線     同じ分類・同じタグの組の端の辺を節点でつないだ連結成分。次数 2 の節点で辺の向きが --turn-deg (既定 45°)
             以上変わる所 (後縁の線と側端の線の角) で分ける (角の節点は両方の線に入る)。2D は端点ごとに 1 本。
  面 (side)  端の要素を含む壁の境界面ごと (knife は上下 2 面、convex は 2 つの壁)。上側・下側はタグで分けて出す
             (同じタグが両面なら壁法線の主成分の符号を付ける)。
  格子線     hex / quad の構造格子を仮定する。端の辺 (t, s) を含む各 hex について t から出る残り 2 本の辺 (2D は quad の t の
             2 辺) を「扇」として環状につなぐ。その面の壁の辺 t→u (u = 境界面の上の t の隣、s でない方) から、
             その境界面を持つセルに入って 1 セル進んだ辺の先 = 第 1 層の点 n1 (壁から法線方向へ 1 接続)、
             2 セル進んだ辺の先 = 下流の点 d (壁の格子線 u→t を端の先へ延ばした格子線)。最初のセルを面の持ち主で
             決めるのは、後縁の線と露出した側端の線の角の節点で上下の面が 3 節点を共有して扇が閉じるため。
             第 k 層 (k = 1..--layers): p_0 = t、p_1 = n1、p_{k+1} = 格子線 p_{k-1}→p_k の直進 (p_k の隣のうち、p_{k-1} との
             共通の隣が p_k だけの点)、u_k・d_k = (u_{k-1}, p_k)・(d_{k-1}, p_k) の共通の隣のうち p_{k-1} でない点。
量 (節点 p_k と面ごと。座標は h5 の単位):
  (i)   kink    = ∠(p_k − u_k, d_k − p_k) [度]。前後の辺の角度差 (0 = 直線)
  (ii)  disp/h  = (d_k − p_k)·n̂ / ((p_{k+1} − p_k)·n̂)。n̂ = その面の壁法線 (境界面の Newell 法線を第 1 層の側へ向けたもの、
                  全層で同じ n̂)。分子 = 下流の辺の壁法線方向の変位、分母 = その節点の壁法線方向の層間隔 (k = 0 で第一層厚)。
                  正 = その面の流体の側へ、負 = 壁の延長面を越えて反対側へ
  (iii) len     = |d_k − p_k| / |p_k − u_k| (下流/上流の辺長)
  (iv)  LSQ     内部隣接の数、M = Σ d̂ d̂ᵀ (d = x_j − x_i。forge の gradLSQ 2 の w = 1/|d|² と同じで辺長は消える) の固有値
                (昇順。2D は xy の 2 つ)、条件数、forge のスペクトル打ち切り (λ < 1e-2 λmax) で落ちる数
  φ (その端の辺の流体角) と扇のセル数も行ごとに出す。
出力: 終端線の一覧、終端線 × 面 × 層の要約 (節点数、kink の最大・中央値、|disp/h| の最大、len の最大・最小、条件数の最大、
  kink 最大の節点と座標)、--sort の上位 N 節点の詳細。CSV は <prefix>_lines.csv (要約) と <prefix>_nodes.csv (全行)。
flag (その行の量が格子線の延長として信用できない理由。空なら問題なし):
  fan<n>          扇のセル数が標準 (knife 4・convex 3・open 2) でない (非構造・特異点)
  face_cell<n>    壁の境界面を持つセルが扇の中に 1 つでない
  u_not_in_face_cell / u_in_<n>_cells / fan_end / fan_branch   扇をたどれない (n1・d が無い)
  F<n>            端の要素を含む元の境界面が 2 枚でない (非多様体・開いた境界)
  not_edge_u / not_edge_d   u・d との辺が内部面 (LSQ の近傍) に無い
  no_p            その層の節点が無い (前の層で止まった)
  amb_<u|d|p>     u_k・d_k・p_{k+1} が一意に決まらない・無い (外部境界に届いた層など)。p だけなら disp/h が出ない
  h<=0            層間隔の n̂ 成分が 0 以下 (層の格子線が壁に沿って寝ている)

限界:
  - hex (3D) / quad (2D) だけの格子 (他のセル型があればエラー)。「2 セル進む」は構造格子の格子線の延長で、
    扇のセル数が標準 (knife 4・convex 3・open 2) でないときは下流の点が格子線の延長とは限らない → flag "fan<n>"。
  - φ は各 hex の二面角を端の両端で、そこから出る 2 本の hex の辺の射影で測った平均 (非平面の面では近似)。
    分類の閾値は大きく離れた値 (90/180/270/360°) を分けるためのもので、細かい角度の判定には使わない。
  - 層間隔は p_{k+1} − p_k の n̂ 成分。層の格子線が壁に垂直でなければ実際の層厚と違う。n̂ は端の節点の境界面の法線で、
    層ごとに測り直さない。
  - 合否の判定はしない (閾値と比較の相手は plan 側)。流れの量は読まない。
  - 端の要素の周りだけを解析するが、PLANES/STRUCT と VIZMESH/CONNE は全体を読む。100 万節点級の 3D は数 GB 要るので
    AWS で回す。
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import h5py
import numpy as np

TOOL = "check_wall_end_geometry"
WALL_KINDS = ("wall", "wall_isothermal")          # ソルバの isWallKind (mesh.cpp) と同じ
LSQ_DROP = 1.0e-2                                 # forge の node LSQ のスペクトル打ち切り (te_grad_diag.py と同じ)
CELL_NNODES = {9: 8, 8: 6, 7: 5, 6: 4, 5: 4, 4: 3}
CELL_NAME = {9: "hex", 8: "prism", 7: "pyramid", 6: "tetra", 5: "quad", 4: "triangle"}
# gmsh hex (底面 0-3・上面 4-7) の辺と、各局所節点の 3 つの隣
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
HEX_NBR = ((1, 3, 4), (0, 2, 5), (1, 3, 6), (2, 0, 7), (5, 7, 0), (4, 6, 1), (5, 7, 2), (6, 4, 3))
EXPECTED_FAN = {"knife": 4, "convex": 3, "open": 2}
CLASS_ORDER = {"knife": 0, "convex": 1, "open": 2}
SORT_KEYS = ("kink", "disp", "len", "cond")


class WallEndError(RuntimeError):
    """入力の不備 (タグ無し・接続なし・未対応のセル型など)。CLI はこれで終了コード 2。"""


# ---------------------------------------------------------------------------
# 読み込み
# ---------------------------------------------------------------------------
def _s(v) -> str:
    return v.decode() if isinstance(v, bytes) else str(v)


def read_internal_pairs(f: h5py.File) -> np.ndarray:
    """内部面 ip < nNormalPlanes の隣接 CV 対 (int64 (nI, 2))。/PLANES/STRUCT は面ごとに [nNodes, nodes..., nCells, cells...]。"""
    nI = int(f["MESH"].attrs["nNormalPlanes"])
    s = np.asarray(f["PLANES/STRUCT"][()], dtype=np.int64)
    # node 変換の内部面は [2, a, b, 2, c0, c1] の 6 語 (先頭から帰納的に成り立つときだけ使う)
    if s.size >= 6 * nI:
        s6 = s[:6 * nI].reshape(nI, 6)
        if nI == 0 or (np.all(s6[:, 0] == 2) and np.all(s6[:, 3] == 2)):
            return np.ascontiguousarray(s6[:, 4:6])
    out = np.empty((nI, 2), dtype=np.int64)
    sl = s.tolist()
    p = 0
    for ip in range(nI):
        nn = sl[p]; p += 1 + nn
        nc = sl[p]; p += 1
        if nc != 2:
            raise WallEndError(f"内部面 {ip} の隣接 CV が {nc} 個 (2 であるべき)")
        out[ip, 0] = sl[p]; out[ip, 1] = sl[p + 1]
        p += nc
    return out


def parse_conne(conne: np.ndarray) -> dict:
    """/VIZMESH/CONNE ([code, nodes..., code, nodes...]) を {code: (n, nnodes) int64}。"""
    conne = np.asarray(conne, dtype=np.int64).reshape(-1)
    for code in (9, 5):
        w = 1 + CELL_NNODES[code]
        if conne.size and conne.size % w == 0 and np.all(conne[0::w] == code):
            return {code: np.ascontiguousarray(conne.reshape(-1, w)[:, 1:])}
    starts = defaultdict(list)
    sl = conne.tolist()
    p = 0
    while p < len(sl):
        c = sl[p]
        if c not in CELL_NNODES:
            raise WallEndError(f"/VIZMESH/CONNE に未知のセル型コード {c} (位置 {p})")
        starts[c].append(p + 1)
        p += 1 + CELL_NNODES[c]
    return {c: conne[np.asarray(st)[:, None] + np.arange(CELL_NNODES[c])[None, :]] for c, st in starts.items()}


class Mesh:
    """解析に使う配列をまとめたもの。"""

    def __init__(self, path: str):
        self.path = path
        with h5py.File(path, "r") as f:
            a = f["MESH"].attrs
            self.n = int(a["nNodes"])
            if "VIZMESH" not in f or int(a["nCells"]) != self.n:
                raise WallEndError(f"{path} は node 変換の h5 でない (/VIZMESH が無いか CV 数 ≠ 節点数)。"
                                   " discretization: node で convertGmshToForge を回した h5 に使う")
            X = np.asarray(f["MESH/COORD"][()], dtype=np.float64).reshape(-1)
            if X.size != 3 * self.n:
                raise WallEndError(f"/MESH/COORD の長さ {X.size} が 3 x nNodes ({3 * self.n}) でない")
            self.X = X.reshape(-1, 3)
            self.pairs = read_internal_pairs(f)
            self.cells = parse_conne(f["VIZMESH/CONNE"][()])
            self.bc = {}
            for k in f["BCONDS"].keys():
                g = f["BCONDS/" + k]
                kind = _s(g.attrs.get("bcondKind", "?"))
                if "vizBfaceSizes" in g and "vizBfaceNodes" in g:
                    sz = np.asarray(g["vizBfaceSizes"][()], dtype=np.int64).reshape(-1)
                    nd = np.asarray(g["vizBfaceNodes"][()], dtype=np.int64).reshape(-1)
                else:
                    sz = nd = None
                self.bc[int(k)] = {"kind": kind, "sizes": sz, "nodes": nd}
        if self.pairs.size and (self.pairs.min() < 0 or self.pairs.max() >= self.n):
            raise WallEndError("内部面の隣接 CV が節点数の範囲外 (ghost を含む面が混ざっている)")
        # 次元: 元の境界面の最大節点数 (2 = 線 → 2D)
        mx = max((int(b["sizes"].max()) for b in self.bc.values() if b["sizes"] is not None and b["sizes"].size), default=0)
        if mx == 0:
            raise WallEndError("どの境界にも元の境界面の接続 (vizBfaceNodes) が無い (古い変換器の h5?)")
        self.dim = 2 if mx == 2 else 3
        want = 9 if self.dim == 3 else 5
        other = {CELL_NAME[c]: int(v.shape[0]) for c, v in self.cells.items() if c != want}
        if other or want not in self.cells:
            raise WallEndError(f"{self.dim}D は {CELL_NAME[want]} だけの格子を前提にしている (他のセル型 {other})")
        self.C = self.cells[want]
        # 内部隣接 (CSR) と辺の存在検査用のキー
        a_, b_ = self.pairs[:, 0], self.pairs[:, 1]
        src = np.concatenate([a_, b_]); dst = np.concatenate([b_, a_])
        order = np.argsort(src, kind="stable")
        self.indices = dst[order]
        self.indptr = np.zeros(self.n + 1, dtype=np.int64)
        np.cumsum(np.bincount(src, minlength=self.n), out=self.indptr[1:])
        self.ekeys = np.unique(np.minimum(a_, b_) * np.int64(self.n) + np.maximum(a_, b_))
        self._nb = {}

    def nbrs(self, i: int) -> frozenset:
        s = self._nb.get(i)
        if s is None:
            s = frozenset(self.indices[self.indptr[i]:self.indptr[i + 1]].tolist())
            self._nb[i] = s
        return s

    def has_edge(self, i: int, j: int) -> bool:
        k = min(i, j) * self.n + max(i, j)
        p = int(np.searchsorted(self.ekeys, k))
        return p < self.ekeys.size and int(self.ekeys[p]) == k


# ---------------------------------------------------------------------------
# タグ → physID と壁の選択
# ---------------------------------------------------------------------------
def resolve_tags(mesh: Mesh, bcond_config, wall_args, walls_csv):
    """(names {pid: 名前}, walls [pid], 注記 [str]) を返す。"""
    notes = []
    names = {pid: f"pid{pid}" for pid in mesh.bc}
    kinds = {pid: b["kind"] for pid, b in mesh.bc.items()}
    cfg_pids = {}
    if bcond_config:
        import yaml
        with open(bcond_config) as fh:
            bc = yaml.safe_load(fh) or {}
        for name, v in bc.items():
            if isinstance(v, dict) and "physID" in v:
                pid = int(v["physID"])
                cfg_pids[name] = pid
                if pid in names:
                    names[pid] = name
                    k = str(v.get("kind", ""))
                    if k and k != kinds[pid]:
                        notes.append(f"注意: {name} (physID {pid}) の kind が bcondConfig '{k}' と h5 '{kinds[pid]}' で違う (bcondConfig を使う)")
                    if k:
                        kinds[pid] = k
    if wall_args:
        walls = []
        for t in wall_args:
            if "=" not in t:
                raise WallEndError(f"--wall は NAME=physID の形 ('{t}')")
            nm, p = t.split("=", 1)
            try:
                pid = int(p)
            except ValueError:
                raise WallEndError(f"--wall の physID が整数でない ('{t}')")
            if pid not in mesh.bc:
                raise WallEndError(f"--wall {nm}={pid}: physID {pid} が h5 の /BCONDS に無い")
            names[pid] = nm.strip()
            walls.append(pid)
    else:
        walls = [pid for pid in mesh.bc if kinds[pid] in WALL_KINDS]
        if not bcond_config:
            notes.append("注意: --bcond-config も --wall も無いので h5 の bcondKind で壁を選び、名前は pid<N> で出す")
    if walls_csv:
        want = [s.strip() for s in walls_csv.split(",") if s.strip()]
        by_name = {v: k for k, v in names.items()}
        bad = [w for w in want if w not in by_name]
        if bad:
            raise WallEndError(f"--walls のタグ {bad} が無い (境界条件名の誤記? 既知: {sorted(by_name)})")
        walls = [by_name[w] for w in want]
    for pid in walls:
        if mesh.bc[pid]["sizes"] is None or mesh.bc[pid]["sizes"].size == 0:
            raise WallEndError(f"壁 {names[pid]} (physID {pid}) に元の境界面の接続 (vizBfaceNodes) が無い・空")
    for pid, b in mesh.bc.items():
        if b["sizes"] is None and pid not in walls:
            notes.append(f"注意: {names[pid]} (physID {pid}) に元の境界面の接続が無い → その面に接する辺の F が欠ける")
    if not walls:
        raise WallEndError("壁のタグが 1 つも無い")
    return names, sorted(set(walls)), notes, cfg_pids


# ---------------------------------------------------------------------------
# 境界面の表
# ---------------------------------------------------------------------------
class BFaces:
    """全タグの元の境界面を通し番号で持つ。"""

    def __init__(self, mesh: Mesh):
        pids, sizes, flat = [], [], []
        for pid in sorted(mesh.bc):
            b = mesh.bc[pid]
            if b["sizes"] is None or b["sizes"].size == 0:
                continue
            if int(b["sizes"].sum()) != b["nodes"].size:
                raise WallEndError(f"physID {pid}: vizBfaceSizes の総和と vizBfaceNodes の長さが違う")
            if b["nodes"].size and (b["nodes"].min() < 0 or b["nodes"].max() >= mesh.n):
                raise WallEndError(f"physID {pid}: 境界面接続の節点 ID が範囲外")
            pids.append(np.full(b["sizes"].size, pid, dtype=np.int64)); sizes.append(b["sizes"]); flat.append(b["nodes"])
        self.pid = np.concatenate(pids)
        self.size = np.concatenate(sizes)
        self.nodes = np.concatenate(flat)
        self.start = np.concatenate([[0], np.cumsum(self.size)[:-1]]).astype(np.int64)

    def face(self, fi: int) -> np.ndarray:
        return self.nodes[self.start[fi]:self.start[fi] + self.size[fi]]

    def edges(self, n: int, sel=None):
        """多角形の辺 (閉じる) のキー min*n+max と面の番号。sel は面の bool マスク。"""
        keys, fids = [], []
        for s in np.unique(self.size):
            m = self.size == s
            if sel is not None:
                m &= sel
            idx = np.nonzero(m)[0]
            if idx.size == 0:
                continue
            F = self.nodes[self.start[idx][:, None] + np.arange(int(s))[None, :]]
            for j in range(int(s)):
                a = F[:, j]; b = F[:, (j + 1) % int(s)]
                keys.append(np.minimum(a, b) * np.int64(n) + np.maximum(a, b)); fids.append(idx)
        if not keys:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        return np.concatenate(keys), np.concatenate(fids)

    def vertices(self, sel=None):
        """2D の線の端点と面の番号。"""
        m = self.size == 2
        if sel is not None:
            m &= sel
        idx = np.nonzero(m)[0]
        F = self.nodes[self.start[idx][:, None] + np.arange(2)[None, :]]
        return np.concatenate([F[:, 0], F[:, 1]]), np.concatenate([idx, idx])


def _group(keys_sorted: np.ndarray, vals_sorted: np.ndarray, cand: np.ndarray) -> list:
    """cand (昇順) の各キーに対応する vals のリスト。keys_sorted は昇順。"""
    lo = np.searchsorted(keys_sorted, cand, side="left")
    hi = np.searchsorted(keys_sorted, cand, side="right")
    return [vals_sorted[a:b].tolist() for a, b in zip(lo, hi)]


# ---------------------------------------------------------------------------
# 流体角 φ
# ---------------------------------------------------------------------------
def _proj_angle(e: np.ndarray, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """e に垂直な面へ v1・v2 を射影した間の角 [rad] (行ごと)。"""
    eh = e / np.linalg.norm(e, axis=1, keepdims=True)
    p1 = v1 - np.sum(v1 * eh, axis=1, keepdims=True) * eh
    p2 = v2 - np.sum(v2 * eh, axis=1, keepdims=True) * eh
    return np.arctan2(np.linalg.norm(np.cross(p1, p2), axis=1), np.sum(p1 * p2, axis=1))


def _angle(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    return np.arctan2(np.linalg.norm(np.cross(v1, v2), axis=-1), np.sum(v1 * v2, axis=-1))


def fluid_angle_3d(mesh: Mesh, Hs: np.ndarray, cand: np.ndarray) -> np.ndarray:
    """候補の辺 (キー、昇順) ごとの流体角 [deg] = その辺を含む hex の二面角の和 (両端の平均)。"""
    X, n = mesh.X, mesh.n
    phi = np.zeros(cand.size)
    if cand.size == 0:
        return phi
    for a, b in HEX_EDGES:
        A = Hs[:, a]; B = Hs[:, b]
        key = np.minimum(A, B) * np.int64(n) + np.maximum(A, B)
        pos_c = np.minimum(np.searchsorted(cand, key), cand.size - 1)
        m = cand[pos_c] == key
        if not np.any(m):
            continue
        H = Hs[m]
        ca = [c for c in HEX_NBR[a] if c != b]
        cb = [c for c in HEX_NBR[b] if c != a]
        e = X[H[:, b]] - X[H[:, a]]
        ang_a = _proj_angle(e, X[H[:, ca[0]]] - X[H[:, a]], X[H[:, ca[1]]] - X[H[:, a]])
        ang_b = _proj_angle(e, X[H[:, cb[0]]] - X[H[:, b]], X[H[:, cb[1]]] - X[H[:, b]])
        np.add.at(phi, pos_c[m], 0.5 * (ang_a + ang_b))
    return np.degrees(phi)


def fluid_angle_2d(mesh: Mesh, Qs: np.ndarray, cand: np.ndarray) -> np.ndarray:
    """候補の節点 (昇順) ごとの流体角 [deg] = その節点を含む quad の内角の和。"""
    X = mesh.X
    phi = np.zeros(cand.size)
    if cand.size == 0:
        return phi
    for i in range(4):
        v = Qs[:, i]
        pos_c = np.minimum(np.searchsorted(cand, v), cand.size - 1)
        m = cand[pos_c] == v
        if not np.any(m):
            continue
        Q = Qs[m]
        ang = _angle(X[Q[:, (i - 1) % 4]] - X[Q[:, i]], X[Q[:, (i + 1) % 4]] - X[Q[:, i]])
        np.add.at(phi, pos_c[m], ang)
    return np.degrees(phi)


def classify(phi: float, has_nonwall: bool, th) -> str | None:
    if phi >= th["knife"]:
        return "knife"
    if phi > th["convex"]:
        return "convex"
    if phi >= th["open"] and has_nonwall:
        return "open"
    return None


# ---------------------------------------------------------------------------
# 端の要素の抽出と終端線
# ---------------------------------------------------------------------------
def find_ends(mesh: Mesh, bf: BFaces, walls: list, th: dict):
    """端の要素の一覧を返す。各要素 = dict(nodes=(t,s) or (t,), cls, phi, faces [面の番号], pids)。
    あわせて全候補の数・分類別の数・壁の hex/quad 部分集合を返す。"""
    n = mesh.n
    wallset = set(walls)
    is_wall_face = np.isin(bf.pid, np.asarray(walls, dtype=np.int64))
    wall_node = np.zeros(n, dtype=bool)
    wall_node[bf.nodes[np.repeat(is_wall_face, bf.size)]] = True
    Cs = mesh.C[wall_node[mesh.C].any(axis=1)]           # 壁の節点を含むセルだけ
    if mesh.dim == 3:
        k_all, f_all = bf.edges(n)
        cand = np.unique(bf.edges(n, sel=is_wall_face)[0])
        phi = fluid_angle_3d(mesh, Cs, cand)
    else:
        k_all, f_all = bf.vertices()
        cand = np.unique(bf.vertices(sel=is_wall_face)[0])
        phi = fluid_angle_2d(mesh, Cs, cand)
    sel = np.isin(k_all, cand)
    order = np.argsort(k_all[sel], kind="stable")
    faces_of = _group(k_all[sel][order], f_all[sel][order], cand)
    ends = []
    for ci, key in enumerate(cand.tolist()):
        fs = faces_of[ci]
        pids = [int(bf.pid[x]) for x in fs]
        has_nonwall = any(p not in wallset for p in pids)
        cls = classify(float(phi[ci]), has_nonwall, th)
        if cls is None:
            continue
        nodes = (key // n, key % n) if mesh.dim == 3 else (key,)
        ends.append({"nodes": nodes, "cls": cls, "phi": float(phi[ci]), "faces": fs, "pids": pids})
    return ends, int(cand.size), Cs


def shared_count(mesh: Mesh, bf: BFaces, pa: int, pb: int):
    """2 タグの面が同じ節点 ID で共有する辺 (2D は端点) のキーの配列。"""
    if mesh.dim == 3:
        ka = np.unique(bf.edges(mesh.n, sel=bf.pid == pa)[0])
        kb = np.unique(bf.edges(mesh.n, sel=bf.pid == pb)[0])
    else:
        ka = np.unique(bf.vertices(sel=bf.pid == pa)[0])
        kb = np.unique(bf.vertices(sel=bf.pid == pb)[0])
    return np.intersect1d(ka, kb)


def build_lines(mesh: Mesh, ends: list, names: dict, turn_deg: float) -> list:
    """端の要素を終端線にまとめる。"""
    X = mesh.X
    groups = defaultdict(list)
    for e in ends:
        tags = tuple(sorted({names[p] for p in e["pids"]}))
        groups[(e["cls"], tags)].append(e)
    lines = []
    for (cls, tags), es in groups.items():
        if mesh.dim == 2:
            for e in es:
                lines.append({"cls": cls, "tags": tags, "elems": [e], "nodes": [e["nodes"][0]]})
            continue
        parent = list(range(len(es)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        at = defaultdict(list)
        for ei, e in enumerate(es):
            for v in e["nodes"]:
                at[v].append(ei)
        for v, eis in at.items():
            if len(eis) != 2:
                continue
            e1, e2 = es[eis[0]]["nodes"], es[eis[1]]["nodes"]
            o1 = e1[0] if e1[1] == v else e1[1]
            o2 = e2[0] if e2[1] == v else e2[1]
            ang = math.degrees(float(_angle(X[v] - X[o1], X[o2] - X[v])))
            if ang < turn_deg:
                parent[find(eis[0])] = find(eis[1])
        comp = defaultdict(list)
        for ei in range(len(es)):
            comp[find(ei)].append(es[ei])
        for elems in comp.values():
            nodes = sorted({v for e in elems for v in e["nodes"]})
            lines.append({"cls": cls, "tags": tags, "elems": elems, "nodes": nodes})
    for L in lines:
        P = X[L["nodes"]]
        c = P.mean(axis=0)
        if len(L["nodes"]) > 1:
            _, _, vt = np.linalg.svd(P - c, full_matrices=False)
            ax = vt[0]
            ax = ax * (1.0 if ax[int(np.argmax(np.abs(ax)))] > 0 else -1.0)
            L["nodes"] = [L["nodes"][i] for i in np.argsort((P - c) @ ax, kind="stable")]
        else:
            ax = np.zeros(3)
        L["dir"] = ax
        L["bbox"] = (P.min(axis=0), P.max(axis=0))
        L["length"] = float(sum(np.linalg.norm(X[e["nodes"][1]] - X[e["nodes"][0]]) for e in L["elems"])) if mesh.dim == 3 else 0.0
        L["phi"] = (min(e["phi"] for e in L["elems"]), max(e["phi"] for e in L["elems"]))
    lines.sort(key=lambda L: (CLASS_ORDER[L["cls"]], L["tags"], float(L["bbox"][0][0]), float(L["bbox"][0][1]), float(L["bbox"][0][2])))
    for i, L in enumerate(lines, start=1):
        L["id"] = f"L{i}"
    return lines


# ---------------------------------------------------------------------------
# 節点ごとの解析
# ---------------------------------------------------------------------------
class CellIndex:
    """節点 → 壁の近くのセル (Cs の行) の CSR。"""

    def __init__(self, Cs: np.ndarray, n: int):
        self.C = Cs
        src = Cs.reshape(-1)
        cid = np.repeat(np.arange(Cs.shape[0], dtype=np.int64), Cs.shape[1])
        order = np.argsort(src, kind="stable")
        self.cells = cid[order]
        self.ptr = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(np.bincount(src, minlength=n), out=self.ptr[1:])

    def of(self, v: int) -> np.ndarray:
        return self.cells[self.ptr[v]:self.ptr[v + 1]]


def fan(mesh: Mesh, ci: CellIndex, elem: tuple):
    """端の要素の周りのセルごとの (ray1, ray2, 角[rad], セルの行番号)。3D は辺 (t, s) を含む hex、2D は t を含む quad。"""
    X = mesh.X
    t = elem[0]
    out = []
    for c in ci.of(t).tolist():
        row = ci.C[c]
        lt = int(np.nonzero(row == t)[0][0])
        if mesh.dim == 3:
            s = elem[1]
            nb = [row[k] for k in HEX_NBR[lt]]
            if s not in nb:
                continue
            r = [int(v) for v in nb if v != s]
            e = (X[s] - X[t])[None, :]
            ang = float(_proj_angle(e, (X[r[0]] - X[t])[None, :], (X[r[1]] - X[t])[None, :])[0])
        else:
            r = [int(row[(lt - 1) % 4]), int(row[(lt + 1) % 4])]
            ang = float(_angle(X[r[0]] - X[t], X[r[1]] - X[t]))
        out.append((r[0], r[1], ang, int(c)))
    return out


def walk_fan(cells: list, u: int, start=None, steps: int = 2):
    """扇を壁の辺 u から steps セル進んだ辺の先の列 [u, r1, r2, ...] (途中で終われば短い)。
    start = 最初に入るセル (扇の中の番号) = その壁の境界面を持つセル。後縁の線と側端の線の角の節点では上下の面が
    3 節点を共有し、扇が閉じて u が 2 セルに入るので、どちらの側へ進むかを面の持ち主のセルで決める。"""
    of = defaultdict(list)
    for k, (a, b, _, _) in enumerate(cells):
        of[a].append(k); of[b].append(k)
    if start is None:
        if len(of.get(u, [])) != 1:
            return [u], f"u_in_{len(of.get(u, []))}_cells"
        start = of[u][0]
    elif u not in cells[start][:2]:
        return [u], "u_not_in_face_cell"
    seq, prev, cur, nxt_cell = [u], None, u, start
    for _ in range(steps):
        if nxt_cell is None:
            cs = [k for k in of[cur] if k != prev]
            if len(cs) != 1:
                return seq, ("fan_end" if not cs else "fan_branch")
            nxt_cell = cs[0]
        a, b = cells[nxt_cell][:2]
        nxt = b if a == cur else a
        seq.append(nxt); prev = nxt_cell; cur = nxt; nxt_cell = None
    return seq, ""


def straight(mesh: Mesh, prev: int, cur: int):
    """格子線 prev→cur を cur の先へ直進した節点 (cur の隣のうち prev との共通の隣が cur だけの点)。一意でなければ None。"""
    np_ = mesh.nbrs(prev)
    c = [q for q in mesh.nbrs(cur) if q != prev and (mesh.nbrs(q) & np_) == {cur}]
    return c[0] if len(c) == 1 else None


def common_other(mesh: Mesh, a: int, b: int, excl: int):
    c = (mesh.nbrs(a) & mesh.nbrs(b)) - {excl}
    return next(iter(c)) if len(c) == 1 else None


def face_normal(P: np.ndarray) -> np.ndarray:
    """多角形の Newell 法線 (単位)。"""
    Q = P - P.mean(axis=0)
    n = np.zeros(3)
    for i in range(Q.shape[0]):
        n += np.cross(Q[i], Q[(i + 1) % Q.shape[0]])
    return n / max(np.linalg.norm(n), 1e-300)


def lsq_stats(mesh: Mesh, i: int):
    """内部隣接の数、M = Σ d̂d̂ᵀ の固有値 (昇順)、条件数、打ち切りで落ちる数。"""
    nb = sorted(mesh.nbrs(i))
    if not nb:
        return 0, [], float("inf"), 0
    d = mesh.X[nb] - mesh.X[i]
    if mesh.dim == 2:
        d = d[:, :2]
    L = np.linalg.norm(d, axis=1)
    d = d[L > 0] / L[L > 0][:, None]
    M = d.T @ d
    lam = np.linalg.eigvalsh(M)
    cond = float(lam[-1] / lam[0]) if lam[0] > 0 else float("inf")
    ndrop = int(np.count_nonzero(lam < LSQ_DROP * lam[-1]))
    return len(nb), lam.tolist(), cond, ndrop


def analyze(mesh: Mesh, bf: BFaces, ci: CellIndex, lines: list, names: dict, walls: list, n_layers: int) -> list:
    """全終端線の全節点 × 面 × 層の行を返す。"""
    X = mesh.X
    wallset = set(walls)
    rows = []
    for L in lines:
        by_node = defaultdict(list)
        for e in L["elems"]:
            for v in e["nodes"]:
                by_node[v].append(e)
        for t in L["nodes"]:
            e = by_node[t][0]
            elem = (t, e["nodes"][1] if e["nodes"][0] == t else e["nodes"][0]) if mesh.dim == 3 else (t,)
            cells = fan(mesh, ci, elem)
            phi_t = math.degrees(sum(c[2] for c in cells))
            side_faces = [fi for fi in e["faces"] if int(bf.pid[fi]) in wallset]
            sides = []
            for fi in side_faces:
                F = bf.face(fi).tolist()
                if mesh.dim == 3:
                    k = F.index(t)
                    cand = [F[k - 1], F[(k + 1) % len(F)]]
                    u = [v for v in cand if v != elem[1]]
                    u = u[0] if len(u) == 1 else None
                else:
                    u = F[1] if F[0] == t else F[0]
                sides.append((fi, int(bf.pid[fi]), u))
            pid_count = defaultdict(int)
            for _, pid, _ in sides:
                pid_count[pid] += 1
            for fi, pid, u in sides:
                flags = []
                if len(cells) != EXPECTED_FAN[L["cls"]]:
                    flags.append(f"fan{len(cells)}")
                fset = set(bf.face(fi).tolist())
                own = [k for k, c in enumerate(cells) if fset <= set(ci.C[c[3]].tolist())]
                if len(own) != 1:
                    flags.append(f"face_cell{len(own)}")
                seq, why = walk_fan(cells, u, own[0] if len(own) == 1 else None) if u is not None else ([None], "no_u")
                if why:
                    flags.append(why)
                n1 = seq[1] if len(seq) > 1 else None
                d = seq[2] if len(seq) > 2 else None
                # 壁法線 (その面の流体の側 = 第 1 層の側が正)
                if mesh.dim == 3:
                    nw = face_normal(X[bf.face(fi)])
                elif u is not None:
                    w = X[u] - X[t]
                    nw = np.array([-w[1], w[0], 0.0]) / max(np.hypot(w[0], w[1]), 1e-300)
                else:
                    nw = np.zeros(3)
                if n1 is not None and float((X[n1] - X[t]) @ nw) < 0.0:
                    nw = -nw
                side = names[pid]
                if pid_count[pid] > 1:
                    ax = int(np.argmax(np.abs(nw)))
                    side += f"({'+' if nw[ax] >= 0 else '-'}{'xyz'[ax]})"
                if len(e["faces"]) != 2:
                    flags.append(f"F{len(e['faces'])}")      # 端の要素を含む境界面が 2 枚でない (非多様体・開いた境界)
                if u is not None and not mesh.has_edge(t, u):
                    flags.append("not_edge_u")
                if d is not None and not mesh.has_edge(t, d):
                    flags.append("not_edge_d")
                # 層をたどる
                P, U, D = [t, n1], [u], [d]
                for k in range(1, n_layers + 1):
                    pk, pprev = P[k], P[k - 1]
                    U.append(common_other(mesh, U[k - 1], pk, pprev) if (pk is not None and U[k - 1] is not None) else None)
                    D.append(common_other(mesh, D[k - 1], pk, pprev) if (pk is not None and D[k - 1] is not None) else None)
                    P.append(straight(mesh, pprev, pk) if (pk is not None and pprev is not None) else None)
                for k in range(0, n_layers + 1):
                    pk, uk, dk, pn = P[k], U[k], D[k], P[k + 1]
                    r = {"line": L["id"], "cls": L["cls"], "tags": "|".join(L["tags"]), "side": side, "layer": k,
                         "node": pk, "up_node": uk, "down_node": dk, "next_node": pn, "phi_deg": phi_t, "fan_cells": len(cells)}
                    fl = list(flags) if k == 0 else []
                    if pk is None:
                        fl.append("no_p")
                    if pk is not None and (uk is None or dk is None or pn is None):
                        fl.append("amb_" + "".join(c for c, v in (("u", uk), ("d", dk), ("p", pn)) if v is None))
                    if pk is not None:
                        r.update(x=X[pk, 0], y=X[pk, 1], z=X[pk, 2])
                        nnb, lam, cond, ndrop = lsq_stats(mesh, pk)
                        r.update(n_lsq=nnb, lam=lam, cond=cond, n_drop=ndrop)
                    if None not in (pk, uk, dk):
                        a = X[pk] - X[uk]; b = X[dk] - X[pk]
                        La, Lb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
                        r.update(L_up=La, L_down=Lb, kink_deg=math.degrees(float(_angle(a, b))),
                                 len_ratio=Lb / La if La > 0 else float("inf"), disp=float(b @ nw))
                        if pn is not None:
                            h = float((X[pn] - X[pk]) @ nw)
                            r["h_layer"] = h
                            if h > 0:
                                r["disp_ratio"] = r["disp"] / h
                            else:
                                fl.append("h<=0")
                    r["flag"] = ",".join(fl)
                    rows.append(r)
    return rows


# ---------------------------------------------------------------------------
# 要約と出力
# ---------------------------------------------------------------------------
def _f(v, fmt):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        w = "".join(ch for ch in fmt.split(".")[0] if ch.isdigit())
        return "-".rjust(int(w)) if w else "-"
    return format(v, fmt)


def summarize(rows: list) -> list:
    groups = defaultdict(list)
    for r in rows:
        groups[(r["line"], r["side"], r["layer"])].append(r)
    out = []
    for (line, side, layer), rs in groups.items():
        ok = [r for r in rs if "kink_deg" in r]
        kinks = np.array([r["kink_deg"] for r in ok]) if ok else np.zeros(0)
        dr = [abs(r["disp_ratio"]) for r in ok if "disp_ratio" in r]
        lr = [r["len_ratio"] for r in ok]
        conds = [r["cond"] for r in rs if "cond" in r]
        lmin = [r["lam"][0] for r in rs if r.get("lam")]
        worst = max(ok, key=lambda r: r["kink_deg"]) if ok else None
        out.append({"line": line, "cls": rs[0]["cls"], "tags": rs[0]["tags"], "side": side, "layer": layer,
                    "n": len(rs), "n_ok": len(ok),
                    "kink_max": float(kinks.max()) if kinks.size else None, "kink_med": float(np.median(kinks)) if kinks.size else None,
                    "absdisp_ratio_max": max(dr) if dr else None,
                    "len_ratio_max": max(lr) if lr else None, "len_ratio_min": min(lr) if lr else None,
                    "cond_max": max(conds) if conds else None, "lam_min": min(lmin) if lmin else None,
                    "n_drop_max": max((r.get("n_drop", 0) for r in rs), default=0),
                    "n_flag": sum(1 for r in rs if r["flag"]),
                    "worst_node": worst["node"] if worst else None,
                    "worst_xyz": (worst["x"], worst["y"], worst["z"]) if worst else None})
    order = {}
    for r in rows:
        order.setdefault(r["line"], len(order))
    out.sort(key=lambda s: (order[s["line"]], s["side"], s["layer"]))
    return out


def sort_value(r: dict, key: str) -> float:
    if key == "kink":
        return r.get("kink_deg", -1.0)
    if key == "disp":
        return abs(r["disp_ratio"]) if "disp_ratio" in r else -1.0
    if key == "len":
        v = r.get("len_ratio")
        return max(v, 1.0 / v) if v and v > 0 and math.isfinite(v) else -1.0
    return r["cond"] if "cond" in r and math.isfinite(r["cond"]) else -1.0


NODE_COLS = ("line", "cls", "tags", "side", "layer", "node", "x", "y", "z", "up_node", "down_node", "next_node",
             "L_up", "L_down", "kink_deg", "disp", "h_layer", "disp_ratio", "len_ratio",
             "n_lsq", "lam1", "lam2", "lam3", "cond", "n_drop", "phi_deg", "fan_cells", "flag")
LINE_COLS = ("line", "cls", "tags", "side", "layer", "n", "n_ok", "kink_max", "kink_med", "absdisp_ratio_max",
             "len_ratio_max", "len_ratio_min", "cond_max", "lam_min", "n_drop_max", "n_flag",
             "worst_node", "worst_x", "worst_y", "worst_z",
             "line_n_nodes", "line_phi_min", "line_phi_max", "line_length", "bbox_min_x", "bbox_min_y", "bbox_min_z",
             "bbox_max_x", "bbox_max_y", "bbox_max_z")


def write_csv(prefix: str, rows: list, summ: list, lines: list):
    pn, pl = prefix + "_nodes.csv", prefix + "_lines.csv"
    with open(pn, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(NODE_COLS)
        for r in rows:
            lam = list(r.get("lam") or []) + [None] * 3
            v = dict(r, lam1=lam[0], lam2=lam[1], lam3=lam[2])
            w.writerow(["" if v.get(c) is None else v.get(c) for c in NODE_COLS])
    lmap = {L["id"]: L for L in lines}
    with open(pl, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(LINE_COLS)
        for s in summ:
            L = lmap[s["line"]]
            xyz = s["worst_xyz"] or (None, None, None)
            v = dict(s, worst_x=xyz[0], worst_y=xyz[1], worst_z=xyz[2], line_n_nodes=len(L["nodes"]),
                     line_phi_min=L["phi"][0], line_phi_max=L["phi"][1], line_length=L["length"],
                     bbox_min_x=L["bbox"][0][0], bbox_min_y=L["bbox"][0][1], bbox_min_z=L["bbox"][0][2],
                     bbox_max_x=L["bbox"][1][0], bbox_max_y=L["bbox"][1][1], bbox_max_z=L["bbox"][1][2])
            w.writerow(["" if v.get(c) is None else v.get(c) for c in LINE_COLS])
    return pl, pn


def _xyz(v):
    return "(" + ", ".join(f"{float(c):.6g}" for c in v) + ")"


def run(args) -> dict:
    """CLI の本体 (試験から呼ぶ)。結果の dict を返す。"""
    th = {"knife": args.knife_deg, "convex": args.convex_deg, "open": args.open_deg}
    if not (th["open"] < th["convex"] < th["knife"] <= 360.0):
        raise WallEndError("閾値は open < convex < knife ≤ 360 で指定する")
    if args.layers < 0:
        raise WallEndError("--layers は 0 以上")
    mesh = Mesh(args.h5)
    names, walls, notes, _ = resolve_tags(mesh, args.bcond_config, args.wall, args.walls)
    bf = BFaces(mesh)
    by_name = {v: k for k, v in names.items()}
    pair_res = []
    for pr in args.pair or []:
        ab = [s.strip() for s in pr.split(",") if s.strip()]
        if len(ab) != 2 or ab[0] == ab[1]:
            raise WallEndError(f"--pair は異なる 2 タグ A,B ('{pr}')")
        for nm in ab:
            if nm not in by_name:
                raise WallEndError(f"--pair のタグ '{nm}' が無い (境界条件名の誤記? 既知: {sorted(by_name)})")
        pair_res.append((ab[0], ab[1], shared_count(mesh, bf, by_name[ab[0]], by_name[ab[1]])))
    ends, n_cand, Cs = find_ends(mesh, bf, walls, th)
    lines = build_lines(mesh, ends, names, args.turn_deg)
    ci = CellIndex(Cs, mesh.n)
    rows = analyze(mesh, bf, ci, lines, names, walls, args.layers)
    summ = summarize(rows)
    end_keys = {}
    for e in ends:
        k = (min(e["nodes"]) * mesh.n + max(e["nodes"])) if mesh.dim == 3 else e["nodes"][0]
        end_keys[k] = e["cls"]
    pair_out = []
    for a, b, keys in pair_res:
        nk = sum(1 for k in keys.tolist() if end_keys.get(k) == "knife")
        pair_out.append({"a": a, "b": b, "n_shared": int(keys.size), "n_knife": int(nk)})

    # ---- 表示 ----
    P = print if not args.quiet else (lambda *a, **k: None)
    unit = "辺" if mesh.dim == 3 else "端点"
    P(f"[{TOOL}] {args.h5}: {mesh.dim}D、節点 {mesh.n}、{'hex' if mesh.dim == 3 else 'quad'} {mesh.C.shape[0]}、内部面 {mesh.pairs.shape[0]}")
    P("  壁: " + ", ".join(f"{names[p]}:{p}({mesh.bc[p]['kind']})" for p in walls))
    for s in notes:
        P("  " + s)
    ncls = defaultdict(int)
    for e in ends:
        ncls[e["cls"]] += 1
    P(f"  壁の{unit}の候補 {n_cand} → 端 knife {ncls['knife']}・convex {ncls['convex']}・open {ncls['open']} "
      f"(閾値 φ: knife ≥ {th['knife']:g}°、convex > {th['convex']:g}°、open ≥ {th['open']:g}° かつ壁でない面に接する)")
    for p in pair_out:
        P(f"  --pair {p['a']},{p['b']}: 同じ節点 ID で共有する{unit} {p['n_shared']} (うち knife {p['n_knife']})")
    if not lines:
        P("  壁の終端は見つからない")
    else:
        P("\n終端線 (φ = 流体角 [deg]、長さ・範囲は h5 の座標の単位):")
        P(f"  {'id':<5} {'分類':<6} {'タグ':<28} {'節点':>5} {'長さ':>10} {'φ min':>7} {'φ max':>7}  {'向き':<22} 範囲")
        for L in lines:
            P(f"  {L['id']:<5} {L['cls']:<6} {'|'.join(L['tags']):<28} {len(L['nodes']):>5} {L['length']:>10.4g} "
              f"{L['phi'][0]:>7.1f} {L['phi'][1]:>7.1f}  {_xyz(np.round(L['dir'], 3) + 0.0):<22} {_xyz(L['bbox'][0])} – {_xyz(L['bbox'][1])}")
        P("\n要約 (終端線 × 面 × 層。kink = 前後の辺の角度差 [deg]、disp/h = 下流の辺の壁法線方向の変位 / 層間隔、"
          "len = 下流/上流の辺長、cond = LSQ の M の条件数):")
        P(f"  {'id':<5} {'面':<18} {'層':>2} {'節点':>5} {'kink最大':>8} {'kink中央':>8} {'|disp/h|最大':>12} {'len最大':>8} "
          f"{'len最小':>8} {'cond最大':>9} {'drop':>4} {'flag':>4}  kink 最大の節点 (座標)")
        for s in summ:
            P(f"  {s['line']:<5} {s['side']:<18} {s['layer']:>2} {s['n']:>5} {_f(s['kink_max'], '8.2f')} {_f(s['kink_med'], '8.2f')} "
              f"{_f(s['absdisp_ratio_max'], '12.4g')} {_f(s['len_ratio_max'], '8.3f')} {_f(s['len_ratio_min'], '8.3f')} "
              f"{_f(s['cond_max'], '9.3g')} {s['n_drop_max']:>4} {s['n_flag']:>4}  "
              + (f"{s['worst_node']} {_xyz(s['worst_xyz'])}" if s["worst_node"] is not None else "-"))
        top = sorted(rows, key=lambda r: sort_value(r, args.sort), reverse=True)[:args.top]
        P(f"\n上位 {len(top)} 節点 (--sort {args.sort}):")
        P(f"  {'id':<5} {'面':<18} {'層':>2} {'節点':>9} {'u→t→d':<22} {'L_up':>10} {'L_down':>10} {'kink':>7} {'disp':>11} "
          f"{'h':>10} {'disp/h':>10} {'len':>7} {'nLSQ':>4} {'λ (昇順)':<24} {'cond':>8} {'φ':>6}  座標 / flag")
        for r in top:
            lam = "(" + ", ".join(f"{x:.3g}" for x in (r.get("lam") or [])) + ")"
            P(f"  {r['line']:<5} {r['side']:<18} {r['layer']:>2} {_f(r['node'], '>9')} "
              f"{(str(r['up_node']) + '→' + str(r['node']) + '→' + str(r['down_node'])):<22} "
              f"{_f(r.get('L_up'), '10.4g')} {_f(r.get('L_down'), '10.4g')} {_f(r.get('kink_deg'), '7.2f')} "
              f"{_f(r.get('disp'), '11.4g')} {_f(r.get('h_layer'), '10.4g')} {_f(r.get('disp_ratio'), '10.4g')} "
              f"{_f(r.get('len_ratio'), '7.3f')} {_f(r.get('n_lsq'), '>4')} {lam:<24} {_f(r.get('cond'), '8.3g')} "
              f"{r['phi_deg']:>6.1f}  " + (_xyz((r['x'], r['y'], r['z'])) if 'x' in r else "-")
              + (f" [{r['flag']}]" if r["flag"] else ""))
    csv_paths = None
    if not args.no_csv:
        prefix = args.out_prefix or (os.path.splitext(args.h5)[0] + "_wall_ends")
        csv_paths = write_csv(prefix, rows, summ, lines)
        P(f"\nCSV: {csv_paths[0]} (要約)、{csv_paths[1]} (全行)")
    for p in pair_out:
        if p["n_shared"] == 0:
            raise WallEndError(f"--pair {p['a']},{p['b']}: 同じ節点 ID で共有する{unit}が無い "
                               "(タグの誤記・有限厚の板・座標一致の別 ID のどれか)")
    return {"mesh": mesh, "names": names, "walls": walls, "ends": ends, "lines": lines, "rows": rows, "summary": summ,
            "pairs": pair_out, "csv": csv_paths, "n_cand": n_cand}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog=__doc__[__doc__.index("定義:"):],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5", help="node 変換済みの入力 h5 (convertGmshToForge の出力、例 run_*/sern.h5)。読むだけ")
    ap.add_argument("--bcond-config", help="run の bcondConfig.yaml (タグ名 → physID と kind)。kind が wall / wall_isothermal のタグを壁にする")
    ap.add_argument("--wall", action="append", help="NAME=physID (複数回)。指定したタグだけを壁にする (bcondConfig より優先)")
    ap.add_argument("--walls", help="壁を名前で絞る (例 cowl_in,cowl_out,ramp)。名前は --bcond-config / --wall のもの")
    ap.add_argument("--pair", action="append", help="A,B: 2 タグが同じ節点 ID で共有する辺 (2D は端点) の数と knife の数を出す。0 ならエラー (複数回可)")
    ap.add_argument("--layers", type=int, default=3, help="端の節点から壁法線方向へ何層まで見るか (既定 3 → 第 0〜3 層)")
    ap.add_argument("--top", type=int, default=20, help="詳細を出す上位の行数 (既定 20)")
    ap.add_argument("--sort", choices=SORT_KEYS, default="kink", help="上位の並べ替え: kink / disp (|disp/h|) / len (max(len,1/len)) / cond")
    ap.add_argument("--knife-deg", type=float, default=300.0, help="knife の流体角の下限 [deg] (既定 300)")
    ap.add_argument("--convex-deg", type=float, default=225.0, help="convex の流体角の下限 (これより大) [deg] (既定 225。45° 未満の壁の折れを除く)")
    ap.add_argument("--open-deg", type=float, default=160.0, help="open の流体角の下限 [deg] (既定 160)")
    ap.add_argument("--turn-deg", type=float, default=45.0, help="終端線を分ける辺の向きの変化 [deg] (既定 45)")
    ap.add_argument("--out-prefix", help="CSV の接頭辞 (既定: <h5 のパスから .h5 を除いたもの>_wall_ends)")
    ap.add_argument("--no-csv", action="store_true", help="CSV を書かない")
    ap.add_argument("--quiet", action="store_true", help="表を出さない (試験用)")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run(args)
        return 0
    except WallEndError as e:
        print(f"[{TOOL}] ERROR: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
