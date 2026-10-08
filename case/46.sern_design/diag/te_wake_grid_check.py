#!/usr/bin/env python3
"""後縁下流の中間線の局所変形 (格子 B) を投入前に確かめる: 旧格子 A と B を並べて表にする (読むだけ)。

plan convection-zero-thickness-edge-reconstruction §4.2「投入前の格子の確認」・§5.1 #2c (codex diagnose 2026-10-08 te-grid-kink)。

  (i)   後縁の列 (i = i_te) の全節点・全 z の**最初の下流の辺**の角度。下側の層 (j < jm)・上側の層 (j > jm)・後縁節点 (j = jm)
        を分け、上流の辺の角度と折れ角 (下流 − 上流) を並べる。後縁節点は下側壁面 (cowl_out) と上側壁面 (cowl_in) の
        両方の上流の辺に対して出す (g3 で約 −4.4° / −5.6° と別物)
  (ii)  i 線 (j, k 固定) に沿った隣接辺の角度変化 (折れ角) の最大値とその位置。後縁 (i = i_te)・復帰区間
        (i_te < i ≤ B の te_wake_i_last + 1)・それ以外に分ける。上側壁面の i 線 (j = jm のカウル上コピー) も別に出す
  (iii) 近壁間隔: 後縁の列と復帰区間の各列で、中間線に近い n 層の y 方向の間隔の B/A と、i 線に垂直な間隔 (Δy·cos θ)
  (iv)  全ヘキサの 8 頂点のヤコビアン (各頂点で 3 方向の前進差分の行列式) が正 (多数派の向き) か
  (v)   節点数・x・z・ヘキサの接続・境界面の接続が A と同一か。動いた節点 (y が違う) の数・範囲・量・境界タグ

使い方:
  python3 case/46.sern_design/diag/te_wake_grid_check.py A B [--n-layers 10] [--info-a JSON] [--info-b JSON]

  A / B = run ディレクトリ (sern.h5、無ければ sern.msh と prepare_info.json を読む)、または格子ファイル
  (.h5 = forge の入力 h5 か res_*.h5、.msh = メッシャの出力) と --info-a/--info-b (prepare_info.json、またはメッシャの info の JSON)。
  格子の構造 (ni, NJ, nz, jm, i_te, k_sw, i_sw, n_dup_cowl) はメッシャの info から取る。節点番号はメッシャの順
  (base(i, j, k) = (i·NJ + j)·nz + k、その後にカウル上コピー・側壁外コピー・外部流ブロック) で、変換器はこの順を保つ
  (`mesh.renumber: rcm` で再番号付けした h5 は /MESH/RENUMBER_PERM でメッシャの番号に戻して読む)。
  run ディレクトリの格子は solverConfig.yaml の `mesh.meshFileName` (無ければ sern.h5、それも無ければ sern.msh)。

品質 (`check_mesh_quality.py`) と双対幾何 (`check_dual_closure.py`) の VERDICT は本道具では出さない。回し方は出力の末尾に書く。

--admission (plan §6.0「投入条件 = 最終 node 格子の検査」の表を判定する。許容差は表のまま、測る前に固定):
  python3 case/46.sern_design/diag/te_wake_grid_check.py A_run B_run --admission \\
      [--quality-a MESH_QUALITY.txt] [--quality-b MESH_QUALITY.txt] --dual-closure-b CLOSURE_B.txt [--dual-closure-a CLOSURE_A.txt]

  表の各行を PASS / FAIL / 判定不能 で出し (行ごとに測った値と閾値)、最後に `ADMISSION VERDICT: PASS|FAIL|UNDECIDABLE` を
  1 行出す (終了コード 0 / 1 / 2)。1 つでも FAIL なら FAIL、FAIL が無く判定不能が 1 つでもあれば UNDECIDABLE
  (AGENTS.md「判定ツールが判定不能を返したら、それは合格ではない」)。行の定義 (第 n 層 = 中間線から n 番目の j、全 span = 全 k):
    前提            A・B とも最終 node 格子 (node 変換の h5)、L_b が A 0・B 1.0 H (事前登録の値)、格子の構造が同じ
    最初の下流の辺  B の後縁の列 (i = i_te) の下側の第 1〜第 10 層 (j = jm−1..jm−10)・全 k で、辺 (i_te → i_te+1) の角度と
                    下側の接線 m0 (B の info の te_wake_slope_te_deg) の差 ≤ 3°
    後縁の折れ      B の i = i_te で上下の第 1〜第 10 層と共有の後縁節点 (j = jm)・全 k の折れ角 (下流の辺 − 同じ i 線の上流の辺)。
                    後縁節点は下側壁面 (cowl_out、板の外は中間線) と上側壁面 (cowl_in) の上流の辺の両方に対して ≤ 3°
    復帰区間の折れ  B の主ブロックの全 i 線 (全 j・全 k) で、節点 i_te < i ≤ i_last + 1 (= 最初の未変形の辺を含む) の折れ角 ≤ 6°
    固定すべき量    節点数・x・z (全節点ビット一致)・ヘキサの接続・境界面の接続 (タグ別)・上流 (x ≤ 後縁の x) の座標・
                    固体 (壁タグ) と外部境界 (sym・side_far 以外) の節点が不変。sym・side_far 上の面内 (y) の移動は数を記録
    層の Δy         i = i_te .. i_last + 1 の列・上下・全 k。第 1 層 |Δy_B − Δy_A| ≤ 1e−8·Δy_A + r、第 1〜第 10 層 ≤ 0.05·Δy_A + r。
                    r = 座標の丸めの上限 = 4·ε₃₂·max|対象の y| (float32 の h5)、0 (倍精度 = 変換前)。msh (10 桁) は判定不能
    primal の品質   B の全ヘキサの 8 頂点のヤコビアン (corner_jacobians、向きは A の多数派) の非正・非有限 0、スケール済みの最小 ≥ 0.65、
                    skew (check_mesh_quality.py の metrics_vectorized、3D の面ごとの equiangle skew) 最大 ≤ 0.9、A からの増加 ≤ 0.10
    品質・双対      (a) MESH_QUALITY.txt の VERDICT が B ≥ A (PASS > SOFT-PASS > FAIL、閾値 AR・skew が同じ、cells が格子のヘキサ数)
                    (b) 変形領域 (動いた節点を含むヘキサ) に新しい AR 閾値超過 (B で超過・A で非超過) が 0 (AR は check_mesh_quality と同じ定義を
                        この道具で計算、閾値は MESH_QUALITY.txt の AR<=、無ければ --ar-max)
                    (c) --dual-closure-b のファイル (check_dual_closure.py の出力) が VERDICT: PASS、tol 1e−5、CV 数 = B の節点数
  入力が欠ける・読めない項目は判定不能 (合格扱いにしない)。

--provenance-out FILE (plan §4.2「来歴」): A・B の実入力格子 (最終 node の h5) から再計算した格子署名
  (`solver_density_cuda/tools/mark_zero_thickness_edges.py` の mesh_signature) と双対幾何のハッシュ (`stage_manifest.dual_geometry_hash`、
  署名に入らない変換器の違いを検出する。plan tooling-sern-te-wake-grid §4)、曲線版 (メッシャの info の te_wake_curve_version)、
  te_wake_blend_H、生成コードの commit (prepare_info に記録があれば。無ければ --code-commit-a/-b の申告を「未検証」として)、
  設定ファイル (run 直下の *.yaml・prepare_info.json・問題 YAML) の sha256、--attach-a/-b の記録 (初期場の転送ログ等) の sha256 を JSON で書く。
  --admission と併用すると判定の行も入れる。

--same-mesh A.h5 B.h5: 2 つの h5 (か run ディレクトリ) の格子署名と双対幾何のハッシュを比べる (継続 run で格子が不変か)。
  両方が同じときだけ YES (署名が同じでも変換器の違いで双対幾何が違えば NO。plan tooling-sern-te-wake-grid §4)。
  `SAME MESH: YES|NO|UNDECIDABLE` (終了コード 0 / 1 / 2)。
"""
import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HEX_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]   # gmsh 型 5 の頂点順
_IDX = {c: n for n, c in enumerate(HEX_CORNERS)}
PHYS_NAME = {1: "inlet_nozzle", 2: "inlet_ext", 3: "outlet", 4: "ramp", 5: "cowl_in", 6: "cowl_out", 7: "bottom", 8: "top_out",
             9: "sym", 10: "side_far", 11: "sidewall_in", 12: "sidewall_out", 14: "vehicle", 15: "vehicle_top",
             16: "underside_far", 17: "vehicle_side", 18: "vehicle_base"}   # mesh_sern3d.PHYS_SERN3D
IN_PLANE_TAGS = (9, 10)   # sym (z = 0)・side_far (z = Z_far): 節点が面内 (y) で動くのは形状の変更ではない
PHYS_ID = {v: k for k, v in PHYS_NAME.items()}
WALL_TAGS = (4, 5, 6, 11, 12, 14, 15, 17, 18)   # 固体形状: ramp・cowl_in・cowl_out・sidewall_*・vehicle*

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "solver_density_cuda" / "tools"
DESIGN = ROOT / "design"

# plan §6.0「投入条件 = 最終 node 格子の検査」の許容差 (事前登録 2026-10-08、codex diagnose te-wake-blend-length)。変えない
ADM_LAYERS = 10                 # 第 1〜第 10 層
ADM_FIRST_EDGE_DEG = 3.0        # 最初の下流の辺: m0 から 3° 以内
ADM_TE_FOLD_DEG = 3.0           # 後縁の折れ ≤ 3°
ADM_RETURN_FOLD_DEG = 6.0       # 復帰区間の折れ ≤ 6°
ADM_DY_FIRST_REL = 1.0e-8       # 第 1 層の Δy の相対差
ADM_DY_LAYERS_REL = 0.05        # 第 1〜第 10 層の各層厚
ADM_ROUND_ULPS = 4.0            # 変換後の座標の丸めの上限 = 4·ε₃₂·max|対象座標|
ADM_SCALED_JAC_MIN = 0.65
ADM_SKEW_MAX = 0.9
ADM_SKEW_INC_MAX = 0.10
ADM_DUAL_TOL = 1.0e-5
ADM_BLEND_H = 1.0               # 事前登録の L_b (§4.2: 1.0 H に固定。結果を見て変えない)
EPS32 = float(np.finfo(np.float32).eps)


# ---------------------------------------------------------------------------------------------- 読み込み
def _read_msh41(path):
    """mesh_sern3d.write_msh41_3d の出力 (gmsh 4.1 ASCII) → (coords, hexes, {physID: (sizes, nodes)})。"""
    with open(path) as f:
        lines = f.read().split("\n")
    p = 0
    ent_phys = {}
    coords = hexes = None
    bf = {}
    while p < len(lines):
        s = lines[p].strip()
        if s == "$Entities":
            npt, ncu, nsu, nvo = map(int, lines[p + 1].split()); q = p + 2 + npt + ncu
            for _ in range(nsu):
                t = lines[q].split(); tag = int(t[0]); nph = int(t[7])
                ent_phys[tag] = int(t[8]) if nph > 0 else None; q += 1
            p = q
        elif s == "$Nodes":
            nb, nn = map(int, lines[p + 1].split()[:2]); q = p + 2
            coords = np.zeros((nn, 3))
            for _ in range(nb):
                _, _, _, cnt = map(int, lines[q].split()); q += 1
                tags = np.array([int(x) for x in lines[q:q + cnt]]) - 1; q += cnt
                coords[tags] = np.array([[float(v) for v in l.split()[:3]] for l in lines[q:q + cnt]]); q += cnt
            p = q
        elif s == "$Elements":
            nb = int(lines[p + 1].split()[0]); q = p + 2
            hx = []
            for _ in range(nb):
                dim, etag, typ, cnt = map(int, lines[q].split()); q += 1
                arr = np.array([[int(v) for v in l.split()[1:]] for l in lines[q:q + cnt]], dtype=np.int64) - 1; q += cnt
                if typ == 5:
                    hx.append(arr)
                elif typ == 3:
                    pid = ent_phys.get(etag)
                    sz, nd = bf.get(pid, (np.zeros(0, np.int64), np.zeros(0, np.int64)))
                    bf[pid] = (np.concatenate([sz, np.full(len(arr), 4, np.int64)]), np.concatenate([nd, arr.reshape(-1)]))
            hexes = np.concatenate(hx) if hx else np.zeros((0, 8), np.int64)
            p = q
        else:
            p += 1
    return coords, hexes, bf


def _read_h5(path):
    """forge の入力 h5 (VIZMESH/CONNE) か res h5 (MESH/CONNE) → (coords, hexes, {physID: (sizes, nodes)} or None)。"""
    import h5py
    with h5py.File(path, "r") as f:
        coords = np.asarray(f["MESH/COORD"][...], dtype=np.float64).reshape(-1, 3)
        hexes = None
        for key in ("VIZMESH/CONNE", "MESH/CONNE"):
            if key in f:
                c = np.asarray(f[key][...], dtype=np.int64)
                if c.size and c.size % 9 == 0 and np.all(c[::9] == 9):
                    hexes = c.reshape(-1, 9)[:, 1:]
                    break
        if hexes is None:
            raise SystemExit(f"{path}: ヘキサの接続 (VIZMESH/CONNE か MESH/CONNE の XDMF 型 9 だけの列) が無い")
        bf = None
        if "BCONDS" in f:
            bf = {}
            for k in f["BCONDS"]:
                g = f["BCONDS/" + k]
                if "vizBfaceSizes" in g and "vizBfaceNodes" in g:
                    bf[int(k)] = (np.asarray(g["vizBfaceSizes"][...], np.int64), np.asarray(g["vizBfaceNodes"][...], np.int64))
        perm = np.asarray(f["MESH/RENUMBER_PERM"][...], np.int64) if "MESH/RENUMBER_PERM" in f else None
    if perm is not None:
        # 変換器の RCM 再番号付け (mesh.renumber: rcm): perm[新] = 旧 (メッシャの番号)。メッシャの番号に戻して構造を当てる
        if perm.shape[0] != coords.shape[0]:
            raise SystemExit(f"{path}: MESH/RENUMBER_PERM の長さ {perm.shape[0]} が節点数 {coords.shape[0]} と違う")
        c0 = np.empty_like(coords); c0[perm] = coords; coords = c0
        hexes = perm[hexes]
        if bf is not None:
            bf = {k: (sz, perm[nd]) for k, (sz, nd) in bf.items()}
        print(f"[te_wake_grid_check] {path}: 再番号付け (RENUMBER_PERM) をメッシャの番号に戻して読んだ")
    return coords, hexes, bf


def mesh_file_of(run_dir):
    """run の solverConfig.yaml の `mesh.meshFileName` (forge が読む格子)。無ければ None。"""
    p = os.path.join(run_dir, "solverConfig.yaml")
    if not os.path.exists(p):
        return None
    m = re.search(r"meshFileName\s*:\s*[\"']?([^\"',}\s]+)", open(p).read())
    return m.group(1) if m else None


def grid_of(path):
    """run ディレクトリなら forge が読む格子 (meshFileName → sern.h5 → sern.msh)、ファイルならそのまま。"""
    if not os.path.isdir(path):
        return path
    for c in (mesh_file_of(path), "sern.h5", "sern.msh"):
        if c and os.path.exists(os.path.join(path, c)):
            return os.path.join(path, c)
    return os.path.join(path, mesh_file_of(path) or "sern.h5")


def _h5_meta(path):
    """h5 の座標の型と、node (median-dual) 変換か (VIZMESH があり CV 数 = 節点数。mark_zero_thickness_edges と同じ判定)。"""
    import h5py
    with h5py.File(path, "r") as f:
        dt = str(f["MESH/COORD"].dtype)
        a = f["MESH"].attrs
        node = "VIZMESH" in f and "nCells" in a and "nNodes" in a and int(a["nCells"]) == int(a["nNodes"])
    return {"coord_dtype": dt, "node": bool(node)}


def load(path, info_path=None):
    """格子と構造情報。戻り = dict(coords [m], hexes, bfaces, info (メッシャ), H, label, 入力の種別)。
    入力の種別: source = h5 / msh、coord_repr = float32 (変換後の h5) / float64 / msh10g (msh の 10 桁)、node_h5 = node 変換の h5。"""
    run_dir = path if os.path.isdir(path) else None
    grid = grid_of(path)
    if run_dir:
        info_path = info_path or os.path.join(path, "prepare_info.json")
    if not info_path or not os.path.exists(info_path):
        raise SystemExit(f"{path}: 構造情報 (prepare_info.json か --info-*) が無い")
    if not os.path.exists(grid):
        raise SystemExit(f"{path}: 格子 (sern.h5 / sern.msh) が無い")
    raw = json.load(open(info_path))
    info = raw["mesh"] if "mesh" in raw else raw
    H = float(raw.get("H_m", info.get("H_m", 1.0)))
    if grid.endswith(".msh"):
        coords, hexes, bf = _read_msh41(grid)
        kind = {"source": "msh", "coord_repr": "msh10g", "node_h5": False}
    else:
        coords, hexes, bf = _read_h5(grid)
        meta = _h5_meta(grid)
        kind = {"source": "h5", "coord_repr": "float32" if meta["coord_dtype"] == "float32" else "float64", "node_h5": meta["node"]}
    if coords.shape[0] != int(info["nodes"]):
        raise SystemExit(f"{grid}: 節点数 {coords.shape[0]} が info の {info['nodes']} と違う (別の格子の info)")
    return {"coords": coords, "hexes": hexes, "bfaces": bf, "info": info, "H": H, "label": grid, **kind,
            "grid_path": grid, "run_dir": run_dir, "info_path": info_path}


def from_mesher(coords, hexes, B, info, H, label="(メモリ上の格子)"):
    """メッシャの出力 (generate_sern_mesh3d の coords [m]・hexes・{タグ名: [(a, b, c, d), ...]}・info) を load() と同じ形にする。
    座標は倍精度のまま (= plan §6.0 の「変換前」)。"""
    bf = {}
    for nm, q in B.items():
        if len(q):
            q = np.asarray(q, np.int64).reshape(-1, 4)
            bf[int(PHYS_ID[nm])] = (np.full(len(q), 4, np.int64), q.reshape(-1))
    return {"coords": np.asarray(coords, np.float64), "hexes": np.asarray(hexes, np.int64), "bfaces": bf, "info": info,
            "H": float(H), "label": label, "source": "memory", "coord_repr": "float64", "node_h5": False,
            "grid_path": None, "run_dir": None, "info_path": None}


# ---------------------------------------------------------------------------------------------- 構造
class Struct:
    """メッシャの節点番号 (mesh_sern3d.generate_sern_mesh3d) を再現する。"""

    def __init__(self, info):
        self.ni, self.NJ, self.nz = int(info["ni"]), int(info["NJ"]), int(info["nz"])
        self.jm, self.i_te, self.k_sw, self.i_sw = int(info["jm"]), int(info["i_te"]), int(info["k_sw"]), int(info["i_sw"])
        self.N_base = self.ni * self.NJ * self.nz
        no_outer = self.nz == self.k_sw + 1
        full = self.i_te * (self.k_sw + 1)
        shared = full - (max(self.i_te - self.i_sw, 0) if not no_outer else 0)
        n = int(info["n_dup_cowl"])
        if n == shared:
            share = True
        elif n == full:
            share = False
        else:
            raise SystemExit(f"カウル上コピーの数 {n} がメッシャの規則 ({shared} / {full}) と合わない")
        self.dup1 = {}
        nid = self.N_base
        for i in range(self.i_te):
            for k in range(self.k_sw + 1):
                if k == self.k_sw and i >= self.i_sw and not no_outer and share:
                    continue
                self.dup1[(i, k)] = nid; nid += 1

    def base(self, i, j, k):
        return (i * self.NJ + j) * self.nz + k

    def grid(self, coords):
        """主ブロックの base 節点を (ni, NJ, nz, 3) に並べる。"""
        return coords[:self.N_base].reshape(self.ni, self.NJ, self.nz, 3)

    def upper_surface_line(self, coords, k):
        """j = jm の上側壁面 (cowl_in) の i 線 (後縁より上流はカウル上コピー)。(ni, 3)"""
        ids = [self.dup1.get((i, k), self.base(i, self.jm, k)) for i in range(self.ni)]
        return coords[np.asarray(ids)]


def _ang(e):
    return np.degrees(np.arctan2(e[..., 1], e[..., 0]))


def _wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def _stats(a):
    a = np.asarray(a, float).ravel()
    return (float(a.min()), float(np.median(a)), float(a.max())) if a.size else (np.nan, np.nan, np.nan)


def _fmt3(t):
    return f"{t[0]:8.3f} {t[1]:8.3f} {t[2]:8.3f}"


# ---------------------------------------------------------------------------------------------- (i)
def te_column(g, S, n_layers):
    """後縁の列の下流・上流の辺の角度。戻り = {行名: (下流 (min,med,max), 上流 (min,med,max), 折れ (max |.|), 節点数)}"""
    G = S.grid(g["coords"])
    it, jm = S.i_te, S.jm
    th_d = _ang(G[it + 1] - G[it])         # (NJ, nz) 最初の下流の辺
    th_u = _ang(G[it] - G[it - 1])         # (NJ, nz) 最後の上流の辺 (j = jm は下側壁面 / 板の外は中間線)
    up_upper = np.array([_ang(G[it, jm, k] - g["coords"][S.dup1.get((it - 1, k), S.base(it - 1, jm, k))]) for k in range(S.nz)])
    kin = np.arange(S.nz) <= S.k_sw
    rows = {}

    def add(name, jsel, ksel, up=None):
        d = th_d[np.ix_(jsel, np.where(ksel)[0])] if np.ndim(jsel) else th_d[jsel, ksel]
        u = (th_u[np.ix_(jsel, np.where(ksel)[0])] if np.ndim(jsel) else th_u[jsel, ksel]) if up is None else up[ksel]
        if np.size(d) == 0:
            return
        rows[name] = (_stats(d), _stats(u), float(np.max(np.abs(_wrap(d - u)))), int(np.size(d)))
    lo_near = np.arange(max(jm - n_layers, 0), jm); up_near = np.arange(jm + 1, min(jm + 1 + n_layers, S.NJ))
    add("後縁節点 j=jm  z≤W/2  上流=下側壁面 (cowl_out)", jm, kin)
    add("後縁節点 j=jm  z≤W/2  上流=上側壁面 (cowl_in)", jm, kin, up=up_upper)
    add("後縁節点 j=jm  z>W/2  上流=中間線 (板なし)", jm, ~kin)
    add(f"下側の層 j=jm-1..jm-{n_layers}  z≤W/2", lo_near, kin)
    add(f"下側の層 j=jm-1..jm-{n_layers}  z>W/2", lo_near, ~kin)
    add(f"上側の層 j=jm+1..jm+{n_layers}  z≤W/2", up_near, kin)
    add(f"上側の層 j=jm+1..jm+{n_layers}  z>W/2", up_near, ~kin)
    add("下側 全層 (j<jm)", np.arange(0, jm), np.ones(S.nz, bool))
    add("上側 全層 (j>jm)", np.arange(jm + 1, S.NJ), np.ones(S.nz, bool))
    return rows


# ---------------------------------------------------------------------------------------------- (ii)
def turn_maxima(g, S, i_ret):
    """i 線の折れ角 |Δθ| の最大とその位置。戻り = {(側, 区間): (max, (i, j, k), xyz)}"""
    G = S.grid(g["coords"])
    th = _ang(np.diff(G, axis=0))                         # (ni-1, NJ, nz)
    turn = np.abs(_wrap(np.diff(th, axis=0)))             # (ni-2, NJ, nz): 節点 i = 1..ni-2
    ii = np.arange(1, S.ni - 1)
    reg = {"後縁 (i=i_te)": ii == S.i_te,
           f"復帰区間 (i_te<i≤{i_ret})": (ii > S.i_te) & (ii <= i_ret),
           "それ以外": (ii != S.i_te) & ~((ii > S.i_te) & (ii <= i_ret))}
    sides = {"下側 (j<jm)": np.arange(S.NJ) < S.jm, "中間線/下側壁面 (j=jm)": np.arange(S.NJ) == S.jm,
             "上側 (j>jm)": np.arange(S.NJ) > S.jm}
    out = {}
    for sn, sm in sides.items():
        for rn, rm in reg.items():
            sub = turn[np.ix_(rm, sm, np.ones(S.nz, bool))]
            if sub.size == 0:
                continue
            a, b, c = np.unravel_index(int(np.argmax(sub)), sub.shape)
            i = int(ii[rm][a]); j = int(np.where(sm)[0][b]); k = int(c)
            out[(sn, rn)] = (float(sub[a, b, c]), (i, j, k), G[i, j, k])
    # 上側壁面 (カウル上コピー) の i 線、z ≤ W/2
    lines = np.stack([S.upper_surface_line(g["coords"], k) for k in range(S.k_sw + 1)], axis=1)   # (ni, k_sw+1, 3)
    thu = _ang(np.diff(lines, axis=0)); tu = np.abs(_wrap(np.diff(thu, axis=0)))                  # (ni-2, k_sw+1)
    for rn, rm in reg.items():
        sub = tu[rm]
        if sub.size == 0:
            continue
        a, c = np.unravel_index(int(np.argmax(sub)), sub.shape)
        i = int(ii[rm][a])
        out[("上側壁面 (cowl_in の i 線)", rn)] = (float(sub[a, c]), (i, S.jm, int(c)), lines[i, c])
    return out


# ---------------------------------------------------------------------------------------------- (iii)
def spacings(g, S, cols, n_layers):
    """列 cols・全 k の中間線から n 層の y 間隔 |Δy| (下側・上側) と、i 線に垂直な間隔 |Δy|·cos θ。
    θ は節点 (i, 層の中間線側の j) の前後の辺の角度の平均 (= 局所の i 線の向き) と、上流の辺だけの角度の 2 通り
    (後縁の列は i 線が折れているので前後の平均に意味が無い。上流の辺 = 壁の向きで測る)。
    戻り = {"lo": (dy, dn_bisector, dn_upstream) 各 (len(cols), n, nz), "up": ...}"""
    G = S.grid(g["coords"])
    jm = S.jm
    out = {}
    for side, js in (("lo", [(jm - n + 1, jm - n) for n in range(1, n_layers + 1) if jm - n >= 0]),
                     ("up", [(jm + n, jm + n - 1) for n in range(1, n_layers + 1) if jm + n < S.NJ])):
        dy = np.abs(np.array([[G[i, a, :, 1] - G[i, b, :, 1] for (a, b) in js] for i in cols]))      # (nc, n, nz)
        th_in = np.array([[_ang(G[i, a] - G[i - 1, a]) for (a, b) in js] for i in cols])
        th_out = np.array([[_ang(G[min(i + 1, S.ni - 1), a] - G[i, a]) for (a, b) in js] for i in cols])
        th_mid = th_in + 0.5 * _wrap(th_out - th_in)
        out[side] = (dy, dy * np.cos(np.radians(th_mid)), dy * np.cos(np.radians(th_in)))
    return out


# ---------------------------------------------------------------------------------------------- (iv)
def corner_jacobians(coords, hexes):
    """各ヘキサの 8 頂点のヤコビアン det[e_ξ, e_η, e_ζ] (各頂点で 3 方向の前進差分) とスケール済み値。(M, 8)"""
    X = coords[hexes]                                    # (M, 8, 3)
    J = np.zeros(hexes.shape); Js = np.zeros(hexes.shape)
    for n, (a, b, c) in enumerate(HEX_CORNERS):
        ex = X[:, _IDX[(1, b, c)]] - X[:, _IDX[(0, b, c)]]
        ey = X[:, _IDX[(a, 1, c)]] - X[:, _IDX[(a, 0, c)]]
        ez = X[:, _IDX[(a, b, 1)]] - X[:, _IDX[(a, b, 0)]]
        d = np.einsum("ij,ij->i", ex, np.cross(ey, ez))
        J[:, n] = d
        Js[:, n] = d / np.maximum(np.linalg.norm(ex, axis=1) * np.linalg.norm(ey, axis=1) * np.linalg.norm(ez, axis=1), 1e-300)
    return J, Js


# ---------------------------------------------------------------------------------------------- (v)
def identity(A, B, S, ib):
    cA, cB = A["coords"], B["coords"]
    r = {"n_nodes": (cA.shape[0], cB.shape[0])}
    same_n = cA.shape[0] == cB.shape[0]
    r["x_equal"] = bool(same_n and np.array_equal(cA[:, 0], cB[:, 0]))
    r["z_equal"] = bool(same_n and np.array_equal(cA[:, 2], cB[:, 2]))
    r["hex_equal"] = bool(A["hexes"].shape == B["hexes"].shape and np.array_equal(A["hexes"], B["hexes"]))
    if A["bfaces"] is not None and B["bfaces"] is not None:
        r["bface_equal"] = bool(set(A["bfaces"]) == set(B["bfaces"]) and all(
            np.array_equal(A["bfaces"][k][0], B["bfaces"][k][0]) and np.array_equal(A["bfaces"][k][1], B["bfaces"][k][1]) for k in A["bfaces"]))
    else:
        r["bface_equal"] = None
    if not same_n:
        return r
    dy = cB[:, 1] - cA[:, 1]
    moved = np.where(dy != 0.0)[0]
    r["moved"] = moved
    r["max_dy_m"] = float(np.max(np.abs(dy))) if moved.size else 0.0
    base = moved[moved < S.N_base]
    r["moved_nonbase"] = int(moved.size - base.size)
    i = base // (S.NJ * S.nz); j = (base // S.nz) % S.NJ
    r["moved_i"] = (int(i.min()), int(i.max())) if base.size else None
    r["moved_j"] = (int(j.min()), int(j.max())) if base.size else None
    r["moved_x_m"] = (float(cA[moved, 0].min()), float(cA[moved, 0].max())) if moved.size else None
    tags = {}
    if A["bfaces"] is not None:
        mset = np.zeros(cA.shape[0], bool); mset[moved] = True
        for pid, (sz, nd) in A["bfaces"].items():
            u = np.unique(nd)
            tags[pid] = int(mset[u].sum())
    r["moved_by_tag"] = tags
    return r


# ---------------------------------------------------------------------------------------------- --admission
PASS, FAIL, UND = "PASS", "FAIL", "判定不能"
VERDICT_WORD = {PASS: "PASS", FAIL: "FAIL", UND: "UNDECIDABLE"}
EXIT_CODE = {PASS: 0, FAIL: 1, UND: 2}
QUALITY_RANK = {"PASS": 2, "SOFT-PASS": 1, "FAIL": 0}


def combine(statuses):
    """1 つでも FAIL → FAIL、FAIL が無く判定不能が 1 つでもあれば 判定不能、全部 PASS → PASS (空は判定不能)。"""
    s = list(statuses)
    if FAIL in s:
        return FAIL
    if UND in s or not s:
        return UND
    return PASS


class Rows:
    """判定の行。row = §6.0 の表の行名、item = その中の確認項目。"""

    def __init__(self):
        self.rows = []

    def add(self, row, item, status, value="", ref="", limit="", note="", num=None):
        """num = 主な測定値 (B) の数値 (JSON・試験用)。"""
        self.rows.append({"row": row, "item": item, "status": status, "value": value, "ref_A": ref, "limit": limit, "note": note,
                          "num": None if num is None else float(num)})

    def status_of(self, row):
        return combine(r["status"] for r in self.rows if r["row"] == row)


def _kind_str(g):
    return {"h5": "h5 (" + ("node 変換" if g.get("node_h5") else "node 変換でない") + f", 座標 {g.get('coord_repr')})",
            "msh": "msh (10 桁の text、変換前)", "memory": "メモリ上 (倍精度、変換前)"}.get(g.get("source"), str(g.get("source")))


def _blend_of(info):
    return float(info.get("te_wake_blend_H", 0.0) or 0.0)


def parse_quality(path):
    """check_mesh_quality.py の出力 (prepare の MESH_QUALITY.txt) → dict(verdict, ar_max, skew_max, cells) か、読めなければ (None, 理由)。
    VERDICT 行が複数あれば `AR<=` を含む行 (成立性の FATAL のときは FATAL 直後の行) を使う。"""
    if not path:
        return None, "ファイルの指定が無い"
    if not os.path.exists(path):
        return None, f"{path} が無い"
    lines = open(path, errors="replace").read().splitlines()
    vlines = [l for l in lines if l.strip().startswith("VERDICT:")]
    cand = [l for l in vlines if "AR<=" in l]
    fatal = any(l.strip().startswith("FATAL:") for l in lines)
    if not cand and fatal:
        idx = next(i for i, l in enumerate(lines) if l.strip().startswith("FATAL:"))
        cand = [l for l in lines[idx + 1:] if l.strip().startswith("VERDICT:")][:1]
    if len(cand) != 1:
        return None, f"{path}: check_mesh_quality の VERDICT 行が特定できない ({len(cand)} 行)"
    m = re.match(r"\s*VERDICT:\s*(SOFT-PASS|PASS|FAIL)", cand[0])
    if not m:
        return None, f"{path}: VERDICT を読めない: {cand[0].strip()}"
    out = {"verdict": m.group(1), "ar_max": None, "skew_max": None, "cells": None, "line": cand[0].strip()}
    t = re.search(r"AR<=\s*([0-9.eE+\-]+),\s*skew<=\s*([0-9.eE+\-]+)", cand[0])
    if t:
        out["ar_max"], out["skew_max"] = float(t.group(1)), float(t.group(2))
    c = [re.match(r"\s*cells:\s*(\d+)", l) for l in lines]
    c = [int(x.group(1)) for x in c if x]
    out["cells"] = c[0] if len(c) == 1 else None
    return out, ""


def parse_dual(path):
    """check_dual_closure.py の出力 → dict(verdict, tol, n_cv) か (None, 理由)。VERDICT は「閉性」行より後の最初の行。"""
    if not path:
        return None, "ファイルの指定が無い (--dual-closure-b)"
    if not os.path.exists(path):
        return None, f"{path} が無い"
    lines = open(path, errors="replace").read().splitlines()
    out = {"verdict": None, "tol": None, "n_cv": None, "undecidable": False}
    i_clo = next((i for i, l in enumerate(lines) if "閉性" in l and "|ΣS|" in l), None)
    v = [l for l in (lines[i_clo + 1:] if i_clo is not None else lines) if l.strip().startswith("VERDICT:")]
    if not v:
        return None, f"{path}: VERDICT 行が無い"
    m = re.match(r"\s*VERDICT:\s*(PASS|FAIL)", v[0])
    if not m:
        return None, f"{path}: VERDICT を読めない: {v[0].strip()}"
    out["verdict"] = m.group(1)
    out["undecidable"] = "判定不能" in v[0]
    if i_clo is not None:
        t = re.search(r"\(>\s*([0-9.eE+\-]+):", lines[i_clo])
        out["tol"] = float(t.group(1)) if t else None
    c = next((re.match(r"\s*CV\s+(\d+)\s*/", l) for l in lines if re.match(r"\s*CV\s+\d+\s*/", l)), None)
    out["n_cv"] = int(c.group(1)) if c else None
    return out, ""


def _cmq():
    """check_mesh_quality.py (skew・AR の定義をそのまま使う)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import check_mesh_quality
    return check_mesh_quality


def hex_quality(coords, hexes, sgn, chunk=400000):
    """全ヘキサの (非正・非有限の頂点ヤコビアンの数, スケール済みの最小, AR (M,), skew (M,))。
    ヤコビアン = corner_jacobians (各頂点で 3 方向の前進差分)、向きは sgn (A の多数派)。AR・skew = check_mesh_quality.metrics_vectorized (3D)。"""
    cmq = _cmq()
    M = hexes.shape[0]
    nbad = 0; jmin = np.inf
    ars = np.empty(M); skews = np.empty(M)
    for s in range(0, M, chunk):
        h = hexes[s:s + chunk]
        J, Js = corner_jacobians(coords, h)
        bad = ~np.isfinite(J) | ~np.isfinite(Js) | (sgn * J <= 0.0)
        nbad += int(bad.sum())
        good = (sgn * Js)[np.isfinite(Js)]
        if good.size:
            jmin = min(jmin, float(good.min()))
        conn = h.reshape(-1); offs = 8 * np.arange(1, h.shape[0] + 1); types = np.full(h.shape[0], 12)
        a, k = cmq.metrics_vectorized(coords, conn, offs, types, True)
        ars[s:s + chunk] = a; skews[s:s + chunk] = k
    return nbad, jmin, ars, skews


def _fmt_loc(G, i, j, k):
    p = G[i, j, k]
    return f"(i,j,k)=({i},{j},{k}) x,y,z=({p[0]:.6f},{p[1]:.6f},{p[2]:.6f}) m"


def _argmax_loc(arr, ii, jj):
    """arr (len(ii), len(jj), nz) の最大とその (i, j, k)。非有限が 1 つでもあれば inf (= 閾値を満たさない) とその位置。"""
    arr = np.asarray(arr, dtype=float)
    bad = ~np.isfinite(arr)
    if bad.any():
        a, b, c = (int(x) for x in np.argwhere(bad)[0])
        return float("inf"), (int(ii[a]), int(jj[b]), c)
    a, b, c = np.unravel_index(int(np.argmax(arr)), arr.shape)
    return float(arr[a, b, c]), (int(ii[a]), int(jj[b]), int(c))


def admission(A, B, quality_a=None, quality_b=None, dual_b=None, dual_a=None, ar_max=None, expect_blend_H=ADM_BLEND_H):
    """plan §6.0 の投入条件の表を判定する。A・B は load() / from_mesher() の戻り。戻り = (総合, Rows)。"""
    R = Rows()
    n = ADM_LAYERS
    # ---- 前提
    for lab, g in (("A", A), ("B", B)):
        ok = g.get("source") == "h5" and bool(g.get("node_h5"))
        R.add("前提", f"{lab} の入力が最終 node 格子", PASS if ok else UND, value=_kind_str(g), limit="node 変換の h5",
              note="" if ok else "§6.0 は最終 node 格子の検査 (変換前の格子では投入条件を判定できない)")
    LbA, LbB = _blend_of(A["info"]), _blend_of(B["info"])
    R.add("前提", "L_b (te_wake_blend_H)", PASS if (LbA == 0.0 and LbB == float(expect_blend_H)) else FAIL,
          value=f"A {LbA:g} H / B {LbB:g} H", limit=f"A 0 / B {expect_blend_H:g} H (事前登録、結果を見て変えない)")
    try:
        SA, SB = Struct(A["info"]), Struct(B["info"])
    except SystemExit as e:
        R.add("前提", "格子の構造 (info)", UND, note=str(e))
        return combine(r["status"] for r in R.rows), R
    keys = ("ni", "NJ", "nz", "jm", "i_te", "k_sw", "i_sw", "N_base")
    diff = {k: (getattr(SA, k), getattr(SB, k)) for k in keys if getattr(SA, k) != getattr(SB, k)}
    if SA.dup1 != SB.dup1:
        diff["カウル上コピー"] = ("", "")
    nA, nB = A["coords"].shape[0], B["coords"].shape[0]
    R.add("前提", "格子の構造 (ni, NJ, nz, jm, i_te, k_sw, i_sw, カウル上コピー) が同じ", PASS if not diff else FAIL,
          value="同じ" if not diff else f"違う {diff}")
    if diff or nA != nB:
        R.add("固定すべき量", "節点数", PASS if nA == nB else FAIL, value=f"A {nA} / B {nB}", limit="同じ")
        for row in ("最初の下流の辺", "後縁の折れ", "復帰区間の折れ", "層の Δy", "primal の品質", "品質・双対の判定"):
            R.add(row, "(構造・節点数が違うので比べられない)", UND)
        return combine(r["status"] for r in R.rows), R
    S = SA
    if S.jm < n or S.NJ - 1 - S.jm < n or S.i_te < 1 or S.i_te + 2 > S.ni - 1:
        R.add("前提", f"上下に {n} 層・後縁の前後に辺がある", UND, value=f"jm {S.jm}, NJ {S.NJ}, i_te {S.i_te}, ni {S.ni}")
        return combine(r["status"] for r in R.rows), R
    GA, GB = S.grid(A["coords"]), S.grid(B["coords"])
    it, jm = S.i_te, S.jm
    ib = B["info"]
    idn = identity(A, B, S, ib)
    moved = idn.get("moved", np.zeros(0, np.int64))
    mv_base = moved[moved < S.N_base]
    mv_i = mv_base // (S.NJ * S.nz)
    i_last = max(int(ib.get("te_wake_i_last", it)) if LbB > 0.0 else it, int(mv_i.max()) if mv_i.size else it)
    i_ret = i_last + 1
    lo = np.arange(jm - n, jm)            # 下側の第 n..第 1 層
    up = np.arange(jm + 1, jm + 1 + n)    # 上側の第 1..第 n 層

    def ang_d(G):
        return _ang(G[it + 1] - G[it])          # (NJ, nz) 最初の下流の辺

    def ang_u(G):
        return _ang(G[it] - G[it - 1])          # (NJ, nz) 同じ i 線の最後の上流の辺 (j = jm は下側壁面 / 板の外は中間線)

    def ang_u_upper(g, G):
        return np.array([_ang(G[it, jm, k] - g["coords"][S.dup1.get((it - 1, k), S.base(it - 1, jm, k))]) for k in range(S.k_sw + 1)])

    # ---- 最初の下流の辺
    row = "最初の下流の辺"
    m0_info = ib.get("te_wake_slope_te_deg")
    m0_meas = float(_ang(GB[it, jm, 0] - GB[it - 1, jm, 0]))       # 下側壁面 (cowl_out) の最後の上流の辺、k = 0
    if m0_info is None:
        m0, src = m0_meas, "B の座標 (下側壁面の最後の上流の辺、k=0。info に te_wake_slope_te_deg が無い)"
    else:
        m0, src = float(m0_info), "B の info te_wake_slope_te_deg"
    if m0_info is not None and abs(_wrap(m0_meas - float(m0_info))) > 0.01:
        R.add(row, "m0 (info) と B の下側壁面の辺の一致", UND, value=f"info {float(m0_info):.4f}° / 座標 {m0_meas:.4f}°", limit="差 ≤ 0.01°",
              note="info が別の格子のもの")
    dB = np.abs(_wrap(ang_d(GB)[lo] - m0)); dA = np.abs(_wrap(ang_d(GA)[lo] - m0))
    vB, lB = _argmax_loc(dB[None], [it], lo)
    vA, _ = _argmax_loc(dA[None], [it], lo)
    dte = float(np.max(np.abs(_wrap(ang_d(GB)[jm] - m0))))
    R.add(row, f"下側の第 1〜第 {n} 層・全 span: max |θ(i_te→i_te+1) − m0|", PASS if vB <= ADM_FIRST_EDGE_DEG else FAIL,
          value=f"{vB:.4f}° @ {_fmt_loc(GB, *lB)}  [B の辺 {float(ang_d(GB)[lo].min()):.4f}〜{float(ang_d(GB)[lo].max()):.4f}°]",
          ref=f"{vA:.4f}°", limit=f"≤ {ADM_FIRST_EDGE_DEG:g}° (m0 = {m0:.4f}°: {src})",
          note=f"参考 (判定外): 後縁節点 j=jm の最初の下流の辺と m0 の差 {dte:.4f}°", num=vB)

    # ---- 後縁の折れ
    row = "後縁の折れ"

    def te_folds(g, G):
        d, u = ang_d(G), ang_u(G)
        f_lo = np.abs(_wrap(d[lo] - u[lo]))
        f_up = np.abs(_wrap(d[up] - u[up]))
        f_te_lo = np.abs(_wrap(d[jm] - u[jm]))                                       # (nz)
        f_te_up = np.abs(_wrap(d[jm, :S.k_sw + 1] - ang_u_upper(g, G)))              # (k_sw+1)
        return f_lo, f_up, f_te_lo, f_te_up
    fB, fA = te_folds(B, GB), te_folds(A, GA)
    labels = (f"下側の第 1〜第 {n} 層・全 span", f"上側の第 1〜第 {n} 層・全 span",
              "共有の後縁節点・全 span (上流 = 下側壁面 cowl_out、板の外 z>W/2 は中間線)", "共有の後縁節点・z≤W/2 (上流 = 上側壁面 cowl_in)")
    for idx, lab in enumerate(labels):
        b, a = fB[idx], fA[idx]
        if idx < 2:
            jj = lo if idx == 0 else up
            v, loc = _argmax_loc(b[None], [it], jj); va = float(a.max())
            where = _fmt_loc(GB, *loc)
        else:
            k = int(np.argmax(b)); v = float(b[k]); va = float(a.max())
            where = _fmt_loc(GB, it, jm, k)
        R.add(row, lab + ": max 折れ角", PASS if v <= ADM_TE_FOLD_DEG else FAIL, value=f"{v:.4f}° @ {where}", ref=f"{va:.4f}°",
              limit=f"≤ {ADM_TE_FOLD_DEG:g}°", num=v)

    # ---- 復帰区間の折れ
    row = "復帰区間の折れ"
    if LbB <= 0.0 and not mv_i.size:
        R.add(row, "復帰区間", UND, note="B に変形区間が無い (te_wake_i_last も動いた節点も無い)")
    elif i_ret > S.ni - 2:
        R.add(row, "復帰区間", UND, note=f"最初の未変形の辺が無い (i_ret {i_ret}、ni {S.ni})")
    else:
        def turns(G):
            th = _ang(np.diff(G, axis=0))                         # (ni-1, NJ, nz)
            return np.abs(_wrap(np.diff(th, axis=0)))             # (ni-2, NJ, nz): 節点 i = 1..ni-2
        ii = np.arange(it + 1, i_ret + 1)
        tB = turns(GB)[ii - 1]; tA = turns(GA)[ii - 1]
        v, loc = _argmax_loc(tB, ii, np.arange(S.NJ))
        nonbase = int(moved.size - mv_base.size)
        st = PASS if v <= ADM_RETURN_FOLD_DEG else FAIL
        if nonbase:
            st = combine([st, UND])
        R.add(row, f"主ブロックの全 i 線・全 span、節点 i = {it + 1}..{i_ret} (最初の未変形の辺を含む): max 折れ角", st,
              value=f"{v:.4f}° @ {_fmt_loc(GB, *loc)}", ref=f"{float(tA.max()):.4f}°", limit=f"≤ {ADM_RETURN_FOLD_DEG:g}°",
              note=(f"主ブロック外で動いた節点 {nonbase} の i 線は測っていない" if nonbase else
                    f"復帰区間 = info の te_wake_i_last {ib.get('te_wake_i_last')} と動いた節点の i の最大の大きい方 + 1"), num=v)

    # ---- 固定すべき量
    row = "固定すべき量"
    R.add(row, "節点数", PASS, value=f"A {nA} / B {nB}", limit="同じ")
    R.add(row, "x 座標 (全節点)", PASS if idn["x_equal"] else FAIL, value="ビット一致" if idn["x_equal"] else
          f"相違 {int(np.sum(A['coords'][:, 0] != B['coords'][:, 0]))} 節点", limit="ビット一致")
    R.add(row, "z 座標 (全節点)", PASS if idn["z_equal"] else FAIL, value="ビット一致" if idn["z_equal"] else
          f"相違 {int(np.sum(A['coords'][:, 2] != B['coords'][:, 2]))} 節点", limit="ビット一致")
    R.add(row, "ヘキサの接続", PASS if idn["hex_equal"] else FAIL, value="同一" if idn["hex_equal"] else "相違", limit="同一")
    if idn["bface_equal"] is None:
        R.add(row, "境界面の接続 (タグ別)", UND, value="境界面の接続が無い入力", limit="同一")
    else:
        R.add(row, "境界面の接続 (タグ別)", PASS if idn["bface_equal"] else FAIL, value="同一" if idn["bface_equal"] else "相違", limit="同一")
    x_te = float(GA[it, jm, 0, 0])
    ups = A["coords"][:, 0] <= x_te
    nchg = int(np.sum(np.any(A["coords"][ups] != B["coords"][ups], axis=1)))
    R.add(row, f"上流の座標 (x ≤ 後縁の x = {x_te:.6f} m の {int(ups.sum())} 節点、3 成分)", PASS if nchg == 0 else FAIL,
          value=f"相違 {nchg} 節点", limit="0")
    if A["bfaces"] is None or B["bfaces"] is None:
        R.add(row, "固体形状・外部境界の形状", UND, value="境界面の接続が無い入力")
    else:
        mt = idn["moved_by_tag"]
        walls = {PHYS_NAME.get(p, p): c for p, c in mt.items() if c and p in WALL_TAGS}
        ext = {PHYS_NAME.get(p, p): c for p, c in mt.items() if c and p not in WALL_TAGS and p not in IN_PLANE_TAGS}
        inpl = {PHYS_NAME.get(p, p): mt.get(p, 0) for p in IN_PLANE_TAGS}
        R.add(row, "固体形状 (壁タグの節点が動かない)", PASS if not walls else FAIL, value=f"動いた {walls or 0}", limit="0")
        R.add(row, "外部境界の形状 (sym・side_far 以外の境界タグの節点が動かない)", PASS if not ext else FAIL,
              value=f"動いた {ext or 0}", limit="0",
              note=f"記録: sym・side_far 上の面内 (y) の移動 {inpl} (z は全節点ビット一致なので面内)、動いた節点の総数 {moved.size}"
                   f" (主ブロック外 {int(moved.size - mv_base.size)})")

    # ---- 層の Δy
    row = "層の Δy"
    rep = {A.get("coord_repr"), B.get("coord_repr")}
    if "msh10g" in rep or len(rep) != 1:
        R.add(row, "第 1 層・第 1〜第 10 層", UND, value=f"座標 A {A.get('coord_repr')} / B {B.get('coord_repr')}",
              note="丸めの上限が事前登録にあるのは倍精度 (変換前) と float32 の h5 (変換後) だけ。A・B は同じ表現で比べる")
    else:
        rnd = ADM_ROUND_ULPS * EPS32 if rep == {"float32"} else 0.0
        cols = np.arange(it, min(i_ret, S.ni - 1) + 1)
        for side, lab in (("lo", "下側"), ("up", "上側")):
            if side == "lo":
                pa, pb = jm - np.arange(0, n), jm - np.arange(1, n + 1)
            else:
                pa, pb = jm + np.arange(1, n + 1), jm + np.arange(0, n)
            yAa, yAb = GA[np.ix_(cols, pa)][..., 1], GA[np.ix_(cols, pb)][..., 1]        # (nc, n, nz)
            yBa, yBb = GB[np.ix_(cols, pa)][..., 1], GB[np.ix_(cols, pb)][..., 1]
            dyA, dyB = np.abs(yAa - yAb), np.abs(yBa - yBb)
            ymax = np.max(np.abs(np.stack([yAa, yAb, yBa, yBb])), axis=0)
            r_ = rnd * ymax
            dd = np.abs(dyB - dyA)
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = dd / dyA
            use1 = dd[:, 0] / (ADM_DY_FIRST_REL * dyA[:, 0] + r_[:, 0])                 # ≤ 1 で合格
            usen = dd / (ADM_DY_LAYERS_REL * dyA + r_)
            ok1 = bool(np.all(dd[:, 0] <= ADM_DY_FIRST_REL * dyA[:, 0] + r_[:, 0]) and np.all(dyA[:, 0] > 0))
            okn = bool(np.all(dd <= ADM_DY_LAYERS_REL * dyA + r_) and np.all(dyA > 0))
            lim1 = f"|Δy_B − Δy_A| ≤ {ADM_DY_FIRST_REL:g}·Δy_A" + (f" + {ADM_ROUND_ULPS:g}·ε₃₂·max|y|" if rnd else "")
            limn = f"|Δy_B − Δy_A| ≤ {ADM_DY_LAYERS_REL:g}·Δy_A" + (f" + {ADM_ROUND_ULPS:g}·ε₃₂·max|y|" if rnd else "")
            R.add(row, f"{lab}の第 1 層 (列 i = {it}..{cols[-1]}・全 span)", PASS if ok1 else FAIL,
                  value=f"max 相対差 {float(np.nanmax(rel[:, 0])):.3e}、許容の使用率 max {float(np.nanmax(use1)):.3g}", limit=lim1,
                  num=float(np.nanmax(rel[:, 0])))
            R.add(row, f"{lab}の第 1〜第 {n} 層の各層厚 (同)", PASS if okn else FAIL,
                  value=f"max 相対差 {float(np.nanmax(rel)):.3e}、許容の使用率 max {float(np.nanmax(usen)):.3g}", limit=limn,
                  num=float(np.nanmax(rel)))

    # ---- primal の品質
    row = "primal の品質"
    JA_s, _ = corner_jacobians(A["coords"], A["hexes"][:min(200000, A["hexes"].shape[0])])
    sgn = 1.0 if np.median(JA_s) > 0 else -1.0
    nbA, jminA, arsA, skA = hex_quality(A["coords"], A["hexes"], sgn)
    nbB, jminB, arsB, skB = hex_quality(B["coords"], B["hexes"], sgn)
    R.add(row, "非正・非有限の頂点ヤコビアン (B、全ヘキサの 8 頂点)", PASS if nbB == 0 else FAIL, value=f"{nbB}", ref=f"{nbA}", limit="0",
          note=f"向き = A の多数派 ({'正' if sgn > 0 else '負'})", num=nbB)
    R.add(row, "スケール済みヤコビアンの最小 (B)", PASS if (np.isfinite(jminB) and jminB >= ADM_SCALED_JAC_MIN) else FAIL,
          value=f"{jminB:.5f}", ref=f"{jminA:.5f}", limit=f"≥ {ADM_SCALED_JAC_MIN:g}", num=jminB)
    skA_max, skB_max = float(np.nanmax(skA)), float(np.nanmax(skB))
    sk_ok = np.all(np.isfinite(skB))
    R.add(row, "skew 最大 (B、check_mesh_quality の 3D 定義)", PASS if (sk_ok and skB_max <= ADM_SKEW_MAX) else FAIL,
          value=f"{skB_max:.5f}", ref=f"{skA_max:.5f}", limit=f"≤ {ADM_SKEW_MAX:g}", num=skB_max)
    R.add(row, "skew 最大の A からの増加", PASS if (sk_ok and skB_max - skA_max <= ADM_SKEW_INC_MAX) else FAIL,
          value=f"{skB_max - skA_max:+.5f}", limit=f"≤ {ADM_SKEW_INC_MAX:g}", num=skB_max - skA_max)

    # ---- 品質・双対の判定
    row = "品質・双対の判定"
    qa, why_a = parse_quality(quality_a)
    qb, why_b = parse_quality(quality_b)
    MA, MB = A["hexes"].shape[0], B["hexes"].shape[0]
    if qa is None or qb is None:
        R.add(row, "(a) check_mesh_quality の VERDICT が A と同等以上", UND, note="; ".join(w for w in (f"A: {why_a}" if qa is None else "",
                                                                                         f"B: {why_b}" if qb is None else "") if w))
    elif qa["ar_max"] is None or qb["ar_max"] is None or (qa["ar_max"], qa["skew_max"]) != (qb["ar_max"], qb["skew_max"]):
        R.add(row, "(a) check_mesh_quality の VERDICT が A と同等以上", UND, value=f"A「{qa['line']}」/ B「{qb['line']}」",
              note="同じ閾値 (AR・skew) で比べられない")
    elif qa["cells"] != MA or qb["cells"] != MB:
        R.add(row, "(a) check_mesh_quality の VERDICT が A と同等以上", UND, value=f"cells A {qa['cells']} / B {qb['cells']}",
              note=f"MESH_QUALITY のセル数が格子のヘキサ数 (A {MA} / B {MB}) と合わない (別の格子の出力)")
    else:
        ok = QUALITY_RANK[qb["verdict"]] >= QUALITY_RANK[qa["verdict"]]
        R.add(row, "(a) check_mesh_quality の VERDICT が A と同等以上", PASS if ok else FAIL, value=f"B {qb['verdict']}", ref=qa["verdict"],
              limit=f"B ≥ A (PASS > SOFT-PASS > FAIL)、閾値 AR<= {qb['ar_max']:g}・skew<= {qb['skew_max']:g}")
    thr, thr_src = None, ""
    if qb is not None and qb["ar_max"] is not None:
        thr, thr_src = qb["ar_max"], "B の MESH_QUALITY"
        if ar_max is not None and float(ar_max) != thr:
            thr = None; thr_src = f"--ar-max {ar_max:g} と MESH_QUALITY の {qb['ar_max']:g} が違う"
    elif ar_max is not None:
        thr, thr_src = float(ar_max), "--ar-max"
    touch = np.zeros(MB, bool)
    if moved.size:
        mset = np.zeros(nB, bool); mset[moved] = True
        touch = mset[B["hexes"]].any(axis=1)
    if thr is None:
        R.add(row, "(b) 変形領域に新しい AR 閾値超過を作らない", UND, note=f"AR 閾値が無い ({thr_src or 'MESH_QUALITY も --ar-max も無い'})")
    else:
        new = touch & (arsB > thr) & ~(arsA > thr)
        R.add(row, "(b) 変形領域に新しい AR 閾値超過を作らない", PASS if int(new.sum()) == 0 else FAIL,
              value=f"新しい超過 {int(new.sum())} (変形領域 {int(touch.sum())} ヘキサ、B の超過 {int((touch & (arsB > thr)).sum())}・"
                    f"A の超過 {int((touch & (arsA > thr)).sum())}、AR 最大 B {float(arsB[touch].max()) if touch.any() else 0.0:.1f}"
                    f" / A {float(arsA[touch].max()) if touch.any() else 0.0:.1f})",
              limit=f"0 (AR > {thr:g}、閾値は {thr_src})", note="AR は check_mesh_quality と同じ定義 (全辺の最長/最短) をこの道具で計算",
              num=int(new.sum()))
    db, why = parse_dual(dual_b)
    if db is None:
        R.add(row, "(c) 最終 node 入力で check_dual_closure.py --tol 1e-5 が PASS", UND, note=why)
    elif db["undecidable"]:
        R.add(row, "(c) 最終 node 入力で check_dual_closure.py --tol 1e-5 が PASS", UND, value="check_dual_closure が判定不能", note=dual_b)
    elif db["tol"] is None or abs(db["tol"] - ADM_DUAL_TOL) > 1e-12 * ADM_DUAL_TOL + 1e-20:
        R.add(row, "(c) 最終 node 入力で check_dual_closure.py --tol 1e-5 が PASS", UND, value=f"tol {db['tol']}", note="--tol 1e-5 の出力でない")
    elif db["n_cv"] != nB:
        R.add(row, "(c) 最終 node 入力で check_dual_closure.py --tol 1e-5 が PASS", UND, value=f"CV {db['n_cv']}",
              note=f"CV 数が B の節点数 {nB} と合わない (別の格子の出力)")
    else:
        R.add(row, "(c) 最終 node 入力で check_dual_closure.py --tol 1e-5 が PASS", PASS if db["verdict"] == "PASS" else FAIL,
              value=f"VERDICT {db['verdict']} (CV {db['n_cv']}, tol {db['tol']:g})", limit="PASS")
    if dual_a:
        da, why = parse_dual(dual_a)
        R.rows[-1]["ref_A"] = (f"{da['verdict']} (tol {da['tol']})" if da else why)
    return combine(r["status"] for r in R.rows), R


def print_admission(verdict, R, A_label, B_label, out=print):
    out(f"A: {A_label}")
    out(f"B: {B_label}")
    out("plan convection-zero-thickness-edge-reconstruction §6.0「投入条件 = 最終 node 格子の検査」(事前登録 2026-10-08)")
    out(f"{'判定':6s} | {'行':14s} | 項目 | 測定 (B) | 参考 (A) | 閾値")
    for r in R.rows:
        out(f"{r['status']:6s} | {r['row']:14s} | {r['item']} | {r['value']} | {r['ref_A'] or '-'} | {r['limit'] or '-'}"
            + (f"\n{'':6s}   └ {r['note']}" if r['note'] else ""))
    out("")
    out("行ごと: " + ", ".join(f"{row} {R.status_of(row)}" for row in dict.fromkeys(r["row"] for r in R.rows)))
    out(f"ADMISSION VERDICT: {VERDICT_WORD[verdict]}")


# ---------------------------------------------------------------------------------------------- 来歴 (--provenance-out・--same-mesh)
def _sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _mark():
    """mark_zero_thickness_edges (格子署名 mesh_signature の正本)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import mark_zero_thickness_edges
    return mark_zero_thickness_edges


def signature_of(path):
    """h5 (か run ディレクトリの格子) の格子署名を再計算する。戻り = (署名, 版, 格子のパス)。読めなければ例外。"""
    grid = grid_of(path)
    if not grid.endswith(".h5"):
        raise ValueError(f"{grid}: h5 でない (格子署名は最終 node の h5 から作る)")
    if not os.path.exists(grid):
        raise FileNotFoundError(grid)
    m = _mark()
    return m.mesh_signature(grid), m.SIG_VERSION, grid


def _stage_manifest():
    """stage_manifest (双対幾何のハッシュ dual_geometry_hash の正本)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import stage_manifest
    return stage_manifest


def dual_of(grid):
    """h5 の双対幾何のハッシュを再計算する (plan tooling-sern-te-wake-grid §4)。戻り = (ハッシュ, 版)。読めなければ例外。"""
    sm = _stage_manifest()
    return sm.dual_geometry_hash(grid), sm.DUAL_HASH_VERSION


def _mesher_curve_version():
    try:
        if str(DESIGN) not in sys.path:
            sys.path.insert(0, str(DESIGN))
        from forge_design.meshing import mesh_sern3d
        return getattr(mesh_sern3d, "TE_WAKE_CURVE_VERSION", None)
    except Exception as e:  # noqa: BLE001
        return f"(取得不能: {type(e).__name__})"


def _git_state():
    def run(*a):
        r = subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True)
        return r.stdout.rstrip() if r.returncode == 0 else None
    head = run("rev-parse", "HEAD")
    dirty = run("status", "--porcelain", "--untracked-files=no", "--", "case/46.sern_design/diag", "design/forge_design/meshing",
                "design/forge_design/evaluate", "solver_density_cuda/tools/mark_zero_thickness_edges.py")
    return {"head": head, "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty_paths": (dirty.splitlines() if dirty else []) if dirty is not None else None,
            "note": "この道具を回した作業ツリーの版 (格子を生成した版とは限らない)"}


def provenance_entry(path, info_path=None, code_commit=None, attach=()):
    """1 つの run (か格子) の来歴。値を転記せず、格子署名は実入力格子から再計算する。"""
    e = {"input": str(path), "run_dir": str(path) if os.path.isdir(path) else None}
    grid = grid_of(path)
    e["grid"] = grid
    e["grid_exists"] = os.path.exists(grid)
    try:
        sig, ver, _ = signature_of(path)
        e["mesh_signature"], e["mesh_signature_version"] = sig, ver
        try:                    # 署名に入らない双対幾何 (変換器の違いを検出する。plan tooling-sern-te-wake-grid §4)
            e["dual_hash"], e["dual_hash_version"] = dual_of(grid)
        except Exception as ex:  # noqa: BLE001
            e["dual_hash"], e["dual_hash_error"] = None, f"{type(ex).__name__}: {ex}"
        meta = _h5_meta(grid)
        e["grid_discretization"] = "node" if meta["node"] else "cell"
        e["grid_coord_dtype"] = meta["coord_dtype"]
    except Exception as ex:  # noqa: BLE001
        e["mesh_signature"] = None
        e.setdefault("dual_hash", None)
        e["mesh_signature_error"] = f"{type(ex).__name__}: {ex}"
    if os.path.isdir(path):
        info_path = info_path or os.path.join(path, "prepare_info.json")
    raw = None
    if info_path and os.path.exists(info_path):
        try:
            raw = json.load(open(info_path))
        except Exception as ex:  # noqa: BLE001
            e["info_error"] = f"{type(ex).__name__}: {ex}"
    e["info_path"] = info_path
    info = (raw.get("mesh", raw) if isinstance(raw, dict) else {}) or {}
    if raw is None:
        e["te_wake_blend_H"] = None
        e["te_wake_blend_H_source"] = "info が無い・読めない"
    elif "te_wake_blend_H" in info:
        e["te_wake_blend_H"] = float(info["te_wake_blend_H"] or 0.0)
        e["te_wake_blend_H_source"] = "prepare_info の mesh"
    else:
        e["te_wake_blend_H"] = 0.0
        e["te_wake_blend_H_source"] = "info にキー無し (オプション導入 [2d280bb1] 前の格子 = 変形なし)"
    if e["te_wake_blend_H"]:
        e["te_wake_curve_version"] = info.get("te_wake_curve_version")
        e["te_wake_curve_version_source"] = "prepare_info の mesh" if "te_wake_curve_version" in info else \
            "info に無い (曲線版の定数の導入前に作った格子。式は当時の mesh_sern3d._te_wake_midline)"
    else:
        e["te_wake_curve_version"] = None
        e["te_wake_curve_version_source"] = "変形なし"
    e["te_wake_info"] = {k: v for k, v in info.items() if str(k).startswith("te_wake_")}
    cc = next(((k, raw[k]) for k in ("code_commit", "git_commit", "commit") if isinstance(raw, dict) and raw.get(k)), None)
    if cc:
        e["code_commit"], e["code_commit_source"] = cc[1], f"prepare_info の {cc[0]}"
    elif code_commit:
        e["code_commit"], e["code_commit_source"] = code_commit, "引数の申告 (未検証。prepare_info に記録が無い)"
    else:
        e["code_commit"], e["code_commit_source"] = None, "記録なし (prepare_info に無く、申告も無い)"
    cfg = {}
    if e["run_dir"]:
        for f in sorted(os.listdir(path)):
            fp = os.path.join(path, f)
            if os.path.isfile(fp) and (f.endswith((".yaml", ".yml")) or f == "prepare_info.json"):
                cfg[f] = _sha256_file(fp)
    e["config_sha256"] = cfg
    prob = raw.get("problem") if isinstance(raw, dict) else None
    if prob:
        cands = [Path(prob)] if os.path.isabs(prob) else [Path.cwd() / prob, ROOT / prob, ROOT / "design" / prob]
        hit = next((c for c in cands if c.exists()), None)
        e["problem_yaml"] = {"recorded": prob, "resolved": str(hit) if hit else None, "sha256": _sha256_file(hit) if hit else None,
                             "note": "" if hit else "見つからない (この作業ツリーに無い)"}
    e["attached"] = [{"path": str(p), "sha256": _sha256_file(p) if os.path.exists(p) else None,
                      "bytes": os.path.getsize(p) if os.path.exists(p) else None} for p in attach]
    return e


def write_provenance(out_path, a_path, b_path, info_a=None, info_b=None, commit_a=None, commit_b=None, attach_a=(), attach_b=(),
                     admission_result=None):
    pa = provenance_entry(a_path, info_a, commit_a, attach_a)
    pb = provenance_entry(b_path, info_b, commit_b, attach_b)
    sa, sb = pa.get("mesh_signature"), pb.get("mesh_signature")
    doc = {"tool": "case/46.sern_design/diag/te_wake_grid_check.py --provenance-out",
           "plan": "plans/active/convection-zero-thickness-edge-reconstruction.md §4.2「来歴」・§6.0",
           "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
           "design_db": "取り込まない (診断の A/B。plan §4.2)",
           "checker_git": _git_state(), "mesher_curve_version_now": _mesher_curve_version(),
           "runs": {"A": pa, "B": pb},
           "A_B_same_mesh_signature": (sa == sb) if (sa and sb) else None,
           "A_B_same_dual_hash": (pa.get("dual_hash") == pb.get("dual_hash")) if (pa.get("dual_hash") and pb.get("dual_hash")) else None}
    if admission_result is not None:
        v, R = admission_result
        doc["admission"] = {"verdict": VERDICT_WORD[v], "rows": R.rows}
    with open(out_path, "w") as f:
        json.dump(doc, f, indent=1, ensure_ascii=False)
    return doc


def same_mesh(a, b, out=print):
    try:
        sa, ver, ga = signature_of(a)
        sb, _, gb = signature_of(b)
        da, dver = dual_of(ga)
        db, _ = dual_of(gb)
    except Exception as e:  # noqa: BLE001
        out(f"格子署名・双対幾何のハッシュを作れない: {type(e).__name__}: {e}")
        out("SAME MESH: UNDECIDABLE")
        return EXIT_CODE[UND]
    out(f"A {ga}: 署名 {sa}  双対 {da}")
    out(f"B {gb}: 署名 {sb}  双対 {db}")
    out(f"(署名の版 {ver}: 離散化・次元・節点順の座標・内部面の接続・境界面の接続。双対の版 {dver}: 面積ベクトル・面積・面重心・"
        f"CV 体積・CV 重心。/VALUE・/AUX は含まない)")
    if sa == sb and da != db:
        out("署名は同じで双対幾何だけが違う (変換器の違いなど)")
    same = sa == sb and da == db
    out("SAME MESH: " + ("YES" if same else "NO"))
    return 0 if same else 1


# ---------------------------------------------------------------------------------------------- 本体
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("A", nargs="?"); ap.add_argument("B", nargs="?")
    ap.add_argument("--info-a"); ap.add_argument("--info-b")
    ap.add_argument("--n-layers", type=int, default=10, help="中間線から数える層の数 ((i)・(iii))。--admission は表の 10 層で固定")
    ap.add_argument("--admission", action="store_true", help="plan §6.0 の投入条件の表を判定する (ADMISSION VERDICT)")
    ap.add_argument("--quality-a", help="A の check_mesh_quality の出力 (既定: A が run ディレクトリなら <A>/MESH_QUALITY.txt)")
    ap.add_argument("--quality-b", help="B の check_mesh_quality の出力 (既定: <B>/MESH_QUALITY.txt)")
    ap.add_argument("--dual-closure-b", help="B の最終 node 入力に対する check_dual_closure.py --tol 1e-5 の出力 (無ければ判定不能)")
    ap.add_argument("--dual-closure-a", help="A の同じ出力 (参考として並べるだけ)")
    ap.add_argument("--ar-max", type=float, help="変形領域の AR 閾値 (MESH_QUALITY.txt に AR<= が無いときだけ。あれば一致を要求)")
    ap.add_argument("--expect-blend-H", type=float, default=ADM_BLEND_H, help="B の L_b の事前登録値 (既定 1.0。plan §4.2)")
    ap.add_argument("--provenance-out", help="来歴 (格子署名・曲線版・L_b・commit・設定の sha256) を JSON で書く")
    ap.add_argument("--code-commit-a", help="A を生成したコードの commit (prepare_info に記録が無いときの申告、未検証として記録)")
    ap.add_argument("--code-commit-b", help="B の同上")
    ap.add_argument("--attach-a", action="append", default=[], help="A の来歴に sha256 を残すファイル (繰り返し可)")
    ap.add_argument("--attach-b", action="append", default=[], help="B の来歴に sha256 を残すファイル (初期場の転送ログ・動いた節点の一覧など)")
    ap.add_argument("--same-mesh", nargs=2, metavar=("A_H5", "B_H5"), help="2 つの h5 (か run) の格子署名を比べる")
    a = ap.parse_args(argv)
    if a.same_mesh:
        return same_mesh(*a.same_mesh)
    if not (a.A and a.B):
        ap.error("A と B (か --same-mesh) を指定する")
    if not a.admission and not a.provenance_out:
        return report(load(a.A, a.info_a), load(a.B, a.info_b), a.n_layers)
    rc, result = 0, None
    if a.admission:
        def qdef(p, q):
            return q if q else (os.path.join(p, "MESH_QUALITY.txt") if os.path.isdir(p) else None)
        try:
            A, B = load(a.A, a.info_a), load(a.B, a.info_b)
        except SystemExit as e:         # 入力が読めない = 判定不能 (合格扱いにしない)
            R = Rows(); R.add("前提", "A・B の格子と info を読む", UND, note=str(e))
            result = (UND, R)
            print_admission(UND, R, a.A, a.B)
        else:
            result = admission(A, B, qdef(a.A, a.quality_a), qdef(a.B, a.quality_b), a.dual_closure_b, a.dual_closure_a,
                               a.ar_max, a.expect_blend_H)
            print_admission(result[0], result[1], A["label"], B["label"])
        rc = EXIT_CODE[result[0]]
    if a.provenance_out:
        doc = write_provenance(a.provenance_out, a.A, a.B, a.info_a, a.info_b, a.code_commit_a, a.code_commit_b, a.attach_a, a.attach_b,
                               result)
        print(f"来歴: {a.provenance_out}  (A 署名 {doc['runs']['A'].get('mesh_signature')}, B 署名 {doc['runs']['B'].get('mesh_signature')})")
    return rc


def report(A, B, n_layers=10):
    """A・B (load() の戻り、またはメッシャの出力から作った同じ形の dict) の表を出す。"""
    SA, SB = Struct(A["info"]), Struct(B["info"])
    for key in ("ni", "NJ", "nz", "jm", "i_te", "k_sw", "i_sw", "N_base"):
        if getattr(SA, key) != getattr(SB, key):
            raise SystemExit(f"A と B で構造が違う ({key}: {getattr(SA, key)} vs {getattr(SB, key)})。同じ格子水準の A/B だけを比べる")
    S = SA
    ib = B["info"]
    Lb = float(ib.get("te_wake_blend_H", 0.0) or 0.0)
    i_ret = int(ib["te_wake_i_last"]) + 1 if Lb > 0.0 and "te_wake_i_last" in ib else S.i_te
    H = A["H"]
    print(f"A: {A['label']}  (te_wake_blend_H {A['info'].get('te_wake_blend_H', 0.0)})")
    print(f"B: {B['label']}  (te_wake_blend_H {Lb}"
          + (f", 後縁の接線 {ib['te_wake_slope_te_deg']:.4f}°, 終端 {ib['te_wake_slope_end_deg']:.4f}°, 区間の station {ib['te_wake_n_stations']}"
             f" [i {ib['te_wake_i_first']}..{ib['te_wake_i_last']}], 中間線の最大移動 {ib['te_wake_max_dy']:.5f} H" if Lb > 0.0 else "") + ")")
    print(f"構造: ni {S.ni}  NJ {S.NJ}  nz {S.nz}  jm {S.jm}  i_te {S.i_te}  k_sw {S.k_sw}  i_sw {S.i_sw}  節点 {A['coords'].shape[0]}  H {H} m")
    G = S.grid(A["coords"])
    print(f"後縁の x = {G[S.i_te, S.jm, 0, 0]:.6f} m、上流の間隔 {(G[S.i_te,0,0,0]-G[S.i_te-1,0,0,0])/H:.5f} H、"
          f"下流の最初の間隔 {(G[S.i_te+1,0,0,0]-G[S.i_te,0,0,0])/H:.5f} H")

    # (v) を先に (同一性が崩れていたら以降の比較は意味が無い)
    idn = identity(A, B, S, ib)
    print("\n(v) 同一性 (B は A と比べて)")
    print(f"  節点数               : A {idn['n_nodes'][0]}  B {idn['n_nodes'][1]}  -> {'同一' if idn['n_nodes'][0] == idn['n_nodes'][1] else '**相違**'}")
    print(f"  x 座標 (全節点)       : {'ビット一致' if idn['x_equal'] else '**相違**'}")
    print(f"  z 座標 (全節点)       : {'ビット一致' if idn['z_equal'] else '**相違**'}")
    print(f"  ヘキサの接続          : {'同一' if idn['hex_equal'] else '**相違**'}  ({A['hexes'].shape[0]} hex)")
    print(f"  境界面の接続 (タグ別) : {'同一' if idn['bface_equal'] else ('**相違**' if idn['bface_equal'] is False else '比較不能 (境界面が無い入力)')}")
    if "moved" in idn:
        mv = idn["moved"]
        print(f"  y が動いた節点        : {mv.size}  (主ブロック外 {idn['moved_nonbase']})"
              + (f"  i {idn['moved_i'][0]}..{idn['moved_i'][1]}  j {idn['moved_j'][0]}..{idn['moved_j'][1]}"
                 f"  x {idn['moved_x_m'][0]:.6f}..{idn['moved_x_m'][1]:.6f} m  最大 |Δy| {idn['max_dy_m']:.4e} m ({idn['max_dy_m']/H:.5f} H)"
                 if mv.size else ""))
        if Lb > 0.0 and mv.size:
            ok_i = idn["moved_i"][0] >= int(ib["te_wake_i_first"]) and idn["moved_i"][1] <= int(ib["te_wake_i_last"])
            ok_j = idn["moved_j"][0] >= 1 and idn["moved_j"][1] <= S.NJ - 2
            print(f"  動いた範囲の検査      : i が変形区間内 {'OK' if ok_i else '**NG**'}、j が下端 (j=0)・上線 (j=NJ-1) を含まない "
                  f"{'OK' if ok_j else '**NG**'}、主ブロック外 0 {'OK' if idn['moved_nonbase'] == 0 else '**NG**'}")
        if idn["moved_by_tag"]:
            bad = {PHYS_NAME.get(p, p): n for p, n in idn["moved_by_tag"].items() if n and p not in IN_PLANE_TAGS}
            inpl = {PHYS_NAME.get(p, p): n for p, n in idn["moved_by_tag"].items() if n and p in IN_PLANE_TAGS}
            print(f"  境界タグの動いた節点  : 面内 (z 一定面) {inpl or 'なし'} / それ以外 {bad or 'なし'}"
                  f"  -> 外部境界・壁 {'不変' if not bad else '**動いている**'}")

    # (i)
    print(f"\n(i) 後縁の列 (i = i_te = {S.i_te}) の辺の角度 [deg]  (x 軸から、負 = 下向き)。各欄 min / median / max")
    rA, rB = te_column(A, S, n_layers), te_column(B, S, n_layers)
    up_same = all(rA[n][1] == rB[n][1] for n in rA)
    print(f"  {'群':44s} {'節点':>6s} | {('上流の辺 (A=B)' if up_same else '上流の辺 (A)'):^26s} | {'A: 最初の下流の辺':^26s} | "
          f"{'B: 最初の下流の辺':^26s} | {'折れ max A':>10s} {'B':>8s}")
    for name in rA:
        dA, uA, tA, n = rA[name]; dB, uB, tB, _ = rB[name]
        print(f"  {name:44s} {n:6d} | {_fmt3(uA)} | {_fmt3(dA)} | {_fmt3(dB)} | {tA:10.3f} {tB:8.3f}")
    if not up_same:
        print("  **上流の辺の角度が A と B で違う** (上流の形状が変わっている)")

    # (ii)
    print(f"\n(ii) i 線に沿った隣接辺の角度変化 |Δθ| の最大 [deg] とその位置 (i, j, k; x, y, z [m])")
    tA, tB = turn_maxima(A, S, i_ret), turn_maxima(B, S, i_ret)
    for key in tA:
        vA, locA, xA = tA[key]; vB, locB, xB = tB.get(key, (np.nan, None, (np.nan,) * 3))
        print(f"  {key[0]:28s} {key[1]:22s}  A {vA:7.3f} @ {locA} ({xA[0]:.5f},{xA[1]:.5f},{xA[2]:.5f})"
              f"   B {vB:7.3f} @ {locB} ({xB[0]:.5f},{xB[1]:.5f},{xB[2]:.5f})")

    # (iii)
    cols_te = [S.i_te]; cols_ret = list(range(S.i_te + 1, i_ret + 1))
    print(f"\n(iii) 近壁間隔 (中間線から {n_layers} 層、全 k)。y 方向の間隔 |Δy| の B/A と、i 線に垂直な間隔 |Δy|·cos θ [m]")
    print("      (後縁の列の垂直間隔は上流の辺 = 壁の向きで、復帰区間の列は前後の辺の平均の向きで測る)")
    for cname, cols, which in (("後縁の列", cols_te, 2), ("復帰区間の列", cols_ret, 1)):
        if not cols:
            continue
        spA, spB = spacings(A, S, cols, n_layers), spacings(B, S, cols, n_layers)
        for side, lab in (("lo", "下側"), ("up", "上側")):
            dyA, dnA = spA[side][0], spA[side][which]; dyB, dnB = spB[side][0], spB[side][which]
            r1 = np.abs(dyB[:, 0] / dyA[:, 0] - 1.0); rn = np.abs(dyB / dyA - 1.0)
            print(f"  {cname:8s} {lab}: 第 1 層 |Δy| の max|B/A−1| {r1.max():.3e}、{n_layers} 層の max|B/A−1| {rn.max():.3e}；"
                  f" 第 1 層の垂直間隔 A {dnA[:, 0].min():.4e}..{dnA[:, 0].max():.4e}  B {dnB[:, 0].min():.4e}..{dnB[:, 0].max():.4e}")

    # (iv)
    print("\n(iv) 全ヘキサの 8 頂点のヤコビアン (各頂点の 3 方向の前進差分の行列式。符号は A の多数派を正とする)")
    JA, JsA = corner_jacobians(A["coords"], A["hexes"]); sgn = 1.0 if np.median(JA) > 0 else -1.0
    JB, JsB = corner_jacobians(B["coords"], B["hexes"])
    touch = np.zeros(B["hexes"].shape[0], bool)
    if "moved" in idn and idn["moved"].size:
        mset = np.zeros(B["coords"].shape[0], bool); mset[idn["moved"]] = True
        touch = mset[B["hexes"]].any(axis=1)
    for lab, J, Js in (("A", JA, JsA), ("B", JB, JsB)):
        npos = int((sgn * J <= 0.0).sum()); nh = int(((sgn * J) <= 0.0).any(axis=1).sum())
        print(f"  {lab}: 非正の頂点 {npos} (ヘキサ {nh} / {J.shape[0]})、スケール済みヤコビアンの最小 {float((sgn * Js).min()):.4f}"
              + (f"、動いた節点を含むヘキサ {int(touch.sum())} の最小 {float((sgn * Js[touch]).min()):.4f}" if touch.any() else ""))
    if touch.any():
        print(f"  A の同じヘキサ {int(touch.sum())} の最小 {float((sgn * JsA[touch]).min()):.4f}")

    print("\n品質・双対幾何は別に回して A と B を並べて貼る:")
    print("  - check_mesh_quality の VERDICT: runner の prepare が <run>/MESH_QUALITY.txt に書く"
          " (`check_mesh_quality.py sern_qc.h5 --mode 3d --ar-max <mesh.ar_max>`、品質確認用の cell 変換に対して)")
    print("  - 双対 CV の閉性・体積: python3 solver_density_cuda/tools/check_dual_closure.py <run>/sern.h5")
    return 0


if __name__ == "__main__":
    sys.exit(main())
