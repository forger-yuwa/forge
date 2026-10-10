#!/usr/bin/env python3
"""g3 → g4 の初期場で、カウル・側壁の双子節点の内外が保たれるかを**時間積分の前 (0 step)** に数える (読むだけ)。

plan tooling-sern-te-wake-grid §5.1 #4「時間積分の前に 0 step の A/B」(codex diagnose 2026-10-10
`notes/reviews/2026-10-10-g4-initial-field-diagnose.md`)。旧メッシャ (`mesh_sern3d.py`) はカウルの上下の面を `dup1`、
側壁の内外を `dup2` という**別 ID の節点**で持ち、厚さ 0 の所では座標が一致する (双子)。座標だけの最近傍
(`interp_field.py`) は双子を区別できないので、g4 の初期場の作り方を実格子で確かめる。

  python3 case/46.sern_design/diag/twin_ic_check.py --g3 G3 --g4 G4_RUN [--src-res G3_res.h5] [--json OUT.json]
          [--examples 8] [--band-h 0.05] [--workdir DIR]

  G3 = g3 の run ディレクトリか node 変換済みの入力 h5 (境界面の接続 `/BCONDS/<physID>/vizBface*` と `VIZMESH/CONNE` を読む)。
  G4_RUN = g4 の run ディレクトリ (`runner_sern3d.prepare` の出力: sern.h5・prepare_info.json・bcondConfig.yaml)。
  --src-res = interp_field に実際に渡す SRC (例 run_1079 の最終 res)。座標が G3 の h5 とビット一致することを確かめてから
              KD 木をそれで作る (無ければ G3 の h5 の座標で作る。interp_field の centroids は node ではどちらも MESH/COORD)。

**正解 (各格子で境界面の接続から決める。B の index 規則・メッシャの番号は使わない)**:
  節点が属する境界面 (元の接続 `vizBfaceNodes`) のタグで、cowl_in・sidewall_in だけ → 排気側 (ノズル内、IN)、
  cowl_out・sidewall_out だけ → 外気側 (OUT)、両方 → 両側で共有 (BOTH)、どれも無い → NONE。physID は run の
  bcondConfig.yaml (無ければ `mesh_sern3d.PHYS_SERN3D`)。
  **双子** = 座標 (h5 の `MESH/COORD`、−0.0 は 0.0 に揃える) が数値として一致する節点の組。各組がちょうど 2 節点で
  IN と OUT を 1 つずつ持てば正解が決まる。そうでない組 (3 節点以上・BOTH・NONE・同じ側どうし) が 1 つでもあれば UNDECIDABLE。
  自由端・後縁で上下が共有する単一節点 (BOTH で組に入らない) は双子でないので除く (件数だけ出す)。
  双子の分類: カウル (cowl_* だけ)・側壁 (sidewall_* だけ)・角 (両方、z = W/2 のカウル側端)。
  「壁節点」= 4 タグのどれかに属し IN か OUT に決まる全節点 (双子 + 板厚のある所の上下の面など。参考の行)。
  **構造上の対** (参考): prepare_info.json のメッシャ情報から `dup1`↔`base(i, jm, k)`・`dup2`↔`base(i, j, k_sw)` を再現した組
  (板厚 `cowl_thickness` > 0 の所では座標が一致しない)。座標一致の組がすべて構造上の対かを照合する。

**A (interp_field と同じ探索)**: `solver_density_cuda/tools/interp_field.py` を import して `is_3d`・`centroids`・`cKDTree` を
  そのまま使う (照合座標の次元も `main()` と同じく両方 3D なら x,y,z)。`tree = cKDTree(centroids(SRC)); tree.query(centroids(DST))`
  を g4 の評価節点について行い (点ごとの照会なので全節点で照会した結果と同じ)、donor の g3 ラベルを g4 へ写す。
  g3 のラベル: 壁節点は上の正解と同じ規則。壁でない donor (境界層の内部節点など) は、g3 の壁節点から**格子の辺に沿った
  距離** (ヘキサの辺、ユークリッド長、IN の壁と OUT の壁からの多始点 Dijkstra) の近い側 (同距離・到達なしは「不定」)。
  双子の内外は格子の上で繋がっていない (カウル・側壁を回り込まないと届かない) ので、近い側 = 板のどちら側か。
  計算は壁節点から max(`--band-h`·H, donor の最近傍距離の 4 倍) 以内の節点に限る。数えるもの: 誤割当 (donor のラベル ≠ 正解:
  反対側・共有節点・不定)、
  双子の 2 節点が同じ donor になった組の数 (座標が同じなら照会結果も同じなので、全組が同じ donor になるのが期待値)。

**B (現行 3D の初期化そのもの)**: `design/forge_design/evaluate/runner_sern3d.py` の `paste_region_ic3d` を**そのまま呼ぶ**。
  prepare と同じ引数 (`minfo` = prepare_info.json の `mesh` を discretization が node のときだけ、`st` = `states`、
  `scale` = `H_m`、`half_W_m`) で、g4 の `CELLS/centCoords` と `VALUE` の該当データセット (同じ名前・型・長さ) だけを持つ一時 h5 に
  書かせ、書かれた値が排気側の状態か外気側の状態かでラベルを読む (状態は `runner_sern.region_ic_arrays` に [True, False] を
  渡して作り、h5 の型に丸めて完全一致で照合。roe は γ に依存するので照合に使わない)。y_mid には呼ぶと止める関数を渡すので、
  座標の分類 (index の分類が使えないときの分岐) に落ちたら UNDECIDABLE にする。prepare_info.json の格子署名
  (`mesh_provenance.mesh_signature`) は g4 の h5 から再計算して照合する (別の格子の情報で分類しない)。
  **B2 (照合)**: g4 の sern.h5 の `VALUE` (prepare が書いた初期場、0 step なら未変更) も同じ方法で読み、B と違えば B_NG。
  読めない (値がどちらの状態でもない = 書き換え済み) ときは B2 は照合不能と書いて B だけで判定する。

判定 (最後の行、終了コード 0 / 1 / 2):
  `TWIN_IC VERDICT: B_OK (A n 件誤り)` … B の誤割当が双子・壁節点とも 0 (B2 も一致か照合不能)。n = A の双子の誤割当数
  `TWIN_IC VERDICT: B_NG` … B (か B2) に誤割当が 1 件以上 → 時間積分へ進まない (plan §5.1 #4)
  `TWIN_IC VERDICT: UNDECIDABLE` … 正解の決まらない双子がある・双子が 0 組・入力が欠ける・B が index の分類を使わない等
  (判定不能は合格ではない)。A の結果は判定に使わない (解釈用。plan §5.1 #4 の分岐は B で決まる)。
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "solver_density_cuda" / "tools"
DESIGN = ROOT / "design"

IN, OUT, BOTH, NONE, UNK = 1, 2, 3, 0, -1          # ラベル (UNK = 内部節点で側が決まらない)
LAB_NAME = {IN: "排気側", OUT: "外気側", BOTH: "共有", NONE: "壁なし", UNK: "不定"}
IN_TAGS, OUT_TAGS = ("cowl_in", "sidewall_in"), ("cowl_out", "sidewall_out")
COWL_TAGS, SIDE_TAGS = ("cowl_in", "cowl_out"), ("sidewall_in", "sidewall_out")
CATS = ("カウル", "側壁", "角 (カウル∩側壁)")
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))   # gmsh 型 5 の辺


class Undecidable(Exception):
    """判定不能 (入力の欠落・正解が決まらない・初期化の経路が想定と違う)。"""


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def interp_field_module():
    """interp_field.py をそのまま読む (探索の実装を写さないため。tools を sys.path に足すのは interp_field 自身も同じ)。"""
    return _load_module(TOOLS / "interp_field.py", "interp_field_for_twin_ic")


def design_modules():
    if str(DESIGN) not in sys.path:
        sys.path.insert(0, str(DESIGN))
    from forge_design.evaluate import runner_sern as R2          # noqa: E402
    from forge_design.evaluate import runner_sern3d as R3        # noqa: E402
    from forge_design.meshing.mesh_sern3d import PHYS_SERN3D     # noqa: E402
    return R2, R3, PHYS_SERN3D


# ---------------------------------------------------------------------------------------------- 読み込み
def _mesh_file_of(run_dir):
    p = os.path.join(run_dir, "solverConfig.yaml")
    if os.path.exists(p):
        m = re.search(r"meshFileName\s*:\s*[\"']?([^\"',}\s]+)", open(p).read())
        if m and os.path.exists(os.path.join(run_dir, m.group(1))):
            return os.path.join(run_dir, m.group(1))
    return os.path.join(run_dir, "sern.h5")


def resolve(path):
    """run ディレクトリか h5 → (h5, run_dir or None)。h5 を直接渡したときも、隣に prepare_info.json があれば run とみなす。"""
    if os.path.isdir(path):
        return _mesh_file_of(path), path
    d = os.path.dirname(os.path.abspath(path))
    return path, (d if os.path.exists(os.path.join(d, "prepare_info.json")) else None)


def phys_ids(run_dir, default):
    """境界名 → physID。run の bcondConfig.yaml を優先 (無ければメッシャの PHYS_SERN3D)。"""
    ids = {k: int(default[k]) for k in IN_TAGS + OUT_TAGS}
    src = "mesh_sern3d.PHYS_SERN3D"
    p = os.path.join(run_dir, "bcondConfig.yaml") if run_dir else None
    if p and os.path.exists(p):
        txt = open(p).read()
        for k in ids:
            m = re.search(r"^\s*" + re.escape(k) + r"\s*:\s*\{\s*physID\s*:\s*(\d+)", txt, re.M)
            if m is None:
                raise Undecidable(f"{p} に {k} の physID が無い")
            ids[k] = int(m.group(1))
        src = p
    return ids, src


def load_grid(h5, run_dir, phys_default, need_hex=False):
    """node 変換済みの入力 h5 → dict (座標・境界面・[ヘキサ]・prepare_info)。"""
    if not os.path.exists(h5):
        raise Undecidable(f"格子 {h5} が無い")
    g = {"h5": h5, "run_dir": run_dir}
    with h5py.File(h5, "r") as f:
        c = np.asarray(f["MESH/COORD"][...])
        g["coord_dtype"] = str(c.dtype)
        c = c.reshape(-1, 3) + c.dtype.type(0.0)          # −0.0 → 0.0 (バイト列で組を作るので数値の一致に揃える)
        g["coords_raw"] = c
        g["coords"] = c.astype(np.float64)
        a = f["MESH"].attrs
        g["node"] = bool("VIZMESH" in f and "nCells" in a and "nNodes" in a and int(a["nCells"]) == int(a["nNodes"]))
        g["n"] = c.shape[0]
        g["perm"] = np.asarray(f["MESH/RENUMBER_PERM"][...], np.int64) if "MESH/RENUMBER_PERM" in f else None
        bf = {}
        if "BCONDS" in f:
            for k in f["BCONDS"]:
                gg = f["BCONDS/" + k]
                if "vizBfaceSizes" in gg and "vizBfaceNodes" in gg:
                    sz = np.asarray(gg["vizBfaceSizes"][...], np.int64).reshape(-1)
                    nd = np.asarray(gg["vizBfaceNodes"][...], np.int64).reshape(-1)
                    if int(sz.sum()) != nd.size:
                        raise Undecidable(f"{h5}: BCONDS/{k} の vizBfaceSizes の総和 {int(sz.sum())} と vizBfaceNodes の長さ {nd.size} が違う")
                    bf[int(k)] = nd
        g["bfaces"] = bf
        if need_hex:
            if "VIZMESH/CONNE" not in f:
                raise Undecidable(f"{h5}: VIZMESH/CONNE (ヘキサの接続) が無い")
            cn = np.asarray(f["VIZMESH/CONNE"][...], np.int64)
            if cn.size == 0 or cn.size % 9 or not np.all(cn[::9] == 9):
                raise Undecidable(f"{h5}: VIZMESH/CONNE がヘキサ (XDMF 型 9) だけの列でない")
            g["hexes"] = cn.reshape(-1, 9)[:, 1:]
    if not g["node"]:
        raise Undecidable(f"{h5}: node 変換の h5 でない (VIZMESH が無いか CV 数 ≠ 節点数)")
    g["phys"], g["phys_src"] = phys_ids(run_dir, phys_default)
    miss = [k for k, v in g["phys"].items() if v not in bf]
    if miss:
        raise Undecidable(f"{h5}: 境界面の接続 (vizBfaceNodes) が無いタグ {miss} → 正解のラベルを決められない")
    ip = os.path.join(run_dir, "prepare_info.json") if run_dir else None
    g["info"] = json.load(open(ip)) if ip and os.path.exists(ip) else None
    return g


# ---------------------------------------------------------------------------------------------- 正解のラベル
def wall_labels(g):
    """境界面のタグから各節点の側 (IN / OUT / BOTH / NONE) と、カウル・側壁のどちらの壁に属するか。"""
    n = g["n"]; P = g["phys"]; bf = g["bfaces"]
    in_m, out_m, cowl_m, side_m = (np.zeros(n, bool) for _ in range(4))
    for k in IN_TAGS:
        in_m[bf[P[k]]] = True
    for k in OUT_TAGS:
        out_m[bf[P[k]]] = True
    for k in COWL_TAGS:
        cowl_m[bf[P[k]]] = True
    for k in SIDE_TAGS:
        side_m[bf[P[k]]] = True
    lab = np.full(n, NONE, np.int8)
    lab[in_m & ~out_m] = IN
    lab[out_m & ~in_m] = OUT
    lab[in_m & out_m] = BOTH
    return lab, cowl_m, side_m


def coincident_groups(coords_raw):
    """座標が一致する節点の組 (2 節点以上)。戻り = (pairs (m, 2) int64, 3 節点以上の組のリスト)。"""
    c = np.ascontiguousarray(coords_raw)
    v = c.view(np.dtype((np.void, c.dtype.itemsize * 3))).ravel()
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    multi = cnt[inv] >= 2
    ids = np.nonzero(multi)[0]
    order = ids[np.argsort(inv[ids], kind="stable")]
    sizes = cnt[np.unique(inv[order])]
    groups = np.split(order, np.cumsum(sizes)[:-1]) if order.size else []
    pairs = np.array([gr for gr in groups if len(gr) == 2], np.int64).reshape(-1, 2)
    big = [gr for gr in groups if len(gr) > 2]
    return pairs, big


def category(cowl_m, side_m, members):
    cw = bool(cowl_m[members].any()); sd = bool(side_m[members].any())
    return CATS[2] if (cw and sd) else (CATS[0] if cw else (CATS[1] if sd else "壁タグなし"))


def truth(g):
    """正解: 双子の組・各組の分類・決まらない組・共有の単一節点・壁節点 (IN/OUT に決まる全節点)。"""
    lab, cowl_m, side_m = wall_labels(g)
    pairs, big = coincident_groups(g["coords_raw"])
    la, lb = lab[pairs[:, 0]], lab[pairs[:, 1]]
    ok = ((la == IN) & (lb == OUT)) | ((la == OUT) & (lb == IN))
    cat = np.array([category(cowl_m, side_m, p) for p in pairs], dtype=object) if len(pairs) else np.zeros(0, object)
    in_group = np.zeros(g["n"], bool)
    in_group[pairs.ravel()] = True
    for gr in big:
        in_group[gr] = True
    undecided = [{"nodes": [int(x) for x in p], "labels": [LAB_NAME[int(lab[x])] for x in p], "category": c,
                  "coord": [float(v) for v in g["coords"][p[0]]]} for p, c in zip(pairs[~ok], cat[~ok])]
    undecided += [{"nodes": [int(x) for x in gr], "labels": [LAB_NAME[int(lab[x])] for x in gr], "category": category(cowl_m, side_m, gr),
                   "coord": [float(v) for v in g["coords"][gr[0]]], "note": f"{len(gr)} 節点の組"} for gr in big]
    shared = np.nonzero((lab == BOTH) & ~in_group)[0]
    wall = np.nonzero((lab == IN) | (lab == OUT))[0]
    return {"lab": lab, "cowl": cowl_m, "side": side_m, "pairs": pairs, "pair_ok": ok, "pair_cat": cat, "big": big,
            "undecided": undecided, "shared": shared, "wall": wall}


def structural_pairs(info, perm):
    """prepare_info.json のメッシャ情報 (mesh) から dup1 (カウル上コピー)・dup2 (側壁外コピー) と元の節点の組を再現する
    (mesh_sern3d.generate_sern_mesh3d の番号付け)。戻り = (pairs (m, 2) = (コピー, 元), 種別の配列)。h5 の番号に写す。"""
    m = info["mesh"]
    ni, NJ, nz, jm = int(m["ni"]), int(m["NJ"]), int(m["nz"]), int(m["jm"])
    i_te, k_sw, i_sw = int(m["i_te"]), int(m["k_sw"]), int(m["i_sw"])
    N_base = ni * NJ * nz
    no_outer = nz == k_sw + 1
    full = i_te * (k_sw + 1)
    shared = full - (max(i_te - i_sw, 0) if not no_outer else 0)
    n1 = int(m["n_dup_cowl"])
    if n1 not in (shared, full):
        raise Undecidable(f"カウル上コピーの数 {n1} がメッシャの規則 ({shared} / {full}) と合わない")
    share = n1 == shared and shared != full
    base = lambda i, j, k: (i * NJ + j) * nz + k
    out, kind = [], []
    nid = N_base
    for i in range(i_te):
        for k in range(k_sw + 1):
            if k == k_sw and i >= i_sw and not no_outer and share:
                continue
            out.append((nid, base(i, jm, k))); kind.append("dup1"); nid += 1
    for i in (range(0) if no_outer else range(i_sw)):
        for j in range(jm + 1, NJ):
            out.append((nid, base(i, j, k_sw))); kind.append("dup2"); nid += 1
    if nid - N_base - n1 != int(m.get("n_dup_side", nid - N_base - n1)):
        raise Undecidable(f"側壁外コピーの数 {nid - N_base - n1} が info の n_dup_side {m.get('n_dup_side')} と違う")
    sp = np.array(out, np.int64).reshape(-1, 2)
    if perm is not None:                       # perm[新] = 旧 (メッシャの番号) → 新 = inv[旧]
        inv = np.empty_like(perm); inv[perm] = np.arange(perm.size)
        sp = inv[sp]
    return sp, np.array(kind)


# ---------------------------------------------------------------------------------------------- A: g3 のラベル (内部節点を含む)
def _multi_source(G, src):
    """多始点 Dijkstra (各節点から最も近い始点までの距離)。min_only の無い古い scipy では、始点へ極小の重みの辺を張った
    仮想の始点 1 つからの最短距離で代える (重み 1e-300 は距離に足しても丸めで消える)。"""
    from scipy.sparse import bmat, coo_matrix
    from scipy.sparse.csgraph import dijkstra
    try:
        return dijkstra(G, directed=False, indices=src, min_only=True)
    except TypeError:
        n = G.shape[0]
        S = coo_matrix((np.full(src.size, 1e-300), (np.zeros(src.size, np.int64), src)), shape=(1, n))
        G2 = bmat([[None, S], [S.T, G]]).tocsr()
        return dijkstra(G2, directed=False, indices=0)[1:]


def side_labels(g, lab, band_m):
    """全節点のラベル: 壁節点は境界面のタグ、壁でない節点は格子の辺に沿った距離で近い側の壁 (band_m [m] 以内の節点だけ)。"""
    from scipy.sparse import coo_matrix
    from scipy.spatial import cKDTree
    out = lab.astype(np.int8).copy()
    out[out == NONE] = UNK
    seeds_all = np.nonzero(lab != NONE)[0]
    if seeds_all.size == 0:
        return out, {"band_nodes": 0}
    X = g["coords"]
    d, _ = cKDTree(X[seeds_all]).query(X, distance_upper_bound=band_m)
    band = np.isfinite(d)
    loc = np.full(g["n"], -1, np.int64); bid = np.nonzero(band)[0]; loc[bid] = np.arange(bid.size)
    H = g["hexes"]
    H = H[band[H].any(axis=1)]
    e = np.concatenate([H[:, [a, b]] for a, b in HEX_EDGES])
    e = e[band[e[:, 0]] & band[e[:, 1]]]
    e = np.unique(np.sort(e, axis=1), axis=0)
    w = np.linalg.norm(X[e[:, 0]] - X[e[:, 1]], axis=1)
    w = np.maximum(w, 1e-300)                 # 長さ 0 の辺 (双子どうしを結ぶ辺は無いはず) も辺として残す
    G = coo_matrix((w, (loc[e[:, 0]], loc[e[:, 1]])), shape=(bid.size, bid.size)).tocsr()
    dist = {}
    for s in (IN, OUT):
        src = loc[np.nonzero(lab == s)[0]]
        src = src[src >= 0]
        dist[s] = _multi_source(G, src) if src.size else np.full(bid.size, np.inf)
    di, do = dist[IN], dist[OUT]
    inner = lab[bid] == NONE
    li = np.full(bid.size, UNK, np.int8)
    li[(di < do)] = IN
    li[(do < di)] = OUT
    out[bid[inner]] = li[inner]
    return out, {"band_nodes": int(bid.size), "band_edges": int(e.shape[0]), "band_m": float(band_m),
                 "d_in": di, "d_out": do, "loc": loc}


# ---------------------------------------------------------------------------------------------- B: paste_region_ic3d
def _no_y_mid(*_a, **_k):
    raise Undecidable("paste_region_ic3d が座標の分類 (y_mid) に落ちた = index の分類が使われない "
                      "(discretization が node でないか、centCoords の長さが info の nodes と違う)")


IC_KEYS_SKIP = ("roe",)       # 照合に使わない (cpg の roe は γ に依存。γ は prepare_info に無い)


def region_states(R2, st, gamma, dtypes):
    """runner_sern.region_ic_arrays に [排気, 外気] を渡して両状態の値を作り、h5 の型に丸める。区別できる量だけ返す。"""
    a = R2.region_ic_arrays(np.array([True, False]), st, gamma)
    ref = {}
    for k, v in a.items():
        if k in IC_KEYS_SKIP or k not in dtypes:
            continue
        v = np.asarray(v).astype(dtypes[k])
        if v[0] != v[1]:
            ref[k] = (v[0], v[1])
    return ref


def decode_region(values, ref, n):
    """書かれた値 → ラベル (IN = 排気の状態、OUT = 外気の状態、UNK = どちらでもない)。"""
    if not ref:
        raise Undecidable("排気と外気の状態を区別できる量が無い")
    is_ex, is_en = np.ones(n, bool), np.ones(n, bool)
    for k, (vex, ven) in ref.items():
        x = values[k]
        is_ex &= x == vex
        is_en &= x == ven
    lab = np.full(n, UNK, np.int8)
    lab[is_ex & ~is_en] = IN
    lab[is_en & ~is_ex] = OUT
    return lab


def paste_labels(g4, workdir=None, gamma=1.4):
    """B: prepare と同じ引数で paste_region_ic3d を一時 h5 に対して呼び、書かれた値からラベルを読む。B2: g4 の h5 の VALUE も読む。"""
    R2, R3, _ = design_modules()
    info = g4["info"]
    if info is None:
        raise Undecidable(f"g4 の prepare_info.json が無い ({g4['run_dir']})")
    for k in ("mesh", "states", "H_m", "half_W_m"):
        if k not in info:
            raise Undecidable(f"g4 の prepare_info.json に {k} が無い")
    disc = info.get("discretization", "node")
    minfo, st = info["mesh"], info["states"]
    # prepare_info の格子情報が g4 の h5 のものか (格子署名を h5 から再計算して照合)
    prov = (info.get("mesh_provenance") or {}).get("mesh_signature")
    sig_note = "prepare_info に格子署名が無い (照合なし)"
    if prov:
        mz = _load_module(TOOLS / "mark_zero_thickness_edges.py", "mzte_for_twin_ic")
        try:
            sig = str(mz.mesh_signature(g4["h5"]))
        except Exception as e:  # noqa: BLE001
            raise Undecidable(f"g4 の h5 から格子署名を再計算できない ({type(e).__name__}: {e}) → prepare_info.json との結び付けを確かめられない")
        if sig != prov:
            raise Undecidable(f"g4 の h5 の格子署名 {sig[:12]}… が prepare_info.json の {prov[:12]}… と違う (別の格子の情報)")
        sig_note = f"格子署名一致 ({sig[:12]}…)"
    tmpd = Path(tempfile.mkdtemp(prefix="twin_ic_", dir=workdir))
    try:
        tmp = tmpd / "paste.h5"
        with h5py.File(g4["h5"], "r") as f, h5py.File(tmp, "w") as t:
            t.create_dataset("CELLS/centCoords", data=f["CELLS/centCoords"][...])
            for k in f["VALUE"]:
                if k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roXi") or re.fullmatch(r"roY\d+", k):
                    ds = f["VALUE/" + k]
                    t.create_dataset("VALUE/" + k, shape=ds.shape, dtype=ds.dtype)   # 0 で初期化 (書かれなければ照合で分かる)
            vals_b2 = {k: np.asarray(f["VALUE/" + k][...]) for k in t["VALUE"]}
        # prepare と同じ呼び方 (runner_sern3d.prepare: minfo=minfo if disc == "node" else None)
        R3.paste_region_ic3d(tmp, _no_y_mid, float(info["H_m"]), float(info["half_W_m"]), st, gamma,
                             minfo=minfo if disc == "node" else None)
        with h5py.File(tmp, "r") as t:
            vals = {k: np.asarray(t["VALUE/" + k][...]) for k in t["VALUE"]}
            dtypes = {k: t["VALUE/" + k].dtype for k in t["VALUE"]}
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    ref = region_states(R2, st, gamma, dtypes)
    n = g4["n"]
    labB = decode_region(vals, ref, n)
    labB2 = decode_region({k: v for k, v in vals_b2.items() if k in ref}, {k: r for k, r in ref.items() if k in vals_b2}, n) \
        if all(k in vals_b2 for k in ref) else np.full(n, UNK, np.int8)
    return labB, labB2, {"keys": sorted(ref), "disc": disc, "sig": sig_note, "gas_model": st.get("gas_model"),
                         "ic_ro": {"exhaust": float(ref["ro"][0]), "ext": float(ref["ro"][1])} if "ro" in ref else None}


# ---------------------------------------------------------------------------------------------- 本体
def run_check(g3_path, g4_path, src_res=None, band_h=0.05, workdir=None, n_examples=8, gamma=1.4):
    """戻り = (verdict 文字列, 終了コード, report dict, 表示用の行)。"""
    R2, R3, PHYS = design_modules()
    L = []
    rep = {"g3": None, "g4": None}
    try:
        h3, d3 = resolve(g3_path)
        h4, d4 = resolve(g4_path)
        g4 = load_grid(h4, d4, PHYS)
        g3 = load_grid(h3, d3, PHYS, need_hex=True)
        rep["g3"], rep["g4"] = h3, h4
        T4, T3 = truth(g4), truth(g3)
        H = float((g4["info"] or {}).get("H_m") or (g3["info"] or {}).get("H_m") or 1.0)
    except Undecidable as e:
        return "UNDECIDABLE", 2, {"reason": str(e)}, [f"判定不能: {e}"]

    L.append(f"g3: {h3}  節点 {g3['n']}  座標 {g3['coord_dtype']}  physID ({g3['phys_src']}) {g3['phys']}")
    L.append(f"g4: {h4}  節点 {g4['n']}  座標 {g4['coord_dtype']}  physID ({g4['phys_src']}) {g4['phys']}")
    reasons = []

    # ---------------------------------------------------------------- 正解
    L.append("")
    L.append("== 正解 (境界面の接続: cowl_in・sidewall_in = 排気側、cowl_out・sidewall_out = 外気側) ==")
    L.append(f"{'格子':4s} {'双子の組':>8s} " + " ".join(f"{c:>14s}" for c in CATS) + f" {'正解なし':>8s} {'共有の単一節点':>14s} {'壁節点':>8s}")
    for nm, T in (("g3", T3), ("g4", T4)):
        cnt = {c: int(np.sum(T["pair_cat"] == c)) for c in CATS}
        L.append(f"{nm:4s} {len(T['pairs']):8d} " + " ".join(f"{cnt[c]:14d}" for c in CATS)
                 + f" {len(T['undecided']):8d} {len(T['shared']):14d} {len(T['wall']):8d}")
        rep[nm + "_truth"] = {"pairs": int(len(T["pairs"])), "by_category": cnt, "undecided": T["undecided"][:50],
                              "n_undecided": len(T["undecided"]), "shared_single": int(len(T["shared"])), "wall_nodes": int(len(T["wall"]))}
        if T["undecided"]:
            reasons.append(f"{nm}: 正解の決まらない双子 {len(T['undecided'])} 組 (例 {T['undecided'][0]})")
    if len(T4["pairs"]) == 0:
        reasons.append("g4 に双子 (座標一致の組) が 0 組 → この A/B の対象が無い")

    # 構造上の対 (参考): 座標一致の組がすべてメッシャの dup1/dup2 か
    struct = {}
    for nm, g, T in (("g3", g3, T3), ("g4", g4, T4)):
        if g["info"] is None or "mesh" not in g["info"]:
            struct[nm] = None
            continue
        try:
            sp, kind = structural_pairs(g["info"], g["perm"])
        except Undecidable as e:
            struct[nm] = None; L.append(f"  {nm}: 構造上の対を再現できない ({e})"); continue
        if sp.size and sp.max() >= g["n"]:
            struct[nm] = None; L.append(f"  {nm}: 構造上の対の番号が節点数を超える (info と h5 が別の格子?)"); continue
        coinc = np.all(g["coords_raw"][sp[:, 0]] == g["coords_raw"][sp[:, 1]], axis=1)
        key = lambda p: (int(min(p)), int(max(p)))
        sset = {key(p) for p in sp}
        non_struct = [key(p) for p in T["pairs"] if key(p) not in sset]
        struct[nm] = {"pairs": sp, "kind": kind, "coinc": coinc}
        L.append(f"  {nm} 構造上の対: dup1 {int(np.sum(kind == 'dup1'))} (座標一致 {int(np.sum(coinc & (kind == 'dup1')))})、"
                 f"dup2 {int(np.sum(kind == 'dup2'))} (座標一致 {int(np.sum(coinc & (kind == 'dup2')))})、"
                 f"座標一致の組のうち構造上の対でないもの {len(non_struct)}")
        rep[nm + "_struct"] = {"dup1": int(np.sum(kind == "dup1")), "dup1_coincident": int(np.sum(coinc & (kind == "dup1"))),
                               "dup2": int(np.sum(kind == "dup2")), "dup2_coincident": int(np.sum(coinc & (kind == "dup2"))),
                               "coincident_not_structural": len(non_struct)}

    # ---------------------------------------------------------------- A
    IF = interp_field_module()
    src = src_res or h3
    with h5py.File(src, "r") as s, h5py.File(h4, "r") as d:
        s3d, d3d = IF.is_3d(s), IF.is_3d(d)
        nd = 3 if (s3d and d3d) else 2                      # interp_field.main と同じ
        cs = IF.centroids(s, nd)
        cd = IF.centroids(d, nd)
        if src_res:
            cres = np.asarray(s["MESH/COORD"][...]).reshape(-1, 3) + np.float32(0.0)
            if cres.shape != g3["coords_raw"].shape or not np.array_equal(cres.astype(np.float64), g3["coords"]):
                return "UNDECIDABLE", 2, {"reason": "--src-res の座標が g3 の h5 と一致しない"}, L + ["判定不能: --src-res の座標が g3 の h5 と一致しない (g3 の正解を当てられない)"]
    if s3d and not d3d:
        reasons.append("interp_field は SRC 3D・DST 2D を拒否する (A を作れない)")
    tree = IF.cKDTree(cs)
    P4, C4 = T4["pairs"][T4["pair_ok"]], T4["pair_cat"][T4["pair_ok"]]     # 正解の決まった組だけを評価する
    ev_twin = P4.ravel()
    ev_wall = T4["wall"]
    ev = np.unique(np.concatenate([ev_twin, ev_wall]))
    if ev.size == 0:
        return "UNDECIDABLE", 2, {"reason": "g4 に評価する節点が無い"}, L + ["判定不能: g4 に双子も壁節点も無い"]
    dist, idx = tree.query(cd[ev])            # interp_field.main の `dist, idx = tree.query(cd)` と同じ (点ごとの照会)
    pos = np.full(g4["n"], -1, np.int64); pos[ev] = np.arange(ev.size)
    rep["_A_eval"], rep["_A_donor"] = ev, idx
    # 帯の幅: 指定 (既定 0.05 H) と、donor の最近傍距離の 4 倍の大きい方 (donor が帯の外に出て「不定」になるのを防ぐ)
    band_m = max(band_h * H, 4.0 * float(dist.max()))
    lab3, side_info = side_labels(g3, T3["lab"], band_m)
    donor_lab = lab3[idx]
    L.append("")
    L.append(f"== A: interp_field の座標最近傍 (照合座標 {'x,y,z' if nd == 3 else 'x,y'}、SRC {src}) ==")
    L.append(f"  g3 の壁でない donor のラベル: 壁節点から {band_m:.4g} m ({band_m / H:.3g} H) 以内 {side_info['band_nodes']} 節点で、"
             f"格子の辺に沿った距離の近い側。donor の最近傍距離 max {float(dist.max()):.4e} m")

    def a_rows(sel_nodes):
        t = T4["lab"][sel_nodes]; dl = donor_lab[pos[sel_nodes]]
        wrong = dl != t
        return {"n": int(sel_nodes.size), "wrong": int(wrong.sum()), "opposite": int(np.sum(wrong & ((dl == IN) | (dl == OUT)))),
                "shared": int(np.sum(dl == BOTH)), "unknown": int(np.sum((dl == UNK) | (dl == NONE))),
                "donor_wall": int(np.sum(T3["lab"][idx[pos[sel_nodes]]] != NONE))}, wrong

    L.append(f"{'対象':22s} {'節点':>7s} {'誤割当':>7s} {'(反対側':>8s} {'共有':>6s} {'不定)':>6s} {'donor が壁節点':>14s} {'同一 donor の組':>16s}")
    A_rep = {}
    for c in CATS:
        P = P4[C4 == c]
        if not len(P):
            L.append(f"{'双子: ' + c:22s} {0:7d}"); A_rep[c] = {"n": 0}; continue
        r, _ = a_rows(P.ravel())
        same = int(np.sum(idx[pos[P[:, 0]]] == idx[pos[P[:, 1]]]))
        r["same_donor_pairs"] = same; r["pairs"] = int(len(P)); A_rep[c] = r
        L.append(f"{'双子: ' + c:22s} {r['n']:7d} {r['wrong']:7d} {r['opposite']:8d} {r['shared']:6d} {r['unknown']:6d} {r['donor_wall']:14d} {same:9d}/{len(P):<6d}")
    rA_twin, wA_twin = a_rows(ev_twin) if ev_twin.size else ({"n": 0, "wrong": 0}, np.zeros(0, bool))
    if ev_twin.size:
        rA_twin["same_donor_pairs"] = int(np.sum(idx[pos[P4[:, 0]]] == idx[pos[P4[:, 1]]]))
        L.append(f"{'双子 (計)':22s} {rA_twin['n']:7d} {rA_twin['wrong']:7d} {rA_twin['opposite']:8d} {rA_twin['shared']:6d} {rA_twin['unknown']:6d} {rA_twin['donor_wall']:14d} {rA_twin['same_donor_pairs']:9d}/{len(P4):<6d}")
    rA_wall, wA_wall = a_rows(ev_wall)
    L.append(f"{'壁節点 (参考)':22s} {rA_wall['n']:7d} {rA_wall['wrong']:7d} {rA_wall['opposite']:8d} {rA_wall['shared']:6d} {rA_wall['unknown']:6d} {rA_wall['donor_wall']:14d}")
    if struct.get("g4") is not None:
        sp = struct["g4"]["pairs"]; nc = ~struct["g4"]["coinc"]
        ok_sp = np.all(pos[sp] >= 0, axis=1)
        spn = sp[nc & ok_sp]
        same_nc = int(np.sum(idx[pos[spn[:, 0]]] == idx[pos[spn[:, 1]]])) if len(spn) else 0
        L.append(f"  構造上の対で座標が一致しない組 (板厚のある所) {len(spn)}: 同一 donor {same_nc}")
        A_rep["struct_noncoincident"] = {"pairs": int(len(spn)), "same_donor_pairs": same_nc}
    A_rep["twins_total"] = rA_twin; A_rep["wall_total"] = rA_wall
    rep["A"] = A_rep
    exA = []
    if ev_twin.size:
        bad = ev_twin[wA_twin]
        for e in bad[:n_examples]:
            k = pos[e]
            exA.append({"g4_node": int(e), "xyz": [float(v) for v in g4["coords"][e]], "truth": LAB_NAME[int(T4['lab'][e])],
                        "g3_donor": int(idx[k]), "donor_xyz": [float(v) for v in g3["coords"][idx[k]]],
                        "donor_label": LAB_NAME[int(donor_lab[k])], "dist_m": float(dist[k])})
    rep["A_examples"] = exA

    # ---------------------------------------------------------------- B
    L.append("")
    L.append("== B: paste_region_ic3d (現行 3D の領域別の一様な初期場、prepare と同じ引数) ==")
    try:
        labB, labB2, binfo = paste_labels(g4, workdir=workdir, gamma=gamma)
    except Undecidable as e:
        reasons.append(f"B: {e}")
        labB = labB2 = None
        L.append(f"  判定不能: {e}")
    B_rep = {}
    nB_bad = nB2_bad = 0
    exB = []
    if labB is not None:
        L.append(f"  照合した量 {binfo['keys']} (gas_model {binfo['gas_model']}、ro 排気 {binfo['ic_ro']['exhaust'] if binfo['ic_ro'] else '-'} / "
                 f"外気 {binfo['ic_ro']['ext'] if binfo['ic_ro'] else '-'})、discretization {binfo['disc']}、{binfo['sig']}")
        undec_B = int(np.sum(labB[ev] == UNK))
        if undec_B:
            reasons.append(f"B: 書かれた値が排気・外気のどちらの状態でもない評価節点 {undec_B} (paste の出力を読めない)")
        b2_ok = bool(np.all(labB2[ev] != UNK))
        L.append(f"{'対象':22s} {'節点':>7s} {'B 誤割当':>9s} {'B2 誤割当':>10s}")
        for c in CATS + ("双子 (計)", "壁節点"):
            if c in CATS:
                sel = P4[C4 == c].ravel()
            elif c == "双子 (計)":
                sel = ev_twin
            else:
                sel = ev_wall
            t = T4["lab"][sel]
            nb = int(np.sum(labB[sel] != t))
            nb2 = int(np.sum(labB2[sel] != t)) if b2_ok else None
            B_rep[c] = {"n": int(sel.size), "wrong": nb, "wrong_B2": nb2}
            L.append(f"{('双子: ' + c) if c in CATS else c:22s} {sel.size:7d} {nb:9d} {('-' if nb2 is None else nb2):>10}")
        nB_bad = B_rep["双子 (計)"]["wrong"] + B_rep["壁節点"]["wrong"]
        if b2_ok:
            nB2_bad = int(np.sum(labB2[ev] != labB[ev]))
            L.append(f"  B2 (g4 の sern.h5 の VALUE = prepare が書いた初期場) と B の不一致: {nB2_bad} 節点")
        else:
            L.append("  B2: g4 の sern.h5 の VALUE が排気・外気の状態として読めない (prepare 後に書き換え済み?) → 照合なし")
        sh = T4["shared"]
        if sh.size:
            L.append(f"  (参考) 共有の単一節点 {sh.size}: B は排気 {int(np.sum(labB[sh] == IN))} / 外気 {int(np.sum(labB[sh] == OUT))}")
        for e in ev[(labB[ev] != T4["lab"][ev])][:n_examples]:
            exB.append({"g4_node": int(e), "xyz": [float(v) for v in g4["coords"][e]], "truth": LAB_NAME[int(T4['lab'][e])],
                        "B": LAB_NAME[int(labB[e])]})
        B_rep["B2_available"] = b2_ok; B_rep["B2_mismatch"] = nB2_bad; B_rep["info"] = binfo
    rep["B"] = B_rep
    rep["B_examples"] = exB

    # ---------------------------------------------------------------- 例
    if exA:
        L.append("")
        L.append("例 (A の誤割当): g4 節点 [x y z] 正解 ← g3 donor [x y z] donor のラベル (距離)")
        for x in exA:
            L.append(f"  {x['g4_node']:9d} [{x['xyz'][0]:.6f} {x['xyz'][1]:.6f} {x['xyz'][2]:.6f}] {x['truth']} ← "
                     f"{x['g3_donor']:9d} [{x['donor_xyz'][0]:.6f} {x['donor_xyz'][1]:.6f} {x['donor_xyz'][2]:.6f}] {x['donor_label']} ({x['dist_m']:.3e} m)")
    if exB:
        L.append("例 (B の誤割当): g4 節点 [x y z] 正解 / B")
        for x in exB:
            L.append(f"  {x['g4_node']:9d} [{x['xyz'][0]:.6f} {x['xyz'][1]:.6f} {x['xyz'][2]:.6f}] {x['truth']} / {x['B']}")
    if len(T4["pairs"]):
        P = T4["pairs"][:min(3, len(T4["pairs"]))]
        L.append("例 (g4 の双子): " + "; ".join(f"{int(a)}/{int(b)} [{g4['coords'][a][0]:.6f} {g4['coords'][a][1]:.6f} {g4['coords'][a][2]:.6f}] "
                                               f"{LAB_NAME[int(T4['lab'][a])]}/{LAB_NAME[int(T4['lab'][b])]} {T4['pair_cat'][i]}"
                                               for i, (a, b) in enumerate(P)))

    # ---------------------------------------------------------------- 判定
    rep["reasons"] = reasons
    if reasons:
        L.append("")
        L.extend("判定不能: " + r for r in reasons)
        return "UNDECIDABLE", 2, rep, L
    if nB_bad > 0 or nB2_bad > 0:
        return "B_NG", 1, rep, L
    return f"B_OK (A {rA_twin['wrong']} 件誤り)", 0, rep, L


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g3", required=True, help="g3 の run ディレクトリか node 変換済みの入力 h5")
    ap.add_argument("--g4", required=True, help="g4 の run ディレクトリ (prepare の出力)")
    ap.add_argument("--src-res", help="interp_field に渡す SRC (g3 の res h5)。座標が g3 の h5 と一致することを確かめて使う")
    ap.add_argument("--band-h", type=float, default=0.05, help="g3 の内部節点のラベルを付ける帯の幅 / H (既定 0.05)")
    ap.add_argument("--gamma", type=float, default=1.4, help="paste_region_ic3d に渡す γ (cpg の roe だけに効く。ラベルには効かない)")
    ap.add_argument("--workdir", help="一時 h5 を置く場所 (既定: システムの一時ディレクトリ)")
    ap.add_argument("--examples", type=int, default=8)
    ap.add_argument("--json", help="結果を JSON で書く")
    a = ap.parse_args(argv)
    try:
        verdict, code, rep, lines = run_check(a.g3, a.g4, src_res=a.src_res, band_h=a.band_h, workdir=a.workdir,
                                              n_examples=a.examples, gamma=a.gamma)
    except Exception as e:  # noqa: BLE001   (例外の終了コード 1 を B_NG と取り違えないよう、判定不能として 2 で返す)
        import traceback
        traceback.print_exc()
        verdict, code, rep, lines = "UNDECIDABLE", 2, {"reason": f"{type(e).__name__}: {e}"}, [f"判定不能: 例外 {type(e).__name__}: {e}"]
    print("\n".join(lines))
    print(f"TWIN_IC VERDICT: {verdict}")
    if a.json:
        rep = {k: v for k, v in rep.items() if not k.startswith("_")}     # 内部の配列 (A の donor) は書かない
        rep["verdict"] = verdict
        with open(a.json, "w") as f:
            json.dump(rep, f, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return code


if __name__ == "__main__":
    sys.exit(main())
