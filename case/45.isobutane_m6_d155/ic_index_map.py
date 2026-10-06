#!/usr/bin/env python3
"""変形メッシュへの初期値写像 (plan tooling-nozzle-throat-monotone-r2 §6 E1、§5.1 #5b; 諮問 2026-10-06 ①)。

腕 B (r″ 単調壁) の G1 格子は、IC の run (run_0114) の格子と節点数・接続が同じで、壁と近傍の節点だけが最大 0.5 µm 動く。
最近傍 (interp_field) で写すと、壁の移動が第 1 セル厚 (0.34 µm) を超える所で 1 層内側の値を拾う (37 節点、局所の密度差 23 %)。
そこで、検査を通したうえで**保存済み保存量を番号でそのまま写す** (index モード)。予備 A/B 用に、同じ検査のうえで
interp_field と同じ最近傍の対応で保存済み保存量を写すモード (nearest) も持つ。

**限定例外**: `restart_field.py` の座標検査は緩めない (本体は変えない)。この写像は「変形メッシュへの初期値写像」の本 plan に限った
例外で、元の座標の場の再現ではない (体積が変わるので領域積分の保存も保証しない)。AGENTS.md への明記はユーザ判断 (plan §5.1 #11)。

検査 (どれか不成立なら**何も書かずに**止める。全項目を調べてから、まとめて理由を出す):
  - 節点数・接続: MESH/CONNE・VIZMESH/CONNE・CELLS/STRUCT・PLANES/STRUCT・CELLS/regionId と MESH の属性が同一
  - 論理位置 (i, j): VIZMESH/CONNE が構造格子の並び (節点番号 = i·nj + j、j が速い) そのものであること、
    両格子で x が i に、r が j に狭義単調増加、j = 0 が軸 (r = 0)
  - 境界種別: BCONDS の physID ごとに bcondKind・節点集合 (vizBfaceNodes)・面の並び・境界値が同一で、
    節点集合が種別どおりの論理辺 (入口 i = 0 / 出口 i = ni−1 / 軸 j = 0 / 壁 j = nj−1) に一致
  - 座標系・単位: solverConfig.yaml の mesh (軸対称・離散化など)・prepare_info.json の scale_m・座標の型・z が同一、
    保存場 (res) の座標がその IC 格子と一致
  - 化学種・エネルギー基準: 保存量のデータセット名・型・形の一致 (写す量はすべて SRC にあり、SRC の保存量はすべて DST にある)、
    solverConfig.yaml の physProp・turbulence・species_meta.yaml の一致、化学種署名 (名前・順序・MW・thermoHrefTemp)、
    SRC の属性 (species_hash と記録の完全性)。書き込み時は restart_field と同じく forge --resolve-species で宛先を解決して
    属性を継承する (--no-species-resolve は乾式確認用で、属性を付けない)
  - 変形後の要素の反転なし: 各四角形の符号付き面積と 4 つの角の向き (隣り合う辺の外積) の符号が SRC と同じで 0 でない
    (面積だけではねじれた四角形を見逃すので角も見る)
  - 座標の移動 ≤ --max-disp-m (既定 1e-6 m = 1 µm。座標は m 単位; r_t 0.0768075 m なら 1.3e-5 r_t)
転送: DST の /VALUE のうち wall_dist 以外 (restart_field が移すのと同じ集合) を、SRC の保存済み値から直接コピーする
  (原始量から組み直さない)。wall_dist と幾何量 (/VALUE 以外) は DST の値を保持する。書き込み後に、転送量が SRC (nearest では
  SRC[対応]) とビット一致すること、幾何量・wall_dist が書き込み前と同一であることを検査する。
記録 (--record、既定 DST の隣の IC_MAP.json): 検査結果、移動量の分布、転送した量、最近傍と番号写像の食い違い (数・位置・j の差) と
  保存量の差 (最大・RMS)。nearest モードではこれが予備 A/B (run_0146) の主記録。

usage: python3 ic_index_map.py SRC_res.h5 DST_input.h5 [--mode index|nearest] [--src-mesh SRC_mesh.h5 (既定 SRC の隣の nozzle.h5)]
         [--src-run DIR (既定 SRC の隣)] [--dst-run DIR (既定 DST の隣)] [--max-disp-m 1e-6] [--record PATH]
         [--no-species-resolve (乾式確認: forge を起動しない。属性は付けない)] [--forge BIN]
終了コード: 0 = VERDICT OK、2 = 検査不成立 (何も書いていない)、1 = 書き込み後の検査で不一致。
"""
import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

C = Path(__file__).resolve().parent
TOOLS = C.parents[1] / "solver_density_cuda/tools"
sys.path.insert(0, str(TOOLS))
import forge_species as fsp  # noqa: E402
from interp_field import centroids, is_3d  # noqa: E402  (最近傍の対応は interp_field と同じ規則で取る)

MAX_DISP_M = 1e-6
KEEP_FROM_DST = ("wall_dist",)                  # メッシュ由来 (restart_field と同じ)
TOPOLOGY = ("MESH/CONNE", "VIZMESH/CONNE", "CELLS/STRUCT", "PLANES/STRUCT", "CELLS/regionId")
BCOND_TOPO = ("iBPlanes", "iCells", "iPlanes", "vizBfaceNodes", "vizBfaceSizes")
# 保存量のデータセット名 (interp_field が扱う種類と同じ)。SRC の res は原始量も持つので、比べるのはこの名前だけ
CONS_RE = re.compile(r"^(ro|roU[xyz]|roe|roK|roOmega|roGamma|roReth|roXi|roY\d+|rog_.+|roQ[012]_.+)$")
# 境界種別 → 論理辺 (構造格子 i = 0..ni−1 (流れ方向)、j = 0..nj−1 (軸 → 壁))
EDGE_OF_KIND = (("inlet", "i0"), ("outlet", "i1"), ("axis", "j0"), ("slip", "j1"), ("wall", "j1"))
PLAN = "plans/accepted/tooling-nozzle-throat-monotone-r2.md §6 E1 (諮問 2026-10-06 ①)"


class MapRefused(Exception):
    """検査不成立 (書き込み前)。failures に全項目の理由を持つ。"""

    def __init__(self, failures, record):
        super().__init__("; ".join(failures))
        self.failures = failures
        self.record = record


def _sha_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


class _NoDupLoader(fsp._StrSafeLoader):
    """fsp.load_yaml_str と同じ解釈 (真偽値を true/false だけに絞る) で、明示キーの重複を拒否する。PyYAML は重複キーの後勝ち、
    ソルバ側の読み方と食い違いうるので、設定の比較には重複の無いものだけを使う (merge キー `<<` で入る値の上書きは重複としない)。"""

    def construct_mapping(self, node, deep=False):
        if isinstance(node, yaml.MappingNode):
            keys = [k.value for k, _ in node.value if k.tag != "tag:yaml.org,2002:merge"]
            dup = sorted({k for k in keys if keys.count(k) > 1})
            if dup:
                raise yaml.constructor.ConstructorError(None, None, f"重複キー {dup}", node.start_mark)
        return super().construct_mapping(node, deep)


def _yaml(p):
    """YAML を読む (無ければ None)。重複キー・構文誤りは ValueError。"""
    if not os.path.exists(p):
        return None
    try:
        with open(p) as f:
            return yaml.load(f, Loader=_NoDupLoader)
    except yaml.YAMLError as e:
        raise ValueError(f"{p} を読めない: {e}") from e


def geometry_digest(h5path) -> str:
    """/VALUE 以外の全データセットと /VALUE/wall_dist (幾何・メッシュ由来の量) の内容ハッシュ (名前・型・形・バイト列)。"""
    h = hashlib.sha256()
    with h5py.File(h5path, "r") as f:
        names = []
        f.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
        for n in sorted(names):
            if n.startswith("VALUE/") and n[len("VALUE/"):] not in KEEP_FROM_DST:
                continue
            a = np.asarray(f[n])
            h.update(n.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()


def mesh_digest(h5path) -> str:
    """格子 (座標・接続・境界の節点と面) の内容ハッシュ。2 つの prep が同じ格子かを比べる用 (/VALUE を含まない)。"""
    h = hashlib.sha256()
    with h5py.File(h5path, "r") as f:
        names = ["MESH/COORD"] + [t for t in TOPOLOGY if t in f]
        names += sorted(f"BCONDS/{b}/{k}" for b in f["BCONDS"] for k in BCOND_TOPO if k in f["BCONDS"][b])
        for n in names:
            a = np.asarray(f[n])
            h.update(n.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()


def structured_shape(viz_conne, n_nodes):
    """VIZMESH/CONNE ([型, n0, n1, n2, n3] の並び) が、節点番号 i·nj + j (j が速い) の構造格子の四角形
    [(i,j), (i+1,j), (i+1,j+1), (i,j+1)] を i 優先・j 速いの順に並べたものと**完全に一致**するかを調べ、(ni, nj) を返す。
    一致しなければ ValueError (論理位置 (i, j) を復元できない)。"""
    a = np.asarray(viz_conne)
    if a.size == 0 or a.size % 5:
        raise ValueError(f"VIZMESH/CONNE の長さ {a.size} が 5 の倍数でない (四角形だけの格子でない)")
    q = a.reshape(-1, 5)
    if np.unique(q[:, 0]).size != 1:
        raise ValueError(f"VIZMESH/CONNE の要素型が混在 ({np.unique(q[:, 0]).tolist()})")
    nj = int(q[0, 2])
    if nj < 2 or n_nodes % nj:
        raise ValueError(f"最初の四角形 {q[0, 1:].tolist()} から nj = {nj} を読めない (節点数 {n_nodes})")
    ni = n_nodes // nj
    i, j = np.meshgrid(np.arange(ni - 1), np.arange(nj - 1), indexing="ij")
    base = (i * nj + j).ravel()
    want = np.stack([base, base + nj, base + nj + 1, base + 1], axis=1)
    if q.shape[0] != want.shape[0] or not np.array_equal(q[:, 1:], want):
        raise ValueError(f"VIZMESH/CONNE が構造格子 (ni {ni}, nj {nj}, 節点番号 i·nj + j) の並びでない")
    return ni, nj


def logical_failures(coord, ni, nj, tag):
    """座標が論理位置どおりか: x は i に、r は j に狭義単調増加、j = 0 が軸 (r = 0)、z は全節点で同じ。"""
    out = []
    X, R = coord[:, 0].reshape(ni, nj), coord[:, 1].reshape(ni, nj)
    if not np.all(np.diff(R, axis=1) > 0):
        k = np.argwhere(~(np.diff(R, axis=1) > 0))[0]
        out.append(f"{tag}: r が j に狭義単調増加でない (最初 i {int(k[0])}, j {int(k[1])}→{int(k[1]) + 1})")
    if not np.all(np.diff(X, axis=0) > 0):
        k = np.argwhere(~(np.diff(X, axis=0) > 0))[0]
        out.append(f"{tag}: x が i に狭義単調増加でない (最初 i {int(k[0])}→{int(k[0]) + 1}, j {int(k[1])})")
    if not np.all(R[:, 0] == 0.0):
        out.append(f"{tag}: j = 0 が軸 (r = 0) でない (max |r| {float(np.abs(R[:, 0]).max()):.3e})")
    if np.unique(coord[:, 2]).size != 1:
        out.append(f"{tag}: z が全節点で同じでない (2D 格子でない)")
    return out


def signed_areas(coord, viz_conne):
    """各四角形の符号付き面積 [m²] と、4 つの角の向き (隣り合う 2 辺の外積) (x, y、float64 で計算)。
    面積の符号だけではねじれた四角形 (蝶形) を見逃す (正負の三角形が打ち消す) ので、角ごとの向きも返す。"""
    q = np.asarray(viz_conne).reshape(-1, 5)[:, 1:]
    x, y = coord[q, 0], coord[q, 1]
    area = 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
    ex1, ey1 = np.roll(x, -1, axis=1) - x, np.roll(y, -1, axis=1) - y       # 頂点 k → k+1
    ex0, ey0 = np.roll(x, 1, axis=1) - x, np.roll(y, 1, axis=1) - y         # 頂点 k → k−1
    corner = ex1 * ey0 - ey1 * ex0
    return area, corner


def _edge(name, ni, nj):
    idx = np.arange(ni * nj).reshape(ni, nj)
    return {"i0": idx[0, :], "i1": idx[-1, :], "j0": idx[:, 0], "j1": idx[:, -1]}[name]


def check_bconds(fs, fd, ni, nj):
    """境界種別: physID ごとの種別・節点集合・面の並び・境界値の一致と、種別どおりの論理辺。戻り値 (failures, 記録)。"""
    out, rec = [], {}
    bs, bd = sorted(fs["BCONDS"]), sorted(fd["BCONDS"])
    if bs != bd:
        return [f"境界の physID が違う (SRC {bs} / DST {bd})"], rec
    for b in bs:
        gs, gd = fs["BCONDS"][b], fd["BCONDS"][b]
        ks, kd = str(gs.attrs.get("bcondKind", "")), str(gd.attrs.get("bcondKind", ""))
        r = {"kind": kd}
        if ks != kd:
            out.append(f"境界 {b}: 種別が違う (SRC {ks} / DST {kd})")
        for k in BCOND_TOPO:
            if (k in gs) != (k in gd) or (k in gs and not np.array_equal(np.asarray(gs[k]), np.asarray(gd[k]))):
                out.append(f"境界 {b} ({kd}): {k} が違う")
        vs, vd = sorted(gs["VALUE"]) if "VALUE" in gs else [], sorted(gd["VALUE"]) if "VALUE" in gd else []
        if vs != vd:
            out.append(f"境界 {b} ({kd}): 境界値の名前が違う")
        else:
            bad = [k for k in vs if np.asarray(gs["VALUE"][k]).tobytes() != np.asarray(gd["VALUE"][k]).tobytes()]
            if bad:
                out.append(f"境界 {b} ({kd}): 境界値が違う {bad}")
        nodes = np.unique(np.asarray(gd["vizBfaceNodes"])) if "vizBfaceNodes" in gd else np.array([], int)
        want = [e for key, e in EDGE_OF_KIND if kd.lower().startswith(key)]
        if not want:
            out.append(f"境界 {b}: 種別 {kd} の論理辺を決められない (入口・出口・軸・壁のどれでもない)")
        elif not np.array_equal(nodes, np.sort(_edge(want[0], ni, nj))):
            out.append(f"境界 {b} ({kd}): 節点集合が論理辺 {want[0]} と一致しない ({nodes.size} 節点)")
        r.update(n_nodes=int(nodes.size), edge=(want[0] if want else None))
        rec[b] = r
    return out, rec


def displacement_record(cs, cd, ni, nj, scale_m):
    """番号ごとの座標の移動 |x_DST − x_SRC| の分布 (最大・位置・節点数・層ごと)。"""
    d = np.linalg.norm(cd - cs, axis=1)
    k = int(np.argmax(d))
    i, j = divmod(k, nj)
    moved = d[d > 0]
    hist_edges = [0.0, 1e-9, 1e-8, 1e-7, 2e-7, 5e-7, 1e-6, math.inf]
    hist = {f"({hist_edges[n]:g}, {hist_edges[n + 1]:g}]": int(np.count_nonzero((d > hist_edges[n]) & (d <= hist_edges[n + 1])))
            for n in range(len(hist_edges) - 1)}
    D = d.reshape(ni, nj)
    layers = {str(L): {"n_moved": int(np.count_nonzero(D[:, nj - 1 - L] > 0)), "max_m": float(D[:, nj - 1 - L].max())}
              for L in range(nj)}
    out = {"max_m": float(d.max()), "max_um": float(d.max() * 1e6),
           "max_rt": (float(d.max() / scale_m) if scale_m else None),
           "argmax": {"index": k, "i": i, "j": j, "x_m": float(cd[k, 0]), "r_m": float(cd[k, 1]),
                      "x_rt": (float(cd[k, 0] / scale_m) if scale_m else None), "layer_from_wall": nj - 1 - j},
           "n_nodes": int(d.size), "n_moved": int(moved.size),
           "moved_quantiles_m": ({q: float(np.quantile(moved, float(q))) for q in ("0.5", "0.9", "0.99")} if moved.size else None),
           "histogram_m": hist, "by_layer_from_wall": layers,
           "x_range_moved_rt": ([float(cd[d > 0, 0].min() / scale_m), float(cd[d > 0, 0].max() / scale_m)]
                                if (moved.size and scale_m) else None)}
    return d, out


def nearest_index(src_res, dst_h5):
    """interp_field と同じ対応規則: SRC の DOF 座標 (centroids) の cKDTree で DST の各 DOF の最近傍を取る。"""
    from scipy.spatial import cKDTree
    with h5py.File(src_res, "r") as s, h5py.File(dst_h5, "r") as d:
        nd = 3 if (is_3d(s) and is_3d(d)) else 2
        cs, cd = centroids(s, nd), centroids(d, nd)
    dist, idx = cKDTree(cs).query(cd)
    return idx.astype(np.int64), dist, nd


def nearest_vs_index(idx, dist, cd, ni, nj, scale_m, src_vals, names, list_max=200):
    """最近傍の対応が番号写像 (自分自身の番号) と食い違う節点の一覧と、2 つの転送で保存量がどれだけ違うか。"""
    n = idx.size
    own = np.arange(n)
    bad = np.flatnonzero(idx != own)
    i, j = np.divmod(own, nj)
    ii, jj = np.divmod(idx, nj)
    rows = [{"index": int(k), "i": int(i[k]), "j": int(j[k]), "i_nn": int(ii[k]), "j_nn": int(jj[k]),
             "dj": int(jj[k] - j[k]), "di": int(ii[k] - i[k]), "x_m": float(cd[k, 0]), "r_m": float(cd[k, 1]),
             "x_rt": (float(cd[k, 0] / scale_m) if scale_m else None), "layer_from_wall": int(nj - 1 - j[k]),
             "nn_dist_m": float(dist[k])} for k in bad[:list_max]]
    dj = (jj - j)[bad]
    diff = {}
    for name in names:
        v = src_vals[name].astype(np.float64)
        dv = v[idx] - v
        sc = float(np.abs(v).max()) or 1.0
        with np.errstate(divide="ignore", invalid="ignore"):
            loc = np.where(v != 0, np.abs(dv) / np.abs(v), np.where(dv != 0, np.inf, 0.0))
        diff[name] = {"max_abs": float(np.abs(dv).max()), "rms": float(np.sqrt(np.mean(dv * dv))),
                      "max_abs_over_max_abs_value": float(np.abs(dv).max() / sc), "n_differs": int(np.count_nonzero(dv)),
                      "max_local_rel": float(loc.max())}
    return {"n_mismatch": int(bad.size), "match_rate": float(1.0 - bad.size / n), "nn_dist_max_m": float(dist.max()),
            "dj_counts": {str(int(v)): int(np.count_nonzero(dj == v)) for v in np.unique(dj)},
            "di_counts": {str(int(v)): int(np.count_nonzero((ii - i)[bad] == v)) for v in np.unique((ii - i)[bad])},
            "x_range_rt": ([float(cd[bad, 0].min() / scale_m), float(cd[bad, 0].max() / scale_m)] if (bad.size and scale_m) else None),
            "j_range": ([int(j[bad].min()), int(j[bad].max())] if bad.size else None),
            "nodes": rows, "nodes_truncated": bool(bad.size > list_max),
            "conserved_diff_nearest_minus_index": diff}


def run_checks(src_res, dst_h5, src_mesh, src_run, dst_run, max_disp_m=MAX_DISP_M, resolve_species=True, forge=None):
    """書き込み前の検査。成立なら (record, ctx)、不成立なら MapRefused (全項目の理由)。"""
    fails, rec = [], {"checks": {}}

    def ck(name, items, detail=None):
        rec["checks"][name] = {"ok": not items, "failures": items, **({"detail": detail} if detail is not None else {})}
        fails.extend(f"[{name}] {x}" for x in items)

    if not (max_disp_m > 0 and math.isfinite(max_disp_m)):
        raise ValueError(f"--max-disp-m {max_disp_m} は正の有限値であること")
    ctx = {}
    with h5py.File(src_mesh, "r") as fs, h5py.File(dst_h5, "r") as fd, h5py.File(src_res, "r") as fr:
        # --- 節点数・接続 ---
        it = []
        ns, nd_ = int(fs["MESH/COORD"].shape[0]) // 3, int(fd["MESH/COORD"].shape[0]) // 3
        if ns != nd_:
            it.append(f"節点数が違う (SRC {ns} / DST {nd_})")
        for t in TOPOLOGY:
            if (t in fs) != (t in fd):
                it.append(f"{t} が片側にしか無い")
            elif t in fs and not np.array_equal(np.asarray(fs[t]), np.asarray(fd[t])):
                it.append(f"{t} が違う")
        if dict(fs["MESH"].attrs) != dict(fd["MESH"].attrs):
            it.append("MESH の属性 (nNodes・nCells・nPlanes など) が違う")
        ck("connectivity", it, {"n_nodes": nd_})
        # --- 論理位置 (i, j) ---
        it, ni, nj = [], None, None
        cs = np.asarray(fs["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        cd = np.asarray(fd["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        try:
            ni, nj = structured_shape(fd["VIZMESH/CONNE"], nd_)
            if "VIZMESH/CONNE" in fs:
                structured_shape(fs["VIZMESH/CONNE"], ns)
        except ValueError as e:
            it.append(str(e))
        if ni is not None and ns == nd_:
            it += logical_failures(cs, ni, nj, "SRC") + logical_failures(cd, ni, nj, "DST")
        for tag, run in (("SRC", src_run), ("DST", dst_run)):
            pi = os.path.join(run, "prepare_info.json")
            if ni is not None and os.path.exists(pi):
                m = (json.load(open(pi)).get("mesh") or {})
                if m.get("ni") is not None and (int(m["ni"]), int(m.get("nj", nj))) != (ni, nj):
                    it.append(f"{tag} の prepare_info の mesh (ni {m.get('ni')}, nj {m.get('nj')}) が接続から復元した (ni {ni}, nj {nj}) と違う")
        ck("logical_ij", it, {"ni": ni, "nj": nj})
        # --- 境界種別 ---
        if ni is not None:
            it, brec = check_bconds(fs, fd, ni, nj)
        else:
            it, brec = ["論理位置を復元できないので境界の論理辺を照合できない"], {}
        ck("boundaries", it, brec)
        # --- 座標系・単位 ---
        it = []
        cfg = {}
        for tag, run in (("SRC", src_run), ("DST", dst_run)):
            try:
                cfg[tag] = _yaml(os.path.join(run, "solverConfig.yaml"))
            except ValueError as e:
                it.append(f"{tag} の solverConfig.yaml を読めない (重複キー・構文誤り): {e}")
                cfg[tag] = None
                continue
            if cfg[tag] is None:
                it.append(f"{tag} の solverConfig.yaml が無い")
            elif not isinstance(cfg[tag], dict):
                it.append(f"{tag} の solverConfig.yaml が対応表でない")
                cfg[tag] = None
        cs_, cd_ = cfg["SRC"], cfg["DST"]
        if cs_ is not None and cd_ is not None:
            ms = {k: v for k, v in (cs_.get("mesh") or {}).items() if k not in ("meshFileName", "valueFileName")}
            md = {k: v for k, v in (cd_.get("mesh") or {}).items() if k not in ("meshFileName", "valueFileName")}
            if ms != md:
                it.append(f"solverConfig.yaml の mesh が違う (SRC {ms} / DST {md})")
        scale = {}
        for tag, run in (("SRC", src_run), ("DST", dst_run)):
            pi = os.path.join(run, "prepare_info.json")
            sc = json.load(open(pi)).get("scale_m") if os.path.exists(pi) else None
            scale[tag] = float(sc) if sc is not None else None
        if scale["SRC"] != scale["DST"]:
            it.append(f"prepare_info.json の scale_m が違う (SRC {scale['SRC']} / DST {scale['DST']})")
        if fs["MESH/COORD"].dtype != fd["MESH/COORD"].dtype:
            it.append(f"座標の型が違う (SRC {fs['MESH/COORD'].dtype} / DST {fd['MESH/COORD'].dtype})")
        if ns == nd_ and not np.array_equal(cs[:, 2], cd[:, 2]):
            it.append("z 座標が違う")
        if "MESH/COORD" not in fr or not np.array_equal(np.asarray(fr["MESH/COORD"]), np.asarray(fs["MESH/COORD"])):
            it.append("保存場 (res) の座標が IC 格子 (--src-mesh) の座標と一致しない (別の格子の場)")
        ck("coordinate_system_units", it, {"scale_m": scale, "coord_dtype": str(fd["MESH/COORD"].dtype)})
        # --- 化学種・エネルギー基準 (保存量のデータセット名・型・属性) ---
        it, srec = [], {}
        dnames = [n for n in sorted(fd["VALUE"]) if n not in KEEP_FROM_DST]
        scons = sorted(n for n in fr["VALUE"] if CONS_RE.match(n))
        miss = [n for n in dnames if n not in fr["VALUE"]]
        if miss:
            it.append(f"写す量が SRC に無い (量欠落): {miss}")
        extra = [n for n in scons if n not in dnames]
        if extra:
            it.append(f"SRC の保存量が DST に無い (黙って落とさない): {extra}")
        notcons = [n for n in dnames if not CONS_RE.match(n)]
        if notcons:
            it.append(f"DST の /VALUE に保存量でない名前がある: {notcons}")
        for n in dnames:
            if n in fr["VALUE"] and (fr["VALUE"][n].shape != fd["VALUE"][n].shape or fr["VALUE"][n].dtype != fd["VALUE"][n].dtype):
                it.append(f"{n}: 形・型が違う (SRC {fr['VALUE'][n].shape} {fr['VALUE'][n].dtype} / DST {fd['VALUE'][n].shape} {fd['VALUE'][n].dtype})")
            elif n in fr["VALUE"] and not np.all(np.isfinite(np.asarray(fr["VALUE"][n]))):
                it.append(f"{n}: SRC に非有限値がある")
        if "ro" in fr["VALUE"] and not np.all(np.asarray(fr["VALUE"]["ro"]) > 0):
            it.append("SRC の ro に 0 以下がある")
        if cs_ is not None and cd_ is not None:
            for key in ("physProp", "turbulence"):
                if (cs_.get(key) or {}) != (cd_.get(key) or {}):
                    it.append(f"solverConfig.yaml の {key} が違う (エネルギー基準・化学種の定義)")
        ms_, md_ = os.path.join(src_run, "species_meta.yaml"), os.path.join(dst_run, "species_meta.yaml")
        if os.path.exists(ms_) != os.path.exists(md_) or (os.path.exists(ms_) and _sha_file(ms_) != _sha_file(md_)):
            it.append("species_meta.yaml が違う")
        try:
            sig_s, sig_d = fsp.species_signature(src_run), fsp.species_signature(dst_run)
            diffs = fsp.compare_signatures(sig_s, sig_d)
            unv = [x for x in diffs if "unverifiable" in x]
            hard = [x for x in diffs if "unverifiable" not in x]
            it += [f"化学種署名: {x}" for x in hard]
            req = fsp.required_conserved(sig_s)
            lack = [n for n in req if n not in dnames]
            if lack:
                it.append(f"化学種署名が要求する保存量が写す量に無い: {lack}")
            srec.update(names=sig_s["names"], thermoHrefTemp=sig_s["thermoHrefTemp"], required=req,
                        config_unverifiable=unv)
        except Exception as e:  # noqa: BLE001 — 署名が作れないことも不成立として記録する
            it.append(f"化学種署名を作れない: {e}")
        st = fsp.source_species_state(src_res)
        srec.update(src_state=st["state"], src_attrs=st["attrs"], src_why=st["why"])
        if "h0" in fr["VALUE"]:
            srec["src_h0_includes_k"] = (int(fr["VALUE"]["h0"].attrs["h0_includes_k"]) if "h0_includes_k" in fr["VALUE"]["h0"].attrs else None)
        tm = str(((cs_ or {}).get("physProp") or {}).get("thermalMethod", 0)).strip() if cs_ else None
        if tm == "2" and st["state"] != "verified":
            it.append(f"SRC の化学種属性を検証できない ({st['state']}: {st['why']})")
        ck("species_energy", it, srec)
        ctx.update(ni=ni, nj=nj, cs=cs, cd=cd, names=dnames, scale=scale["DST"])
    # --- 変形後の要素の反転 ---
    it = []
    if ni is not None and ns == nd_:
        with h5py.File(dst_h5, "r") as fd:
            vc = np.asarray(fd["VIZMESH/CONNE"])
        (As, Cs), (Ad, Cd) = signed_areas(cs, vc), signed_areas(cd, vc)
        bad_c = np.any((np.sign(Cs) != np.sign(Cd)) | (Cd == 0), axis=1)          # 角の向きが SRC と違う・退化 (ねじれを含む)
        flip = np.flatnonzero((np.sign(As) != np.sign(Ad)) | (Ad == 0) | bad_c)
        if flip.size:
            k = int(flip[0])
            it.append(f"反転・ねじれ・退化した要素が {flip.size} 個 (最初 要素 {k}: i {k // (nj - 1)}, j {k % (nj - 1)})")
        ratio = Ad / As
        ck("no_inversion", it, {"n_elements": int(As.size), "n_flipped_or_zero": int(flip.size),
                                "n_corner_orientation_changed": int(np.count_nonzero(bad_c)),
                                "area_ratio_min": float(ratio.min()), "area_ratio_max": float(ratio.max())})
    else:
        ck("no_inversion", ["節点数・論理位置が不成立なので要素の反転を調べられない"])
    # --- 座標の移動 ---
    if ni is not None and ns == nd_:
        d, drec = displacement_record(cs, cd, ni, nj, ctx["scale"])
        it = [] if drec["max_m"] <= max_disp_m else [f"最大移動 {drec['max_m']:.4e} m > {max_disp_m:g} m (i {drec['argmax']['i']}, j {drec['argmax']['j']})"]
        drec["limit_m"] = max_disp_m
        ck("displacement", it, drec)
    else:
        ck("displacement", ["節点数・論理位置が不成立なので移動量を測れない"])
    # --- 化学種の属性の継承判定 (restart_field と同じ; 書き込み前)。ほかの検査が不成立なら forge を起動せずに止める ---
    species_plan = None
    if fails:
        raise MapRefused(fails, rec)
    if resolve_species:
        try:
            species_plan = fsp.plan_inherit(str(src_res), str(dst_run), forge=forge, tool="ic_index_map", inplace=True)
            rec["checks"]["species_energy"]["detail"]["inherit"] = (
                {"species_hash": species_plan["species_hash"], "species_input_unverified": species_plan["species_input_unverified"]}
                if species_plan else None)
        except fsp.SpeciesCheckError as e:
            ck("species_inherit", [str(e)])
    else:
        rec["checks"]["species_energy"]["detail"]["inherit"] = "未解決 (--no-species-resolve: 乾式確認。属性は付けない)"
    ctx["species_plan"] = species_plan
    if fails:
        raise MapRefused(fails, rec)
    return rec, ctx


def apply_map(src_res, dst_h5, src_mesh=None, src_run=None, dst_run=None, mode="index", max_disp_m=MAX_DISP_M,
              record=None, resolve_species=True, forge=None) -> dict:
    """検査 → 転送 → 書き込み後の検査 → 記録。検査不成立は MapRefused (何も書かない; 記録だけ書く)。"""
    src_res, dst_h5 = Path(src_res).resolve(), Path(dst_h5).resolve()
    src_mesh = Path(src_mesh).resolve() if src_mesh else src_res.parent / "nozzle.h5"
    src_run = str(Path(src_run).resolve() if src_run else src_res.parent)
    dst_run = str(Path(dst_run).resolve() if dst_run else dst_h5.parent)
    record = Path(record) if record else dst_h5.parent / "IC_MAP.json"
    if mode not in ("index", "nearest"):
        raise ValueError(f"mode は index / nearest: {mode}")
    head = {"tool": "ic_index_map.py", "plan": PLAN, "mode": mode, "src_res": str(src_res), "src_mesh": str(src_mesh),
            "src_run": src_run, "dst": str(dst_h5), "dst_run": dst_run, "species_resolved": bool(resolve_species),
            "note": ("番号写像は変形した論理格子への初期値の移送で、元の座標の場の再現ではない (体積が変わるので領域積分の保存も保証しない)。"
                     "restart_field の座標検査は緩めない (本 plan に限った例外)")}
    try:
        rec, ctx = run_checks(src_res, dst_h5, src_mesh, src_run, dst_run, max_disp_m, resolve_species, forge)
    except MapRefused as e:
        out = {**head, **e.record, "VERDICT": "REFUSED", "failures": e.failures}
        record.write_text(json.dumps(out, indent=1, ensure_ascii=False))
        raise
    names, ni, nj = ctx["names"], ctx["ni"], ctx["nj"]
    with h5py.File(src_res, "r") as s:
        sv = {n: np.asarray(s["VALUE"][n]) for n in names}
    idx_nn, dist, nd = nearest_index(src_res, dst_h5)
    cmp_ = nearest_vs_index(idx_nn, dist, ctx["cd"], ni, nj, ctx["scale"], sv, names)
    cmp_["nd"] = nd
    idx = np.arange(ctx["cd"].shape[0]) if mode == "index" else idx_nn
    geo0, sha0 = geometry_digest(dst_h5), _sha_file(dst_h5)
    fsp.write_species_attrs(str(dst_h5), None)      # 書き込み途中で失敗しても古い属性が残らないように先に消す (restart_field と同じ)
    with h5py.File(dst_h5, "r+") as d:
        for n in names:
            d["VALUE"][n][...] = sv[n][idx]
        if resolve_species:
            fsp.commit_inherit(d, ctx["species_plan"])
    # --- 書き込み後の検査: 転送量のビット一致・幾何量と wall_dist の保持 ---
    bad = []
    with h5py.File(dst_h5, "r") as d:
        for n in names:
            if np.asarray(d["VALUE"][n]).tobytes() != sv[n][idx].tobytes():
                bad.append(n)
    geo1 = geometry_digest(dst_h5)
    post = {"bit_exact_vs_src": not bad, "not_bit_exact": bad, "geometry_and_wall_dist_unchanged": geo0 == geo1}
    out = {**head, **rec, "transferred": names, "kept_from_dst": [n for n in KEEP_FROM_DST], "post_write": post,
           "nearest_vs_index": cmp_, "dst_sha256_before": sha0, "dst_sha256_after": _sha_file(dst_h5),
           "dst_mesh_digest": mesh_digest(dst_h5),
           "species_attrs_written": (fsp.field_species_attrs(str(dst_h5)) if resolve_species else None)}
    ok = (not bad) and geo0 == geo1
    out["VERDICT"] = "OK" if ok else "FAIL (書き込み後の検査で不一致)"
    record.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    if not ok:
        raise RuntimeError(f"書き込み後の検査で不一致: ビット不一致 {bad}, 幾何量・wall_dist 保持 {geo0 == geo1}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", help="IC の保存場 res_*.h5")
    ap.add_argument("dst", help="宛先の forge 入力 h5 (/VALUE を上書きする)")
    ap.add_argument("--mode", default="index", choices=("index", "nearest"))
    ap.add_argument("--src-mesh", help="IC の格子 (既定: SRC の隣の nozzle.h5)")
    ap.add_argument("--src-run", help="IC の run ディレクトリ (既定: SRC の隣)")
    ap.add_argument("--dst-run", help="宛先 run ディレクトリ (既定: DST の隣)")
    ap.add_argument("--max-disp-m", type=float, default=MAX_DISP_M, help="座標の移動の上限 [m] (既定 1e-6 = 1 µm)")
    ap.add_argument("--record", help="記録 JSON (既定: DST の隣の IC_MAP.json)")
    ap.add_argument("--no-species-resolve", action="store_true", help="乾式確認: forge --resolve-species を起動しない (属性は付けない)")
    ap.add_argument("--forge", help="--resolve-species を持つ forge (既定: FORGE_BIN)")
    a = ap.parse_args(argv)
    try:
        out = apply_map(a.src, a.dst, a.src_mesh, a.src_run, a.dst_run, a.mode, a.max_disp_m, a.record,
                        resolve_species=not a.no_species_resolve, forge=a.forge)
    except MapRefused as e:
        print("[ic_index_map] REFUSED (何も書いていない):\n  " + "\n  ".join(e.failures))
        return 2
    except RuntimeError as e:
        print(f"[ic_index_map] FAIL: {e}")
        return 1
    dr = out["checks"]["displacement"]["detail"]
    nv = out["nearest_vs_index"]
    print(f"[ic_index_map] mode {out['mode']}: 移した量 {out['transferred']} (wall_dist は DST を保持)")
    print(f"[ic_index_map] 移動 最大 {dr['max_um']:.4f} µm (i {dr['argmax']['i']}, j {dr['argmax']['j']}, x/r_t {dr['argmax']['x_rt']}) "
          f"/ 動いた節点 {dr['n_moved']} / {dr['n_nodes']}")
    print(f"[ic_index_map] 最近傍が番号と食い違う節点 {nv['n_mismatch']} (j の差 {nv['dj_counts']}, x/r_t {nv['x_range_rt']}, j {nv['j_range']})")
    print(f"VERDICT: {out['VERDICT']} (転送量 SRC{'' if out['mode'] == 'index' else '[最近傍]'} とビット一致、幾何量・wall_dist 保持)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
