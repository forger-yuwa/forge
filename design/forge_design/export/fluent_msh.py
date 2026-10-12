"""Fluent の ASCII メッシュ (2D 軸対称) の書き出しと読み直し検査、gmsh 経由の CGNS、境界条件表
(plans/accepted/tooling-nozzle-wt-deliverables.md §4.4・§6 V4)。

- **面の向き**: 2D の Fluent は面 $(n_0, n_1, c_0, c_1)$ で **$n_0\\to n_1$ の進行方向の左側が $c_0$**。
  TGrid 2D の実物 (OpenFOAM-10 の elbow.msh、1,454 面) で全面がこの向きだった (plan §3)。3D の右手系の規則を $z$ に
  押し出した推定は逆になるので使わない。ファイル冒頭の注記 `(node-0 node-1 right-cell left-cell)` とも逆なので、
  注記も根拠にしない。
- 向きの測定 (`measure_orientation_2d`) は書き出し器の規則を使わず、セル重心と面の幾何から $c_0$ が左か右かを数える
  (plan §6 V4(a))。ParaView の `vtkFLUENTReader` は $c_0$ の面の節点順でセルを組むだけで向きの誤りを止めないので、
  向きの検証には使わない。
- 番号: 10・12・13 節の見出しと面の並びは 16 進。45 節 (ゾーン名) のゾーン ID は 10 進 (Fluent の書式。ParaView の
  `vtkFLUENTReader` も ID 10〜15 を 10 進で対応付けることを確認, 2026-10-13)。既定のゾーン数 (6) では同じ表記になる。
"""
from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path

import numpy as np

# 面ゾーンの型 (Fluent の境界条件の型番号、10 進)。ファイルには 16 進で書く (37 → 25)
ZONE_TYPE_NAMES = {2: "interior", 3: "wall", 4: "pressure-inlet", 5: "pressure-outlet", 7: "symmetry",
                   9: "pressure-far-field", 10: "velocity-inlet", 20: "mass-flow-inlet", 36: "outflow", 37: "axis"}
DEFAULT_ZONE_TYPES = {"wall": 3, "inlet": 4, "outlet": 5, "axis": 37}
# 境界ゾーンの ID は名前順でなくこの順に 3 から振る (無い名前は飛ばし、ほかの名前は後ろに名前順で足す)
ZONE_ORDER = ("wall", "inlet", "outlet", "axis")
CELL_ZONE_NAME = "fluid"
INTERIOR_ZONE_NAME = "interior-fluid"
# 動作圧 (Fluent の Operating Pressure)。境界条件表の圧力の入力欄は「絶対値 − 動作圧」
OPERATING_PRESSURE_PA = 0.0

_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_\-]*$")


# ---------------------------------------------------------------------------
# 面の組み立て
# ---------------------------------------------------------------------------

def _planar_xy(coords) -> np.ndarray:
    """(N,2) または (N,3) の座標から (x, y) を取り出す。(N,3) なら z が全て 0 であること (2D 平面メッシュ)。"""
    c = np.asarray(coords, dtype=float)
    if c.ndim != 2 or c.shape[1] not in (2, 3):
        raise ValueError(f"coords の形 {c.shape} は (N,2) か (N,3)")
    if c.shape[1] == 3 and np.any(c[:, 2] != 0.0):
        raise ValueError("coords の z が 0 でない節点がある (2D 平面メッシュだけを書く)")
    if not np.all(np.isfinite(c[:, :2])):
        raise ValueError("coords に有限でない値がある")
    return np.ascontiguousarray(c[:, :2])


def _quad_areas(xy: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """四角形の符号付き面積 (反時計回りで正)。桁落ちを避けるため各セルの節点の平均を原点にして靴紐公式を取る。"""
    p = xy[quads]                                   # (M,4,2)
    p = p - p.mean(axis=1, keepdims=True)
    x, y = p[..., 0], p[..., 1]
    return 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)


def _as_edges(e, name) -> np.ndarray:
    a = np.asarray(e, dtype=np.int64)
    if a.size == 0:
        return a.reshape(0, 2)
    if a.ndim != 2 or a.shape[1] != 2:
        raise ValueError(f"境界 '{name}' の辺の形 {a.shape} は (K,2)")
    return a


def build_faces_2d(coords, quads, bedges) -> dict:
    """四角形メッシュ (構造化に限らない) から Fluent の面を作る。

    各面は (n0, n1, c0, c1) (0 始まり、境界の c1 = -1)。**c0 が n0→n1 の左側**: 四角形を反時計回りにたどる辺は
    そのセルが左なので、その向きのまま c0 = そのセル。内部面は先に現れた側のセルの辺の向きを使い、他方を c1 にする。
    境界辺は `bedges` の並びのまま、向きだけをセルの辺から決める。

    例外: 四角形が時計回り・面積 0、節点の重複、3 枚以上のセルが共有する辺、隣り合うセルの向きの不整合、
    境界辺がどの四角形の辺とも一致しない (または内部の辺)、境界辺の重複、どの境界にも属さない外周辺。
    戻り: {"xy", "n_cells", "interior" (F,4), "boundary" {名前: (K,4)}}。
    """
    xy = _planar_xy(coords)
    N = len(xy)
    q = np.asarray(quads, dtype=np.int64)
    if q.ndim != 2 or q.shape[1] != 4 or len(q) == 0:
        raise ValueError(f"quads の形 {q.shape} は (M,4) (M ≥ 1)")
    if q.min() < 0 or q.max() >= N:
        raise ValueError(f"quads の節点番号が範囲外 [{q.min()}, {q.max()}] (節点 {N})")
    qs = np.sort(q, axis=1)
    if np.any(qs[:, 1:] == qs[:, :-1]):
        raise ValueError("同じ節点を 2 回以上使う四角形がある")
    area = _quad_areas(xy, q)
    bad = np.flatnonzero(~(area > 0.0))
    if bad.size:
        raise ValueError(f"時計回り (または面積 0) の四角形が {bad.size} 個 (最初のセル {int(bad[0])}、面積 {area[bad[0]]:.3e})。"
                         "四角形は (x, r) 平面で反時計回りに並べる")
    M = len(q)
    a = q.ravel()
    b = np.roll(q, -1, axis=1).ravel()
    cell = np.repeat(np.arange(M, dtype=np.int64), 4)
    key = np.minimum(a, b) * N + np.maximum(a, b)
    order = np.argsort(key, kind="stable")          # 同じ辺の中では半辺の番号順 (= セル番号順)
    ks = key[order]
    start = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    cnt = np.diff(np.r_[start, len(ks)])
    if cnt.max() > 2:
        k = int(ks[start[np.argmax(cnt)]])
        raise ValueError(f"3 枚以上の四角形が共有する辺がある (節点 {k // N}, {k % N})")
    h0 = order[start[cnt == 2]]
    h1 = order[start[cnt == 2] + 1]
    if not np.all((a[h0] == b[h1]) & (b[h0] == a[h1])):
        raise ValueError("隣り合う四角形が共有する辺を同じ向きにたどっている (向きが不整合)")
    sel = np.argsort(h0, kind="stable")
    h0, h1 = h0[sel], h1[sel]
    interior = np.stack([a[h0], b[h0], cell[h0], cell[h1]], axis=1)
    interior_keys = np.sort(key[h0])
    hb = order[start[cnt == 1]]                       # 外周の半辺 (キー昇順)
    bkeys = key[hb]
    used = np.zeros(len(hb), dtype=bool)
    boundary = {}
    for name in bedges:
        E = _as_edges(bedges[name], name)
        if E.size and (E.min() < 0 or E.max() >= N):
            raise ValueError(f"境界 '{name}' の節点番号が範囲外")
        ke = np.minimum(E[:, 0], E[:, 1]) * N + np.maximum(E[:, 0], E[:, 1])
        pos = np.searchsorted(bkeys, ke)
        posc = np.minimum(pos, max(len(bkeys) - 1, 0))
        hit = (pos < len(bkeys)) & (bkeys[posc] == ke) if len(bkeys) else np.zeros(len(ke), bool)
        if not np.all(hit):
            k = int(np.flatnonzero(~hit)[0])
            j = np.searchsorted(interior_keys, ke[k])
            inner = j < len(interior_keys) and interior_keys[j] == ke[k]
            raise ValueError(f"境界 '{name}' の辺 {k} (節点 {E[k, 0]}, {E[k, 1]}) は"
                             + ("内部の辺 (2 枚の四角形が共有)" if inner else "どの四角形の辺とも一致しない"))
        if np.any(used[pos]) or len(np.unique(pos)) != len(pos):
            raise ValueError(f"境界 '{name}' の辺が重複している (同じ辺が 2 回、または別の境界と重複)")
        used[pos] = True
        h = hb[pos]
        boundary[name] = np.stack([a[h], b[h], cell[h], np.full(len(h), -1, np.int64)], axis=1)
    if not np.all(used):
        k = int(np.flatnonzero(~used)[0])
        h = hb[k]
        raise ValueError(f"どの境界にも属さない外周辺が {int((~used).sum())} 本ある (最初: 節点 {a[h]}, {b[h]}、セル {cell[h]})")
    return {"xy": xy, "n_cells": M, "interior": interior, "boundary": boundary}


def _zone_plan(bedges, zone_types) -> list:
    """境界ゾーンの (名前, 型番号) を ID の順に並べる。"""
    types = dict(DEFAULT_ZONE_TYPES)
    if zone_types:
        types.update({str(k): int(v) for k, v in zone_types.items()})
    names = [n for n in ZONE_ORDER if n in bedges] + sorted(n for n in bedges if n not in ZONE_ORDER)
    out = []
    for n in names:
        if not _NAME_RE.match(n) or n in (CELL_ZONE_NAME, INTERIOR_ZONE_NAME):
            raise ValueError(f"境界の名前 '{n}' は Fluent のゾーン名にできない (英字で始まる英数字・'_'・'-'、"
                             f"'{CELL_ZONE_NAME}'・'{INTERIOR_ZONE_NAME}' は予約)")
        if n not in types:
            raise ValueError(f"境界 '{n}' の Fluent の型が無い (zone_types で与える)")
        if types[n] not in ZONE_TYPE_NAMES or types[n] == 2:
            raise ValueError(f"境界 '{n}' の型 {types[n]} は未対応 (使える型: "
                             f"{ {k: v for k, v in ZONE_TYPE_NAMES.items() if k != 2} })")
        out.append((n, types[n]))
    return out


# ---------------------------------------------------------------------------
# 書き出し
# ---------------------------------------------------------------------------

def write_fluent_msh_2d(path, coords, quads, bedges, zone_types=None,
                        title: str = "forge nozzle 2D axisymmetric mesh (x = axis, y = radius, unit m)") -> dict:
    """Fluent の ASCII メッシュ (2D) を書く。座標 [m] は 17 有効桁、番号は 16 進・1 始まり、境界の c1 = 0。

    ゾーン: 1 = fluid (セル、型 1・要素型 3 = 四角形)、2 = interior-fluid、3 以降は境界 (`wall, inlet, outlet, axis` の順、
    型は既定 `DEFAULT_ZONE_TYPES`、`zone_types` で上書き)。面の向きは `build_faces_2d` (c0 が n0→n1 の左側)。
    戻り: {"nodes", "cells", "faces", "zones": {名前: {"id", "type", "faces"}}}。
    """
    if not title.isascii() or '"' in title:
        raise ValueError("title は ASCII で '\"' を含まないこと (Fluent の注記節)")
    F = build_faces_2d(coords, quads, bedges)
    plan = _zone_plan(bedges, zone_types)
    xy = F["xy"]
    N, M = len(xy), F["n_cells"]
    zones = [(2, INTERIOR_ZONE_NAME, 2, F["interior"])]
    for zid, (name, zt) in enumerate(plan, start=3):
        if len(F["boundary"][name]) == 0:
            raise ValueError(f"境界 '{name}' に辺が無い (空の面ゾーンは書けない)")
        zones.append((zid, name, zt, F["boundary"][name]))
    nf = sum(len(z[3]) for z in zones)
    out = [f'(0 "{title}")', "(2 2)", f"(10 (0 1 {N:x} 0 2))", f"(10 (1 1 {N:x} 1 2)("]
    out.append("\n".join(f"{x:.16e} {y:.16e}" for x, y in xy.tolist()))
    out += ["))", f"(12 (0 1 {M:x} 0))", f"(12 (1 1 {M:x} 1 3))", f"(13 (0 1 {nf:x} 0))"]
    first = 1
    summary = {}
    for zid, name, zt, fc in zones:
        last = first + len(fc) - 1
        out.append(f"(13 ({zid:x} {first:x} {last:x} {zt:x} 2)(")
        f1 = fc + 1                                   # 1 始まり; 境界の c1 = -1 → 0
        out.append("\n".join(f"{n0:x} {n1:x} {c0:x} {c1:x}" for n0, n1, c0, c1 in f1.tolist()))
        out.append("))")
        summary[name] = {"id": zid, "type": ZONE_TYPE_NAMES[zt], "faces": int(len(fc))}
        first = last + 1
    out.append(f"(45 (1 fluid {CELL_ZONE_NAME})())")
    for zid, name, zt, _ in zones:
        out.append(f"(45 ({zid:d} {ZONE_TYPE_NAMES[zt]} {name})())")
    Path(path).write_text("\n".join(out) + "\n", encoding="ascii")
    return {"nodes": N, "cells": M, "faces": nf, "zones": summary}


# ---------------------------------------------------------------------------
# 読み込み (ASCII の 2・10・12・13・39/45 節だけ)
# ---------------------------------------------------------------------------

def _hex_ints(s: str) -> np.ndarray:
    return np.fromiter((int(t, 16) for t in s.split()), dtype=np.int64)


def read_fluent_ascii(path) -> dict:
    """Fluent の ASCII メッシュを読む (検査用)。最上位の節を括弧の対応で切り出し、2・10・12・13・39/45 節を読む。
    バイナリの節 (2010 など) は例外。注記 (0・1 節) とほかの節は読み飛ばす。

    戻り: {"dim", "decl" {"nodes","cells","faces"}, "xy" (最大節点番号+1, nd; 未定義は NaN), "node_defined" (bool),
    "cell_zones" [{"id","first","last","type","etype"}], "face_zones" [{"id","first","last","bctype","ftype","faces" (F,4)}],
    "zone_names" {id: (型名, 名前)}}。面の配列は (n0, n1, c0, c1) で 1 始まりのまま (境界の c1 = 0)。
    """
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    secs = []
    depth, inq, s0 = 0, False, -1
    for m in re.finditer(r'[()"]', txt):            # 数値の並びには括弧が無いので一致は節の数ほど
        ch = m.group()
        if inq:
            inq = ch != '"'
            continue
        if ch == '"':
            inq = True
        elif ch == "(":
            if depth == 0:
                s0 = m.start()
            depth += 1
        else:
            depth -= 1
            if depth < 0:
                raise ValueError(f"{path}: 括弧の対応が取れない (位置 {m.start()})")
            if depth == 0:
                secs.append((s0, m.start()))
    if depth != 0:
        raise ValueError(f"{path}: 括弧が閉じていない")
    out = {"dim": None, "decl": {}, "node_zones": [], "cell_zones": [], "face_zones": [], "zone_names": {}}
    node_blocks = []
    body_re = re.compile(r"\s*(\d+)\s*\(([^()]*)\)\s*(?:\(([^()]*)\))?\s*$", re.S)
    for s, e in secs:
        body = txt[s + 1:e]
        mi = re.match(r"\s*(\d+)", body)
        if not mi:
            continue
        idx = int(mi.group(1))
        if idx in (0, 1):
            continue
        if idx == 2:
            out["dim"] = int(body.split()[1])
            continue
        if idx in (39, 45):
            mh = re.match(r"\s*\d+\s*\(([^()]*)\)", body)
            f = mh.group(1).split()
            out["zone_names"][int(f[0])] = (f[1], f[2])
            continue
        if idx >= 2000 and idx % 1000 in (10, 12, 13):
            raise ValueError(f"{path}: バイナリの節 ({idx}) は未対応 (ASCII だけを読む)")
        if idx not in (10, 12, 13):
            continue
        mb = body_re.match(body)
        if not mb:
            raise ValueError(f"{path}: 節 {idx} の書式を読めない: {body[:80]!r}")
        hdr = [int(t, 16) for t in mb.group(2).split()]
        data = mb.group(3)
        zid, first, last = hdr[0], hdr[1], hdr[2]
        if zid == 0:
            out["decl"][{10: "nodes", 12: "cells", 13: "faces"}[idx]] = last - first + 1
            continue
        n = last - first + 1
        if idx == 10:
            nd = hdr[4] if len(hdr) > 4 else (out["dim"] or 2)
            v = np.array(data.split(), dtype=float).reshape(-1, nd) if data is not None else np.zeros((0, nd))
            if len(v) != n:
                raise ValueError(f"{path}: 節点ゾーン {zid} の節点数 {len(v)} ≠ 見出し {n}")
            node_blocks.append((first, v))
            out["node_zones"].append({"id": zid, "first": first, "last": last, "type": hdr[3], "nd": nd})
        elif idx == 12:
            out["cell_zones"].append({"id": zid, "first": first, "last": last, "type": hdr[3],
                                      "etype": hdr[4] if len(hdr) > 4 else None})
        else:
            ft = hdr[4] if len(hdr) > 4 else 0
            tok = _hex_ints(data or "")
            if ft == 0:                               # 混在: 各面が節点数から始まる
                rows, i = [], 0
                while i < len(tok):
                    k = int(tok[i])
                    if k != 2:
                        raise ValueError(f"{path}: 面ゾーン {zid} に節点数 {k} の面 (2D は 2 だけ)")
                    rows.append(tok[i + 1:i + 5])
                    i += 5
                faces = np.array(rows, dtype=np.int64).reshape(-1, 4)
            elif ft == 2:
                faces = tok.reshape(-1, 4)
            else:
                raise ValueError(f"{path}: 面ゾーン {zid} の面型 {ft} は 2D では未対応 (0 か 2)")
            if len(faces) != n:
                raise ValueError(f"{path}: 面ゾーン {zid} の面数 {len(faces)} ≠ 見出し {n}")
            out["face_zones"].append({"id": zid, "first": first, "last": last, "bctype": hdr[3], "ftype": ft,
                                      "faces": faces})
    if out["dim"] is not None and out["dim"] != 2:
        raise ValueError(f"{path}: 次元 {out['dim']} (2D だけを読む)")
    nmax = max((f + len(v) - 1 for f, v in node_blocks), default=0)
    nd = max((v.shape[1] for _, v in node_blocks), default=2)
    xy = np.full((nmax + 1, nd), np.nan)
    defined = np.zeros(nmax + 1, dtype=bool)
    for f, v in node_blocks:
        if np.any(defined[f:f + len(v)]):
            raise ValueError(f"{path}: 節点番号が重複して定義されている")
        xy[f:f + len(v)] = v
        defined[f:f + len(v)] = True
    out["xy"] = xy
    out["node_defined"] = defined
    return out


# ---------------------------------------------------------------------------
# 向きの測定 (書き出し器の規則を使わず幾何から測る)
# ---------------------------------------------------------------------------

def _all_faces(rd: dict) -> tuple:
    if not rd["face_zones"]:
        raise ValueError("面ゾーンが無い")
    faces = np.concatenate([z["faces"] for z in rd["face_zones"]])
    zid = np.concatenate([np.full(len(z["faces"]), z["id"]) for z in rd["face_zones"]])
    return faces, zid


def _cell_centroids(xy: np.ndarray, faces: np.ndarray) -> tuple:
    """セルの重心 = そのセルの面に現れる節点 (重複を除く) の平均。戻り (重心 (最大セル番号+1, 2), 節点数)。"""
    n0, n1, c0, c1 = faces.T
    c = np.concatenate([c0, c0, c1, c1])
    n = np.concatenate([n0, n1, n0, n1])
    m = c > 0
    c, n = c[m], n[m]
    nn = len(xy)
    u = np.unique(c * nn + n)
    uc, un = u // nn, u % nn
    nc = int(uc.max()) + 1
    cnt = np.bincount(uc, minlength=nc)
    cen = np.full((nc, 2), np.nan)
    ok = cnt > 0
    for k in range(2):
        s = np.bincount(uc, weights=xy[un, k], minlength=nc)
        cen[ok, k] = s[ok] / cnt[ok]
    return cen, cnt


def measure_orientation_2d(path) -> dict:
    """Fluent の 2D ASCII メッシュで、各面の c0 が n0→n1 の進行方向の左か右かを幾何から数える (plan §6 V4(a))。

    セル重心は面の節点の平均。面の中点から c0 の重心へのベクトルと、進行方向の右向きの法線 (t_y, −t_x) の内積の符号で
    左右を決める。内部面は c1 が c0 と反対側にあるかも数える。面型 0 (混在) も読む。
    戻り: {"nodes", "cells", "faces", "interior" {"left","right","zero","c1_opposite"}, "boundary" {"left","right","zero"},
    "c0_missing" (c0 = 0 の面), "left_fraction" (全面に対する左の割合)}。
    """
    rd = read_fluent_ascii(path)
    xy = rd["xy"][:, :2]
    faces, _ = _all_faces(rd)
    n0, n1, c0, c1 = faces.T
    if np.any(~rd["node_defined"][np.concatenate([n0, n1])]):
        raise ValueError(f"{path}: 面が未定義の節点を参照している")
    cen, _ = _cell_centroids(xy, faces)
    p0, p1 = xy[n0], xy[n1]
    t = p1 - p0
    nr = np.stack([t[:, 1], -t[:, 0]], axis=1)           # 進行方向の右向き
    mid = 0.5 * (p0 + p1)
    has0 = c0 > 0
    s0 = np.full(len(faces), np.nan)
    s0[has0] = np.einsum("ij,ij->i", cen[c0[has0]] - mid[has0], nr[has0])
    inner = has0 & (c1 > 0)
    bnd = has0 & (c1 == 0)
    s1 = np.einsum("ij,ij->i", cen[c1[inner]] - mid[inner], nr[inner])

    def tally(s):
        return {"left": int(np.sum(s < 0)), "right": int(np.sum(s > 0)), "zero": int(np.sum(s == 0))}

    ti = tally(s0[inner])
    ti["c1_opposite"] = int(np.sum(s0[inner] * s1 < 0))
    tb = tally(s0[bnd])
    nf = len(faces)
    return {"nodes": int(rd["node_defined"].sum()), "cells": int(np.sum(~np.isnan(cen[:, 0]))), "faces": int(nf),
            "interior": ti, "boundary": tb, "c0_missing": int(np.sum(~has0)),
            "left_fraction": (ti["left"] + tb["left"]) / nf if nf else float("nan")}


# ---------------------------------------------------------------------------
# 書いたファイルの読み直し検査
# ---------------------------------------------------------------------------

def _directed_cell_edges(faces: np.ndarray) -> tuple:
    """セルごとの有向辺 (セルが左になる向き): c0 の面は n0→n1、c1 の面は n1→n0。0 始まりで返す。"""
    n0, n1, c0, c1 = (faces - 1).T
    m1 = c1 >= 0
    cell = np.concatenate([c0, c1[m1]])
    a = np.concatenate([n0, n1[m1]])
    b = np.concatenate([n1, n0[m1]])
    return cell, a, b


def check_fluent_msh_2d(path, coords, quads, bedges, zone_types=None) -> dict:
    """書いた Fluent メッシュを読み直して元のメッシュと照合する (plan §6 V4(a))。

    (1) 節点数・セル数・ゾーンごとの面数 (内部面は四角形から独立に数えた共有辺の数、境界は `bedges` の辺数) と
    宣言 (0 番ゾーン) の数、ゾーン名と型 (45 節と 13 節の見出しの両方)。(2) 座標が元と相対 1e-15 以内。
    (3) `measure_orientation_2d` で左が 100 %。(4) 各セルを c0/c1 の面から組んだ有向の多角形が元の四角形の反時計回りの辺と
    一致し、面積が全て正、総和が元の四角形の面積の総和と相対 1e-12。(5) 境界面が 1 本の閉路で、囲む面積が同じ総和と相対 1e-12。
    戻り: 各量と "failures" (不合格の説明の並び)・"ok"。
    """
    fails = []
    res = {"path": str(path)}
    xy0 = _planar_xy(coords)
    q = np.asarray(quads, dtype=np.int64)
    N, M = len(xy0), len(q)
    rd = read_fluent_ascii(path)
    # (1) 数
    nn = int(rd["node_defined"].sum())
    nc = sum(z["last"] - z["first"] + 1 for z in rd["cell_zones"])
    res.update(nodes=nn, cells=nc)
    if nn != N or len(rd["node_defined"]) != N + 1 or not np.all(rd["node_defined"][1:]):
        fails.append(f"節点数 {nn} ≠ 元 {N} (または番号が 1..N でない)")
    if nc != M:
        fails.append(f"セル数 {nc} ≠ 元 {M}")
    faces, _ = _all_faces(rd)
    nf = len(faces)
    for k, v in (("nodes", N), ("cells", M), ("faces", nf)):
        if rd["decl"].get(k) != v:
            fails.append(f"宣言の {k} {rd['decl'].get(k)} ≠ {v}")
    a = q.ravel()
    b = np.roll(q, -1, axis=1).ravel()
    ek = np.minimum(a, b) * N + np.maximum(a, b)
    _, ecnt = np.unique(ek, return_counts=True)
    expect = {INTERIOR_ZONE_NAME: (int(np.sum(ecnt == 2)), "interior")}
    types = dict(DEFAULT_ZONE_TYPES)
    if zone_types:
        types.update(zone_types)
    for name in bedges:
        expect[name] = (len(_as_edges(bedges[name], name)), ZONE_TYPE_NAMES.get(types.get(name), "?"))
    names = rd["zone_names"]
    zone_counts = {}
    for z in rd["face_zones"]:
        tname, zname = names.get(z["id"], ("?", f"(名前なし {z['id']})"))
        zone_counts[zname] = len(z["faces"])
        if ZONE_TYPE_NAMES.get(z["bctype"]) != tname:
            fails.append(f"ゾーン {zname}: 13 節の型 {z['bctype']} と 45 節の型 '{tname}' が食い違う")
        if zname not in expect:
            fails.append(f"予期しない面ゾーン '{zname}'")
        elif expect[zname] != (len(z["faces"]), tname):
            fails.append(f"ゾーン {zname}: (面数, 型) {(len(z['faces']), tname)} ≠ 期待 {expect[zname]}")
    for zname in expect:
        if zname not in zone_counts:
            fails.append(f"面ゾーン '{zname}' が無い")
    cz = [names.get(z["id"]) for z in rd["cell_zones"]]
    if cz != [("fluid", CELL_ZONE_NAME)] or [z["etype"] for z in rd["cell_zones"]] != [3]:
        fails.append(f"セルゾーン {cz} (要素型 {[z['etype'] for z in rd['cell_zones']]}) ≠ 期待 [('fluid', '{CELL_ZONE_NAME}')] 型 3")
    res["zone_faces"] = zone_counts
    if fails:                                       # 数が合わなければ以降の照合は意味がない
        res.update(failures=fails, ok=False)
        return res
    # (2) 座標
    xy = rd["xy"][1:, :2]
    den = np.maximum(np.abs(xy0), np.finfo(float).tiny)
    rel = float(np.max(np.abs(xy - xy0) / den))
    res["coord_max_rel"] = rel
    if not rel <= 1e-15:
        fails.append(f"座標の相対差 {rel:.3e} > 1e-15")
    # (3) 向き (幾何から)
    ori = measure_orientation_2d(path)
    res["orientation"] = ori
    if not (ori["left_fraction"] == 1.0 and ori["interior"]["c1_opposite"] == ori["interior"]["left"] + ori["interior"]["right"]
            + ori["interior"]["zero"] and ori["c0_missing"] == 0):
        fails.append(f"c0 が左の割合 {ori['left_fraction']:.6f} (内部 {ori['interior']}、境界 {ori['boundary']})")
    # (4) セルの多角形
    cell, ea, eb = _directed_cell_edges(faces)
    if cell.min() < 0 or cell.max() >= M:
        fails.append("面のセル番号が範囲外")
        res.update(failures=fails, ok=False)
        return res
    got = np.stack([cell, ea, eb], axis=1)
    got = got[np.lexsort((got[:, 2], got[:, 1], got[:, 0]))]
    ref = np.stack([np.repeat(np.arange(M), 4), a, b], axis=1)
    ref = ref[np.lexsort((ref[:, 2], ref[:, 1], ref[:, 0]))]
    same_poly = got.shape == ref.shape and bool(np.all(got == ref))
    res["cell_polygons_match"] = same_poly
    if not same_poly:
        fails.append("面から組んだセルの有向の辺が元の四角形の反時計回りの辺と一致しない")
    cnt = np.bincount(cell, minlength=M)
    ref_pt = np.stack([np.bincount(cell, weights=xy[ea, k], minlength=M) / np.maximum(cnt, 1) for k in range(2)], axis=1)
    pa, pb = xy[ea] - ref_pt[cell], xy[eb] - ref_pt[cell]
    cell_area = 0.5 * np.bincount(cell, weights=pa[:, 0] * pb[:, 1] - pb[:, 0] * pa[:, 1], minlength=M)
    total0 = float(np.sum(_quad_areas(xy0, q)))
    total = float(np.sum(cell_area))
    res.update(cell_area_min=float(cell_area.min()), n_cell_area_nonpos=int(np.sum(~(cell_area > 0))),
               area_total=total, area_total_ref=total0, area_total_rel=abs(total - total0) / abs(total0))
    if res["n_cell_area_nonpos"]:
        fails.append(f"面積が正でないセル {res['n_cell_area_nonpos']} 個")
    if not res["area_total_rel"] <= 1e-12:
        fails.append(f"セルの面積の総和の相対差 {res['area_total_rel']:.3e} > 1e-12")
    # (5) 境界の閉路 (c0 が左 → 領域を反時計回り)
    bf = faces[faces[:, 3] == 0]
    s, t = bf[:, 0], bf[:, 1]
    loop_ok = len(np.unique(s)) == len(s) and len(np.unique(t)) == len(t) and set(s.tolist()) == set(t.tolist())
    n_loop = 0
    if loop_ok:
        nxt = dict(zip(s.tolist(), t.tolist()))
        p = int(s[0])
        while True:
            p = nxt[p]
            n_loop += 1
            if p == int(s[0]) or n_loop > len(s):
                break
        loop_ok = n_loop == len(s)
    res.update(boundary_faces=int(len(bf)), boundary_loop_len=n_loop, boundary_single_loop=bool(loop_ok))
    if not loop_ok:
        fails.append(f"境界面が 1 本の閉路でない (面 {len(bf)}、たどれた長さ {n_loop})")
    else:
        o = xy[np.unique(s - 1)].mean(axis=0)
        u, v = xy[s - 1] - o, xy[t - 1] - o
        la = 0.5 * float(np.sum(u[:, 0] * v[:, 1] - v[:, 0] * u[:, 1]))
        res.update(boundary_loop_area=la, boundary_loop_area_rel=abs(la - total0) / abs(total0))
        if not res["boundary_loop_area_rel"] <= 1e-12:
            fails.append(f"境界の閉路が囲む面積の相対差 {res['boundary_loop_area_rel']:.3e} > 1e-12")
    res.update(failures=fails, ok=not fails)
    return res


# ---------------------------------------------------------------------------
# gmsh 経由の CGNS
# ---------------------------------------------------------------------------

_GMSH_CGNS_SCRIPT = r'''
import json, sys, gmsh
msh, cgns = sys.argv[1], sys.argv[2]
gmsh.initialize(["", "-noenv"])
gmsh.option.setNumber("General.Terminal", 0)

def stats():
    nt, _, _ = gmsh.model.mesh.getNodes()
    et, tags, _ = gmsh.model.mesh.getElements()
    elems = {gmsh.model.mesh.getElementProperties(int(t))[0]: int(len(g)) for t, g in zip(et, tags)}
    groups = []
    for d, t in gmsh.model.getPhysicalGroups():
        n = 0
        for e in gmsh.model.getEntitiesForPhysicalGroup(d, t):
            _, tg, _ = gmsh.model.mesh.getElements(d, int(e))
            n += sum(len(g) for g in tg)
        groups.append([int(d), int(t), gmsh.model.getPhysicalName(d, t), int(n)])
    return {"nodes": int(len(nt)), "elements": elems, "physical": sorted(groups, key=lambda g: (g[0], g[2]))}

gmsh.open(msh)
a = stats()
gmsh.write(cgns)
gmsh.clear()
gmsh.open(cgns)
b = stats()
gmsh.finalize()
print("FORGE_JSON " + json.dumps({"msh": a, "cgns": b}))
'''


def write_cgns_via_gmsh(msh_path, cgns_path, gmsh_python, timeout: float = 600.0) -> dict:
    """gmsh (`.venv-mesh` の python を `gmsh_python` で渡す) で msh を開いて CGNS を書き、CGNS を読み直す。
    戻り: {"msh": {"nodes","elements" {型名: 数},"physical" [[次元, タグ, 名前, 要素数]]}, "cgns": 同じ, "match" (節点数・
    要素数・物理名と要素数が一致), "tags_preserved" (物理タグの番号まで一致したか), "seconds"}。
    gmsh は CGNS を読み直すと物理タグの番号を振り直す (名前と要素数は保たれる) ので、照合は名前で行う。
    境界の Fluent の型は CGNS では付かない (Fluent 側で付け直す; plan §4.4)。"""
    t0 = time.perf_counter()
    p = subprocess.run([str(gmsh_python), "-I", "-c", _GMSH_CGNS_SCRIPT, str(Path(msh_path).resolve()),
                        str(Path(cgns_path).resolve())], capture_output=True, text=True, timeout=timeout)
    line = [s for s in p.stdout.splitlines() if s.startswith("FORGE_JSON ")]
    if p.returncode != 0 or not line:
        raise RuntimeError(f"gmsh の CGNS 書き出しが失敗 (終了 {p.returncode}): {p.stderr.strip()[-2000:]}")
    out = json.loads(line[-1][len("FORGE_JSON "):])
    a, b = out["msh"], out["cgns"]
    out["match"] = (a["nodes"] == b["nodes"] and a["elements"] == b["elements"]
                    and [g[0::2] + [g[3]] for g in a["physical"]] == [g[0::2] + [g[3]] for g in b["physical"]])
    out["tags_preserved"] = a["physical"] == b["physical"]
    out["seconds"] = time.perf_counter() - t0
    out["cgns_path"] = str(cgns_path)
    return out


# ---------------------------------------------------------------------------
# 境界条件表
# ---------------------------------------------------------------------------

def _load(obj, kind):
    if obj is None:
        return {}
    if isinstance(obj, dict):
        return obj
    p = Path(obj)
    if kind == "json":
        return json.loads(p.read_text())
    import yaml
    return yaml.safe_load(p.read_text()) or {}


def _get(d, *keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _num(v) -> str:
    return repr(float(v))


def _isentropic_subsonic_p_ratio(area_ratio: float, gamma: float) -> tuple:
    """亜音速の等エントロピー流で A/A* から (M, p/Pt) を解く (一定比熱)。"""
    g = gamma

    def f(m):
        return (1.0 / m) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * m * m)) ** ((g + 1.0) / (2.0 * (g - 1.0))) - area_ratio

    lo, hi = 1e-12, 1.0                               # f(lo) > 0, f(1) = 1 − A/A* < 0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    m = 0.5 * (lo + hi)
    return m, (1.0 + 0.5 * (g - 1.0) * m * m) ** (-g / (g - 1.0))


def _inlet_static_guess(problem: dict, prep: dict) -> tuple:
    """入口の等エントロピー静圧の目安 (p/Pt, 説明) か (None, 理由)。入口半径は問題 YAML の geometry.r_inlet、
    スロートは prepare_info の物理スロートの半径 (どちらも r_t 単位)、γ は問題 YAML の gas.gamma (一定比熱の近似)。"""
    r_in = _get(problem, "geometry", "r_inlet")
    r_th = _get(prep, "physical_wall", "throat", "r")
    if r_th is None:
        r_th = _get(prep, "throat_physical", "r")
    gamma = _get(problem, "gas", "gamma")
    miss = [n for n, v in (("問題 YAML の geometry.r_inlet", r_in), ("prepare_info の物理スロートの半径", r_th),
                           ("問題 YAML の gas.gamma", gamma)) if v is None]
    if miss:
        return None, "目安なし (" + "・".join(miss) + " が無い)"
    ar = (float(r_in) / float(r_th)) ** 2
    g = float(gamma)
    if not (ar > 1.0 and g > 1.0):
        return None, f"目安なし (A/A* = {ar:.6g}、γ = {g:.6g} では解けない)"
    m, pr = _isentropic_subsonic_p_ratio(ar, g)
    return pr, (f"等エントロピーの目安: A/A* = (r_inlet / r_throat,phys)² = ({float(r_in):.6g} / {float(r_th):.8g})² = {ar:.6g}、"
                f"γ = {g:.6g} 一定 (gas.gamma の参照値) で M = {m:.4g}、p/Pt = {pr:.8f}")


def write_fluent_bc_table(path, problem_yaml_or_dict, bcond_yaml, prepare_info=None, zone_types=None,
                          species_names=None, msh_physical=None) -> dict:
    """Fluent 用の境界条件表 (Markdown、日本語) を書く (plan §4.4・§6 V4(d))。

    値は `bcondConfig.yaml` (圧力・温度・k/ω・化学種) と問題 YAML・`prepare_info.json` (ガス・物性・寸法) から写し、
    無いものは「不明」と書く (推測で補わない)。Operating Pressure は 0 Pa とし、各圧力は「絶対値」と「Fluent の入力欄の値
    (= 絶対値 − 動作圧)」を併記する (codex plan レビュー M5)。`species_names` を渡すと bcond の Y0.. に種名を付ける。
    `msh_physical` (msh の物理名 → タグ、`read_msh41_mesh2d` の 4 番目の戻り値) を渡すと、forge が境界を対応付ける physID から
    ゾーン名 (= msh の物理名 = Fluent のゾーン名) を引く (bcond のキー名と違えば注記する)。
    戻り: {"path", "rows" [{"zone","fluent_type","item","abs","input","unit","source"}], "operating_pressure_pa"}。
    """
    problem = _load(problem_yaml_or_dict, "yaml")
    bc = _load(bcond_yaml, "yaml")
    prep = _load(prepare_info, "json")
    types = dict(DEFAULT_ZONE_TYPES)
    if zone_types:
        types.update(zone_types)
    op = OPERATING_PRESSURE_PA
    bsrc = "`bcondConfig.yaml`"
    rows = []
    notes = []

    def row(zone, item, value, unit, source, pressure=False):
        if value is None:
            rows.append({"zone": zone, "fluent_type": ftype, "item": item, "abs": "不明", "input": "不明", "unit": unit,
                         "source": source})
            return
        v = float(value)
        rows.append({"zone": zone, "fluent_type": ftype, "item": item, "abs": _num(v),
                     "input": _num(v - op) if pressure else _num(v), "unit": unit, "source": source})

    order = [n for n in ZONE_ORDER if n in bc] + sorted(n for n in bc if n not in ZONE_ORDER)
    wt = _get(problem, "spec", "wall_thermal") or _get(prep, "wall_thermal")
    tag2name = {int(v): str(k) for k, v in msh_physical.items()} if msh_physical else None
    for key in order:
        ent = bc[key] or {}
        zone = key
        if tag2name is not None:
            pid = ent.get("physID")
            zone = tag2name.get(int(pid)) if pid is not None else None
            if zone is None:
                zone = f"不明 (physID {pid} が msh に無い)"
            elif zone != key:
                notes.append(f"bcond の `{key}` (physID {pid}) は msh のゾーン `{zone}`。")
        kind = str(ent.get("kind"))
        fl = ent.get("floats") or {}
        ftype = ZONE_TYPE_NAMES.get(types.get(zone), "不明")
        src = f"{bsrc} {key}.floats"
        if kind in ("inlet_Pressure", "inlet_Pressure_dir"):
            row(zone, "Gauge Total Pressure (全圧 Pt)", fl.get("Pt"), "Pa", f"{src}.Pt", pressure=True)
            pr, why = _inlet_static_guess(problem, prep)
            if pr is not None and fl.get("Pt") is not None:
                row(zone, "Supersonic/Initial Gauge Pressure (初期ゲージ静圧の目安)", pr * float(fl["Pt"]), "Pa", why, pressure=True)
            else:
                rows.append({"zone": zone, "fluent_type": ftype, "item": "Supersonic/Initial Gauge Pressure (初期ゲージ静圧の目安)",
                             "abs": "目安なし", "input": "目安なし", "unit": "Pa",
                             "source": why if pr is None else "目安なし (bcond の Pt が無い)"})
            row(zone, "Total Temperature (全温 Tt)", fl.get("Tt"), "K", f"{src}.Tt")
            if "k" in fl:
                row(zone, "Turbulent Kinetic Energy (k)", fl["k"], "m²/s²", f"{src}.k")
            if "omega" in fl:
                row(zone, "Specific Dissipation Rate (ω)", fl["omega"], "1/s", f"{src}.omega")
            ys = sorted((k for k in fl if re.fullmatch(r"Y\d+", str(k))), key=lambda k: int(k[1:]))
            for k in ys:
                nm = species_names[int(k[1:])] if species_names and int(k[1:]) < len(species_names) else None
                row(zone, f"Species Mass Fraction {k}" + (f" ({nm})" if nm else " (種名は run の solverConfig.yaml の並び)"),
                    fl[k], "-", f"{src}.{k}")
            if kind == "inlet_Pressure":
                notes.append(f"`{zone}`: forge の `inlet_Pressure` は速度の向きを内部から外挿する (`methods/boundary.md`)。"
                             "Fluent の pressure-inlet の向きの指定 (既定は Normal to Boundary) との違いに注意。")
        elif kind == "outlet_statPress":
            row(zone, "Gauge Pressure (背圧 Ps)", fl.get("Ps"), "Pa", f"{src}.Ps", pressure=True)
            if "Pt" in fl:
                row(zone, "Backflow Total Pressure (逆流時の全圧; forge の逆流条件の Pt)", fl["Pt"], "Pa", f"{src}.Pt", pressure=True)
            if "Tt" in fl:
                row(zone, "Backflow Total Temperature (逆流時の全温; forge の逆流条件の Tt)", fl["Tt"], "K", f"{src}.Tt")
        elif kind == "wall":
            rows.append({"zone": zone, "fluent_type": ftype, "item": "Thermal: Heat Flux (断熱壁、滑りなし)", "abs": "0.0",
                         "input": "0.0", "unit": "W/m²", "source": f"{bsrc} {key}.kind = wall (断熱)"})
        elif kind == "wall_isothermal":
            row(zone, "Thermal: Temperature (等温壁 Ts、滑りなし)", fl.get("Ts"), "K", f"{src}.Ts")
        elif kind == "axis":
            rows.append({"zone": zone, "fluent_type": ftype, "item": "(入力なし)", "abs": "—", "input": "—", "unit": "—",
                         "source": f"{bsrc} {key}.kind = axis"})
        elif kind == "slip":
            rows.append({"zone": zone, "fluent_type": ftype, "item": "Shear Condition: Specified Shear 0 (滑り壁)", "abs": "0.0",
                         "input": "0.0", "unit": "Pa", "source": f"{bsrc} {key}.kind = slip"})
        else:
            rows.append({"zone": zone, "fluent_type": ftype, "item": f"不明 (forge の kind `{kind}` の対応はこの表では決めない)",
                         "abs": "不明", "input": "不明", "unit": "—", "source": f"{bsrc} {key}.kind"})
        if kind in ("wall", "wall_isothermal") and isinstance(wt, dict):
            want = "isothermal" if kind == "wall_isothermal" else "adiabatic"
            if str(wt.get("mode")) != want:
                notes.append(f"**注意**: `{key}` の bcond は `{kind}` だが、問題/prepare の wall_thermal は {wt} (食い違い)。")

    L = ["# Fluent 用の境界条件表 (forge の設定から写したもの)", ""]
    L += ["この表は forge の `bcondConfig.yaml`・問題 YAML・`prepare_info.json` の値を写したもの。値が無いものは「不明」とし、推測で補っていない。",
          "**Fluent でこの表と `.msh` を読んだ確認はまだしていない** (plan §5.1 #10、ユーザ確認待ち)。", ""]
    L += ["## 解析の設定", "", "| 項目 | 値 |", "| --- | --- |",
          "| 空間 | 2D 軸対称 (General > 2D Space = Axisymmetric) |",
          "| 対称軸 | $x$ 軸 (半径方向 = $y$、$y \\ge 0$ の上半分) |",
          "| 長さの単位 | m (メッシュは m で書いてある。Mesh > Scale は不要) |",
          f"| Operating Pressure | {_num(op)} Pa (Operating Conditions で 0 にする) |",
          "| 圧力の入力欄 | Fluent の圧力は動作圧に対するゲージ値。下表の「Fluent の入力欄の値」= 絶対値 − 動作圧 (動作圧 0 Pa なので絶対値と同じ) |", ""]
    L += ["## 境界ごとの値", "", "| ゾーン名 | Fluent の型 | 項目 (Fluent の入力欄) | 絶対値 | Fluent の入力欄の値 | 単位 | 出典 |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        L.append(f"| {r['zone']} | {r['fluent_type']} | {r['item']} | {r['abs']} | {r['input']} | {r['unit']} | {r['source']} |")
    L.append("")
    if notes:
        L += ["### 境界の注記", ""] + [f"- {n}" for n in notes] + [""]
    # ガスと物性
    L += ["## ガスと物性 (forge 側)", ""]
    species = _get(problem, "gas", "species")
    basis = str(_get(problem, "gas", "composition_basis") or "mass")
    L += ["| 項目 | 値 | 出典 |", "| --- | --- | --- |"]
    if isinstance(species, dict) and species:
        lab = "質量分率" if basis.lower() == "mass" else f"{basis} 分率 (入力のまま)"
        for k, v in species.items():
            L.append(f"| 組成 {k} ({lab}) | {_num(v)} | 問題 YAML gas.species |")
    else:
        L.append("| 組成 | 不明 | 問題 YAML gas.species が無い |")
    gm = _get(problem, "gas", "model")
    L.append(f"| 気体モデル (forge) | {gm if gm is not None else (_get(prep, 'initializer', 'gas') or '不明')} | "
             + ("問題 YAML gas.model" if gm is not None else "prepare_info initializer.gas") + " |")
    tr = _get(problem, "gas", "transport")
    if tr is None:
        tr = _get(prep, "transport", "models")
    trs = ", ".join(f"{k}: {v}" for k, v in tr.items()) if isinstance(tr, dict) else (str(tr) if tr is not None else "不明")
    L.append(f"| 輸送物性 (forge) | {trs} | 問題 YAML gas.transport / prepare_info transport.models |")
    vm = _get(prep, "transport", "viscMethod")
    L.append(f"| viscMethod (forge) | {vm if vm is not None else '不明'} | prepare_info transport.viscMethod |")
    for key, lab in (("gamma", "γ (参照値)"), ("cp", "cp (参照値) [J/(kg·K)]")):
        v = _get(problem, "gas", key)
        L.append(f"| {lab} | {_num(v) if v is not None else '不明'} | 問題 YAML gas.{key} |")
    model = gm if gm is not None else _get(prep, "initializer", "gas")
    L.append("")
    if str(model) == "semiperfect":
        L.append("- forge の semiperfect は NASA-9 (CEA) 係数の熱的完全気体・凍結組成。輸送物性は上の種ごとのモデル。")
    L += ["- **Fluent 側の物性の選び方 (比熱の温度依存・粘性・熱伝導率・混合則) はこの表では決めない**。forge と同じにする方法を断定できない"
          "ため、Fluent との差の要因になる。", ""]
    # メッシュの注記
    scale = _get(prep, "scale_m")
    L += ["## メッシュの注記", "",
          "- forge は節点中心 (節点 = 未知数)、Fluent はセル中心。同じ一次メッシュでも Fluent の第 1 未知数 (第 1 セルの中心) の"
          "壁距離は forge の第 1 内部節点の約半分になり、$y_1^+$ も約半分になる。",
          f"- 設計スロート半径 r_t = {_num(scale) + ' m' if scale is not None else '不明'} (prepare_info scale_m)。", ""]
    Path(path).write_text("\n".join(L), encoding="utf-8")
    return {"path": str(path), "rows": rows, "operating_pressure_pa": op}
