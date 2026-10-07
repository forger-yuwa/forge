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

品質 (`check_mesh_quality.py`) と双対幾何 (`check_dual_closure.py`) の VERDICT は本道具では出さない。回し方は出力の末尾に書く。
"""
import argparse
import json
import os
import sys

import numpy as np

HEX_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]   # gmsh 型 5 の頂点順
_IDX = {c: n for n, c in enumerate(HEX_CORNERS)}
PHYS_NAME = {1: "inlet_nozzle", 2: "inlet_ext", 3: "outlet", 4: "ramp", 5: "cowl_in", 6: "cowl_out", 7: "bottom", 8: "top_out",
             9: "sym", 10: "side_far", 11: "sidewall_in", 12: "sidewall_out", 14: "vehicle", 15: "vehicle_top",
             16: "underside_far", 17: "vehicle_side", 18: "vehicle_base"}   # mesh_sern3d.PHYS_SERN3D
IN_PLANE_TAGS = (9, 10)   # sym (z = 0)・side_far (z = Z_far): 節点が面内 (y) で動くのは形状の変更ではない


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


def load(path, info_path=None):
    """格子と構造情報。戻り = dict(coords [m], hexes, bfaces, info (メッシャ), H, label)。"""
    if os.path.isdir(path):
        grid = os.path.join(path, "sern.h5")
        if not os.path.exists(grid):
            grid = os.path.join(path, "sern.msh")
        info_path = info_path or os.path.join(path, "prepare_info.json")
    else:
        grid = path
    if not info_path or not os.path.exists(info_path):
        raise SystemExit(f"{path}: 構造情報 (prepare_info.json か --info-*) が無い")
    if not os.path.exists(grid):
        raise SystemExit(f"{path}: 格子 (sern.h5 / sern.msh) が無い")
    raw = json.load(open(info_path))
    info = raw["mesh"] if "mesh" in raw else raw
    H = float(raw.get("H_m", info.get("H_m", 1.0)))
    coords, hexes, bf = _read_msh41(grid) if grid.endswith(".msh") else _read_h5(grid)
    if coords.shape[0] != int(info["nodes"]):
        raise SystemExit(f"{grid}: 節点数 {coords.shape[0]} が info の {info['nodes']} と違う (別の格子の info)")
    return {"coords": coords, "hexes": hexes, "bfaces": bf, "info": info, "H": H, "label": grid}


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


# ---------------------------------------------------------------------------------------------- 本体
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--info-a"); ap.add_argument("--info-b")
    ap.add_argument("--n-layers", type=int, default=10, help="中間線から数える層の数 ((i)・(iii))")
    a = ap.parse_args(argv)
    return report(load(a.A, a.info_a), load(a.B, a.info_b), a.n_layers)


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
