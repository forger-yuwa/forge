r"""STEP の壁曲線から gmsh で「同じ構造」の 2D 軸対称メッシュを作り直す単体スクリプト
(plans/active/tooling-nozzle-wt-deliverables.md §4.3・§6 V3)。

forge にも forge_design にも依存しない (標準ライブラリ・numpy・gmsh だけ)。成果物フォルダにそのまま写して配る。

    python rebuild_from_step.py STEP STRUCTURE_JSON OUT_PREFIX [--compare ORIG.msh] [--n-sample 20001]

- STEP: 単位 mm、$x$ 軸 = 流れ方向 (回転軸)、$y$ = 半径 ($r \ge 0$ の上半分)、$z=0$ の平面曲線 1 本
  (forge の `export/wall_step.py` が書く形。CAD で編集した曲線でもよい)。
- STRUCTURE_JSON: `mesh_structure.json` (`forge.mesh_structure/v1`、接続の規則 `mesh2d_ij_v1`)。station $x_i$ [mm]、
  半径方向の比 $s_{ij}$、物理名とタグ。
- 出力 (座標は m):
  - `OUT_PREFIX.msh`: **forge の変換器 (`convertGmshToForge`) に入れる形**。forge の `meshing/mesh2d.py` の `write_msh41_2d` と
    同じ書式 (msh4.1 ASCII、節点ブロック 1 つ・面 1・曲線 4、節点・要素番号は元の msh と同じでファイル順に 1..N) で、座標だけ
    17 有効桁。forge の変換器は節点のタグ行を読まず「ファイル内の順番 = 節点番号」とみなすので、この形でないと接続が壊れる。
  - `OUT_PREFIX_cad.msh`: gmsh が書く msh4.1 (節点を CAD 要素ごとのブロックに分け、壁の節点は曲線に分類)。gmsh・Fluent 側で
    使う。番号は元と同じだがファイル内の順番は番号順でない (**forge の変換器には入れない**)。助変数は書かない。
  - `OUT_PREFIX_v22.msh` (msh2.2 ASCII)・`OUT_PREFIX.cgns`: gmsh が書く。節点・要素の番号 (CGNS は物理タグも) は gmsh が
    振り直す (四角形の並びと座標・物理名は同じ)。
  - `--compare` のとき `OUT_PREFIX_compare.json`。

手順:
1. STEP を `gmsh.model.occ.importShapes` で読む (単位 mm に換算)。曲線がちょうど 1 本 (面・体なし) であることを確かめる。
2. 曲線を助変数について密に (既定 20001 点) 評価し、平面 ($|z|\le10^{-9}$ mm)・$r>0$・$x(u)$ の狭義単調 (標本の差分と
   $dx/du$ の符号; 折り返す = 多価な曲線は拒否)・全 station を含む (端は $10^{-6}$ mm の余裕) を検査する。
3. 各 station で **$x(u)=x_i$ を $u$ について解く** (標本で挟んでから二分法と割線法 (はさみうち) を交互に、
   $|x(u_i)-x_i|\le10^{-9}$ mm まで)。CAD で編集した曲線では $x(u)=u$ が成り立たないため ($u$ をそのまま $x$ に使わない)。
   曲線の端と端の station の差が $10^{-6}$ mm 以内ならその端の助変数を使い、曲線が station より長ければ $u_0$・$u_{n_i-1}$ で切る。
4. 壁半径 $r_w(x_i)=r(u_i)$、節点 $(x_i,\ s_{ij}\,r_w(x_i))$。節点の $x$ は station のまま (縦線を保つ)。
5. 軸・入口・出口の直線と面を OCC で作り (壁は STEP の曲線そのもの)、全体を mm → m に縮める。
   節点を角の点・各曲線の内点 (壁は助変数 $u_i$ つき)・面の内点に、要素を線分 (型 1)・四角形 (型 3) で載せる。
   節点番号は元の msh と同じ $i\,n_j+j+1$、要素番号も元と同じ並び (四角形が先、境界辺は構造ファイルの順)。
6. 物理グループ (名前・タグは構造ファイルどおり) を付けて書く (gmsh で `_cad.msh`・`_v22.msh`・`.cgns`、自前の書き出し器で
   forge の形の `.msh`)。
7. `--compare ORIG.msh`: 元の msh と、書いた `.msh` (forge の形) と `_cad.msh` をそれぞれ gmsh で読み直し、(a) 四角形の接続が
   全要素で一致、(b) 節点の距離の最大 [mm] と各 station の第 1 層厚に対する比の最大、(c) 物理グループの名前・タグ・要素数の
   一致を JSON に書く。加えて (d) `.msh` を直接読んで forge の形 (節点ブロック 1 つ・節点番号と要素番号がファイル順に 1..N) か。
   判定 (plan §6 V3): 両方の msh で最大距離 ≤ 1e-5 mm かつ比 ≤ 0.01、(a)・(c) が一致、(d) が成り立てば `verdict: PASS`。

終了コード: 0 = 成功 (比較ありなら PASS)、1 = 比較が FAIL、2 = 入力を拒否 (理由を標準エラーに書く)、
3 = 形状の組み立て中のその他のエラー (gmsh・OCC)。
曲線の検査は標本 (既定 20001 点) とその点の $dx/du$ で行うので、標本の間隔より短い折り返しは見落としうる
(`--n-sample` で増やせる)。
CAD で曲線を修正したときは、作り直したメッシュの品質を forge の `check_mesh_quality.py` で別に確かめること。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

FORMAT = "forge.mesh_structure/v1"
CONNECTIVITY = "mesh2d_ij_v1"
SIDES = ("inlet", "outlet", "wall", "axis")
PLANE_TOL_MM = 1e-9      # 曲線の |z| の許容
END_TOL_MM = 1e-6        # 曲線の端と端の station の差の許容 (これ以内は曲線の端を使う)
SOLVE_TOL_MM = 1e-9      # |x(u_i) − x_i| の許容
MAX_ITER = 200
MM_TO_M = 1e-3
DIST_TOL_MM = 1e-5       # plan §6 V3 (b)
RATIO_TOL = 0.01         # 第 1 層厚に対する比


class RebuildError(ValueError):
    """入力 (STEP・構造ファイル) を拒否する理由。"""


# --------------------------------------------------------------------------------------------
# 構造ファイル
# --------------------------------------------------------------------------------------------

def load_structure(path) -> dict:
    st = json.loads(Path(path).read_text(encoding="utf-8"))
    if st.get("format") != FORMAT or st.get("connectivity") != CONNECTIVITY:
        raise RebuildError(f"構造ファイルの format/connectivity が {FORMAT}/{CONNECTIVITY} でない "
                           f"({st.get('format')}/{st.get('connectivity')})")
    ni, nj = int(st["ni"]), int(st["nj"])
    x = np.asarray(st["x_mm"], dtype=np.float64)
    if x.shape != (ni,) or not np.all(np.diff(x) > 0):
        raise RebuildError("構造ファイルの x_mm が長さ ni の狭義増加列でない")
    if "s_common" in st:
        S = np.repeat(np.asarray(st["s_common"], dtype=np.float64)[None, :], ni, axis=0)
    else:
        S = np.asarray(st["s"], dtype=np.float64)
    if S.shape != (ni, nj) or not (np.all(S[:, 0] == 0.0) and np.all(S[:, -1] == 1.0) and np.all(np.diff(S, axis=1) > 0)):
        raise RebuildError("構造ファイルの比 s が (ni, nj) で 0 → 1 の狭義増加でない")
    phys = st["physical"]
    surf = [k for k in phys if k not in SIDES]
    if not set(SIDES) <= set(phys) or len(surf) != 1 or phys[surf[0]]["dim"] != 2:
        raise RebuildError(f"構造ファイルの物理名が {SIDES} + 面 1 つでない ({sorted(phys)})")
    order = list(st.get("element_tag_order") or [surf[0], *SIDES])
    if order[0] != surf[0] or sorted(order[1:]) != sorted(SIDES):
        raise RebuildError(f"構造ファイルの element_tag_order が不正 ({order})")
    return {"raw": st, "ni": ni, "nj": nj, "x": x, "S": S, "phys": phys, "surf": surf[0], "order": order}


def connectivity(ni: int, nj: int):
    """規則 mesh2d_ij_v1 の四角形と境界辺 (0 始まり)。"""
    nid = np.arange(ni * nj, dtype=np.int64).reshape(ni, nj)
    quads = np.stack([nid[:-1, :-1], nid[1:, :-1], nid[1:, 1:], nid[:-1, 1:]], axis=-1).reshape(-1, 4)
    bedges = {
        "inlet": np.stack([nid[0, :-1], nid[0, 1:]], axis=-1),
        "outlet": np.stack([nid[-1, :-1], nid[-1, 1:]], axis=-1),
        "wall": np.stack([nid[:-1, -1], nid[1:, -1]], axis=-1),
        "axis": np.stack([nid[:-1, 0], nid[1:, 0]], axis=-1),
    }
    return quads, bedges


# --------------------------------------------------------------------------------------------
# STEP の曲線
# --------------------------------------------------------------------------------------------

def _xyz(gmsh, tag: int, u) -> np.ndarray:
    return np.asarray(gmsh.model.getValue(1, tag, list(map(float, u))), dtype=np.float64).reshape(-1, 3)


def import_curve(gmsh, step: str) -> int:
    gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
    gmsh.model.occ.importShapes(str(step))
    gmsh.model.occ.synchronize()
    c1, c2, c3 = (gmsh.model.getEntities(d) for d in (1, 2, 3))
    if len(c1) != 1 or c2 or c3:
        raise RebuildError(f"STEP の中身が曲線 1 本でない (曲線 {len(c1)}・面 {len(c2)}・体 {len(c3)})")
    return c1[0][1]


def curve_checks(gmsh, tag: int, x_st: np.ndarray, n_sample: int) -> dict:
    """密な評価で平面・r>0・x(u) の狭義単調・station を含むことを検査する。拒否なら RebuildError。"""
    lo, hi = (float(v[0]) for v in gmsh.model.getParametrizationBounds(1, tag))
    u = np.linspace(lo, hi, int(n_sample))
    P = _xyz(gmsh, tag, u)
    D = np.asarray(gmsh.model.getDerivative(1, tag, u.tolist()), dtype=np.float64).reshape(-1, 3)
    zmax = float(np.abs(P[:, 2]).max())
    if not np.isfinite(P).all() or zmax > PLANE_TOL_MM:
        raise RebuildError(f"曲線が z = 0 の平面にない (|z| 最大 {zmax:.3e} mm、許容 {PLANE_TOL_MM:g} mm)")
    if not np.all(P[:, 1] > 0):
        k = int(np.argmin(P[:, 1]))
        raise RebuildError(f"曲線の r が 0 以下になる (u = {u[k]:.9g}・x = {P[k, 0]:.9g} mm で r = {P[k, 1]:.3e} mm)")
    sgn = 1.0 if P[-1, 0] > P[0, 0] else -1.0
    dx = sgn * np.diff(P[:, 0])
    dxdu = sgn * D[:, 0]
    if not (np.all(dx > 0) and np.all(dxdu > 0)):
        bad_dx = np.flatnonzero(~(dx > 0))
        bad_d = np.flatnonzero(~(dxdu > 0))
        k = int(bad_dx[0]) if bad_dx.size else int(bad_d[0])
        raise RebuildError(
            f"x(u) が狭義単調でない (折り返す曲線; x が r の多価関数になる): u = {u[k]:.9g} (x = {P[k, 0]:.9g} mm) 付近で "
            f"向きが反転 (標本の差分の反転 {bad_dx.size} 箇所・dx/du の反転 {bad_d.size} 箇所)")
    xmin, xmax = float(P[:, 0].min()), float(P[:, 0].max())
    short_in, short_out = xmin - x_st[0], x_st[-1] - xmax
    if short_in > END_TOL_MM or short_out > END_TOL_MM:
        raise RebuildError(
            f"曲線の x の範囲 [{xmin:.9g}, {xmax:.9g}] mm が station の範囲 [{x_st[0]:.9g}, {x_st[-1]:.9g}] mm を含まない "
            f"(入口側の不足 {max(short_in, 0):.3e} mm・出口側 {max(short_out, 0):.3e} mm、許容 {END_TOL_MM:g} mm)")
    return {"u": u, "P": P, "sgn": sgn, "param_range": [lo, hi], "x_range_mm": [xmin, xmax],
            "max_abs_z_mm": zmax, "min_r_mm": float(P[:, 1].min()),
            "x_minus_u_max_abs_mm": float(np.abs(P[:, 0] - u).max()), "n_sample": int(n_sample),
            "direction": "x 増加" if sgn > 0 else "x 減少"}


def solve_stations(gmsh, tag: int, x_st: np.ndarray, cc: dict) -> dict:
    """x(u) = x_i を各 station で解く。端の station は曲線の端から END_TOL_MM 以内なら端の助変数。"""
    u, P, sgn = cc["u"], cc["P"], cc["sgn"]
    xs = sgn * P[:, 0]                      # u について狭義増加
    tgt = sgn * x_st
    n = len(x_st)
    u_sol = np.full(n, np.nan)
    end_info = []
    # 曲線の入口側の端 (x 最小) と出口側の端 (x 最大)
    k_in, k_out = (0, len(u) - 1) if sgn > 0 else (len(u) - 1, 0)
    trim = {"inlet": False, "outlet": False}
    for side, i, k in (("inlet", 0, k_in), ("outlet", n - 1, k_out)):
        gap = float(P[k, 0] - x_st[i])
        if abs(gap) <= END_TOL_MM:
            u_sol[i] = u[k]
            end_info.append({"side": side, "station": i, "use": "曲線の端", "curve_end_minus_station_mm": gap})
        else:
            trim[side] = True               # 曲線が station より長い (短い側は curve_checks で拒否済み)
            end_info.append({"side": side, "station": i, "use": "解いて切る", "curve_end_minus_station_mm": gap})
    act = np.flatnonzero(np.isnan(u_sol))
    k = np.clip(np.searchsorted(xs, tgt[act], side="left"), 1, len(u) - 1)
    a, b = u[k - 1].copy(), u[k].copy()
    fa, fb = xs[k - 1] - tgt[act], xs[k] - tgt[act]
    if not (np.all(fa <= 0) and np.all(fb >= 0)):
        raise RebuildError("station を標本で挟めない (内部の不整合)")
    sol = np.full(act.size, np.nan)
    res = np.full(act.size, np.nan)
    m = np.abs(fb) <= SOLVE_TOL_MM
    sol[m], res[m] = b[m], fb[m]
    m = np.isnan(sol) & (np.abs(fa) <= SOLVE_TOL_MM)
    sol[m], res[m] = a[m], fa[m]
    n_iter = 0
    for it in range(MAX_ITER):
        live = np.flatnonzero(np.isnan(sol))
        if live.size == 0:
            break
        n_iter = it + 1
        al, bl, fal, fbl = a[live], b[live], fa[live], fb[live]
        mid = 0.5 * (al + bl)
        if it % 2 == 0:                     # はさみうち (割線) と二分法を交互に: 収束を保証しつつ速い
            den = fbl - fal
            c = np.where(den > 0, al - fal * (bl - al) / np.where(den > 0, den, 1.0), mid)
            c = np.where((c > al) & (c < bl), c, mid)
        else:
            c = mid
        if np.any(~((c > al) & (c < bl))):
            raise RebuildError(f"x(u) = x_i が {SOLVE_TOL_MM:g} mm まで解けない (区間が浮動小数の幅まで縮んだ)")
        fc = sgn * _xyz(gmsh, tag, c)[:, 0] - tgt[act][live]
        done = np.abs(fc) <= SOLVE_TOL_MM
        sol[live[done]], res[live[done]] = c[done], fc[done]
        left = fc < 0
        a[live] = np.where(left, c, al)
        fa[live] = np.where(left, fc, fal)
        b[live] = np.where(left, bl, c)
        fb[live] = np.where(left, fbl, fc)
    if np.isnan(sol).any():
        raise RebuildError(f"x(u) = x_i が {MAX_ITER} 回で {SOLVE_TOL_MM:g} mm まで解けない ({int(np.isnan(sol).sum())} station)")
    u_sol[act] = sol
    Pw = _xyz(gmsh, tag, u_sol)
    resid = np.abs(Pw[:, 0] - x_st)
    resid_solved = float(resid[act].max()) if act.size else 0.0
    if resid_solved > SOLVE_TOL_MM:
        raise RebuildError(f"x(u_i) − x_i が {resid_solved:.3e} mm (許容 {SOLVE_TOL_MM:g})")
    if not np.all(sgn * np.diff(u_sol) > 0):
        raise RebuildError("解いた助変数が station の順に単調でない")
    if not np.all(Pw[:, 1] > 0):
        raise RebuildError("station の壁半径が 0 以下")
    return {"u": u_sol, "r_wall_mm": Pw[:, 1].copy(), "x_wall_mm": Pw[:, 0].copy(), "trim": trim,
            "ends": end_info, "max_solve_residual_mm": resid_solved,
            "max_end_residual_mm": float(max(resid[0], resid[-1])), "n_iter": n_iter}


def trim_curve(gmsh, tag: int, cuts: list[float], keep_range: tuple[float, float]) -> int:
    """曲線を助変数 cuts で切り、keep_range (station の助変数の範囲) の外側を捨てた曲線のタグを返す
    (OCC の分割は元の助変数を保つ)。"""
    if not cuts:
        return tag
    P = _xyz(gmsh, tag, cuts)
    pts = [gmsh.model.occ.addPoint(*p) for p in P]
    gmsh.model.occ.fragment([(1, tag)], [(0, p) for p in pts])
    gmsh.model.occ.synchronize()
    keep, drop = [], []
    for _, t in gmsh.model.getEntities(1):
        a, b = (float(v[0]) for v in gmsh.model.getParametrizationBounds(1, t))
        (keep if keep_range[0] < 0.5 * (a + b) < keep_range[1] else drop).append(t)
    if len(keep) != 1:
        raise RebuildError(f"曲線を station の端で切れない (残る区間 {len(keep)} 本)")
    bnd = {abs(p[1]) for p in gmsh.model.getBoundary([(1, keep[0])], oriented=False)}
    gmsh.model.occ.remove([(1, t) for t in drop])
    gmsh.model.occ.synchronize()
    orphans = [(0, t) for _, t in gmsh.model.getEntities(0) if t not in bnd]
    if orphans:
        gmsh.model.occ.remove(orphans)
        gmsh.model.occ.synchronize()
    return keep[0]


# --------------------------------------------------------------------------------------------
# 形状と節点・要素
# --------------------------------------------------------------------------------------------

def build(gmsh, stc: dict, step: str, n_sample: int) -> dict:
    ni, nj, x_st, S = stc["ni"], stc["nj"], stc["x"], stc["S"]
    occ = gmsh.model.occ
    wall = import_curve(gmsh, step)
    cc = curve_checks(gmsh, wall, x_st, n_sample)
    sv = solve_stations(gmsh, wall, x_st, cc)
    u = sv["u"]
    cuts = []
    if sv["trim"]["inlet"]:
        cuts.append(float(u[0]))
    if sv["trim"]["outlet"]:
        cuts.append(float(u[-1]))
    wall = trim_curve(gmsh, wall, cuts, (float(u.min()), float(u.max())))
    # 曲線の端点 (入口側 = x が小さい方)
    bpts = [abs(p[1]) for p in gmsh.model.getBoundary([(1, wall)], oriented=False)]
    if len(bpts) != 2:
        raise RebuildError(f"壁曲線の端点が 2 つでない ({bpts})")
    bxyz = [np.asarray(gmsh.model.getValue(0, p, []), dtype=np.float64) for p in bpts]
    w0, w1 = (bpts[0], bpts[1]) if bxyz[0][0] < bxyz[1][0] else (bpts[1], bpts[0])
    xw0 = float(min(b[0] for b in bxyz))
    xw1 = float(max(b[0] for b in bxyz))
    a0 = occ.addPoint(xw0, 0.0, 0.0)
    a1 = occ.addPoint(xw1, 0.0, 0.0)
    l_axis = occ.addLine(a0, a1)
    l_out = occ.addLine(a1, w1)
    l_in = occ.addLine(a0, w0)
    loop = occ.addCurveLoop([l_axis, l_out, wall, l_in])
    surf = occ.addPlaneSurface([loop])
    occ.synchronize()
    sb = {abs(c[1]) for c in gmsh.model.getBoundary([(2, surf)], oriented=False)}
    if sb != {l_axis, l_out, wall, l_in}:
        raise RebuildError(f"面の境界が軸・出口・壁・入口の 4 本でない ({sorted(sb)})")
    occ.dilate(gmsh.model.getEntities(), 0.0, 0.0, 0.0, MM_TO_M, MM_TO_M, MM_TO_M)
    occ.synchronize()
    # 節点 (m)。番号 i·nj + j + 1
    R = S * sv["r_wall_mm"][:, None]
    X = np.repeat(x_st[:, None], nj, axis=1)
    XYZ = np.zeros((ni, nj, 3))
    XYZ[..., 0] = X * MM_TO_M
    XYZ[..., 1] = R * MM_TO_M
    tag = np.arange(1, ni * nj + 1, dtype=np.int64).reshape(ni, nj)
    msh = gmsh.model.mesh

    def put(dim, ent, ii, jj, param=None):
        t = tag[ii, jj].ravel()
        c = XYZ[ii, jj].reshape(-1, 3).ravel()
        if t.size:
            msh.addNodes(dim, ent, t.tolist(), c.tolist(), [] if param is None else list(map(float, param)))

    put(0, a0, [0], [0])
    put(0, w0, [0], [nj - 1])
    put(0, a1, [ni - 1], [0])
    put(0, w1, [ni - 1], [nj - 1])
    inner_i = np.arange(1, ni - 1)
    inner_j = np.arange(1, nj - 1)
    put(1, l_axis, inner_i, np.zeros_like(inner_i))
    put(1, wall, inner_i, np.full_like(inner_i, nj - 1), param=u[1:-1])
    put(1, l_in, np.zeros_like(inner_j), inner_j)
    put(1, l_out, np.full_like(inner_j, ni - 1), inner_j)
    II, JJ = np.meshgrid(inner_i, inner_j, indexing="ij")
    put(2, surf, II.ravel(), JJ.ravel())
    # 要素 (元と同じ番号の並び)
    quads, bedges = connectivity(ni, nj)
    ent = {"inlet": l_in, "outlet": l_out, "wall": wall, "axis": l_axis}
    et = 1
    for nm in stc["order"]:
        if nm == stc["surf"]:
            conn, dim, e, typ = quads, 2, surf, 3
        else:
            conn, dim, e, typ = bedges[nm], 1, ent[nm], 1
        tags = np.arange(et, et + len(conn), dtype=np.int64)
        msh.addElements(dim, e, [typ], [tags.tolist()], [(conn + 1).ravel().tolist()])
        et += len(conn)
    for nm, g in stc["phys"].items():
        e = surf if nm == stc["surf"] else ent[nm]
        gmsh.model.addPhysicalGroup(int(g["dim"]), [e], int(g["tag"]))
        gmsh.model.setPhysicalName(int(g["dim"]), int(g["tag"]), nm)
    info = {k: cc[k] for k in ("param_range", "x_range_mm", "max_abs_z_mm", "min_r_mm", "x_minus_u_max_abs_mm",
                               "n_sample", "direction")}
    info["x_u_identity_note"] = "x_minus_u_max_abs_mm は情報 (wall_step.py の曲線なら ~0、CAD で編集した曲線では大きくてよい)"
    stations = {k: sv[k] for k in ("ends", "max_solve_residual_mm", "max_end_residual_mm", "n_iter")}
    stations["trimmed"] = [s for s, v in sv["trim"].items() if v]
    stations["wall_u_range"] = [float(u.min()), float(u.max())]
    return {"curve": info, "stations": stations, "r_wall_mm": sv["r_wall_mm"], "u_wall": u,
            "coords": XYZ.reshape(-1, 3)}


def write_msh41_forge(path, coords, quads, bedges, phys: dict, curve_order, surf_name: str, digits: int = 17) -> None:
    """forge の変換器に入れる msh4.1 を書く。forge の `meshing/mesh2d.py` の `write_msh41_2d` と同じ並び・書式
    (物理名・実体 (曲線 4 + 面 1、bbox は 9 桁)・節点ブロック 1 つ (番号 1..N をファイル順)・要素は四角形が先で
    続いて曲線ごとの線分、番号はファイル順に 1..M)。違いは座標の有効桁 `digits` だけ (write_msh41_2d は 10)。
    phys: 物理名 → タグ、curve_order: 曲線 (= 境界辺) の名前の並び、surf_name: 面の物理名。"""
    coords = np.asarray(coords, dtype=np.float64)
    n_nodes = coords.shape[0]
    n_elems = len(quads) + sum(len(bedges[g]) for g in curve_order)
    mn = coords.min(0)
    mx = coords.max(0)
    L = []
    ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap(f"$PhysicalNames\n{len(curve_order) + 1}")
    for nm in curve_order:
        ap(f'1 {phys[nm]} "{nm}"')
    ap(f'2 {phys[surf_name]} "{surf_name}"')
    ap("$EndPhysicalNames")
    ap("$Entities")
    ap(f"0 {len(curve_order)} 1 0")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for ci, nm in enumerate(curve_order, start=1):
        ap(f"{ci} {bb} 1 {phys[nm]} 0")
    ap(f"1 {bb} 1 {phys[surf_name]} 0")
    ap("$EndEntities")
    ap("$Nodes")
    ap(f"1 {n_nodes} 1 {n_nodes}")
    ap(f"2 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    fmt = f"{{:.{int(digits)}g}} {{:.{int(digits)}g}} {{:.{int(digits)}g}}"
    ap("\n".join(fmt.format(*c) for c in coords.tolist()))
    ap("$EndNodes")
    ap("$Elements")
    ap(f"{1 + len(curve_order)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"2 1 3 {len(quads)}")
    buf = []
    for q in np.asarray(quads).tolist():
        buf.append(f"{et} {q[0] + 1} {q[1] + 1} {q[2] + 1} {q[3] + 1}")
        et += 1
    ap("\n".join(buf))
    for ci, nm in enumerate(curve_order, start=1):
        edges = np.asarray(bedges[nm]).tolist()
        ap(f"1 {ci} 1 {len(edges)}")
        buf = []
        for a, b in edges:
            buf.append(f"{et} {a + 1} {b + 1}")
            et += 1
        ap("\n".join(buf))
    ap("$EndElements")
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(L) + "\n")


def forge_format_check(path) -> dict:
    """msh4.1 が forge の変換器の読める形か: 節点ブロックが 1 つ、節点番号がファイル順に 1..N、要素番号がファイル順に 1..M
    (変換器は節点のタグ行を読まずファイル内の順番を番号とみなす)。"""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    out = {"path": str(path)}
    try:
        k = lines.index("$Nodes")
        nblk, nnode, tmin, tmax = (int(v) for v in lines[k + 1].split())
        out["n_node_blocks"] = nblk
        tags = []
        p = k + 2
        for _ in range(nblk):
            nn = int(lines[p].split()[3])
            tags.extend(int(v) for v in lines[p + 1:p + 1 + nn])
            p += 1 + 2 * nn
        out["node_tags_in_file_order"] = bool(tags == list(range(1, nnode + 1)) and (tmin, tmax) == (1, nnode))
        k = lines.index("$Elements")
        nblk_e, nel, _, _ = (int(v) for v in lines[k + 1].split())
        etags = []
        p = k + 2
        for _ in range(nblk_e):
            ne = int(lines[p].split()[3])
            etags.extend(int(r.split(maxsplit=1)[0]) for r in lines[p + 1:p + 1 + ne])
            p += 1 + ne
        out["element_tags_in_file_order"] = bool(etags == list(range(1, nel + 1)))
    except (ValueError, IndexError) as e:
        out["error"] = f"{type(e).__name__}: {e}"
    out["ok"] = bool(out.get("n_node_blocks") == 1 and out.get("node_tags_in_file_order")
                     and out.get("element_tags_in_file_order"))
    return out


def write_outputs(gmsh, prefix: str, stc: dict, coords) -> dict:
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    out = {"msh41": f"{prefix}.msh", "msh41_cad": f"{prefix}_cad.msh", "msh22": f"{prefix}_v22.msh",
           "cgns": f"{prefix}.cgns"}
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(out["msh41_cad"])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out["msh22"])
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(out["cgns"])
    quads, bedges = connectivity(stc["ni"], stc["nj"])
    write_msh41_forge(out["msh41"], coords, quads, bedges, {k: int(v["tag"]) for k, v in stc["phys"].items()},
                      stc["order"][1:], stc["surf"])
    return out


# --------------------------------------------------------------------------------------------
# 比較
# --------------------------------------------------------------------------------------------

def read_mesh(gmsh, path) -> dict:
    """gmsh で msh を読み、節点 (番号順)・四角形・線分・物理グループを返す。"""
    gmsh.clear()
    gmsh.open(str(path))
    t, c, _ = gmsh.model.mesh.getNodes()
    t = np.asarray(t, dtype=np.int64)
    c = np.asarray(c, dtype=np.float64).reshape(-1, 3)
    o = np.argsort(t)
    qt, qn = gmsh.model.mesh.getElementsByType(3)
    qt = np.asarray(qt, dtype=np.int64)
    qn = np.asarray(qn, dtype=np.int64).reshape(-1, 4)
    oq = np.argsort(qt)
    groups = {}
    for dim, ptag in gmsh.model.getPhysicalGroups():
        name = gmsh.model.getPhysicalName(dim, ptag)
        n_el = 0
        conn = []
        for e in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
            types, etags, enodes = gmsh.model.mesh.getElements(dim, int(e))
            for ty, tg, nd in zip(types, etags, enodes):
                n_el += len(tg)
                npe = len(nd) // max(len(tg), 1)
                tg = np.asarray(tg, dtype=np.int64)
                nd = np.asarray(nd, dtype=np.int64).reshape(-1, npe)
                conn.append((tg, nd))
        groups[name] = {"dim": int(dim), "tag": int(ptag), "n_elements": int(n_el), "conn": conn}
    return {"tags": t[o], "xyz": c[o], "quad_tags": qt[oq], "quads": qn[oq], "groups": groups}


def _edges_sorted(conn) -> np.ndarray:
    if not conn:
        return np.zeros((0, 3), dtype=np.int64)
    tg = np.concatenate([a for a, _ in conn])
    nd = np.concatenate([b for _, b in conn])
    o = np.argsort(tg)
    return np.c_[tg[o], nd[o]]


def _compare_mesh(A: dict, B: dict, ni: int, nj: int) -> tuple[dict, dict]:
    """元 A と作り直し B (read_mesh の戻り値) の (a) 四角形・(b) 節点の距離・(c) 物理グループ。戻り値 (結果, 合否)。"""
    out: dict = {}
    nodes_same = A["tags"].shape == B["tags"].shape and np.array_equal(A["tags"], B["tags"])
    out["n_nodes"] = {"orig": int(A["tags"].size), "rebuilt": int(B["tags"].size), "tags_match": bool(nodes_same)}
    quads_ok = (A["quads"].shape == B["quads"].shape and np.array_equal(A["quad_tags"], B["quad_tags"])
                and np.array_equal(A["quads"], B["quads"]))
    out["a_quads"] = {"n_orig": int(len(A["quads"])), "n_rebuilt": int(len(B["quads"])), "match": bool(quads_ok)}
    if nodes_same and A["tags"].size == ni * nj and np.array_equal(A["tags"], np.arange(1, ni * nj + 1)):
        d = np.linalg.norm(A["xyz"] - B["xyz"], axis=1) / MM_TO_M          # mm
        Ra = A["xyz"][:, 1].reshape(ni, nj) / MM_TO_M
        first = Ra[:, -1] - Ra[:, -2]
        ratio = d.reshape(ni, nj) / first[:, None]
        k = int(np.argmax(d))
        kr = int(np.argmax(ratio))
        dw = d.reshape(ni, nj)[:, -1]
        out["b_distance"] = {
            "max_dist_mm": float(d.max()), "at_node": {"i": k // nj, "j": k % nj, "x_mm": float(A["xyz"][k, 0] / MM_TO_M)},
            "max_dist_over_first_layer": float(ratio.max()),
            "ratio_at_node": {"i": kr // nj, "j": kr % nj, "first_layer_mm": float(first[kr // nj])},
            "max_wall_node_dist_mm": float(dw.max()), "rms_dist_mm": float(np.sqrt(np.mean(d ** 2))),
            "first_layer_min_mm": float(first.min()),
            "tol": {"max_dist_mm": DIST_TOL_MM, "max_dist_over_first_layer": RATIO_TOL},
        }
        dist_ok = bool(d.max() <= DIST_TOL_MM and ratio.max() <= RATIO_TOL)
    else:
        out["b_distance"] = {"error": "節点の番号の集合が元と違うか 1..ni·nj でない"}
        dist_ok = False
    grp = {}
    groups_ok = set(A["groups"]) == set(B["groups"])
    edges_ok = True
    for nm in sorted(set(A["groups"]) | set(B["groups"])):
        ga, gb = A["groups"].get(nm), B["groups"].get(nm)
        row = {"orig": None if ga is None else {k: ga[k] for k in ("dim", "tag", "n_elements")},
               "rebuilt": None if gb is None else {k: gb[k] for k in ("dim", "tag", "n_elements")}}
        row["match"] = bool(ga is not None and gb is not None and row["orig"] == row["rebuilt"])
        if ga is not None and gb is not None:
            ea, eb = _edges_sorted(ga["conn"]), _edges_sorted(gb["conn"])
            row["elements_identical"] = bool(ea.shape == eb.shape and np.array_equal(ea, eb))
            edges_ok &= row["elements_identical"]
        groups_ok &= row["match"]
        grp[nm] = row
    out["c_physical_groups"] = {"groups": grp, "match": bool(groups_ok)}
    out["info_group_elements_identical"] = bool(edges_ok)
    return out, {"a": bool(quads_ok), "b": dist_ok, "c": bool(groups_ok)}


def compare(gmsh, stc: dict, orig_path: str, files: dict) -> dict:
    """元の msh と作り直した `.msh` (forge の形)・`_cad.msh` を比べる (plan §6 V3)。"""
    ni, nj = stc["ni"], stc["nj"]
    A = read_mesh(gmsh, orig_path)
    out: dict = {"orig": str(orig_path), "rebuilt": str(files["msh41"]), "rebuilt_cad": str(files["msh41_cad"])}
    sha = hashlib.sha256(Path(orig_path).read_bytes()).hexdigest()
    src_sha = stc["raw"].get("source", {}).get("sha256")
    out["orig_sha256"] = sha
    out["orig_sha256_matches_structure"] = (src_sha == sha) if src_sha else None
    # forge の形の .msh (主): (a)(b)(c) は最上位に置く
    B = read_mesh(gmsh, files["msh41"])
    res_f, ok_f = _compare_mesh(A, B, ni, nj)
    out.update(res_f)
    out["d_forge_format"] = forge_format_check(files["msh41"])
    out["d_forge_format"]["note"] = ("forge の変換器は節点のタグ行を読まずファイル内の順番を番号とみなす。"
                                     "節点ブロック 1 つ・番号がファイル順に 1..N でないと接続が壊れる")
    # CAD 要素つきの _cad.msh (gmsh が書く)
    C = read_mesh(gmsh, files["msh41_cad"])
    res_c, ok_c = _compare_mesh(A, C, ni, nj)
    out["cad_msh"] = {k: res_c[k] for k in ("n_nodes", "a_quads", "b_distance", "c_physical_groups",
                                            "info_group_elements_identical")}
    out["cad_msh"]["info_forge_format"] = forge_format_check(files["msh41_cad"])
    # msh2.2 の読み直し (情報)。gmsh の msh2 の書き出しは節点・要素を実体の順に番号を振り直すので、番号ではなく
    # 四角形の角の座標 (要素の並びは同じ) と物理グループの要素数で同じメッシュであることを見る
    D = read_mesh(gmsh, files["msh22"])
    same_corners = bool(D["quads"].shape == C["quads"].shape and np.array_equal(
        D["xyz"][np.searchsorted(D["tags"], D["quads"])], C["xyz"][np.searchsorted(C["tags"], C["quads"])]))
    out["info_v22_readback"] = {
        "numbering_kept": bool(np.array_equal(D["quads"], C["quads"])),
        "quad_corner_coords_identical": same_corners,
        "group_counts_match": {k: v["n_elements"] for k, v in D["groups"].items()} == {k: v["n_elements"] for k, v in C["groups"].items()},
        "note": "gmsh の msh2.2 は番号を振り直す (元と同じ番号は .msh と _cad.msh だけ)",
    }
    reasons = []
    for tag, ok in (("", ok_f), ("_cad.msh の ", ok_c)):
        if not ok["a"]:
            reasons.append(f"(a) {tag}四角形の接続が元と違う")
        if not ok["b"]:
            reasons.append(f"(b) {tag}節点の距離が許容を超える")
        if not ok["c"]:
            reasons.append(f"(c) {tag}物理グループの名前・タグ・要素数が違う")
    if not out["d_forge_format"]["ok"]:
        reasons.append("(d) .msh が forge の変換器の読める形でない (節点ブロック 1 つ・番号がファイル順に 1..N)")
    out["criteria"] = ("plan §6 V3: .msh と _cad.msh の両方で (a) 四角形の接続が全要素で一致、(b) 最大距離 ≤ 1e-5 mm かつ"
                       "第 1 層厚比 ≤ 0.01、(c) 物理グループの名前・タグ・要素数が一致。加えて (d) .msh が forge の形")
    out["verdict"] = "PASS" if not reasons else "FAIL"
    out["fail_reasons"] = reasons
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="STEP の壁曲線 + mesh_structure.json → 同じ構造の gmsh メッシュ")
    ap.add_argument("step")
    ap.add_argument("structure")
    ap.add_argument("out_prefix")
    ap.add_argument("--compare", default=None, help="元の msh (節点の距離・接続・物理グループを比べる)")
    ap.add_argument("--n-sample", type=int, default=20001, help="曲線の検査の標本数")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):        # Windows のリダイレクト (cp932) でも記号で落ちないように
        try:
            stream.reconfigure(errors="backslashreplace")
        except (AttributeError, ValueError):
            pass
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("rebuild")
        try:
            stc = load_structure(a.structure)
            rep = build(gmsh, stc, a.step, a.n_sample)
        except RebuildError as e:
            print(f"REJECT: {e}", file=sys.stderr)
            return 2
        except Exception as e:                     # gmsh・OCC の失敗など
            print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
            return 3
        Path(a.out_prefix).parent.mkdir(parents=True, exist_ok=True)
        files = write_outputs(gmsh, a.out_prefix, stc, rep["coords"])
        print(f"rebuild: ni {stc['ni']} nj {stc['nj']} 曲線 {rep['curve']['direction']} 助変数 {rep['curve']['param_range']} "
              f"x(u)−u 最大 {rep['curve']['x_minus_u_max_abs_mm']:.3e} mm、解の残差 {rep['stations']['max_solve_residual_mm']:.2e} mm "
              f"(反復 {rep['stations']['n_iter']})、切った端 {rep['stations']['trimmed'] or 'なし'}")
        label = {"msh41": "forge の変換器に入れる形 (元と同じ番号)", "msh41_cad": "CAD 要素つき (gmsh・Fluent 用)",
                 "msh22": "msh2.2 (番号は gmsh が振り直す)", "cgns": "CGNS (番号・物理タグは gmsh が振り直す)"}
        for k, v in files.items():
            print(f"  {v}  … {label[k]}")
        if a.compare:
            cmp_ = compare(gmsh, stc, a.compare, files)
            doc = {"format": "forge.mesh_rebuild_compare/v1", "step": str(a.step), "structure": str(a.structure),
                   "outputs": files, "curve": rep["curve"], "stations": rep["stations"], **cmp_}
            path = f"{a.out_prefix}_compare.json"
            Path(path).write_text(json.dumps(_jsonable(doc), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
            bd = cmp_["b_distance"]
            print(f"compare: {cmp_['verdict']} 最大距離 {bd.get('max_dist_mm', float('nan')):.3e} mm・第 1 層厚比 "
                  f"{bd.get('max_dist_over_first_layer', float('nan')):.3e}・四角形 {cmp_['a_quads']['match']}・"
                  f"物理グループ {cmp_['c_physical_groups']['match']}・forge の形 {cmp_['d_forge_format']['ok']}・"
                  f"_cad.msh の最大距離 {cmp_['cad_msh']['b_distance'].get('max_dist_mm', float('nan')):.3e} mm → {path}")
            return 0 if cmp_["verdict"] == "PASS" else 1
        return 0
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
