"""軸対称ノズルの構造化 2D quad メッシュ (トポロジ固定・決定的)。

- x 方向: スロート近傍を滑らかに細分する間隔関数の逆積分で station を配置。
- r 方向: 全 station 共通の正規化分布 s_j (壁側幾何級数クラスタリング)。
  r_j(i) = r_w(x_i) * s_j — 同一トポロジで壁だけ動かせる (親計画 §4.3 R1)。
- 出力: gmsh msh4.1 テキスト (物理タグ inlet=1 outlet=2 wall=3 axis=4 fluid=5)
  → 既存 convertGmshToForge で forge h5 化 (wall_dist・幾何量は変換器が計算)。

軸対称の座標規約: x=軸, y=半径, z=0 平面, 上半分のみ (methods/axisymmetric)。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PHYS = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4, "fluid": 5}


@dataclass
class Mesh2DParams:
    ni: int = 241              # x 方向ノード数
    nj: int = 81               # r 方向ノード数
    wall_first_frac: float = 2.0e-3  # 壁第一セル厚 / 局所半径
    throat_refine: float = 3.0       # スロート近傍の間隔比 (far/throat)
    throat_width: float = 1.5        # 細分の e^-x^2 幅 (r* 単位)
    scale: float = 1.0               # r* [m] (無次元 → m)
    # --- 接合部近傍の局所細分 (B10: x_reach のアンカー安定性・こぶの格子収束用) ---
    local_center: float = 0.0        # 細分中心 [r*] (0 で無効)
    local_refine: float = 1.0        # 中心での間隔比 (1.0 で無効)
    local_width: float = 0.75        # 細分の e^-x^2 幅 [r*]
    # --- 冷却壁用: 第一セル厚の x 依存 (2026-09-12, plan tooling-nozzle-isothermal-wall-chain §4.2) ---
    # 冷却壁の y+ はスロート近傍 (T_e 高・ρ_w/ρ_e 大) で最も高いので、そこだけ第一セルを詰める。
    # None なら従来どおり一様 wall_first_frac。x <= blend_x0 で wall_first_frac_throat、blend_x1 以降で wall_first_frac、
    # 間は smoothstep で繋ぐ (r* 単位)。
    wall_first_frac_throat: float | None = None
    wall_first_blend_x0: float = 0.5
    wall_first_blend_x1: float = 6.0
    # 上流側 (収縮部) も同様に戻す: x <= up_x0 で wall_first_frac、up_x1 以上で wall_first_frac_throat (None = 上流は全域 throat 値)
    wall_first_up_x0: float | None = None
    wall_first_up_x1: float | None = None
    # --- 軸側の最大間隔の上限 (2026-10-04, plan verification-m6-axis-wave-mesh-su2 §4.1) ---
    # None なら従来どおり nj 点の等比。指定すると、壁側は nj から決まる同じ第一セル・同じ比 q で伸ばし、
    # 間隔 (/ r_w) がこの値に達したら軸まで一様にする。nj は導出値になる (壁・境界層の格子は不変)。
    axis_gap_frac: float | None = None
    # --- 軸側の間隔の上限・nj 固定版 (2026-10-06, plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11h) ---
    # None なら従来どおり。指定すると各断面で、壁の第一セル (wall_first_frac / 断面ごとの値) から比 q_i の等比で伸ばし、
    # 間隔 (/ r_w) がこの値 c に達したら軸まで一様にする。nj は全断面で同じ (Σ_j min(fr_i q^j, c) = 1 を q_i について解く)。
    # axis_gap_frac (nj が導出値) とは併用しない。
    axis_cap_frac: float | None = None
    # --- 表で与える第一セル厚と x 方向の密度 (2026-10-08, plan tooling-nozzle-isothermal-wall-chain §5.1 #13) ---
    # 冷却壁は y1+ が断熱の 8〜9 倍になり、必要な第一セル厚が x で 1 桁以上変わる (入口直管・スロートで最も薄く、試験部で厚い)。
    # throat/blend の smoothstep 1 本ではこの形に合わず、AR ≤ 5000 と y1+ ≤ 1 を同時に満たせない。
    # wall_first_frac_table: [[x, 第一セル厚/局所半径], ...] (x 昇順、r* 単位)。log 線形で補間。指定すると wall_first_frac・
    #   wall_first_frac_throat・前後ブレンドより優先する (併用は例外)。表が x の範囲を覆わなければ例外 (端値で黙って外挿しない)。
    # x_density_table: [[x, 相対密度 > 0], ...]。線形補間した密度をスロート細分・局所細分の密度に掛ける。範囲の扱いは同上。
    # どちらも None なら従来どおり (既存の格子はビット同一)。
    wall_first_frac_table: tuple | list | None = None
    x_density_table: tuple | list | None = None
    # --- 壁法線に沿った近壁層 (2026-10-08, plan tooling-nozzle-isothermal-wall-chain §5.1 #15) ---
    # None なら従来どおり (j 方向の線は半径方向)。[d_n, d_b] (局所半径比) を指定すると、壁からの距離 d (= r_w − 半径方向の配置の r) が
    # d ≤ d_n の節点を壁の法線上 (壁点 + d·内向き法線) に置き、d_n〜d_b で半径方向の配置へ smoothstep で戻す (d ≥ d_b と軸は従来の位置)。
    # 傾斜壁 (縮流部で最大 40°) では半径方向の線が壁と直交せず、AR の大きい近壁セルにスキューが入る (AGENTS.md の AR ≤ 5000 の例外の外)。
    # 壁の第一セル厚は法線距離 = d になる (半径方向の配置では法線距離 = d·cos θ_w)。
    wall_normal_layer: tuple | list | None = None


def _check_table(name: str, tbl, x0: float, x1: float, positive: bool = True) -> np.ndarray:
    """[[x, v], ...] を (n, 2) 配列にして検査する: x は狭義単調増加・有限、v は有限 (positive なら > 0)、[x0, x1] を覆う。"""
    t = np.asarray(tbl, dtype=float)
    if t.ndim != 2 or t.shape[1] != 2 or len(t) < 2:
        raise ValueError(f"{name}: [[x, 値], ...] の 2 点以上の表が要る (shape {t.shape})")
    if not np.all(np.isfinite(t)):
        raise ValueError(f"{name}: 非有限値がある")
    if not np.all(np.diff(t[:, 0]) > 0):
        raise ValueError(f"{name}: x は狭義単調増加にする")
    if positive and not np.all(t[:, 1] > 0):
        raise ValueError(f"{name}: 値は正にする")
    tol = 1e-9 * max(1.0, abs(x0), abs(x1))
    if t[0, 0] > x0 + tol or t[-1, 0] < x1 - tol:
        raise ValueError(f"{name}: 表の x 範囲 [{t[0, 0]}, {t[-1, 0]}] が格子の範囲 [{x0}, {x1}] を覆わない (端値で外挿しない)")
    return t


def _first_frac_profile(xs: np.ndarray, prm: "Mesh2DParams") -> np.ndarray:
    """各 station の第一セル厚 / 局所半径。表 → throat/blend → 一様、の順。"""
    if prm.wall_first_frac_table is not None:
        if prm.wall_first_frac_throat is not None or prm.wall_first_up_x0 is not None or prm.wall_first_up_x1 is not None:
            raise ValueError("wall_first_frac_table と wall_first_frac_throat・wall_first_up_* は同時に指定できない")
        t = _check_table("wall_first_frac_table", prm.wall_first_frac_table, float(xs[0]), float(xs[-1]))
        return np.exp(np.interp(xs, t[:, 0], np.log(t[:, 1])))
    if prm.wall_first_frac_throat is None:
        return np.full(len(xs), float(prm.wall_first_frac))
    t = np.clip((xs - prm.wall_first_blend_x0) / max(prm.wall_first_blend_x1 - prm.wall_first_blend_x0, 1e-9), 0.0, 1.0)
    t = t * t * (3.0 - 2.0 * t)
    if prm.wall_first_up_x0 is not None and prm.wall_first_up_x1 is not None:
        tu = np.clip((prm.wall_first_up_x1 - xs) / max(prm.wall_first_up_x1 - prm.wall_first_up_x0, 1e-9), 0.0, 1.0)
        tu = tu * tu * (3.0 - 2.0 * tu)
        t = np.maximum(t, tu)
    return prm.wall_first_frac_throat + (prm.wall_first_frac - prm.wall_first_frac_throat) * t


def _x_stations(x0: float, x1: float, ni: int, refine: float, width: float,
                local_center: float = 0.0, local_refine: float = 1.0,
                local_width: float = 0.75, density_table=None) -> np.ndarray:
    """間隔 h(x) ∝ 1/(密度) の逆積分で station を置く。

    密度 = スロート細分 (中心 x=0) × 局所細分 (中心 `local_center`)。後者は
    **接合部近傍だけ**を細かくして、$x_{reach}$ の局所微分アンカーと接合こぶの
    格子収束を、全域解像度を上げずに調べるためのもの (B10)。
    """
    xs = np.linspace(x0, x1, 8001)
    dens = 1.0 + (refine - 1.0) * np.exp(-((xs / width) ** 2))
    if local_refine > 1.0:
        dens = dens * (1.0 + (local_refine - 1.0)
                       * np.exp(-(((xs - local_center) / local_width) ** 2)))
    if density_table is not None:
        t = _check_table("x_density_table", density_table, x0, x1)
        dens = dens * np.interp(xs, t[:, 0], t[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xs))])
    cum /= cum[-1]
    return np.interp(np.linspace(0.0, 1.0, ni), cum, xs)


def _radial_fracs(nj: int, first_frac: float) -> np.ndarray:
    """s_0=0 (軸) → s_{nj-1}=1 (壁)。壁側の最終区間 = first_frac。幾何級数。"""
    n = nj - 1
    if first_frac * n >= 1.0:
        return np.linspace(0.0, 1.0, nj)
    # 比 q>1 を解く: first*(q^n - 1)/(q - 1) = 1
    q = 1.5
    for _ in range(100):
        f = first_frac * (q ** n - 1.0) / (q - 1.0) - 1.0
        df = first_frac * ((n * q ** (n - 1)) * (q - 1.0) - (q ** n - 1.0)) / (q - 1.0) ** 2
        qn = q - f / df
        if abs(qn - q) < 1e-14:
            q = qn
            break
        q = qn
    gaps = first_frac * q ** np.arange(n)          # 壁側から軸側へ拡大
    s = np.concatenate([[0.0], np.cumsum(gaps[::-1])])
    return s / s[-1]


def _radial_fracs_capped(nj_base: int, first_frac: float, cap: float) -> np.ndarray:
    """壁側は `_radial_fracs(nj_base, first_frac)` と同じ第一セル・比 q の等比、間隔が `cap` に達したら軸まで一様。"""
    s0 = _radial_fracs(nj_base, first_frac)
    g0 = np.diff(s0)[::-1]                      # 壁側から
    q = g0[1] / g0[0]
    gaps = []
    g = g0[0]
    while g < cap and sum(gaps) + g < 1.0:
        gaps.append(g)
        g *= q
    rest = 1.0 - sum(gaps)
    m = int(np.ceil(rest / cap - 1e-12))
    gaps = gaps + [rest / m] * m
    s = np.concatenate([[0.0], np.cumsum(gaps[::-1])])
    return s / s[-1]


def _radial_fracs_capfixed(nj: int, first_frac: float, cap: float) -> tuple[np.ndarray, float]:
    """s_0=0 (軸) → s_{nj-1}=1 (壁)。壁側から first_frac·q^j、上限 cap で頭打ち、総数 nj−1 で和 1 になる q を二分法で解く。
    戻り値 (s, 隣接比の最大)。実現不能 (first > cap、n·first > 1、first + (n−1)·cap < 1) は例外。"""
    n = nj - 1
    fr, c = float(first_frac), float(cap)
    if not (np.isfinite(fr) and np.isfinite(c) and fr > 0 and c > 0):
        raise ValueError(f"axis_cap_frac: 不正な値 first {fr} cap {c}")
    if fr > c:
        raise ValueError(f"axis_cap_frac: 第一セル {fr} が上限 {c} より大きい")
    if n * fr > 1.0 or fr + (n - 1) * c < 1.0:
        raise ValueError(f"axis_cap_frac: nj {nj}・第一セル {fr}・上限 {c} では和が 1 にならない (n·fr ≤ 1 ≤ fr + (n−1)c が必要)")
    j = np.arange(n)
    total = lambda q: float(np.minimum(fr * q ** j, c).sum())  # noqa: E731
    lo, hi = 1.0, 2.0
    while total(hi) < 1.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < 1.0:
            lo = mid
        else:
            hi = mid
    gaps = np.minimum(fr * hi ** j, c)          # 壁側から
    gaps = gaps / gaps.sum()
    if not (np.all(np.isfinite(gaps)) and np.all(gaps > 0)):
        raise ValueError("axis_cap_frac: 間隔が不正")
    ratio = float(np.max(np.maximum(gaps[1:] / gaps[:-1], gaps[:-1] / gaps[1:])))
    s = np.concatenate([[0.0], np.cumsum(gaps[::-1])])
    s = s / s[-1]
    if not np.all(np.diff(s) > 0):
        raise ValueError("axis_cap_frac: 節点が単調でない")
    return s, ratio


def generate_axisym_mesh(wall, prm: Mesh2DParams):
    """wall: NozzleWall。戻り値 (coords (N,3) [m], quads (M,4), 境界辺 dict)。"""
    xs = _x_stations(wall.x_in, wall.x_e, prm.ni, prm.throat_refine, prm.throat_width,
                     prm.local_center, prm.local_refine, prm.local_width, prm.x_density_table)
    rw = wall.r(xs)
    ni, nj = prm.ni, prm.nj
    X = np.repeat(xs[:, None], nj, axis=1)
    if prm.axis_cap_frac is not None:
        if prm.axis_gap_frac is not None:
            raise ValueError("axis_cap_frac と axis_gap_frac は同時に指定できない")
        fr = _first_frac_profile(xs, prm)
        R = np.empty((ni, nj))
        for i in range(ni):
            si, _ = _radial_fracs_capfixed(nj, float(fr[i]), float(prm.axis_cap_frac))
            R[i, :] = rw[i] * si
    elif prm.axis_gap_frac is not None:
        if prm.wall_first_frac_throat is not None or prm.wall_first_frac_table is not None:
            raise ValueError("axis_gap_frac と wall_first_frac_throat・wall_first_frac_table の併用は未対応")
        s = _radial_fracs_capped(nj, prm.wall_first_frac, prm.axis_gap_frac)
        nj = len(s)
        X = np.repeat(xs[:, None], nj, axis=1)
        R = rw[:, None] * s[None, :]
    elif prm.wall_first_frac_throat is None and prm.wall_first_frac_table is None:
        s = _radial_fracs(nj, prm.wall_first_frac)
        R = rw[:, None] * s[None, :]
    else:
        # station ごとに第一セル比を変える (構造は同じ nj、壁側クラスタリングだけ x で滑らかに変化)
        fr = _first_frac_profile(xs, prm)
        R = np.empty((ni, nj))
        for i in range(ni):
            R[i, :] = rw[i] * _radial_fracs(nj, float(fr[i]))
    if prm.wall_normal_layer is not None:
        dn, db = (float(v) for v in prm.wall_normal_layer)
        if not (0.0 < dn < db < 1.0):
            raise ValueError(f"wall_normal_layer: 0 < d_n < d_b < 1 にする ({dn}, {db})")
        try:
            drw = np.asarray(wall.r(xs, 1), dtype=float)       # 壁の解析的な傾き (物理壁は導関数を返す)
        except TypeError:
            drw = None
        if drw is None or drw.shape != xs.shape:
            drw = np.gradient(rw, xs)
        nrm = np.sqrt(1.0 + drw ** 2)
        d = rw[:, None] - R                                  # 壁からの (半径方向の配置での) 距離 [r*]
        t = np.clip((d / rw[:, None] - dn) / (db - dn), 0.0, 1.0)
        beta = 1.0 - t * t * (3.0 - 2.0 * t)                 # 1 (d ≤ d_n) → 0 (d ≥ d_b)
        Xn = xs[:, None] + d * (drw / nrm)[:, None]          # 壁点 + d·内向き法線 (r′, −1)/√(1 + r′²)
        Rn = rw[:, None] - d / nrm[:, None]
        X = X + beta * (Xn - X)
        R = R + beta * (Rn - R)
        # 格子の線が交差していないこと (各 j で x が単調、各 i で r が単調)
        if not (np.all(np.diff(X, axis=0) > 0) and np.all(np.diff(R, axis=1) > 0)):
            raise ValueError("wall_normal_layer: 格子の線が交差する (d_n・d_b を小さくする)")
    coords = np.zeros((ni * nj, 3))
    coords[:, 0] = X.ravel() * prm.scale
    coords[:, 1] = R.ravel() * prm.scale

    def nid(i, j):
        return i * nj + j

    quads = np.empty(((ni - 1) * (nj - 1), 4), dtype=np.int64)
    k = 0
    for i in range(ni - 1):
        for j in range(nj - 1):
            quads[k] = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            k += 1
    bedges = {
        "axis": [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)],
        "wall": [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)],
        "inlet": [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)],
        "outlet": [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)],
    }
    return coords, quads, bedges


def write_msh41_2d(path, coords, quads, bedges, digits: int = 10) -> None:
    """2D 平面メッシュを gmsh msh4.1 テキストで書く (決定的)。"""
    curve_order = ["inlet", "outlet", "wall", "axis"]
    n_nodes = coords.shape[0]
    n_elems = quads.shape[0] + sum(len(bedges[g]) for g in curve_order)
    mn = coords.min(0)
    mx = coords.max(0)
    L = []
    ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap("$PhysicalNames\n5")
    for nm in curve_order:
        ap(f'1 {PHYS[nm]} "{nm}"')
    ap(f'2 {PHYS["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities")
    ap(f"0 {len(curve_order)} 1 0")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for ci, nm in enumerate(curve_order, start=1):
        ap(f"{ci} {bb} 1 {PHYS[nm]} 0")
    ap(f"1 {bb} 1 {PHYS['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes")
    ap(f"1 {n_nodes} 1 {n_nodes}")
    ap(f"2 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    # digits: 座標の有効桁数 (既定 10 = 従来どおり)。冷却壁の格子 (第一層厚 / 半径 3e-7) は 10 桁だと第一層厚が 0.2 % 動くので 17 にする
    ap("\n".join(f"{c[0]:.{digits}g} {c[1]:.{digits}g} {c[2]:.{digits}g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements")
    ap(f"{1 + len(curve_order)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"2 1 3 {quads.shape[0]}")  # dim2, surface ent1, type3=quad
    buf = []
    for q in quads:
        buf.append(f"{et} {q[0]+1} {q[1]+1} {q[2]+1} {q[3]+1}")
        et += 1
    ap("\n".join(buf))
    for ci, nm in enumerate(curve_order, start=1):
        edges = bedges[nm]
        ap(f"1 {ci} 1 {len(edges)}")  # dim1, curve ent, type1=line
        buf = []
        for a, b in edges:
            buf.append(f"{et} {a+1} {b+1}")
            et += 1
        ap("\n".join(buf))
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
