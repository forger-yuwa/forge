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
    # --- 区間の局所の並べ直し (2026-10-13, plan tooling-nozzle-isothermal-wall-chain §4.8-1・§5.1 #29) ---
    # None なら従来どおり (既存の格子はビット同一)。[{x0, x1, blend, factor, ni_ref}] (区間は 1 つだけ) を指定すると、
    # x の station を ni_ref 列で従来どおり (_x_stations) 置いてから、区間 [x0, x1] の中だけ流れ方向の密度を max(1, d/factor) 側へ
    # 寄せて置き直す。区間の外の station は元の配列そのもの (半径方向・壁法線の層は列ごとに x だけで決まるので、区間の外の列はビット一致)。
    # 列の数は計算で決まるので ni はその値を書く (違えば generate_axisym_mesh が止める)。式は _x_stations_respaced を参照。
    x_respace: tuple | list | None = None


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
    t = None if density_table is None else _check_table("x_density_table", density_table, x0, x1)
    dens = _x_density(xs, refine, width, local_center, local_refine, local_width, t)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xs))])
    cum /= cum[-1]
    return np.interp(np.linspace(0.0, 1.0, ni), cum, xs)


def _x_density(xs: np.ndarray, refine: float, width: float, local_center: float = 0.0, local_refine: float = 1.0,
               local_width: float = 0.75, table: np.ndarray | None = None) -> np.ndarray:
    """station の相対密度 = スロート細分 × 局所細分 × 表 (_x_stations と _x_stations_respaced が共有する式)。
    table は _check_table 済みの (n, 2) 配列か None。"""
    dens = 1.0 + (refine - 1.0) * np.exp(-((xs / width) ** 2))
    if local_refine > 1.0:
        dens = dens * (1.0 + (local_refine - 1.0)
                       * np.exp(-(((xs - local_center) / local_width) ** 2)))
    if table is not None:
        dens = dens * np.interp(xs, table[:, 0], table[:, 1])
    return dens


X_RESPACE_KEYS = ("x0", "x1", "blend", "factor", "ni_ref")


def _x_stations_respaced(x_in: float, x_e: float, prm: "Mesh2DParams") -> tuple[np.ndarray, dict]:
    """区間の局所の並べ直し (`Mesh2DParams.x_respace`、plan tooling-nozzle-isothermal-wall-chain §4.8-1・§5.1 #29)。

    式・区間の端・積分の点数・丸め・端の上書きは forge-25 の G0c の作り方 (commit 96197821 の
    case/45.isobutane_m6_d155/core_grid_mesh.py:451-473 `grefc_stations`) と同じ:
    1. 元の station xs を ni_ref 列で _x_stations のとおり作る。ni_ref は「並べ直す前の列の数」(全域の配列の列の数)。
    2. 区間の端 xa・xb = xs のうち x0・x1 に最も近い列 (固定)。ブレンドの内側の端は**指定の** x0 + blend・x1 − blend のまま
       (端の station から blend を取らない。取ると case/45 の G0c で列が 1 つ増える: codex diagnose 2026-10-13)。
    3. [xa, xb] を 20,001 点に刻み、密度 d (_x_density = _x_stations と同じ式) から
       d_n = d·(1 − w·(1 − max(1, d/factor)/d))、w は両端の幅 blend の smoothstep (0 → 1 → 0)。
    4. 区間の列の数 n_new = round(n_old·∫d_n/∫d) (n_old = ib − ia、台形の累積)、d_n の累積の逆写像で置き直し、両端を xa・xb で上書き。
    5. 区間の外の station は xs そのもの。
    forge-25 版は d を x_density_table の線形補間だけで取る。こちらは _x_stations と同じ式なので、throat_refine 1・local_refine 1
    (case/45 の冷却壁の格子) では同じ値になる。ガウスの細分 (throat_refine ≠ 1・local_refine > 1) との組み合わせは旧 G0c の検証の範囲外。
    区間は 1 つだけ (複数の区間の重なり・順序・ni_ref の意味は決めていないので受けない)。戻り値 (station, 情報)。"""
    segs = prm.x_respace
    if not isinstance(segs, (list, tuple)) or len(segs) != 1 or not isinstance(segs[0], dict):
        raise ValueError(f"x_respace: [{{x0, x1, blend, factor, ni_ref}}] の 1 区間だけを受ける (受け取った値: {segs!r})")
    s = segs[0]
    bad = sorted(set(map(str, s)) ^ set(X_RESPACE_KEYS))
    if bad:
        raise ValueError(f"x_respace: キーは {list(X_RESPACE_KEYS)} のちょうど 5 つ (過不足 {bad})")
    x0, x1, blend, factor = (float(s[k]) for k in ("x0", "x1", "blend", "factor"))
    ni_ref = s["ni_ref"]
    if isinstance(ni_ref, bool) or int(ni_ref) != ni_ref or int(ni_ref) < 3:
        raise ValueError(f"x_respace.ni_ref: 3 以上の整数にする ({ni_ref!r})")
    ni_ref = int(ni_ref)
    if not all(np.isfinite(v) for v in (x0, x1, blend, factor)) or blend <= 0.0 or factor <= 0.0:
        raise ValueError(f"x_respace: x0・x1・blend・factor は有限、blend・factor は正にする ({x0}, {x1}, {blend}, {factor})")
    if not (x_in < x0 < x1 < x_e):
        raise ValueError(f"x_respace: x_in {x_in} < x0 {x0} < x1 {x1} < x_e {x_e} にする")
    a1, b1 = x0 + blend, x1 - blend                   # ブレンドの内側の端 (指定の x0・x1 から取る)
    if a1 > b1:
        raise ValueError(f"x_respace: 2·blend ({2 * blend}) が区間の幅 ({x1 - x0}) を超える")
    xs = _x_stations(x_in, x_e, ni_ref, prm.throat_refine, prm.throat_width,
                     prm.local_center, prm.local_refine, prm.local_width, prm.x_density_table)
    ia = int(np.argmin(np.abs(xs - x0)))
    ib = int(np.argmin(np.abs(xs - x1)))
    xa, xb = float(xs[ia]), float(xs[ib])
    if not (0 < ia < ib < ni_ref - 1 and xa < a1 and b1 < xb):
        raise ValueError(f"x_respace: 区間の端の列 (ia {ia} x {xa}, ib {ib} x {xb}) が両端の列の内側で、ブレンドの内側の端 "
                         f"({a1}, {b1}) より外にない (ni_ref {ni_ref} の格子が粗すぎるか blend が狭すぎる)")
    xf = np.linspace(xa, xb, 20001)
    t = None if prm.x_density_table is None else _check_table("x_density_table", prm.x_density_table, x_in, x_e)
    d = _x_density(xf, prm.throat_refine, prm.throat_width, prm.local_center, prm.local_refine, prm.local_width, t)
    w = np.clip((xf - xa) / (a1 - xa), 0, 1) * np.clip((xb - xf) / (xb - b1), 0, 1)       # 0 → 1 → 0 (端でなめらかに)
    w = w * w * (3 - 2 * w)                                                                 # smoothstep
    dn = d * (1 - w * (1 - np.maximum(1.0, d / factor) / d))
    cr = np.concatenate([[0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(xf))])
    cn = np.concatenate([[0], np.cumsum(0.5 * (dn[1:] + dn[:-1]) * np.diff(xf))])
    nseg = int(round((ib - ia) * cn[-1] / cr[-1]))
    if nseg < 1:
        raise ValueError(f"x_respace: 区間の列の数が 1 未満になる (n_old {ib - ia}・∫d_n/∫d {cn[-1] / cr[-1]})")
    seg = np.interp(np.linspace(0, cn[-1], nseg + 1), cn, xf)
    seg[0], seg[-1] = xa, xb
    out = np.concatenate([xs[:ia], seg, xs[ib + 1:]])
    if not np.all(np.diff(out) > 0):
        raise ValueError("x_respace: 並べ直した station が狭義単調増加でない")
    return out, {"ni_ref": ni_ref, "ia": ia, "ib": ib, "xa": xa, "xb": xb, "n_old": ib - ia, "n_new": nseg, "ni": len(out)}


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


def _capped_layout(nj_base: int, first_frac: float, cap: float, max_nj: int | None = None) -> tuple:
    """`_radial_fracs_capped` の配点の内訳 (壁側の等比の間隔の list (壁側から)、残り rest、軸側の一様の区間の数 m)。径方向の点数は
    len(等比) + m + 1。一様の部分はリストを作らない (m は式で決まる; 極小の cap でも巨大な確保をしない)。
    max_nj を渡すと、点数の下限 len(等比) + 2 が max_nj を超えた時点で等比の数え上げをやめ (等比, None, None) を返す
    (plan tooling-nozzle-mesh-size-guard §4.2)。max_nj なしの算術は従来の `_radial_fracs_capped` と同じ (生成座標は不変)。"""
    s0 = _radial_fracs(nj_base, first_frac)
    g0 = np.diff(s0)[::-1]                      # 壁側から
    q = g0[1] / g0[0]
    gaps = []
    g = g0[0]
    while g < cap and sum(gaps) + g < 1.0:
        gaps.append(g)
        g *= q
        if max_nj is not None and len(gaps) + 2 > max_nj:
            return gaps, None, None
    rest = 1.0 - sum(gaps)
    m = int(np.ceil(rest / cap - 1e-12))
    return gaps, rest, m


def _radial_fracs_capped(nj_base: int, first_frac: float, cap: float) -> np.ndarray:
    """壁側は `_radial_fracs(nj_base, first_frac)` と同じ第一セル・比 q の等比、間隔が `cap` に達したら軸まで一様。"""
    gaps, rest, m = _capped_layout(nj_base, first_frac, cap)
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


def radial_node_count(prm: Mesh2DParams, max_nj: int | None = None) -> tuple:
    """`generate_axisym_mesh` が作る径方向の点数 (実効の nj) を、格子を作らずに返す。戻り値 (nj, 正確か)。
    axis_gap_frac のときは導出値 (壁側の等比 + 軸側の一様; `_capped_layout`)、それ以外は prm.nj。
    max_nj を渡すと、nj がそれを超えると分かった時点で数えるのをやめ (下限 > max_nj, False) を返す。
    併用できない指定の組は ValueError (generate_axisym_mesh もこの関数で nj を決めるので、判定は 1 か所)。"""
    if prm.axis_cap_frac is not None:
        if prm.axis_gap_frac is not None:
            raise ValueError("axis_cap_frac と axis_gap_frac は同時に指定できない")
        return prm.nj, True
    if prm.axis_gap_frac is not None:
        if prm.wall_first_frac_throat is not None or prm.wall_first_frac_table is not None:
            raise ValueError("axis_gap_frac と wall_first_frac_throat・wall_first_frac_table の併用は未対応")
        if max_nj is not None and int(prm.nj) - 1 > max_nj:
            # 配列を作る前の打ち切り (2026-10-14 result 段レビュー M1: `_capped_layout` は指定 nj の `_radial_fracs` を作るので、nj が極端に大きいと
            # 数える前に巨大な配列を確保していた)。実効の nj ≥ 指定 nj − 1: 等比の部分は元の配点の先頭の区間そのままで、残りの区間はどれも cap 以上
            # (q ≥ 1 なので打ち切った後の区間は大きくなる一方) だから一様の区間の数 m = ceil(rest / cap) は残りの区間の数以上になる。
            # 前提は `_radial_fracs` の等比の配点が正常なこと (nj が大きく第 1 層が薄いと Newton が退化して第 1 層 0 になる既存の欠陥がある; その格子は品質検査で落ちる)。
            return int(prm.nj) - 1, False
        gaps, _, m = _capped_layout(prm.nj, prm.wall_first_frac, prm.axis_gap_frac, max_nj)
        if m is None:
            return len(gaps) + 2, False
        return len(gaps) + m + 1, True
    return prm.nj, True


def node_count(prm: Mesh2DParams, limit: int | None = None) -> int:
    """`generate_axisym_mesh` が作る節点数 (ni × 実効の nj) を、格子を作らずに返す (plan tooling-nozzle-mesh-size-guard §4.2)。
    limit を渡すと、超えると分かった時点で数えるのをやめて limit + 1 以上の値を返す (極小の axis_gap_frac でも巨大な確保をしない)。"""
    ni = int(prm.ni)
    max_nj = None if limit is None else int(limit) // max(ni, 1)
    nj, _ = radial_node_count(prm, max_nj)
    return ni * int(nj)


def generate_axisym_mesh(wall, prm: Mesh2DParams):
    """wall: NozzleWall。戻り値 (coords (N,3) [m], quads (M,4), 境界辺 dict)。"""
    if prm.x_respace is None:
        xs = _x_stations(wall.x_in, wall.x_e, prm.ni, prm.throat_refine, prm.throat_width,
                         prm.local_center, prm.local_refine, prm.local_width, prm.x_density_table)
    else:
        xs, info = _x_stations_respaced(float(wall.x_in), float(wall.x_e), prm)
        if len(xs) != prm.ni:          # 列の数は並べ直しで決まる。ni と違えば黙って作らない
            raise ValueError(f"x_respace: 並べ直した列の数 {len(xs)} が mesh.ni {prm.ni} と違う (ni に {len(xs)} を書く; {info})")
    rw = wall.r(xs)
    ni = prm.ni
    nj, _ = radial_node_count(prm)      # 実効の nj (axis_gap_frac では導出値)。事前の判定 (node_count) と同じ決まり方
    X = np.repeat(xs[:, None], nj, axis=1)
    if prm.axis_cap_frac is not None:
        fr = _first_frac_profile(xs, prm)
        R = np.empty((ni, nj))
        for i in range(ni):
            si, _ = _radial_fracs_capfixed(nj, float(fr[i]), float(prm.axis_cap_frac))
            R[i, :] = rw[i] * si
    elif prm.axis_gap_frac is not None:
        s = _radial_fracs_capped(prm.nj, prm.wall_first_frac, prm.axis_gap_frac)
        if len(s) != nj:
            raise RuntimeError(f"axis_gap_frac: 径方向の点数 {len(s)} が radial_node_count の {nj} と一致しない")
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
