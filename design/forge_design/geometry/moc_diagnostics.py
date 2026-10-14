r"""軸則 A/B/C 比較の診断一式 (計画: plans/active/tooling-nozzle-axislaw-smoothness.md)。

- `axis_smoothness_diagnostics`: 軸則の $M',M'',M'''$・$J_{\rm axis}$・knot/節点ジャンプ
  (`law.deriv(x, order)` を持つ任意のオブジェクトに対して動く — A/B/C 共通)。
- `characteristic_topology_diagnostics`: `fill_levels` が返す三角充填網の位相健全性
  (signed area・orientation flip・退化セル・隣接特性線間隔・点堆積・非隣接 C±
  polyline 交差)。単なる $\theta<\mu$ 判定の代用ではなく、網そのものを検査する。
- `wall_curvature_diagnostics`: 壁曲率 $\kappa_w(x)$・$d\kappa_w/ds$ を **`AxisMachCFDWall`
  の補間 spline の解析微分**で評価する (生の非一様点への数値差分ではない — 壁 spline
  自体が既に「平滑化しない補間 spline」)。
"""
from __future__ import annotations

import numpy as np


# --- 軸則の滑らかさ -----------------------------------------------------------
def axis_smoothness_diagnostics(law, x_A: float, x_E: float, n: int = 4001,
                                breakpoints: list | None = None) -> dict:
    r"""$M',M'',M'''$ の統計と $J_{\rm axis}=\int(M''')^2dx$ (Gauss 求積)、
    節点 (`breakpoints`, 例: knot 則の $x_K$) 前後の $M'''$ ジャンプ。

    `breakpoints=None` (B/C など明示的な低次節点を持たない、または報告しない場合) は
    ジャンプ診断を省く。"""
    x = np.linspace(x_A, x_E, n)
    Mp = np.asarray(law.deriv(x, 1), dtype=float)
    Mpp = np.asarray(law.deriv(x, 2), dtype=float)
    Mppp = np.asarray(law.deriv(x, 3), dtype=float)
    from numpy.polynomial.legendre import leggauss
    gx, gw = leggauss(8)
    # J_axis: breakpoints があれば区間ごと、なければ全域を 40 分割して Gauss 求積
    edges = sorted(set([x_A, x_E] + (list(breakpoints) if breakpoints else [])))
    segs = np.linspace(0.0, 1.0, 41)
    fine_edges = []
    for a, b in zip(edges[:-1], edges[1:]):
        fine_edges += list(a + segs * (b - a))
    fine_edges = sorted(set(fine_edges))
    J = 0.0
    for a, b in zip(fine_edges[:-1], fine_edges[1:]):
        if b - a < 1e-14:
            continue
        xm, xr = 0.5 * (a + b), 0.5 * (b - a)
        xs = xm + xr * gx
        v = np.asarray(law.deriv(xs, 3), dtype=float)
        J += xr * float(np.sum(gw * v * v))
    out = {
        "max_abs_Mp": float(np.max(np.abs(Mp))), "min_Mp": float(np.min(Mp)),
        "max_abs_Mpp": float(np.max(np.abs(Mpp))),
        "max_abs_Mppp": float(np.max(np.abs(Mppp))),
        "J_axis": float(J),
        "mpp_sign_changes": _sign_changes(Mpp),
    }
    if breakpoints:
        jumps = []
        for xb in breakpoints:
            e = 1e-6
            m_l = float(law.deriv(np.array([xb - e]), 3)[0])
            m_r = float(law.deriv(np.array([xb + e]), 3)[0])
            jumps.append({"x": float(xb), "Mppp_minus": m_l, "Mppp_plus": m_r,
                          "jump": m_r - m_l})
        out["breakpoint_jumps"] = jumps
        out["max_breakpoint_jump"] = float(max(abs(j["jump"]) for j in jumps))
    return out


def _sign_changes(v: np.ndarray) -> int:
    sg = np.sign(v)
    sg = sg[np.abs(v) > 1e-10 * max(1.0, float(np.max(np.abs(v))))]
    return int(np.sum(np.diff(sg) != 0.0)) if len(sg) > 1 else 0


# --- 特性線網トポロジ ----------------------------------------------------------
def characteristic_topology_diagnostics(levels: np.ndarray, n_ax: int,
                                        min_spacing_frac: float = 1e-4,
                                        n_crossing_sample: int = 40,
                                        seam_i_pad: int = 45,
                                        seam_k_max: int = 45,
                                        seam_x_max: float | None = None,
                                        wall: np.ndarray | None = None,
                                        wall_tol: float = 0.02) -> dict:
    r"""`fill_levels` の三角充填網 (`levels[k,i]=[x,r,θ,ν,M]`, NaN=欠損) を直接検査する。

    - **signed area**: 三角形 $(L_{k-1}[i], L_{k-1}[i+1], L_k[i])$ の符号付き面積
      (x-r 平面の外積)。全三角形が同符号なら向きが揃っている (fold なし)。
    - **orientation flip**: 隣接三角形間で符号が変わる箇所の数。
    - **退化セル**: 面積が (基準スケール比で) ほぼ 0 の三角形。
    - **隣接特性線間隔のゼロ化**: 各レベル内の隣接点間距離が場のスケールに対し
      異常に小さくなる箇所 (特性線が収束/交差する兆候)。
    - **点堆積**: 同一レベル内で非隣接な点対が近接 (多価の兆候)。
    - **非隣接 C±線の交差**: サンプルした C⁺ 線と C⁻ 線 (`levels` の対角線 / 行) の
      折れ線同士の交差を、隣接しない組についてだけ調べる (計算量を抑えるため間引き)。

    `θ<μ` の代用診断ではなく、網そのものの構造を見る。

    **既知の良性境界 (seam)**: 初期値線は「軸目標配列 (index 0..n_ax-1)」+
    「スロート特性線配列 (index n_ax..)」の L 字構成 (`inverse_design`)。この 2 つの
    配列種別が切り替わる添字境界 $i\approx n_{\rm ax}$ 付近・レベル $k$ が小さい
    (充填の初期数段) 領域は、担体の種類が変わる幾何的な継ぎ目であり、健全な設計
    (A, 2026-08-16 実測 M6/R3/L_c45) でも±数百個の符号反転が**そこだけに**局在する
    (fold ではない)。既定でこの帯 ($|i-n_{\rm ax}|\le$`seam_i_pad`+初期値線点数、
    $k\le$`seam_k_max`) を除いた「内部反転数」も別途報告する — 除外は診断結果を
    甘くするためではなく、既知の非物理境界と真の fold を区別するため。"""
    n_lev, n_pt, _ = levels.shape
    xs, rs = np.ascontiguousarray(levels[..., 0]), np.ascontiguousarray(levels[..., 1])
    finite = np.isfinite(xs) & np.isfinite(rs)
    scale = float(np.nanmax(rs)) if np.any(finite) else 1.0
    scale = max(scale, 1e-9)

    # 三角形 (L[k-1,i], L[k-1,i+1], L[k,i]) の符号付き面積 (完全ベクトル化 — n_lev x n_pt
    # 規模のループを Python で回すと O(n^2) 点数で数分かかる、実測)。網の半分 (上三角) は NaN なので、3 頂点が有限の三角形だけを
    # 行優先の順 (以前の `area2[valid]` と同じ順) で詰めて計算する (2026-10-15、plan tooling-nozzle-mknot-fast §4.6 (2))。
    # 以下の配列 (Ax・area2・interior …) は全部この詰めた並び。
    tk, ti = np.nonzero(finite[:-1, :-1] & finite[:-1, 1:] & finite[1:, :-1])
    fa = tk * n_pt + ti                                       # 頂点 A の平坦な添字 (B = +1、P = +n_pt)
    xf, rf = xs.ravel(), rs.ravel()
    Ax, Ar = xf.take(fa), rf.take(fa)
    Bx, Br = xf.take(fa + 1), rf.take(fa + 1)
    Px, Pr = xf.take(fa + n_pt), rf.take(fa + n_pt)
    area2 = (Bx - Ax) * (Pr - Ar) - (Br - Ar) * (Px - Ax)
    # 内部 (seam 帯を除く) マスク: i が軸/初期値線境界から十分離れているか k が十分大 (k = 新点 P の段)
    interior = (np.abs(ti - n_ax) > seam_i_pad) | (tk + 1 > seam_k_max)
    if wall is not None:
        # 三角充填網は壁の外側 (壁流線より上の r) まで伸びる (逆 MOC は壁を後から閉包する)。
        # 壁外の三角形はノズル流に無関係なので、`wall` (n,>=2)[x,r] を渡した場合は
        # 頂点 A が壁の内側 (r_A <= r_wall(x_A)(1+wall_tol)) の三角形だけを interior 判定と
        # 余裕指標に使う (D 案計画で追加 — 余裕の最小値が壁外の r≈18 の三角形で決まっていた)。
        w = np.asarray(wall, dtype=float)
        rw = np.interp(Ax, w[:, 0], w[:, 1], left=w[0, 1], right=w[-1, 1])
        inside = np.isfinite(Ax) & (Ar <= rw * (1.0 + wall_tol))
        interior &= inside
    if seam_x_max is not None:
        # index ベースの seam 窓は解像度が上がるほど物理的に同じ近傍を覆うのに必要な
        # index 数が増える (実測: M6/R3/L_c45 で n_axis 1200→2400 でも x<0.22 r_t の
        # 同じ近傍に留まる反転が index 窓の外へはみ出す)。物理座標 (x < seam_x_max)
        # にいる三角形は index 窓に関わらず一律 seam 扱いとする (AND ではなく上書き)。
        interior &= (Ax >= seam_x_max) & (Bx >= seam_x_max) & (Px >= seam_x_max)
    areas = area2
    areas_int = area2[interior]
    if len(areas) == 0:
        return {"n_triangles": 0, "violations": ["特性線網が空 (三角形 0 個)"]}
    sign = np.sign(areas)
    nz = sign[np.abs(areas) > 1e-12 * scale * scale]
    n_flip = int(np.sum(np.diff(nz) != 0.0)) if len(nz) > 1 else 0
    sign_int = np.sign(areas_int)
    nz_int = sign_int[np.abs(areas_int) > 1e-12 * scale * scale]
    n_flip_interior = int(np.sum(np.diff(nz_int) != 0.0)) if len(nz_int) > 1 else 0
    degenerate_frac = float(np.mean(np.abs(areas) < 1e-6 * scale * scale))
    dominant_sign = float(np.sign(np.median(areas[np.abs(areas) > 1e-9 * scale * scale]))) \
        if np.any(np.abs(areas) > 1e-9 * scale * scale) else 0.0
    wrong_sign_frac = float(np.mean((sign != 0) & (sign != dominant_sign)
                                    & (np.abs(areas) > 1e-9 * scale * scale)))
    wrong_sign_frac_interior = float(np.mean((sign_int != 0) & (sign_int != dominant_sign)
                                             & (np.abs(areas_int) > 1e-9 * scale * scale))) \
        if len(areas_int) else 0.0

    # 隣接特性線間隔 (同一レベル内の隣接点距離) の異常な縮小と、点堆積 (同一レベル内で非隣接点対が近接しているか、隣接除外)。
    # 点堆積は raw と seam 除外 (interior) の両方を数える — 粗い n_axis (600) では seam 内 (throat 直後 x<0.3, index 境界) にだけ
    # 点堆積が出て健全な A でも 180 対に達する (実測)。段ごとの走査は numba の 1 関数 (`moc_fill_nb.level_scan_nb`; 以前は段ごとの
    # cKDTree で 1 評価 約 0.5 s、2026-10-15 plan tooling-nozzle-mknot-fast §4.6 (2))。判定は同じ、位置の一覧は段ごとに対の添字の順。
    from .moc_fill_nb import level_scan_nb
    min_gap, gk, gi, pileup_pairs, pileup_pairs_interior, locs = level_scan_nb(
        xs, rs, int(n_ax), int(seam_i_pad), int(seam_k_max),
        float(seam_x_max) if seam_x_max is not None else 0.0, seam_x_max is not None, 20)
    min_gap = float(min_gap)
    min_gap_loc = (int(gk), int(gi)) if np.isfinite(min_gap) else None
    pileup_pairs, pileup_pairs_interior = int(pileup_pairs), int(pileup_pairs_interior)
    pileup_locs = [(int(a), float(b), float(c)) for a, b, c in locs]

    # 非隣接 C+/C- 折れ線の交差 (サンプル)
    n_crossings = _sample_characteristic_crossings(levels, n_ax, n_crossing_sample)

    # トポロジ余裕 (D 案計画): 内部三角形の |signed area| を代表値 (中央値) で規格化した最小値
    # = 反転 (面積 0 を跨ぐ) までの余裕。1 に近いほど一様な網、0 に近いほど退化寸前。
    # ただし面積比は格子伸縮 (axis_dx0 の等比細分・軸近傍の r→0) に支配されるので、
    # スケール不変な **signed Jacobian 相当量** = signed area / (|AB|·|AP|) = 頂点 A での
    # sin(角) も併記する (これが 0 → 退化、負 → 反転。反転までの角度余裕として読める)。
    med_int = float(np.median(np.abs(areas_int))) if len(areas_int) else 0.0
    if med_int > 0 and len(areas_int):
        signed_rel = dominant_sign * areas_int / med_int
        min_signed_area_rel = float(np.min(signed_rel))
    else:
        min_signed_area_rel = float("nan")
    q_int = np.flatnonzero(interior)
    Axi, Ari = Ax[q_int], Ar[q_int]
    lenAB = np.hypot(Bx[q_int] - Axi, Br[q_int] - Ari)
    lenAP = np.hypot(Px[q_int] - Axi, Pr[q_int] - Ari)
    with np.errstate(divide="ignore", invalid="ignore"):
        sinA_int = dominant_sign * areas_int / np.maximum(lenAB * lenAP, 1e-300)
    if len(sinA_int):
        j = int(np.argmin(sinA_int))
        min_sin_angle_int = float(sinA_int[j])
        q = int(q_int[j])
        min_sin_angle_loc = (int(tk[q]), int(ti[q]), float(Ax[q]), float(Ar[q]))
    else:
        min_sin_angle_int = float("nan"); min_sin_angle_loc = None
    # seam 内の反転三角形の物理座標分布 (raw と seam 除外の差分を後から検証できるように)
    seam_wrong = ~interior & (np.sign(area2) != 0) & (np.sign(area2) != dominant_sign) \
        & (np.abs(area2) > 1e-9 * scale * scale)
    if np.any(seam_wrong):
        seam_x = Ax[seam_wrong]; seam_r = Ar[seam_wrong]
        seam_flip_dist = {"n": int(seam_wrong.sum()),
                          "x_min": float(np.nanmin(seam_x)), "x_max": float(np.nanmax(seam_x)),
                          "r_min": float(np.nanmin(seam_r)), "r_max": float(np.nanmax(seam_r))}
    else:
        seam_flip_dist = {"n": 0}

    v = []
    if n_flip_interior > 0:
        v.append(f"三角形の向き反転 (seam 除外・内部) {n_flip_interior} 箇所 (fold の兆候)")
    if wrong_sign_frac_interior > 1e-4:
        v.append(f"支配符号と逆の三角形 (内部) が {100*wrong_sign_frac_interior:.3f}% ある")
    if pileup_pairs_interior > 0:
        v.append(f"非隣接点の近接 (点堆積, seam 除外) {pileup_pairs_interior} 組")
    if n_crossings > 0:
        v.append(f"非隣接 C±線の交差 {n_crossings} 箇所")
    return {
        "n_triangles": int(len(areas)), "n_orientation_flip": n_flip,
        "n_orientation_flip_interior": n_flip_interior,
        "degenerate_fraction": degenerate_frac, "wrong_sign_fraction": wrong_sign_frac,
        "wrong_sign_fraction_interior": wrong_sign_frac_interior,
        "min_adjacent_spacing": float(min_gap) if np.isfinite(min_gap) else None,
        "min_adjacent_spacing_rel": (float(min_gap) / scale) if np.isfinite(min_gap) else None,
        "min_adjacent_spacing_loc": min_gap_loc,
        "pileup_pairs": pileup_pairs, "pileup_pairs_interior": pileup_pairs_interior,
        "pileup_locs_interior": pileup_locs, "n_crossings_sampled": n_crossings,
        "min_signed_area_rel_interior": min_signed_area_rel,
        "min_sin_angle_interior": min_sin_angle_int,
        "min_sin_angle_interior_loc": min_sin_angle_loc,
        "seam_flip_distribution": seam_flip_dist,
        "seam_params": {"seam_i_pad": int(seam_i_pad), "seam_k_max": int(seam_k_max),
                        "seam_x_max": (None if seam_x_max is None else float(seam_x_max))},
        "violations": v,
    }


def _sample_characteristic_crossings(levels: np.ndarray, n_ax: int, n_sample: int) -> int:
    """C⁺ 線 (levels の反対角線, `moc_inverse.cplus_lines` と同じ添字規約) を
    間引いてサンプルし、非隣接な線同士の折れ線交差を数える。バウンディングボックスで
    まず粗く篩い、通った組だけ全対 segment を numpy でベクトル化して判定する
    (Python 二重ループで直接判定すると n_axis~1000 規模で数分かかる、実測)。"""
    from .moc_inverse import cplus_lines
    n = levels.shape[1]
    ms = np.unique(np.linspace(0, n - 1, min(n_sample, n)).astype(int))
    lines = []
    for m in ms:
        ln = cplus_lines(levels, int(m))
        if len(ln) >= 2:
            lines.append(ln[:, :2])
    boxes = [(l[:, 0].min(), l[:, 0].max(), l[:, 1].min(), l[:, 1].max()) for l in lines]
    count = 0
    for a in range(len(lines)):
        xa0, xa1, ra0, ra1 = boxes[a]
        for b in range(a + 2, len(lines)):        # 隣接 (b=a+1) は共有点で必ず接するので除外
            xb0, xb1, rb0, rb1 = boxes[b]
            if xa1 < xb0 or xb1 < xa0 or ra1 < rb0 or rb1 < ra0:
                continue                            # バウンディングボックス非重複
            if _polylines_cross_aabb(lines[a], lines[b]):
                count += 1
    return count


def _polylines_cross_aabb(p: np.ndarray, q: np.ndarray, chunk_size: int = 256) -> bool:
    """折れ線 p, q の交差判定。AABB で候補線分対を絞ってから従来と同じ外積を評価する。

    どちらか一方が x 単調なら、その線分 x 区間が順序付けられることを使い、もう一方の各線分と
    x 区間が重なる添字範囲を `searchsorted` で求める。両方が x 非単調なら、AABB 行列だけを
    チャンク生成する汎用経路へ落とす。`_polylines_cross_vec` と交差の定義は同一。"""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if len(p) < 2 or len(q) < 2:
        return False

    def monotone_x(v: np.ndarray) -> bool:
        dx = np.diff(v[:, 0])
        return bool(np.all(dx >= 0.0) or np.all(dx <= 0.0))

    p_mono = monotone_x(p)
    q_mono = monotone_x(q)
    if p_mono or q_mono:
        # q を x 単調側にする。両方単調なら p を短い側にして searchsorted の照会数を減らす。
        if (not q_mono) or (p_mono and len(p) > len(q)):
            p, q = q, p
        if q[-1, 0] < q[0, 0]:
            q = q[::-1]
        return _polylines_cross_monotone_q(p, q, chunk_size)
    return _polylines_cross_aabb_chunks(p, q, chunk_size)


def _polylines_cross_monotone_q(p: np.ndarray, q: np.ndarray, chunk_size: int) -> bool:
    """q の x 座標が非減少である場合の AABB sweep。"""
    p1, p2 = p[:-1], p[1:]
    q1, q2 = q[:-1], q[1:]
    qx_lo = np.minimum(q1[:, 0], q2[:, 0])
    qx_hi = np.maximum(q1[:, 0], q2[:, 0])
    qr_lo = np.minimum(q1[:, 1], q2[:, 1])
    qr_hi = np.maximum(q1[:, 1], q2[:, 1])

    for i0 in range(0, len(p1), chunk_size):
        i1 = min(i0 + chunk_size, len(p1))
        pa, pb = p1[i0:i1], p2[i0:i1]
        px_lo = np.minimum(pa[:, 0], pb[:, 0])
        px_hi = np.maximum(pa[:, 0], pb[:, 0])
        # qx_hi >= px_lo かつ qx_lo <= px_hi となる閉区間の範囲。
        j0 = np.searchsorted(qx_hi, px_lo, side="left")
        j1 = np.searchsorted(qx_lo, px_hi, side="right")
        counts = np.maximum(j1 - j0, 0)
        n_pairs = int(np.sum(counts))
        if n_pairs == 0:
            continue
        ip = np.repeat(np.arange(len(pa), dtype=np.intp), counts)
        group_start = np.repeat(np.cumsum(counts) - counts, counts)
        jq = np.repeat(j0, counts) + np.arange(n_pairs, dtype=np.intp) - group_start
        pr_lo = np.minimum(pa[ip, 1], pb[ip, 1])
        pr_hi = np.maximum(pa[ip, 1], pb[ip, 1])
        keep = (pr_hi >= qr_lo[jq]) & (qr_hi[jq] >= pr_lo)
        if np.any(keep) and np.any(_segments_cross_pairs(
                pa[ip[keep]], pb[ip[keep]], q1[jq[keep]], q2[jq[keep]])):
            return True
    return False


def _polylines_cross_aabb_chunks(p: np.ndarray, q: np.ndarray, chunk_size: int) -> bool:
    """両折れ線が x 非単調な場合の汎用 AABB チャンク判定。"""
    p1, p2 = p[:-1], p[1:]
    q1, q2 = q[:-1], q[1:]
    qx_lo = np.minimum(q1[:, 0], q2[:, 0])
    qx_hi = np.maximum(q1[:, 0], q2[:, 0])
    qr_lo = np.minimum(q1[:, 1], q2[:, 1])
    qr_hi = np.maximum(q1[:, 1], q2[:, 1])
    for i0 in range(0, len(p1), chunk_size):
        i1 = min(i0 + chunk_size, len(p1))
        pa, pb = p1[i0:i1], p2[i0:i1]
        px_lo = np.minimum(pa[:, 0], pb[:, 0])[:, None]
        px_hi = np.maximum(pa[:, 0], pb[:, 0])[:, None]
        pr_lo = np.minimum(pa[:, 1], pb[:, 1])[:, None]
        pr_hi = np.maximum(pa[:, 1], pb[:, 1])[:, None]
        overlap = (px_hi >= qx_lo[None, :]) & (qx_hi[None, :] >= px_lo) \
            & (pr_hi >= qr_lo[None, :]) & (qr_hi[None, :] >= pr_lo)
        ip, jq = np.nonzero(overlap)
        if len(ip) and np.any(_segments_cross_pairs(pa[ip], pb[ip], q1[jq], q2[jq])):
            return True
    return False


def _polylines_cross_vec(p: np.ndarray, q: np.ndarray) -> bool:
    """参照実装: 折れ線 p, q の全 segment 対を一括して交差判定。"""
    p1, p2 = p[:-1], p[1:]
    q1, q2 = q[:-1], q[1:]
    return bool(_segments_cross_matrix(p1, p2, q1, q2).any())


def _segments_cross_pairs(p1: np.ndarray, p2: np.ndarray, q1: np.ndarray,
                          q2: np.ndarray) -> np.ndarray:
    """同じ長さの線分配列について要素ごとの交差判定を返す。"""
    def cross2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]

    d1 = cross2(q2 - q1, p1 - q1)
    d2 = cross2(q2 - q1, p2 - q1)
    d3 = cross2(p2 - p1, q1 - p1)
    d4 = cross2(p2 - p1, q2 - p1)
    return ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))


def _segments_cross_matrix(p1: np.ndarray, p2: np.ndarray, q1: np.ndarray,
                           q2: np.ndarray) -> np.ndarray:
    """(Np,2)x2 と (Nq,2)x2 の全組合せ交差判定 → (Np,Nq) 真偽値行列。"""
    def cross2(ax, ay, bx, by):
        return ax * by - ay * bx

    d1 = cross2(q2[None, :, 0] - q1[None, :, 0], q2[None, :, 1] - q1[None, :, 1],
               p1[:, None, 0] - q1[None, :, 0], p1[:, None, 1] - q1[None, :, 1])
    d2 = cross2(q2[None, :, 0] - q1[None, :, 0], q2[None, :, 1] - q1[None, :, 1],
               p2[:, None, 0] - q1[None, :, 0], p2[:, None, 1] - q1[None, :, 1])
    d3 = cross2(p2[:, None, 0] - p1[:, None, 0], p2[:, None, 1] - p1[:, None, 1],
               q1[None, :, 0] - p1[:, None, 0], q1[None, :, 1] - p1[:, None, 1])
    d4 = cross2(p2[:, None, 0] - p1[:, None, 0], p2[:, None, 1] - p1[:, None, 1],
               q2[None, :, 0] - p1[:, None, 0], q2[None, :, 1] - p1[:, None, 1])
    return ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))


# --- 壁曲率 (spline 解析微分、生の数値差分でない) -------------------------------
def wall_curvature_diagnostics(wall_obj, x_lo: float, x_hi: float, n: int = 4001,
                               x_split: float | None = 0.25) -> dict:
    r"""$\kappa_w(x)=r''/(1+r'^2)^{3/2}$、$d\kappa_w/ds=(d\kappa_w/dx)/\sqrt{1+r'^2}$
    ($ds=\sqrt{1+r'^2}\,dx$) を `wall_obj.r(x, deriv)` (`AxisMachCFDWall` の補間 5 次
    B-spline、既に生データへの当て直し = 数値差分でない) の解析微分から評価する。

    `x_split` (既定 0.25 $r_t$): **throat 窓 $[x_{lo}, x_{split}]$ と下流 $(x_{split}, x_{hi}]$ に
    分けて**も返す (D 案計画 2026-08-16)。A/B 比較で $\max|d\kappa/ds|$ が $x=0$ に出ており、
    スロート端のクランプ補間に支配されて内部軸則の差を隠す可能性があるため、軸則比較の
    目的関数には原則 throat 窓**外** (`downstream_*`) を使う。"""
    x = np.linspace(x_lo, x_hi, n)
    rp = wall_obj.r(x, 1)
    rpp = wall_obj.r(x, 2)
    rppp = wall_obj.r(x, 3)
    kappa = rpp / (1.0 + rp * rp) ** 1.5
    dkappa_dx = rppp / (1.0 + rp * rp) ** 1.5 - 3.0 * rp * rpp * rpp / (1.0 + rp * rp) ** 2.5
    ds_dx = np.sqrt(1.0 + rp * rp)
    dkappa_ds = dkappa_dx / ds_dx
    J = float(np.trapezoid(dkappa_ds ** 2 * ds_dx, x))
    out = {
        "max_abs_kappa": float(np.max(np.abs(kappa))),
        "max_abs_dkappa_ds": float(np.max(np.abs(dkappa_ds))),
        "J_dkappa_ds2": J,
        "x_at_max_dkappa_ds": float(x[int(np.argmax(np.abs(dkappa_ds)))]),
    }
    if x_split is not None:
        for tag, m in (("throat", x <= x_split), ("downstream", x > x_split)):
            if int(m.sum()) < 3:
                continue
            out[f"{tag}_max_abs_kappa"] = float(np.max(np.abs(kappa[m])))
            out[f"{tag}_max_abs_dkappa_ds"] = float(np.max(np.abs(dkappa_ds[m])))
            out[f"{tag}_J_dkappa_ds2"] = float(np.trapezoid(dkappa_ds[m] ** 2 * ds_dx[m], x[m]))
            out[f"{tag}_x_at_max_dkappa_ds"] = float(x[m][int(np.argmax(np.abs(dkappa_ds[m])))])
        out["x_split"] = float(x_split)
    return out


def wall_margin_diagnostics(wall_tbl: np.ndarray) -> dict:
    r"""$\min(\mu_w-\theta_w)$ (壁マッハ角 − 壁角の最小値、逆 MOC 位相条件への
    補助的な余裕指標。健全性の一次判定には使わない — §topology 診断を参照)。"""
    th = wall_tbl[:, 2]
    M = wall_tbl[:, 3]
    mu = np.arcsin(1.0 / np.maximum(M, 1.0 + 1e-12))
    margin = mu - th
    i = int(np.argmin(margin))
    return {"min_mu_minus_theta_deg": float(np.degrees(margin[i])),
            "x_at_min_margin": float(wall_tbl[i, 0]),
            "theta_deg_at_min_margin": float(np.degrees(th[i])),
            "M_at_min_margin": float(M[i])}
