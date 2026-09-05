"""⑤ SERN の平面 2D 構造 quad メッシュ (2 バンド + スリットカウル、決定的)。

plan: plans/active/tooling-nozzle-sern-chain.md §4.7 / S2。無次元 (H = 1) → `scale` [m]。

配置 (x 右、y 上):
    上境界 = ランプ (x<0 は直管 y=1、0..L_ramp は設計輪郭、以降は角 top_ext_angle の直線 = outlet 扱い)
    中間線 = カウル (x ≤ L_cowl: x<0 は y=0 の直管、0..L_cowl は −x tanθ_c0) → TE 以降は角 interface_angle の
             直線 (せん断層に沿わせる格子線。境界条件なし)
    下境界 = 遠方 (y = y_mid(x_out) − bot_depth)
バンド: 下 (下境界→中間線, nj_bot 点) / 上 (中間線→ランプ, nj_top 点)。x station は共通。
カウル (x ≤ L_cowl) は**スリット**: 中間線ノードを上下 2 重にし、上側 = cowl_in、下側 = cowl_out。
TE 以降の中間線ノードは共有 (内部線)。
**カウル板厚 `cowl_thickness`** (2026-09-04): 厚さ 0 だと上下の壁ノードが同一座標になり、node モードでは双子の
壁ノードが同じ状態 (片側は排気 15 kPa、反対側は外部流 2 kPa なのに同一 p) になって 2 次で発散した
(case/46 run_0008/0009)。厚さ t(x) を「入口から t、TE 手前 20 % で 0 に絞る」形で入れ、t>0 の station だけ 2 重ノードにする (TE は共有。
入口側を 0 に絞ると 4 境界を持つ 1 ノードができて発散するので入口では厚さを保つ)。壁第一セルは上下とも first_wall_frac·H の絶対値。
物理タグ: PHYS_SERN。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mesh2d import _radial_fracs

PHYS_SERN = {"inlet_nozzle": 1, "inlet_ext": 2, "outlet": 3, "ramp": 4, "cowl_in": 5,
             "cowl_out": 6, "bottom": 7, "top_out": 8, "fluid": 9, "vehicle": 10}


@dataclass
class SernMeshParams:
    ni_up: int = 16
    ni_noz: int = 120
    ni_plume: int = 220
    nj_top: int = 101
    nj_bot: int = 61
    L_up: float = 0.5
    x_out_extra: float = 2.0
    bot_depth: float = 3.0
    first_wall_frac: float = 2.0e-3      # 壁第一セル / H (上下バンドとも絶対値で同じ)
    first_bot_frac: float = 0.02         # (未使用)
    cowl_thickness: float = 0.0          # カウル板厚 / H (0 = 厚さ 0 のスリット)。node では 2e-3 を推奨 (docstring)
    interface_angle: float = 0.0
    top_ext_angle: float = 0.0
    scale: float = 1.0
    x_cluster_w: float = 0.15
    x_cluster_a: float = 3.0
    # --- ランプ側外部流ブロック (plan §4.11) ---
    ext_top: bool = False                # 機体上面・base・後流 + 自由流バンドを付ける
    top_depth: float = 2.0               # 機体上面線から上境界までの高さ / H
    nj_ext_top: int = 41
    nj_wake: int = 9                     # base 高さ分の後流ブロックの j 点数
    vehicle_clearance: float = 0.02      # 機体上面 = max(ランプ y) + これ (/H)
    first_top_frac: float = 0.02         # top バンドの第一セル / H (slip 壁なので粗くてよい)
    ramp_fillet: float = 0.0             # ランプ膨張角部 (x=0) の丸め半径 / H。0 = 鋭角 (従来)。
                                         # 鋭角だと SST が θ_r0 ≳ 18° で `roOmega` 発散する (case/46 run_0055-0057)。
                                         # NASA TM X-71972 も "linear segments joined by small radii" と記す
    vehicle_taper: float = 0.0           # >0: 機体後端のテーパ長。**L_ramp に対する比**で与える (0.2 = 後縁手前 20 %)。
                                         # 絶対長 (H) にしていた 2026-09-05 版は短ランプ設計で破綻した: L_ramp 4.29 に対し
                                         # taper 2.0 H が 47 % を占め、θ_e>0 でランプが後縁まで上がるため機体厚が
                                         # 入口から clearance と同オーダーに潰れ、first_top_frac と同じ大きさのセルが並んで
                                         # ω が発散した (run_0071 inf_00/inf_01 の m10_on)。
                                         # 0: 鉛直 base + wake ブロック。node では base の 90° 二重 slip 角で step 5 発散 (run_0034)


def _cluster_stations(x0, x1, n, ends=(True, True), w=0.15, a=3.0):
    xi = np.linspace(0.0, 1.0, 4001)
    dens = np.ones_like(xi)
    L = x1 - x0
    wn = max(w / L, 1e-3)
    if ends[0]:
        dens += (a - 1.0) * np.exp(-(xi / wn) ** 2)
    if ends[1]:
        dens += (a - 1.0) * np.exp(-((1.0 - xi) / wn) ** 2)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xi))])
    cum /= cum[-1]
    return x0 + L * np.interp(np.linspace(0.0, 1.0, n), cum, xi)


def _geom_start(n, first):
    """[0,1] を n 点、始端の第一間隔が first (比) になる幾何級数で切る (first ≥ 1/(n−1) なら一様)。"""
    if first >= 1.0 / (n - 1):
        return np.linspace(0.0, 1.0, n)
    lo, hi = 1.0, 3.0
    for _ in range(100):                       # first·(r^(n−1) − 1)/(r − 1) = 1 の r を二分法で
        r = 0.5 * (lo + hi)
        tot = first * (r ** (n - 1) - 1.0) / (r - 1.0)
        lo, hi = (r, hi) if tot < 1.0 else (lo, r)
    r = 0.5 * (lo + hi)
    s = np.concatenate([[0.0], np.cumsum(first * r ** np.arange(n - 1))])
    s /= s[-1]
    return s


def _tanh_two_sided(n, first):
    xi = np.linspace(0.0, 1.0, n)

    def s_of(d):
        return 0.5 * (1.0 + np.tanh(d * (xi - 0.5)) / np.tanh(0.5 * d))
    lo, hi = 1e-3, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if s_of(mid)[1] > first:
            lo = mid
        else:
            hi = mid
    s = s_of(0.5 * (lo + hi))
    s[0], s[-1] = 0.0, 1.0
    return s


def generate_sern_mesh(design, prm: SernMeshParams):
    """design: SernDesign (無次元)。戻り値 (coords [m], quads, bedges{name: edges}, info)。"""
    L_cowl = float(design.cowl_xy[-1, 0])
    y_te = float(design.cowl_xy[-1, 1])
    tan_c = -y_te / L_cowl if L_cowl > 0 else 0.0
    L_ramp = float(design.L_ramp)
    y_e = float(design.ramp_xy[-1, 1])
    x_out = L_ramp + prm.x_out_extra
    # 丸め区間 [xf1, xf2] は station を「追加」せず**ブロック境界**にする (追加すると既存点と 1e-9 まで接近して AR が飛ぶ)
    Rf0 = float(prm.ramp_fillet)
    if Rf0 > 0.0 and len(design.ramp_xy) > 1:
        _th = float(np.arctan2(design.ramp_xy[1, 1] - design.ramp_xy[0, 1], design.ramp_xy[1, 0] - design.ramp_xy[0, 0]))
        _t = Rf0 * np.tan(0.5 * _th)
        f1, f2, nf = -_t, _t * np.cos(_th), 7
        xs = np.concatenate([
            _cluster_stations(-prm.L_up, f1, max(prm.ni_up - 6, 4), (False, True), prm.x_cluster_w, prm.x_cluster_a),
            np.linspace(f1, 0.0, nf)[1:],
            np.linspace(0.0, f2, nf)[1:],
            _cluster_stations(f2, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
            _cluster_stations(L_cowl, x_out, prm.ni_plume, (True, False), prm.x_cluster_w, prm.x_cluster_a)[1:],
        ])
    else:
        xs = np.concatenate([
            _cluster_stations(-prm.L_up, 0.0, prm.ni_up, (False, True), prm.x_cluster_w, prm.x_cluster_a),
            _cluster_stations(0.0, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
            _cluster_stations(L_cowl, x_out, prm.ni_plume, (True, False), prm.x_cluster_w, prm.x_cluster_a)[1:],
        ])
    # ランプ後縁に station を置く (最寄りを置換)
    k = int(np.argmin(np.abs(xs - L_ramp)))
    xs[k] = L_ramp
    i_te = int(np.argmin(np.abs(xs - L_cowl)))
    assert abs(xs[i_te] - L_cowl) < 1e-12
    rx, ry = design.ramp_xy[:, 0], design.ramp_xy[:, 1]
    # ランプ膨張角部 (x=0, y=1) の丸め: 直管 (勾配 0) と設計輪郭の初期勾配 tanθ を半径 R の円弧で接続する。
    # 接点は x1 = −t (直管側), x2 = t cosθ (輪郭側), t = R tan(θ/2)、中心 (−t, 1+R)。壁は流体と反対側へ膨らむ。
    R_f = float(prm.ramp_fillet)
    if R_f > 0.0 and len(rx) > 1:
        th_r0 = float(np.arctan2(ry[1] - ry[0], rx[1] - rx[0]))
        t_f = R_f * np.tan(0.5 * th_r0)
        xf1, xf2 = -t_f, t_f * np.cos(th_r0)
    else:
        R_f = 0.0; xf1 = xf2 = 0.0

    def y_top(x):
        x = np.asarray(x, dtype=float)
        y = np.where(x < 0.0, 1.0,
                     np.where(x <= L_ramp, np.interp(x, rx, ry), y_e + (x - L_ramp) * np.tan(prm.top_ext_angle)))
        if R_f > 0.0:
            u = x + t_f
            arc = 1.0 + R_f - np.sqrt(np.clip(R_f * R_f - u * u, 0.0, None))
            y = np.where((x >= xf1) & (x <= xf2), np.maximum(y, arc), y)
        return y

    def y_mid(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 0.0,
                        np.where(x <= L_cowl, -x * tan_c, y_te + (x - L_cowl) * np.tan(prm.interface_angle)))

    y_bot = float(y_mid(x_out)) - prm.bot_depth
    yt, ym = y_top(xs), y_mid(xs)
    # カウル板厚 t(x): 入口 (−L_up) から 0.8 L_cowl まで t、TE で 0 に絞る。TE 以降 0。入口側で 0 に絞ると i=0 が
    # 4 種の境界 (inlet×2 + wall×2) を持つ 1 ノードになり soft 段 step 4 で発散した (run_0014) ので入口では厚さを保つ
    t_c = float(prm.cowl_thickness)
    tk = np.zeros_like(xs)
    if t_c > 0.0:
        tk = np.interp(xs, [-prm.L_up, 0.8 * L_cowl, L_cowl], [t_c, t_c, 0.0], left=t_c, right=0.0)
    dup = tk > 1e-15                      # 2 重ノードにする station (厚さ 0 なら i < i_te 全部)
    if t_c <= 0.0:
        dup = np.arange(len(xs)) < i_te
    ym_up, ym_lo = ym + 0.5 * tk, ym - 0.5 * tk
    ni, njt, njb = len(xs), prm.nj_top, prm.nj_bot
    dup_idx = {int(i): k for k, i in enumerate(np.where(dup)[0])}
    N_low = ni * njb
    N_up = ni * (njt - 1)
    N_dup = len(dup_idx)
    coords = np.zeros((N_low + N_up + N_dup, 3))

    def low(i, j):
        return i * njb + j

    def up(i, j):
        if j == 0:
            return N_low + N_up + dup_idx[i] if i in dup_idx else low(i, njb - 1)
        return N_low + i * (njt - 1) + (j - 1)

    for i in range(ni):
        # 壁第一セルは上下とも first_wall_frac·H (絶対値)。バンド高さで割った比に直す
        h_lo = max(ym_lo[i] - y_bot, 1e-12); h_up = max(yt[i] - ym_up[i], 1e-12)
        s_bot = _radial_fracs(njb, min(prm.first_wall_frac / h_lo, 0.5 / (njb - 1)))
        s_top = _tanh_two_sided(njt, min(prm.first_wall_frac / h_up, 0.5 / (njt - 1)))
        yb = y_bot + s_bot * (ym_lo[i] - y_bot)
        coords[low(i, 0):low(i, njb), 0] = xs[i]
        coords[low(i, 0):low(i, njb), 1] = yb
        yu = ym_up[i] + s_top * (yt[i] - ym_up[i])
        coords[up(i, 1):up(i, njt - 1) + 1, 0] = xs[i]
        coords[up(i, 1):up(i, njt - 1) + 1, 1] = yu[1:]
        if i in dup_idx:
            coords[up(i, 0)] = (xs[i], ym_up[i], 0.0)
    quads = []
    for i in range(ni - 1):
        for j in range(njb - 1):
            quads.append((low(i, j), low(i + 1, j), low(i + 1, j + 1), low(i, j + 1)))
        for j in range(njt - 1):
            quads.append((up(i, j), up(i + 1, j), up(i + 1, j + 1), up(i, j + 1)))
    quads = np.asarray(quads, dtype=np.int64)
    bedges = {n: [] for n in PHYS_SERN if n != "fluid"}
    for i in range(ni - 1):
        bedges["bottom"].append((low(i, 0), low(i + 1, 0)))
        xm = 0.5 * (xs[i] + xs[i + 1])
        if xm <= L_ramp:
            bedges["ramp"].append((up(i, njt - 1), up(i + 1, njt - 1)))
        elif not prm.ext_top:                     # ext_top ではプルーム上線は wake/top との内部線
            bedges["top_out"].append((up(i, njt - 1), up(i + 1, njt - 1)))
        if i + 1 <= i_te:
            bedges["cowl_in"].append((up(i, 0), up(i + 1, 0)))
            bedges["cowl_out"].append((low(i, njb - 1), low(i + 1, njb - 1)))
    for j in range(njb - 1):
        bedges["inlet_ext"].append((low(0, j), low(0, j + 1)))
        bedges["outlet"].append((low(ni - 1, j), low(ni - 1, j + 1)))
    for j in range(njt - 1):
        bedges["inlet_nozzle"].append((up(0, j), up(0, j + 1)))
        bedges["outlet"].append((up(ni - 1, j), up(ni - 1, j + 1)))
    ext = {"ext_top": bool(prm.ext_top)}
    if prm.ext_top:
        coords, quads = _add_ext_top(coords, quads, bedges, xs, yt, k, L_ramp, y_e, up, njt, prm, ext)
    coords *= prm.scale
    info = {"ni": ni, "nj_top": njt, "nj_bot": njb, "cells": int(quads.shape[0]), "nodes": int(coords.shape[0]),
            "x_out": x_out, "y_bot": y_bot, "i_te": i_te, "L_cowl": L_cowl, "L_ramp": L_ramp,
            "cowl_thickness": t_c, "ramp_fillet": float(prm.ramp_fillet),
            "dup_stations": [int(i) for i in np.where(dup)[0]], **ext}
    return coords, quads, bedges, info, y_mid, y_top


def _add_ext_top(coords, quads, bedges, xs, yt, k, L_ramp, y_e, up, njt, prm, ext):
    """ランプ側外部流ブロック (top バンド + wake ブロック) を既存ノード・セルの後ろに足す (plan §4.11)。
    k = ランプ後縁の station。yt = noz バンド上境界 (ランプ/プルーム上線)。"""
    ni = len(xs); njT, njW = int(prm.nj_ext_top), int(prm.nj_wake)
    # 機体上面はランプ最大 y の上に置く。クリアランスは絶対値 (H) だが、壁第一セル (first_top_frac) より
    # 十分厚くないとテーパ区間で潰れるので下限を課す (run_0071 の発散対策)
    clr = max(float(prm.vehicle_clearance), 3.0 * float(prm.first_top_frac))
    y_veh = float(yt[:k + 1].max()) + clr
    ii = np.arange(ni)
    if prm.vehicle_taper > 0.0:
        # 機体厚 t_v(x) = (y_veh − ランプ y) × 係数。taper は **L_ramp 比**で、短ランプでも全長を食い潰さない。
        # 係数は **smoothstep** (C¹): 線形 clip だとテーパ開始点で上面の勾配が水平から不連続に折れ、
        # そこが M∞10 のスリップ壁上の凸角になって ω が発散した (run_0071 inf_00/inf_01 の m10_on,
        # 折れ角 −11.8°、テーパを短くすると −25.8° と悪化する)。smoothstep は両端で fac' = 0 なので
        # 上面はテーパ開始点でも TE でも勾配が連続になる。
        taper_len = float(prm.vehicle_taper) * L_ramp
        _t = np.clip((L_ramp - xs) / taper_len, 0.0, 1.0)
        fac = _t * _t * (3.0 - 2.0 * _t)
        y3_veh = yt + (y_veh - yt) * fac
        y3_veh[k] = y_e
        h_base = 0.0
    else:
        y3_veh = np.full(ni, y_veh)
        h_base = y_veh - y_e
    has_wake = h_base > 1e-9
    if has_wake and njW < 3:
        raise ValueError("nj_wake は 3 以上")
    y3 = np.where(ii <= k, y3_veh, (y_e + h_base) + (xs - L_ramp) * np.tan(prm.top_ext_angle))   # top の下線 (TE で連続)
    N0 = coords.shape[0]
    top0_own = list(range(ni)) if has_wake else list(range(k))     # j=0 を自前ノードで持つ station
    top0_idx = {i: n for n, i in enumerate(top0_own)}
    N_top0, N_topj = len(top0_own), ni * (njT - 1)

    def top(i, j):
        if j == 0:
            return N0 + top0_idx[i] if i in top0_idx else up(i, njt - 1)
        return N0 + N_top0 + i * (njT - 1) + (j - 1)
    N1 = N0 + N_top0 + N_topj

    def wake(i, j):
        if j == 0:
            return up(i, njt - 1)
        if j == njW - 1:
            return top(i, 0)
        return N1 + (i - k) * (njW - 2) + (j - 1)
    N2 = N1 + ((ni - k) * (njW - 2) if has_wake else 0)
    coords = np.vstack([coords, np.zeros((N2 - N0, 3))])
    s_t = _geom_start(njT, min(prm.first_top_frac / prm.top_depth, 0.5))
    s_w = _geom_start(njW, min(prm.first_wall_frac / h_base, 0.5)) if has_wake else None
    for i in range(ni):
        yy = y3[i] + s_t * prm.top_depth
        if i in top0_idx:
            coords[top(i, 0)] = (xs[i], y3[i], 0.0)
        coords[top(i, 1):top(i, njT - 1) + 1, 0] = xs[i]
        coords[top(i, 1):top(i, njT - 1) + 1, 1] = yy[1:]
        if has_wake and i >= k:
            yw = yt[i] + s_w * h_base
            for j in range(1, njW - 1):
                coords[wake(i, j)] = (xs[i], yw[j], 0.0)
    extra = []
    for i in range(ni - 1):
        for j in range(njT - 1):
            extra.append((top(i, j), top(i + 1, j), top(i + 1, j + 1), top(i, j + 1)))
    if has_wake:
        for i in range(k, ni - 1):
            for j in range(njW - 1):
                extra.append((wake(i, j), wake(i + 1, j), wake(i + 1, j + 1), wake(i, j + 1)))
    quads = np.vstack([quads, np.asarray(extra, dtype=np.int64)])
    for i in range(k):                                    # 機体上面 (x ≤ L_ramp)
        bedges["vehicle"].append((top(i, 0), top(i + 1, 0)))
    if has_wake:                                          # base (後縁の鉛直面)
        for j in range(njW - 1):
            bedges["vehicle"].append((wake(k, j), wake(k, j + 1)))
    for i in range(ni - 1):
        bedges["top_out"].append((top(i, njT - 1), top(i + 1, njT - 1)))
    for j in range(njT - 1):
        bedges["inlet_ext"].append((top(0, j), top(0, j + 1)))
        bedges["outlet"].append((top(ni - 1, j), top(ni - 1, j + 1)))
    if has_wake:
        for j in range(njW - 1):
            bedges["outlet"].append((wake(ni - 1, j), wake(ni - 1, j + 1)))
    ext.update({"nj_ext_top": njT, "nj_wake": njW if has_wake else 0, "y_veh": y_veh, "h_base": h_base,
                "has_wake": has_wake, "i_ramp_te": int(k), "vehicle_taper": float(prm.vehicle_taper)})
    return coords, quads


def write_msh41_named(path, coords, quads, bedges: dict, phys: dict) -> None:
    """任意の境界名セットで gmsh msh4.1 を書く (mesh2d.write_msh41_2d の一般化)。"""
    names = [n for n in phys if n != "fluid" and len(bedges.get(n, [])) > 0]
    n_nodes = coords.shape[0]
    n_elems = quads.shape[0] + sum(len(bedges[g]) for g in names)
    mn, mx = coords.min(0), coords.max(0)
    L = []
    ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap(f"$PhysicalNames\n{len(names) + 1}")
    for nm in names:
        ap(f'1 {phys[nm]} "{nm}"')
    ap(f'2 {phys["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities")
    ap(f"0 {len(names)} 1 0")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for ci, nm in enumerate(names, start=1):
        ap(f"{ci} {bb} 1 {phys[nm]} 0")
    ap(f"1 {bb} 1 {phys['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes")
    ap(f"1 {n_nodes} 1 {n_nodes}")
    ap(f"2 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    ap("\n".join(f"{c[0]:.10g} {c[1]:.10g} {c[2]:.10g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements")
    ap(f"{1 + len(names)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"2 1 3 {quads.shape[0]}")
    ap("\n".join(f"{et + k} {q[0]+1} {q[1]+1} {q[2]+1} {q[3]+1}" for k, q in enumerate(quads)))
    et += quads.shape[0]
    for ci, nm in enumerate(names, start=1):
        edges = bedges[nm]
        ap(f"1 {ci} 1 {len(edges)}")
        ap("\n".join(f"{et + k} {a+1} {b+1}" for k, (a, b) in enumerate(edges)))
        et += len(edges)
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
