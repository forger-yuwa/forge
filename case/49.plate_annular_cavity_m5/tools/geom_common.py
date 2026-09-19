#!/usr/bin/env python3
"""case/49 の共通幾何関数 (manifest から解決した値だけを使う)。

plan §4.8 / §4.11 (codex 2026-09-19 plan-2 M5):
内側円柱は x 方向に x_off 偏心しうるので、キャビティ流体断面は
    x^2 + y^2 < Ro^2   かつ   (x - x_off)^2 + y^2 > Ri^2
であって `Ri < r < Ro` **ではない**。同心マスクを偏心形状に当てると
(Ro=25, Ri=22.5, x_off=1 で) 真の半環断面の 12 % を取りこぼし、同量の固体側を誤って含む。
CV マスク・すきま中央線・周方向平均の重みは必ずここを経由する。

方位 θ は**外筒軸まわり・上流基準**で測る (2026-09-19 ユーザ指定):
  **θ=0° が −x (上流よどみ側)**、θ=180° が +x (下流側)、θ=90° が +y (横)。
  すなわち方位 θ の向きは d(θ) = (−cosθ, sinθ)。半割 (y>=0) では θ∈[0, 180°]。
"""
import json
from pathlib import Path

import numpy as np


def load_manifest(path=None, run=None):
    """manifest を読む。優先順は  明示 path > run/manifest.json > CASE49_MANIFEST > manifest.json。

    **run ディレクトリに manifest.json があればそれを最優先**にする。偏心スイープのように
    幾何が run ごとに違う場合、共有の `manifest.json` を読むと**別の偏心の CV マスク**で
    評価してしまう (2026-09-19 に実際に起きた: すきま中央 ΔT 11.75 K が 10.53 K、
    侵入深さ 46.0 mm が 4.5 mm になった)。run 自身が持つ幾何を使えばこれが起きない。
    """
    import os
    if path:
        p = Path(path)
    elif run is not None and (Path(run) / "manifest.json").exists():
        p = Path(run) / "manifest.json"
    else:
        p = (Path(__file__).resolve().parent.parent
             / os.environ.get("CASE49_MANIFEST", "manifest.json"))
    if not p.exists():
        raise SystemExit("manifest が無い: %s  (`python3 setup.py --resolve` を先に実行)" % p)
    return json.loads(p.read_text())


def geom(man):
    return man["geometry"]


# ---------------------------------------------------------------- マスク
def in_gap(x, y, man, tol=0.0):
    """点 (x,y) が環状すきま内 (内円の外・外円の内) にあるか。中央面などの自己検査用。"""
    G = man["geometry"]
    ro = np.hypot(x, y)
    ri = np.hypot(x - G["x_off"], y)
    return (ro < G["Ro"] - tol) & (ri > G["Ri"] + tol)


def cavity_mask(x, y, z, man, shrink=0.0):
    """キャビティ流体 (環状スリット) の中か。shrink>0 で両壁から shrink [m] 内側に絞る
    (壁ピン値を避けて「ガスの」温度を採るときに使う)。"""
    g = geom(man)
    Ro, Ri, off, d = g["Ro"], g["Ri"], g["x_off"], g["depth"]
    r_o = np.hypot(x, y)
    r_i = np.hypot(x - off, y)
    return (r_o < Ro - shrink) & (r_i > Ri + shrink) & (z < -shrink) & (z > -d + shrink)


def in_domain_plate(x, y, z, man, tol=1e-9):
    """平板面 (z=0) のうちキャビティ開口の外か (= 固体壁に接する面)。"""
    g = geom(man)
    return (np.abs(z) < tol) & ~((np.hypot(x, y) < g["Ro"]) & (np.hypot(x - g["x_off"], y) > g["Ri"]))


# ---------------------------------------------------------------- すきま形状
def ray_dir(theta):
    """方位 θ の単位ベクトル (x, y)。**θ=0 が上流 (−x)**。"""
    return -np.cos(theta), np.sin(theta)


def inner_radius_at(theta, man):
    """外筒軸から見た方位 θ [rad] のレイが内円柱と交わる半径 r_i(θ) [m]。
    偏心 2 円 (内円中心 (x_off,0)) とレイ d=(−cosθ, sinθ) の交点:
      r_i = −x_off cosθ + sqrt(Ri^2 − x_off^2 sin^2 θ)
    (θ=0 が上流基準なので、下流ずらし x_off>0 では θ=0 側のすきまが**広く**なる)。"""
    g = geom(man)
    off, Ri = g["x_off"], g["Ri"]
    c, s = np.cos(theta), np.sin(theta)
    disc = Ri ** 2 - (off * s) ** 2
    if np.any(disc < 0):
        raise SystemExit("|x_off| が Ri を超えている (レイが内円柱と交わらない)")
    return -off * c + np.sqrt(disc)


def gap_at(theta, man):
    """方位 θ の**半径方向**すきま幅 [m] (レイに沿った定義。厳密な最短距離ではない)。"""
    return geom(man)["Ro"] - inner_radius_at(theta, man)


def gap_center_radius(theta, man):
    """方位 θ のすきま中央半径 [m] (外筒軸まわり)。"""
    return 0.5 * (geom(man)["Ro"] + inner_radius_at(theta, man))


def gap_center_points(theta, z, man):
    """すきま中央・深さ z の (x, y, z) 点列。theta, z はブロードキャスト可能。"""
    rc = gap_center_radius(theta, man)
    dx, dy = ray_dir(theta)
    return rc * dx, rc * dy, np.broadcast_to(z, np.shape(rc) * 1 or (1,))


def azimuth_weights(theta, man):
    """周方向平均の重み: すきま断面積 ∝ (Ro^2 - r_i^2)/2 の dθ 重み。
    偏心時は下流側が狭く上流側が広いので単純平均にしない。"""
    g = geom(man)
    ri = inner_radius_at(theta, man)
    w = 0.5 * (g["Ro"] ** 2 - ri ** 2)
    return w / np.sum(w)


def azimuth_of(x, y, man):
    """点の方位 θ [rad] (外筒軸まわり・**上流基準**, 0..π が半割の範囲)。"""
    return np.arctan2(y, -x)


# ---------------------------------------------------------------- 測点
def probe_points(man, depth_frac, n_theta=None, theta=None):
    """すきま中央・指定深さ (depth_frac: 0=開口, 1=床) の測点。周方向平均用。"""
    g = geom(man)
    if theta is None:
        n = n_theta or man["eval"]["n_theta"]
        theta = np.linspace(0.0, np.pi, n)
    z = -depth_frac * g["depth"]
    rc = gap_center_radius(theta, man)
    dx, dy = ray_dir(theta)
    return np.stack([rc * dx, rc * dy, np.full_like(rc, z)], axis=1), theta


def opening_area(man, half=True):
    """開口 (= 床) の環状面積 [m^2]。偏心に依らず π(Ro^2-Ri^2)(/2)。"""
    g = geom(man)
    a = np.pi * (g["Ro"] ** 2 - g["Ri"] ** 2)
    return 0.5 * a if half else a


def wall_area(man, name, half=True):
    """壁グループの解析面積 [m^2] (manifest の group_area と同じ値)。"""
    a = man["geometry"]["group_area_m2"][name]
    return a if half else 2.0 * a


def self_test(man):
    """同心/偏心で面積とマスクの整合を確認する (plan §4.8 の 12 % 取りこぼし検出)。"""
    g = geom(man)
    rng = np.random.default_rng(0)
    N = 2_000_000
    Ro, Ri, off = g["Ro"], g["Ri"], g["x_off"]
    x = rng.uniform(-Ro, Ro, N); y = rng.uniform(0.0, Ro, N)
    z = np.full(N, -0.5 * g["depth"])
    m_true = cavity_mask(x, y, z, man)
    a_mc = m_true.mean() * (2 * Ro) * Ro
    a_ex = 0.5 * np.pi * (Ro ** 2 - Ri ** 2)
    # 同心マスク (誤り) との差
    r = np.hypot(x, y)
    m_conc = (r > Ri) & (r < Ro) & (y >= 0)
    miss = (m_true & ~m_conc).mean() * (2 * Ro) * Ro
    print("cavity 断面 MC %.6e / 解析 %.6e (差 %.3f %%)" % (a_mc, a_ex, 100 * abs(a_mc / a_ex - 1)))
    print("同心マスクの取りこぼし %.6e m^2 = %.2f %%" % (miss, 100 * miss / a_ex))
    th = np.linspace(0, np.pi, 5)
    print("gap(θ=0,45,90,135,180deg) [mm]  (θ=0 が上流):", np.round(gap_at(th, man) * 1e3, 4))


if __name__ == "__main__":
    self_test(load_manifest())
