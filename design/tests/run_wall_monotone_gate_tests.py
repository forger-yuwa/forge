#!/usr/bin/env python3
"""設計壁の単調性の複合ゲート (`wall_axismach.monotone_gate`) と `AxisMachCFDWall.validate` の置き換えの試験。
plan: plans/active/tooling-nozzle-deltastar-band-root-zone.md §5.1 #4b (諮問 2026-10-13 design-wall-monotone-gate)。

ゲート: 区間 [a, b] の区間多項式で (i) 最小の r′ ≥ −1.74533e-6 (= −tan(1e-4°))、(ii) 最大の累積減少
D = max_x{max_{u≤x} r(u) − r(x)} ≤ 1e-7 r_t。非有限・評価不能は不合格、BSpline / PPoly でなければ対象外 (従来の標本の検査)。

1. 定数が plan の値。
2. 解析解のある区間多項式で min r′・D が厳密 (相対 1e-9)。
3. 負例 (事前固定の形): 格子点の間の収縮 (min θ = −2e-4°・D < 1e-7) は傾きだけで、緩い収縮 (min θ ≥ −1e-4°・D = 2e-7) は累積だけで不合格。
   区間の途中と終わりの直前 (その後ろは平行) の両方。従来の標本の検査 (4000 点・隣接差分 < −1e-9) が見逃す収縮を拒否する。
4. 対照: 厳密単調・平行 (r′ ≡ 0)・r″ ≡ 0 の区間・重ノットは合格。
5. 境界の値: 閾値ちょうどは合格、外は不合格 (既定の閾値と 2 進で厳密な閾値の両方)。
6. 評価不能: 非有限の係数・定義域の外・多次元・不正な区間は不合格。BSpline / PPoly 以外は対象外。
7. 密な標本との照合 (実際の 5 次スプラインに局所の多項式を足したもの): ゲートの min r′ ≤ 標本の min r′、D ≥ 標本の D、差は 1e-12 以内。
8. 設計壁の型 (interp・lsq・joint) の `_spl` はすべて BSpline でゲートの対象。validate は負例で「設計壁区間で半径が非単調」を出し、
   対象外の型 (BSpline / PPoly でない呼び出し可能なもの) では従来の標本の検査を使う。

usage: python3 design/tests/run_wall_monotone_gate_tests.py   (実行時間 数秒、run データ不要)
"""
import copy
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline, CubicSpline, PPoly, make_interp_spline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.wall_axismach import (MONO_D_MAX, MONO_R1_MIN, AxisMachCFDWall,  # noqa: E402
                                                 JointFitCFDWall, LSQBsplineCFDWall, monotone_gate)

FAIL = 0
DEG = np.pi / 180.0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def pp_from_pieces(breaks, polys):
    """区間ごとの多項式 (局所座標 s = x − 区間の始点、高次から) → PPoly。"""
    deg = max(len(p) for p in polys) - 1
    c = np.zeros((deg + 1, len(polys)))
    for j, p in enumerate(polys):
        c[deg + 1 - len(p):, j] = p
    return PPoly(c, np.asarray(breaks, dtype=float))


def dip_pieces(a_amp, L, rise=0.01, pre=1.0, tail=None, tail_slope=None):
    """[0, pre] で傾き rise の上り → [pre, pre+L] で r′ = −a_amp (1 − u²)² (u ∈ [−1, 1]) → 後ろ (tail) は傾き tail_slope。
    min r′ = −a_amp、累積減少 D = a_amp·L·8/15 (a_amp > 0 のとき)。r は連続 (r′ は継ぎ目で跳んでよい)。"""
    tail = pre if tail is None else tail
    tail_slope = rise if tail_slope is None else tail_slope
    p1 = np.array([rise, 1.0])
    r_pre = float(np.polyval(p1, pre))
    Q = np.poly1d([1.0, 0.0, -2.0, 0.0, 1.0])(np.poly1d([2.0 / L, -1.0]))     # (1 − u²)², u = 2s/L − 1
    p2 = np.polyint((-a_amp * Q).coeffs, k=r_pre)
    r_mid = float(np.polyval(p2, L))
    p3 = np.array([tail_slope, r_mid])
    return pp_from_pieces([0.0, pre, pre + L, pre + L + tail], [p1, p2, p3])


def sampled_old(pp, a, b, n=4000, tol=1e-9):
    """従来の検査の式 ([a, b] の n 点の等間隔標本の隣接差分 < −tol) を区間多項式に当てたもの。"""
    xs = np.linspace(a, b, n)
    return not bool(np.any(np.diff(pp(xs)) < -tol))


# --- 1. 定数 -------------------------------------------------------------------------------------------------------
check(f"MONO_R1_MIN = −1.74533e-6 ({MONO_R1_MIN!r})", MONO_R1_MIN == -1.74533e-6)
check(f"MONO_R1_MIN ≈ −tan(1e-4°) (差 {abs(MONO_R1_MIN + np.tan(1e-4 * DEG)):.1e} ≤ 1e-11)",
      abs(MONO_R1_MIN + np.tan(1e-4 * DEG)) <= 1e-11)
check(f"MONO_D_MAX = 1e-7 ({MONO_D_MAX!r})", MONO_D_MAX == 1e-7)

# --- 2. 解析解: r = (x − c)³/3 − ε (x − c) on [0, 2c] (r″ = 0 の根 x = c で min r′ = −ε、D = (4/3) ε^{3/2}) -----------
eps, cc = 1e-6, 1.0
cub = np.poly1d([1.0 / 3.0, 0.0, -eps, 0.0])(np.poly1d([1.0, -cc])).coeffs
g = monotone_gate(pp_from_pieces([0.0, 2 * cc], [cub]), 0.0, 2 * cc)
check(f"解析解: min r′ = −ε (相対差 {abs(g['min_r1'] / -eps - 1):.1e} ≤ 1e-9)、x = c ({g['x_min_theta']:.12f})",
      abs(g["min_r1"] / -eps - 1) <= 1e-9 and abs(g["x_min_theta"] - cc) <= 1e-9)
D_ex = 4.0 / 3.0 * eps ** 1.5
check(f"解析解: D = (4/3) ε^1.5 = {D_ex:.6e} (ゲート {g['D']:.6e}、相対差 {abs(g['D'] / D_ex - 1):.1e} ≤ 1e-6)、"
      f"x_D_peak = c − √ε・x_D = c + √ε",
      abs(g["D"] / D_ex - 1) <= 1e-6 and abs(g["x_D_peak"] - (cc - eps ** 0.5)) <= 1e-9
      and abs(g["x_D"] - (cc + eps ** 0.5)) <= 1e-9)
g = monotone_gate(dip_pieces(1.5e-6, 0.25), 0.0, 2.25)
check(f"解析解 (区間 3 つ): min r′ = −1.5e-6 ({g['min_r1']:.9e})・D = 2e-7 ({g['D']:.9e})",
      abs(g["min_r1"] / -1.5e-6 - 1) <= 1e-9 and abs(g["D"] / 2e-7 - 1) <= 1e-9)

# --- 3. 負例 (事前固定の形) ------------------------------------------------------------------------------------------
a_slope = float(np.tan(2e-4 * DEG))
for where, kw in (("途中", {}), ("終わりの直前 (後ろは平行)", {"tail": 0.05, "tail_slope": 0.0})):
    pp = dip_pieces(a_slope, 0.01, **kw)
    b = float(pp.x[-1])
    g = monotone_gate(pp, 0.0, b)
    check(f"負例 傾き・{where}: 不合格は傾きだけ (min θ {g['min_theta_deg']:.4e}° < −1e-4°、D {g['D']:.3e} ≤ 1e-7) {g['reasons']}",
          (not g["pass"]) and g["min_r1"] < MONO_R1_MIN and g["D"] <= MONO_D_MAX
          and abs(g["min_theta_deg"] / -2e-4 - 1) <= 1e-6 and len(g["reasons"]) == 1 and "壁角" in g["reasons"][0])
    pp = dip_pieces(1.5e-6, 0.25, **kw)
    b = float(pp.x[-1])
    g = monotone_gate(pp, 0.0, b)
    check(f"負例 累積・{where}: 不合格は累積だけ (min θ {g['min_theta_deg']:.4e}° ≥ −1e-4°、D {g['D']:.4e} = 2e-7) {g['reasons']}",
          (not g["pass"]) and g["min_r1"] >= MONO_R1_MIN and abs(g["D"] / 2e-7 - 1) <= 1e-9
          and len(g["reasons"]) == 1 and "累積" in g["reasons"][0])
# 標本の間の収縮: r = (x − 20)³/3 − ε (x − 20) on [0, 39.99] (4000 点の間隔 0.01 r_t > 負の区間の幅 2√ε = 3.7e-3)
cub = np.poly1d([1.0 / 3.0, 0.0, -a_slope, 0.0])(np.poly1d([1.0, -20.0])).coeffs
pp = pp_from_pieces([0.0, 39.99], [cub])
g = monotone_gate(pp, 0.0, 39.99)
check(f"標本の間の収縮: 従来の標本の検査は見逃す ({sampled_old(pp, 0.0, 39.99)}) がゲートは傾きで拒否 "
      f"(min θ {g['min_theta_deg']:.4e}°、D {g['D']:.2e})",
      sampled_old(pp, 0.0, 39.99) and not g["pass"] and g["min_r1"] < MONO_R1_MIN)

# --- 4. 対照 -------------------------------------------------------------------------------------------------------
for name, pp in (("厳密単調 (平坦部 q > 0)", dip_pieces(-1.5e-6, 0.25)),
                 ("厳密単調 (平坦部 q > 0、後ろは平行)", dip_pieces(-1.5e-6, 0.25, tail=0.05, tail_slope=0.0)),
                 ("平行 (r′ ≡ 0、r ≡ 1)", pp_from_pieces([0.0, 1.0, 2.0], [np.array([1.0]), np.array([0.0, 0.0, 1.0])])),
                 ("平行の区間 (r″ ≡ 0) を挟む上り", dip_pieces(0.0, 0.25))):
    g = monotone_gate(pp, float(pp.x[0]), float(pp.x[-1]))
    check(f"対照 {name}: 合格 (min r′ {g['min_r1']:.2e}・D {g['D']:.2e})", g["pass"] and g["D"] <= 1e-15)
# 重ノット (2 次の BSpline の内部の 2 重ノットで r′ が不連続、r は連続): 継ぎ目の左右の極限を両方評価する
# 左の極限 r′ = 2、右の極限 r′ = 4 (c3 − c2) = −3e-6 (その後は上り)
t = np.r_[[0.0] * 3, 0.5, 0.5, [1.0] * 3]
spl = BSpline(t, np.array([0.0, 0.5, 1.0, 1.0 - 7.5e-7, 2.0]), 2)
g = monotone_gate(spl, 0.0, 1.0)
check(f"重ノットの BSpline (k=2): 継ぎ目の右の極限 r′ = −3e-6 で不合格 (min r′ {g['min_r1']:.6e} @ {g['x_min_theta']})",
      (not g["pass"]) and abs(g["min_r1"] / -3e-6 - 1) <= 1e-9 and g["x_min_theta"] == 0.5 and g["D"] < 1e-12)

# --- 5. 境界の値 ---------------------------------------------------------------------------------------------------
lin = lambda s, L: pp_from_pieces([0.0, L], [np.array([s, 0.0])])  # noqa: E731
g = monotone_gate(lin(MONO_R1_MIN, 0.05), 0.0, 0.05)
check(f"境界: r′ = MONO_R1_MIN ちょうど (D {g['D']:.3e}) は合格", g["pass"])
g = monotone_gate(lin(MONO_R1_MIN * (1 + 1e-6), 0.05), 0.0, 0.05)
check(f"境界: r′ = MONO_R1_MIN·(1 + 1e-6) は傾きで不合格 {g['reasons']}", (not g["pass"]) and len(g["reasons"]) == 1)
g = monotone_gate(lin(-1e-6, 0.1 * (1 - 1e-6)), 0.0, 0.1 * (1 - 1e-6))
check(f"境界: D = 1e-7·(1 − 1e-6) ({g['D']:.9e}) は合格", g["pass"])
g = monotone_gate(lin(-1e-6, 0.1 * (1 + 1e-6)), 0.0, 0.1 * (1 + 1e-6))
check(f"境界: D = 1e-7·(1 + 1e-6) ({g['D']:.9e}) は累積で不合格 {g['reasons']}",
      (not g["pass"]) and len(g["reasons"]) == 1 and "累積" in g["reasons"][0])
s2, L2 = -2.0 ** -20, 2.0 ** -3                                    # 2 進で厳密: min r′ = −2^-20、D = 2^-23
g = monotone_gate(lin(s2, L2), 0.0, L2, r1_min=s2, d_max=2.0 ** -23)
check(f"境界 (厳密): min r′ = r1_min・D = d_max ちょうどは合格 (D {g['D']!r})", g["pass"] and g["D"] == 2.0 ** -23)
g = monotone_gate(lin(s2, L2), 0.0, L2, r1_min=s2 * (1 - 2.0 ** -10), d_max=2.0 ** -23)
check("境界 (厳密): r1_min を 1 段きつくすると傾きで不合格", (not g["pass"]) and g["min_r1"] < g["r1_min"])
g = monotone_gate(lin(s2, L2), 0.0, L2, r1_min=s2, d_max=2.0 ** -23 * (1 - 2.0 ** -10))
check("境界 (厳密): d_max を 1 段きつくすると累積で不合格", (not g["pass"]) and g["D"] > g["d_max"])

# --- 6. 評価不能・対象外 --------------------------------------------------------------------------------------------
base = dip_pieces(-1.5e-6, 0.25)
c_nan = base.c.copy(); c_nan[0, 1] = np.nan
g = monotone_gate(PPoly(c_nan, base.x), 0.0, 2.25)
check(f"評価不能: NaN の係数は不合格 {g['reasons']}", (not g["pass"]) and g["applicable"] and "非有限" in g["reasons"][0])
spl_inf = BSpline(np.r_[[0.0] * 6, [1.0] * 6], np.r_[0.0, 0.1, 0.2, np.inf, 0.4, 0.5], 5)
g = monotone_gate(spl_inf, 0.0, 1.0)
check(f"評価不能: inf の係数の BSpline は不合格 {g['reasons']}", (not g["pass"]) and g["applicable"])
g = monotone_gate(base, -0.1, 2.25)
check(f"評価不能: 区間が定義域の外は不合格 {g['reasons']}", (not g["pass"]) and "定義域" in g["reasons"][0])
g = monotone_gate(base, 1.0, 1.0)
check(f"評価不能: 長さ 0 の区間は不合格 {g['reasons']}", (not g["pass"]) and "不正" in g["reasons"][0])
g = monotone_gate(base, 0.0, np.nan)
check(f"評価不能: NaN の区間端は不合格 {g['reasons']}", (not g["pass"]) and "不正" in g["reasons"][0])
g = monotone_gate(PPoly(np.zeros((2, 1, 2)), [0.0, 1.0]), 0.0, 1.0)
check(f"評価不能: 多次元の値は不合格 {g['reasons']}", (not g["pass"]) and "多次元" in g["reasons"][0])
g = monotone_gate(lambda x, nu=0: x, 0.0, 1.0)
check(f"対象外: BSpline / PPoly でない呼び出し可能なもの ({g.get('reason')})", g["applicable"] is False and g["pass"] is None)
g = monotone_gate(CubicSpline([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]), 0.0, 3.0)
check(f"CubicSpline (PPoly の派生) は対象で合格 ({g['kind']})", g["applicable"] and g["pass"])

# --- 7. 合成の設計壁の表 (MOC 風: r = 1 + A (1 − exp(−x²/(2 A R))), r′(0) = 0・r″(0) = 1/R、単調増加) ---------------
R, A_ = 2.0, 2.0
xt = np.r_[np.linspace(0.0, 0.5, 41)[:-1], np.linspace(0.5, 8.0, 160)]
rt_ = 1.0 + A_ * (1.0 - np.exp(-xt ** 2 / (2.0 * A_ * R)))
tht = np.arctan(xt / R * np.exp(-xt ** 2 / (2.0 * A_ * R)))
tbl = np.c_[xt, rt_, tht, np.full_like(xt, 3.0)]
walls = {"interp": AxisMachCFDWall(tbl, R=R), "lsq": LSQBsplineCFDWall(tbl, R=R, n_cp=24, tol_dr=1e-2, tol_dtheta_deg=1.5),
         "joint": JointFitCFDWall(tbl, R=R)}
for name, w in walls.items():
    g = monotone_gate(w._spl, w.x0, w.x_e)
    check(f"設計壁の型 {name} ({type(w).__name__}): _spl は {type(w._spl).__name__} でゲートの対象・合格 "
          f"(min θ {g['min_theta_deg']:.2e}°、D {g['D']:.1e})",
          isinstance(w._spl, BSpline) and g["applicable"] and g["pass"])
    check(f"設計壁の型 {name}: validate に非単調なし", not any("非単調" in m for m in w.validate()))

# 実際の 5 次スプラインに局所の C² の多項式 −α s³ (h − s)³ を 1 区間だけ足す (端で値・1 階・2 階が 0)
w0 = walls["interp"]
Spp = PPoly.from_spline(w0._spl)
j = int(np.searchsorted(Spp.x, 7.0, side="right") - 1)
h = float(Spp.x[j + 1] - Spp.x[j])
bump = np.poly1d([1.0, 0.0, 0.0, 0.0])(np.poly1d([1.0, 0.0])) * np.poly1d([-1.0, h]) ** 3       # s³ (h − s)³
r1_here = float(np.max(w0._spl(np.linspace(Spp.x[j], Spp.x[j + 1], 101), 1)))
sg = np.linspace(0.0, 0.5, 100001)
f_max = float(np.max(3.0 * sg ** 2 * (1.0 - sg) ** 2 * (1.0 - 2.0 * sg)))   # r′ の減り = α h⁵ · 3σ²(1−σ)²(1−2σ) の最大


def with_bump(alpha):
    c = np.zeros((7, Spp.c.shape[1]))
    c[1:, :] = Spp.c
    c[:, j] += np.r_[np.zeros(7 - len(bump.coeffs)), -alpha * bump.coeffs]
    return PPoly(c, Spp.x)


# 区間の前半で r′ を負にする α (区間の中の S′ の最大の 1.5 倍を減らす)
alpha = 1.5 * r1_here / (f_max * h ** 5)
ppb = with_bump(alpha)
g = monotone_gate(ppb, w0.x0, w0.x_e)
xs = np.unique(np.r_[np.linspace(w0.x0, w0.x_e, 400001), np.linspace(Spp.x[j], Spp.x[j + 1], 400001)])
r1s, rs = ppb(xs, 1), ppb(xs)
Ds = float((np.maximum.accumulate(rs) - rs).max())
check(f"密な標本と照合: ゲートの min r′ {g['min_r1']:.12e} ≤ 標本 {r1s.min():.12e} (差 ≤ 1e-12)",
      g["min_r1"] <= r1s.min() + 1e-15 and r1s.min() - g["min_r1"] <= 1e-12)
check(f"密な標本と照合: ゲートの D {g['D']:.12e} ≥ 標本 {Ds:.12e} (差 ≤ 1e-12)、不合格",
      g["D"] >= Ds - 1e-15 and g["D"] - Ds <= 1e-12 and not g["pass"])

# --- 8. validate への組み込み --------------------------------------------------------------------------------------
wn = copy.copy(w0)
wn._spl = ppb
msgs = wn.validate()
check(f"validate: 負例の _spl (PPoly) で「設計壁区間で半径が非単調」とゲートの理由 ({[m[:70] for m in msgs if '非単調' in m]})",
      any(m.startswith("設計壁区間で半径が非単調 (") for m in msgs))
wf = copy.copy(w0)
wf._spl = lambda x, nu=0: w0._spl(x, nu)                          # 対象外の型 (呼び出し可能だが BSpline / PPoly でない)
check("validate: 対象外の型で単調な壁は非単調なし (従来の標本の検査)", not any("非単調" in m for m in wf.validate()))
dip = with_bump(alpha * 40.0)                                       # 標本の検査でも見える深い収縮
wf2 = copy.copy(w0)
wf2._spl = lambda x, nu=0: dip(x, nu)
m2 = [m for m in wf2.validate() if "非単調" in m]
check(f"validate: 対象外の型は従来の標本の検査で拒否 (メッセージ {m2})", m2 == ["設計壁区間で半径が非単調"])

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
