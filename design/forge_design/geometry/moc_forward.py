r"""壁を与えた前進 MOC (直接問題) と、その受理ゲート (前進 MOC ゲート)。

計画: plans/active/tooling-nozzle-bell-two-mode.md §4.5 (第一の判定量)・§5.1 #3b (共用モジュール、所有はこの系列)・#3c (ガード)。
③ベルのゲート、④デュアルベルの基部受理、剥離 plan の TC1A の全流れ壁圧が同じこのモジュールを使う。

**問い**: 与えた壁 $r_w(x)$ ($x\in[x_0,L]$) の下で、初期線 (スロートの C⁻ 特性線) から出口 $x=L$ まで、
衝撃波 (同族特性線の交差)・網の反転・単位過程の失敗・亜音速なしに等エントロピーの MOC 場を構成できるか。

単位過程:

- **内点** は①風洞の生産で検証済みの配列版 `moc_kernel.interior_vec` (`axis_limit="analytic"`・`corrector="converge"`)。
  対ごとの分類 `stats["status"]` (`moc_kernel.py` の `PAIR_*`) をそのまま不受理の理由にする。スカラー版 `KernelMOC` は
  使わない (亜音速入力を音速にクリップし、失敗点を `continue` で捨てるため)。
- **壁点** (本モジュール): 下の点 $A$ からの C⁺ と壁の交点。$\theta_W=\arctan r_w'(x_W)$、C⁺ の適合式
  $\theta_W-\nu_W=(\theta_A-\nu_A)+S_+$。係数平均・源項の端点評価・収束判定は `interior_vec` と同じ形
  (更新量 ≤ `tol` まで、上限 `max_corr` 回)。
- **軸点** (本モジュール): 上の点 $B$ ($G_k[1]$) からの C⁻ と軸 ($r=0,\ \theta=0$)。軸点の源項 $\sin\theta/r$ の極限は
  $B$ の値で代用する (`axis_thr_neighbor`)。前進問題には風洞の逆 MOC の軸則 (`target_dM`) が無いので、解析極限
  $\theta_r$ が厳密に分かるのは**初期線の軸端だけ** (Hall の軸上アンカー)。`analytic` ではその値と、計算した軸点に持たせた
  代用値を `interior_vec` の `thrA` に渡す (軸端点は対の相手から借りない)。`legacy` は `interior_vec` の従来の代用
  (対の相手から借りる)。高次の推定が不安定だった記録は `axis_thr_neighbor` の docstring。

**網の組み方** (逆 MOC `moc_inverse.InverseMOC.fill_levels` を手本に、依存のない点を段ごとに一括処理):
前線 $k$ = 壁点 $G_k[n]$ から軸点 $G_k[0]$ までの C⁻ ($k=0$ が初期線、$n+1$ 点で一定)。
$G_k[n]$ = $G_{k-1}[n-1]$ からの C⁺ と壁、$G_k[j]$ ($1\le j\le n-1$) = $A=G_{k-1}[j-1]$ (C⁺ 担体) と
$B=G_k[j+1]$ (C⁻ 担体)、$G_k[0]$ = $G_k[1]$ からの C⁻ と軸。C⁺ 線は $k-j$ 一定の対角線。
依存は $(k-1,j-1)$ と $(k,j+1)$ だけなので、段 $\ell=2k-j$ が同じ点は互いに独立 (1 段に内点 ~$n/2$ 個・壁点と軸点が
高々 1 個ずつ)。$x>L$ の点は出口の外として計算に使わない (依存は下流へしか伝わらないので $x\le L$ の場に影響しない)。
C⁺ が出口より下流で壁に当たったら前線の追加をやめる (= 出口まで届いた)。

**ゲート** (`forward_moc(...)["gate"]`): 不受理の理由を符号 (`code`) と文 (`reasons`) で返す。
1. 初期線: 非有限・**亜音速 (M ≤ 1、`pm_nu` のクリップ前に検査)**・壁点が壁に乗っていない
2. $x\le L$ の単位過程の失敗: `interior_vec` の分類 (幾何的棄却・反復の失敗・軸より下)、壁点・軸点の失敗、
   **亜音速 ($\nu<0$ — `mach_of_nu` が音速にクリップする前の値で判定)**・表の外 ($\nu>\nu_{max}$)・壁の外の点・
   C⁺ が壁に届かない
3. 収束した点の最終残差 > `RESID_TOL`
4. **網の反転**: 新点が担体より特性線の向きで手前 (C⁺: $(x_P-x_A)/\cos\alpha_+\le0$、C⁻ も同様)、壁点・軸点の $x$ が単調でない
5. **同族特性線の交差**: `moc_diagnostics` の非隣接折れ線交差の検査 (`_polylines_cross_aabb`) を C⁺ 族 (対角線)・
   C⁻ 族 (前線) の両方に、間引いた線の組で
6. 出口まで届かない (失敗で前線が途切れた・前線数の上限)

**判定は 3 値** (`gate["verdict"]`、plan §5.1 #3c・diagnostician 2026-10-11): `"pass"` (通過) / `"reject"` (壁の不受理) /
`"undecidable"` (判定不能 — 壁の性質として数えない)。`gate["pass"]` は `verdict == "pass"` と同じ。判定不能は次の 3 つ
(`gate["undecidable"]["kind"]`):

- `"initial_line"`: **前線 1 (k = 1) に網の反転がある** (C⁺ 沿い・C⁻ 沿いの反転、壁点・軸点の x の逆行が前線 0→1 の間)。
  初期線と前進の適合式の不整合のアーチファクト (Hall R_t 1・n 641 で全点が不受理になり分類 B と誤判定された件)。
  初期線の C⁻ 沿いの反転 (k = 0) も同じ扱い。他の理由があってもこれを優先する。
- `"near_axis"`: 不受理の根拠となった点 (単位過程の失敗・網の反転・同族特性線の交差・最終残差の最大点・出口に届かない) の
  **軸からの添字 j がすべて j ≤ 2** (`NEAR_AXIS_J`)。前進の軸点の θ_r は 1 次の代用 (`axis_thr_neighbor`) なので、
  軸の極近傍だけに出る不受理は壁の性質と言えない。位置の分からない根拠 (前線数の上限など) が 1 つでもあれば不受理のまま。
- `"initial_line_check"`: 初期線の検査 (亜音速・非有限・壁から外れる・順序) で止まった (マーチしない)。
根拠の点は `gate["evidence"]` (種類・k・j・x・r) と、その j の要約 `gate["evidence_j"]` に記録する。

**内部衝撃波の指標** (連続量、合否とは別に返す): 同族特性線の間隔の比。C⁻ 族は前線 $k-1,k$ の間隔を C⁺ に沿って測り、
壁での間隔で割る ($d_k[j]/d_k[n]$)。C⁺ 族は隣り合う C⁺ の間隔を C⁻ に沿って測り、その 2 本が最初にそろった前線での間隔で割る。
1 から 0 へ近づくほど詰まっている (0 で包絡線 = 衝撃波、負で交差)。最小値と位置を返す。
"""
from __future__ import annotations

import time

import numpy as np

from .moc_inverse import _mass_flux_density
from .moc_kernel import (AXIS_R_EPS, CORR_MAX, CORR_TOL, PAIR_BELOW_AXIS, PAIR_CONVERGED, PAIR_GEOM,
                         PAIR_MAXITER, PAIR_MISSING, PAIR_NONFINITE, RESID_TOL, _check_corr_params, _is_gas,
                         _src_vec, axis_theta_r, interior_vec, mu_vec, pm_mach_vec, pm_nu)

# --- 点の状態 (net["status"]) ---------------------------------------------------------------
# 0..6 は moc_kernel の PAIR_* と同じ意味 (内点は interior_vec の分類をそのまま入れる)。
ST_NONE = -1                 # 計算していない (最後の前線より後)
ST_OK = PAIR_CONVERGED       # 0 収束
ST_MISSING = PAIR_MISSING    # 1 担体が欠けている (上流の失敗の伝播)
ST_GEOM = PAIR_GEOM          # 2 幾何的棄却: 特性線が平行
ST_NONFINITE = PAIR_NONFINITE  # 3 反復中に NaN・Inf
ST_MAXITER = PAIR_MAXITER    # 4 修正子の上限到達
ST_BELOW_AXIS = PAIR_BELOW_AXIS  # 6 新点が軸より下
ST_SUBSONIC = 10             # 収束後の ν < 0 (M < 1。ガスの表・Newton は音速にクリップするのでその前の値で判定)
ST_ABOVE_TABLE = 11          # ν > ν_max (ガスの M↔ν 表の外 — 表は端でクリップする)
ST_WALL_MISS = 12            # 壁点: C⁺ が壁に届かない (交点なし・特性線が後ろ向き)
ST_OUTSIDE_WALL = 13         # 新点 (または壁点の担体) が壁より外 (r > r_w(x))
ST_BEYOND_EXIT = 20          # x > L (出口の外。失敗ではない — 以後の計算に使わない)

FAIL_STATES = (ST_GEOM, ST_NONFINITE, ST_MAXITER, ST_BELOW_AXIS, ST_SUBSONIC, ST_ABOVE_TABLE,
               ST_WALL_MISS, ST_OUTSIDE_WALL)
STATE_NAMES = {ST_NONE: "none", ST_OK: "converged", ST_MISSING: "missing", ST_GEOM: "geom_parallel",
               ST_NONFINITE: "iter_nonfinite", ST_MAXITER: "iter_maxiter", ST_BELOW_AXIS: "below_axis",
               ST_SUBSONIC: "subsonic", ST_ABOVE_TABLE: "above_table", ST_WALL_MISS: "wall_miss",
               ST_OUTSIDE_WALL: "outside_wall", ST_BEYOND_EXIT: "beyond_exit"}

WALL_ON_TOL = 1e-9           # 初期線の壁点が壁に乗っているかの許容差 [r_t]
OUTSIDE_WALL_TOL = 1e-9      # 新点が壁の外とみなす超過量 [r_t]
N_CROSS_SAMPLE = 40          # 非隣接交差の検査で間引く線の本数 (族ごと; moc_diagnostics と同じ既定)
NEAR_AXIS_J = 2              # 不受理の根拠がすべて j ≤ この値なら「判定不能 (軸近傍)」(j = 0 が軸、§5.1 #3c)
HALL_N_MAX_RT_LE1 = 321      # Hall の初期線 (R_t ≤ 1) で許す点数の上限 (§5.1 #3c; 321 可・641 不可の実測)


# ------------------------------------------------------------------------------------------------
# ガスまわり
def nu_max_of(g) -> float:
    """ν の上限 (ガスモデルは M↔ν 表の端、CPG は M→∞ の極限)。"""
    if _is_gas(g):
        tab = getattr(g, "_nu", None)
        if tab is not None:
            return float(tab[-1])
        gam = float(getattr(g, "gamma_ref"))
    else:
        gam = float(g)
    return float(0.5 * np.pi * (np.sqrt((gam + 1.0) / (gam - 1.0)) - 1.0))


def _mach(nu, g):
    return np.asarray(pm_mach_vec(np.asarray(nu, dtype=float), g), dtype=float)


# ------------------------------------------------------------------------------------------------
# 初期線
def hall_initial_line(R: float, g, n: int) -> dict:
    r"""スロートの C⁻ 特性線 (`HallThroat.throat_characteristic`、壁点 $(0,1)$ → 軸) を $n$ 点 ($r$ 等間隔)、軸→壁の順で。

    Hall の級数は定数 γ なので、①風洞の semi-perfect と同じく**スロートの局所 γ\*** (`gas.gamma_throat()`) を渡し、
    $\nu$ は呼び出し側がガスモデルで作る (`evaluate/runner_axismach.py` の `g_hall`)。軸端の $\theta_r$ は Hall の
    軸上アンカー $M'$ (`axis_anchor`) から `axis_theta_r(M, M', g)` (風洞の `axis_theta_r_init` のアンカー経路と同じ)。

    **既知の制約 (2026-10-11 実測、plan 未決)**: Hall の級数の $(\theta, M)$ は C⁻ の適合式
    $d(\theta+\nu)=\sin\mu\,(\sin\theta/r)\,d\sigma$ を厳密には満たさない。線に沿った $\theta+\nu-\int$源項 の壁と軸端の差は
    $R$ 1 で −1.09e-2 rad (③の EXH、CPG γ 1.4 でも −1.07e-2)、$R$ 2 で −2.1e-3 rad。前進マーチは前線 1 をこの線から
    C⁺ 区間 1 つだけ下流に作るので、前線の間隔がこの不整合の効果より細かくなると (③の R 1・n ≥ 641 で実測)
    前線 1 がスロート直後で初期線と交差し、壁によらず不受理になる。n ≤ 321 では起きなかった。

    **ガード** (plan §5.1 #3c、diagnostician 2026-10-11): $R\le1$ では $n>321$ (`HALL_N_MAX_RT_LE1`) を `ValueError` で拒否する。
    実測 (③の EXH、R 1): n 161・321 は (D) と TOP の 13 点で合否が整合初期線 n 321/641/1281 と一致、
    **n 641 は前線 1 がスロート直後 (x ≈ 0.30) で初期線と交差して全点が不受理** (分類 B と誤判定)。$R>1$ は制限しない
    (R 2 の不整合は −2.1e-3 rad と小さい)。不整合の大きさを測る診断には制限のない `_hall_line_raw` を使う。"""
    if float(R) <= 1.0 and int(n) > HALL_N_MAX_RT_LE1:
        raise ValueError(f"hall_initial_line: R = {float(R):g} ≤ 1 で n = {int(n)} > {HALL_N_MAX_RT_LE1} は使えない "
                         f"(Hall の初期線と C⁻ の適合式の不整合 ~1e-2 rad で前線 1 が初期線と交差する — 321 可・641 不可の実測、"
                         f"plan §5.1 #3c)")
    return _hall_line_raw(R, g, n)


def _hall_line_raw(R: float, g, n: int) -> dict:
    """`hall_initial_line` の本体 (点数の制限なし)。初期線の不整合を測る診断と試験だけが直接使う。"""
    from .transonic import HallThroat
    g_hall = float(g.gamma_throat()) if _is_gas(g) else float(g)
    ht = HallThroat(R=float(R), gamma=g_hall)
    x, r, M, th = ht.throat_characteristic(n=int(n))
    M_A, Mp_A, _ = (float(v) for v in ht.axis_anchor(float(x[0])))
    thr0 = float(axis_theta_r(M_A, Mp_A, g))
    return {"x": np.asarray(x, dtype=float), "r": np.asarray(r, dtype=float), "M": np.asarray(M, dtype=float),
            "theta": np.asarray(th, dtype=float), "thr_axis": thr0, "gamma_hall": g_hall, "R": float(R),
            "source": "hall_throat_characteristic", "M_axis_anchor": M_A, "Mp_axis_anchor": Mp_A}


# ------------------------------------------------------------------------------------------------
# 壁 (与えられた r_w(x) を出口から先は接線で延長 — 出口の外の交点の判定用)
class CallableWall:
    """`r(x)`・`drdx(x)` の 2 関数 (例: `geometry.wall.NozzleWall` の TOP) を `r(x, deriv)` の形にする。"""

    def __init__(self, r, drdx):
        self._r, self._rp = r, drdx

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if deriv == 0:
            return np.asarray(self._r(x), dtype=float)
        if deriv == 1:
            return np.asarray(self._rp(x), dtype=float)
        raise ValueError("CallableWall: deriv は 0 か 1")


class _Wall:
    def __init__(self, wall, L: float):
        self.w = wall
        self.L = float(L)
        self.rL = float(self._raw(self.L, 0))
        self.rpL = float(self._raw(self.L, 1))

    def _raw(self, x, d):
        return float(np.asarray(self.w.r(np.array([float(x)]), d), dtype=float).ravel()[0])

    def r(self, x):
        return self._raw(x, 0) if x <= self.L else self.rL + self.rpL * (x - self.L)

    def rp(self, x):
        return self._raw(x, 1) if x <= self.L else self.rpL

    def intersect(self, xA: float, rA: float, m: float):
        r"""直線 $r=r_A+m(x-x_A)$ ($x>x_A$) と壁の最初の交点。戻り (x, 状態) — 状態は "ok" / "above" (担体が壁の上) /
        "miss" (交わらない)。出口 $L$ より先は接線の延長で解く (交点が $x>L$ なら呼び出し側が出口の外と判定する)。"""
        f = lambda x: self.r(x) - rA - m * (x - xA)
        f0 = f(xA)
        if not np.isfinite(f0):
            return np.nan, "miss"
        if f0 <= 0.0:
            return np.nan, "above"
        if xA < self.L:
            fL = f(self.L)
            if fL > 0.0:                         # 出口までに交わらない → 延長 (直線) 上で解く
                return (self.L + fL / (m - self.rpL), "ok") if m > self.rpL else (np.nan, "miss")
            # [xA, L] に根がある: 最初の符号変化を細かく探してから Brent
            from scipy.optimize import brentq
            sl = m - self.rp(xA)
            step = f0 / sl if sl > 0.0 else 1e-3 * max(self.L - xA, 1e-6)
            step = min(max(step, 1e-12), self.L - xA)
            lo, hi = xA, min(xA + step, self.L)
            while f(hi) > 0.0:
                lo = hi
                step *= 2.0
                hi = min(lo + step, self.L)
            if f(lo) == 0.0:
                return lo, "ok"
            return brentq(f, lo, hi, xtol=1e-15, rtol=4.0 * np.finfo(float).eps, maxiter=200), "ok"
        return (xA + f0 / (m - self.rpL), "ok") if m > self.rpL else (np.nan, "miss")


# ------------------------------------------------------------------------------------------------
# 壁点・軸点の単位過程 (スカラー。1 段に高々 1 個ずつ)
def _wall_unit(A, wall: _Wall, g, delta, analytic, thrA, tol, max_corr):
    """A = (x, r, θ, ν, M) (壁の下の点、前の前線の j = n−1)。戻り dict: state, x, r, th, nu, n_iter
    (収束時は resid_geom・resid_comp も)。出口の外に着地したら state = ST_BEYOND_EXIT。"""
    xA, rA, thA, nuA, MA = (float(v) for v in A)
    muA = float(mu_vec(np.array([MA]))[0])
    thrA_a = None if thrA is None else np.array([thrA])
    m0 = np.tan(thA + muA)
    if not (np.cos(thA + muA) > 0.0):
        return {"state": ST_WALL_MISS, "why": "C⁺ が後ろ向き (θ+μ ≥ 90°)", "n_iter": 0}
    xW, how = wall.intersect(xA, rA, m0)
    if how == "above":
        return {"state": ST_OUTSIDE_WALL, "why": "担体が壁の外", "n_iter": 0}
    if how == "miss":
        return {"state": ST_WALL_MISS, "why": "C⁺ が壁に届かない", "n_iter": 0}
    if xW > wall.L:
        return {"state": ST_BEYOND_EXIT, "x": xW, "n_iter": 0}
    thW = float(np.arctan(wall.rp(xW)))
    nuW = nuA + (thW - thA)                     # 単純波の初期値 (C⁺ 沿いに θ−ν 一定)
    state = None
    for k in range(int(max_corr) + 1):
        MW = float(_mach(np.array([nuW]), g)[0])
        muW = float(mu_vec(np.array([MW]))[0])
        ang = 0.5 * (thA + thW) + 0.5 * (muA + muW)
        if not (np.cos(ang) > 0.0):
            return {"state": ST_WALL_MISS, "why": "C⁺ が後ろ向き (θ+μ ≥ 90°)", "n_iter": k}
        mp = np.tan(ang)
        xn, how = wall.intersect(xA, rA, mp)
        if how == "above":
            return {"state": ST_OUTSIDE_WALL, "why": "担体が壁の外", "n_iter": k}
        if how == "miss" or not np.isfinite(xn):
            return {"state": ST_WALL_MISS, "why": "C⁺ が壁に届かない", "n_iter": k}
        if xn > wall.L:
            return {"state": ST_BEYOND_EXIT, "x": xn, "n_iter": k}
        rn = wall.r(xn)
        thn = float(np.arctan(wall.rp(xn)))
        srcA = float(_src_vec(np.array([rA]), np.array([thA]), np.array([rn]), np.array([thn]), thrA_a, analytic)[0])
        srcW = np.sin(thn) / rn
        fA = 0.5 * (np.sin(muA) * srcA + np.sin(muW) * srcW)
        Sp = -delta * fA / np.cos(ang) * (xn - xA)
        nun = thn - (thA - nuA) - Sp
        if not (np.isfinite(xn) and np.isfinite(nun) and np.isfinite(thn)):
            state = ST_NONFINITE
            xW, thW, nuW = xn, thn, nun
            break
        upd = max(abs(thn - thW), abs(nun - nuW))
        xW, thW, nuW = xn, thn, nun
        if k >= 1 and upd <= tol:
            state = ST_OK
            break
    if state is None:
        state = ST_MAXITER
    out = {"state": state, "x": xW, "r": wall.r(xW), "th": thW, "nu": nuW, "n_iter": k}
    if state == ST_OK:
        MW = float(_mach(np.array([nuW]), g)[0])
        muW = float(mu_vec(np.array([MW]))[0])
        ang = 0.5 * (thA + thW) + 0.5 * (muA + muW)
        rW = wall.r(xW)
        srcA = float(_src_vec(np.array([rA]), np.array([thA]), np.array([rW]), np.array([thW]), thrA_a, analytic)[0])
        fA = 0.5 * (np.sin(muA) * srcA + np.sin(muW) * np.sin(thW) / rW)
        Sp = -delta * fA / np.cos(ang) * (xW - xA)
        out["resid_geom"] = abs(rW - (rA + np.tan(ang) * (xW - xA)))
        out["resid_comp"] = max(abs((thW - nuW) - ((thA - nuA) + Sp)), abs(thW - float(np.arctan(wall.rp(xW)))))
    return out


def axis_thr_neighbor(r1: float, th1: float) -> float:
    r"""計算した軸点に持たせる解析極限 $\theta_r$ の推定 = 同じ前線 (C⁻) の最も軸に近い点 $(r_1,\theta_1)$ の $\sin\theta/r$
    (L'Hôpital を隣の点で代用、誤差 $O(h)$)。

    前進問題には風洞の逆 MOC の軸則 (`target_dM`) が無いので、軸点の $\theta_r$ は計算した場から作るしかない。
    高次の推定は**不安定だった** (2026-10-11、放射源流 γ 1.4・半角 10°、初期線 = x 1.2 の壁点からの厳密な C⁻、
    初期線の点数 n 25/49/97/193、場の誤差は x ≤ 2.4):
    (a) 軸上の $\nu$ 列の 3 点後退差分から $\theta_r=\tfrac12\sqrt{M^2-1}\,d\nu/dx$ → $\theta_r$ の相対誤差が符号を変えながら
    ±0.3 % → ±2.4 % に増え、場の誤差が解像度で減らない (1.2e-4 → 6.9e-5 → 5.0e-5 → 1.0e-4)、
    (b) 近軸 2 点の $\sin\theta/r$ の $r\to0$ への線形外挿 → 前線ごとに約 −1.4 倍で発散 (n 97 で網が反転)、
    (c) $r^2$ の外挿 (偶関数を仮定) → 細かい解像度で次数が 0.75 に落ちる。
    どれも、軸点の $\theta_r$ が次の前線の最下点の $\theta$ を通じて次の軸点の $\theta_r$ に $1/r_1$ で増幅されて戻る結合が原因。
    この代用 (d) は安定で、場の $\nu$ の最大誤差 (軸上) は 2.2e-4 → 6.5e-5 → 2.0e-5 → 6.3e-6 (次数 1.7)、壁は 2 次。"""
    return float(np.sin(th1) / r1)


def _axis_unit(B, g, delta, srcP, tol, max_corr):
    """B = (x, r, θ, ν, M) (軸の上の点、前線の j = 1)。srcP = 軸点の源項 sinθ/r の極限 (`axis_thr_neighbor`)。
    新点の ν と μ の結合だけを反復する。戻り dict: state, x, nu, n_iter (収束時は thr・resid_comp も)。"""
    xB, rB, thB, nuB, MB = (float(v) for v in B)
    if not rB > AXIS_R_EPS:
        return {"state": ST_GEOM, "why": "C⁻ 担体が軸上", "n_iter": 0}
    if not np.isfinite(srcP):
        return {"state": ST_NONFINITE, "why": "軸点の θ_r が非有限", "n_iter": 0}
    muB = float(mu_vec(np.array([MB]))[0])
    srcB = np.sin(thB) / rB
    nuP = thB + nuB
    state = None
    xP = np.nan
    for k in range(int(max_corr) + 1):
        MP = float(_mach(np.array([nuP]), g)[0])
        muP = float(mu_vec(np.array([MP]))[0])
        ang = 0.5 * thB - 0.5 * (muB + muP)
        mm = np.tan(ang)
        if not (mm < 0.0 and np.cos(ang) > 0.0):
            state = ST_GEOM
            break
        xP = xB - rB / mm
        fB = 0.5 * (np.sin(muB) * srcB + np.sin(muP) * srcP)
        Sm = delta * fB / np.cos(ang) * (xP - xB)
        nun = thB + nuB + Sm
        if not (np.isfinite(nun) and np.isfinite(xP)):
            state = ST_NONFINITE
            nuP = nun
            break
        upd = abs(nun - nuP)
        nuP = nun
        if k >= 1 and upd <= tol:
            state = ST_OK
            break
    if state is None:
        state = ST_MAXITER
    out = {"state": state, "x": xP, "nu": nuP, "n_iter": k}
    if state == ST_OK:
        MP = float(_mach(np.array([nuP]), g)[0])
        muP = float(mu_vec(np.array([MP]))[0])
        ang = 0.5 * thB - 0.5 * (muB + muP)
        xP = xB - rB / np.tan(ang)
        fB = 0.5 * (np.sin(muB) * srcB + np.sin(muP) * srcP)
        Sm = delta * fB / np.cos(ang) * (xP - xB)
        out.update(x=xP, thr=float(srcP), resid_comp=abs(nuP - (thB + nuB + Sm)), resid_geom=0.0)
    return out


# ------------------------------------------------------------------------------------------------
# 前進マーチ
class _Net:
    """前線 k × 点 j の配列 (必要に応じて前線方向に伸ばす)。"""
    FIELDS = ("x", "r", "th", "nu", "M")

    def __init__(self, n_pt: int, k_alloc: int):
        self.n_pt = n_pt
        self.a = {f: np.full((k_alloc, n_pt), np.nan) for f in self.FIELDS}
        self.status = np.full((k_alloc, n_pt), ST_NONE, dtype=np.int8)
        self.n_iter = np.full((k_alloc, n_pt), -1, dtype=np.int16)
        self.resid = np.full((k_alloc, n_pt), np.nan)
        self.thr_axis = np.full(k_alloc, np.nan)

    @property
    def k_alloc(self):
        return self.status.shape[0]

    def grow(self, k_need: int):
        if k_need < self.k_alloc:
            return
        new = max(2 * self.k_alloc, k_need + 1)
        pad = new - self.k_alloc
        for f in self.FIELDS:
            self.a[f] = np.vstack([self.a[f], np.full((pad, self.n_pt), np.nan)])
        self.status = np.vstack([self.status, np.full((pad, self.n_pt), ST_NONE, dtype=np.int8)])
        self.n_iter = np.vstack([self.n_iter, np.full((pad, self.n_pt), -1, dtype=np.int16)])
        self.resid = np.vstack([self.resid, np.full((pad, self.n_pt), np.nan)])
        self.thr_axis = np.concatenate([self.thr_axis, np.full(pad, np.nan)])

    def put(self, k, j, x, r, th, nu, M, st):
        self.a["x"][k, j], self.a["r"][k, j], self.a["th"][k, j] = x, r, th
        self.a["nu"][k, j], self.a["M"][k, j] = nu, M
        self.status[k, j] = st

    def trim(self, K: int):
        for f in self.FIELDS:
            self.a[f] = self.a[f][:K + 1]
        self.status, self.n_iter, self.resid = self.status[:K + 1], self.n_iter[:K + 1], self.resid[:K + 1]
        self.thr_axis = self.thr_axis[:K + 1]


def _check_line(line: dict, wall: _Wall, g) -> tuple:
    """初期線の検査 (クリップ前)。戻り (codes, reasons, info)。"""
    codes, reasons = [], []
    x, r, M, th = (np.asarray(line[k], dtype=float) for k in ("x", "r", "M", "theta"))
    info = {"n": int(len(x))}
    if not (len(x) == len(r) == len(M) == len(th)) or len(x) < 4:
        return ["init_shape"], [f"初期線の長さが不正 (x {len(x)}・r {len(r)}・M {len(M)}・θ {len(th)}; 4 点以上)"], info
    fin = np.isfinite(x) & np.isfinite(r) & np.isfinite(M) & np.isfinite(th)
    if not np.all(fin):
        codes.append("init_nonfinite")
        reasons.append(f"初期線に非有限の値 {int(np.sum(~fin))} 点")
    sub = fin & ~(M > 1.0)
    info["M_min"] = float(np.nanmin(M)) if np.any(np.isfinite(M)) else None
    if np.any(sub):
        idx = np.flatnonzero(sub)
        codes.append("init_subsonic")
        reasons.append(f"亜音速入力: 初期線の {len(idx)} 点で M ≤ 1 (最小 {float(np.min(M[sub])):.6g} @ "
                       f"(x, r) = ({x[idx[0]]:.4g}, {r[idx[0]]:.4g}))")
    if fin[-1]:
        try:
            dr = r[-1] - wall.r(float(x[-1]))
            dth = th[-1] - float(np.arctan(wall.rp(float(x[-1]))))
        except ValueError as e:
            dr, dth = np.nan, np.nan
            reasons.append(f"初期線の壁点で壁が評価できない ({e})")
            codes.append("init_wall_eval")
        info.update(wall_dr=float(dr), wall_dtheta=float(dth))
        if not (abs(dr) <= WALL_ON_TOL):
            codes.append("init_off_wall")
            reasons.append(f"初期線の壁点が壁に乗っていない (Δr = {dr:.3e})")
    if fin[0] and not abs(r[0]) <= AXIS_R_EPS:
        codes.append("init_axis")
        reasons.append(f"初期線の始点が軸上にない (r = {r[0]:.3e})")
    if np.any(fin[1:] & fin[:-1] & ~(np.diff(r) > 0.0)):
        codes.append("init_order")
        reasons.append("初期線が軸→壁の順に r 単調増加でない")
    return codes, reasons, info


def forward_moc(wall, L: float, g, line: dict, *, axis_limit: str = "analytic", corrector: str = "converge",
                tol: float = CORR_TOL, max_corr: int = CORR_MAX, delta: float = 1.0,
                max_fronts: int | None = None, n_cross_sample: int = N_CROSS_SAMPLE) -> dict:
    r"""壁 `wall` (`wall.r(x, deriv)` を $x\in[x_{line,wall},L]$ で評価できること、deriv 0・1) の下で、初期線 `line`
    (dict: `x`・`r`・`M`・`theta` を**軸→壁**の順、`thr_axis` = 軸端の $\theta_r$) から前進 MOC で場を作り、ゲートを返す。

    `g`: γ (float) またはガスモデル (`moc_kernel` の規約)。`axis_limit`: `"analytic"` (既定) / `"legacy"`。
    `corrector`: `"converge"` のみ (前進ゲートは収束の分類を理由に使うため)。`max_fronts`: 前線数の上限 (既定 400 n)。
    戻り dict: `gate` (`pass`・`verdict` [pass / reject / undecidable]・`undecidable` [種類と説明]・`codes`・`reasons`・
    `indicators`・`evidence` [不受理の根拠の点と j]・`evidence_j`・`counts` [点の状態の内訳]・`failures_inside`・
    `failure_locations`)、`net` (`_Net`: 前線 × 点の配列 `a[x|r|th|nu|M]`・`status`・`n_iter`・`resid`・`thr_axis`)、
    `end` (マーチの終わり方)、`n`・`n_fronts`・`mdot_star`・`resid_max`・`wall_info`・`axis_info`・`timing`。
    壁点・軸点は `wall_points(res)`・`axis_points(res)`、CSV 1 行の要約は `gate_record(res)`。

    初期線に問題があれば (`gate["stage"] == "initial_line"`) マーチせずに返す (`net` は None)。"""
    t0 = time.time()
    if axis_limit not in ("analytic", "legacy"):
        raise ValueError(f"axis_limit は 'analytic' か 'legacy' (受け取った値: {axis_limit!r})")
    if corrector != "converge":
        raise ValueError(f"前進 MOC ゲートは corrector='converge' のみ (受け取った値: {corrector!r})")
    _check_corr_params(tol, max_corr)
    analytic = axis_limit == "analytic"
    L = float(L)
    W = _Wall(wall, L)
    nu_max = nu_max_of(g)
    codes, reasons, line_info = _check_line(line, W, g)
    base = {"L": L, "axis_limit": axis_limit, "corrector": corrector, "tol": float(tol), "max_corr": int(max_corr),
            "nu_max": nu_max, "line": line_info}
    if codes:
        gate = {"pass": False, "verdict": "undecidable", "codes": codes, "reasons": reasons, "indicators": {},
                "stage": "initial_line",
                "undecidable": {"kind": "initial_line_check", "detail": "初期線の検査で止まった (マーチしない)"},
                "evidence": [], "evidence_j": _evidence_j_summary([])}
        return {**base, "gate": gate, "net": None, "timing": {"total_s": time.time() - t0}}
    x0, r0, M0, th0 = (np.asarray(line[k], dtype=float) for k in ("x", "r", "M", "theta"))
    nu0 = np.asarray(pm_nu(M0, g), dtype=float)
    n = len(x0) - 1                               # 前線の点 j = 0 (軸) .. n (壁)
    thr0 = line.get("thr_axis")
    if analytic and not (thr0 is not None and np.isfinite(thr0)):
        raise ValueError("axis_limit='analytic' には初期線の軸端の θ_r (line['thr_axis']) が必要")
    k_cap_hard = int(max_fronts) if max_fronts is not None else 400 * n
    net = _Net(n + 1, max(64, 2 * n))              # 前線数は実測で ~1.1 n (足りなければ倍々に伸ばす)
    net.put(0, slice(None), x0, r0, th0, nu0, M0, ST_OK)
    net.status[0, :] = ST_OK
    net.thr_axis[0] = float(thr0) if thr0 is not None else np.nan
    k_cap = k_cap_hard                            # 前線の上限 (壁点が出口の外・失敗で確定する)
    end = {"reason": None, "front": None}
    wall_info = {}                                # k → 壁点の単位過程の情報
    axis_info = {}
    res_max = {"geom": 0.0, "comp": 0.0}
    lev = 2 - n
    while True:
        k_lo = max(1, -(-lev // 2))               # ceil(lev/2)
        k_hi = min((lev + n) // 2, k_cap)
        if k_lo > k_cap or (k_hi < k_lo and lev > 2 * k_cap):
            break
        if k_hi >= k_lo:
            net.grow(k_hi + 1)
            ks = np.arange(k_lo, k_hi + 1)
            js = 2 * ks - lev
            # --- 壁点 (j = n) ---
            mw = js == n
            if np.any(mw):
                k = int(ks[mw][0])
                A_st = net.status[k - 1, n - 1]
                if A_st != ST_OK:
                    if A_st == ST_BEYOND_EXIT:
                        end = {"reason": "exit_reached", "front": k - 1, "detail": "壁点の担体が出口の外"}
                    else:
                        end = {"reason": "carrier_lost", "front": k - 1,
                               "detail": f"壁点の担体 G[{k - 1},{n - 1}] が {STATE_NAMES.get(int(A_st), A_st)}"}
                    k_cap = k - 1
                else:
                    A = tuple(net.a[f][k - 1, n - 1] for f in _Net.FIELDS)
                    thrA = net.thr_axis[k - 1] if (n - 1 == 0 and analytic) else None
                    wu = _wall_unit(A, W, g, delta, analytic, thrA, tol, max_corr)
                    wall_info[k] = {kk: vv for kk, vv in wu.items() if kk in ("state", "why", "n_iter")}
                    st = wu["state"]
                    if st == ST_BEYOND_EXIT:
                        end = {"reason": "exit_reached", "front": k - 1, "detail": f"C⁺ の着地 x = {wu['x']:.6g} > L"}
                        k_cap = k - 1
                    elif st == ST_OK:
                        nuW = wu["nu"]
                        stf = ST_SUBSONIC if nuW < 0.0 else (ST_ABOVE_TABLE if nuW > nu_max else ST_OK)
                        MW = float(_mach(np.array([nuW]), g)[0])
                        net.put(k, n, wu["x"], wu["r"], wu["th"], nuW, MW, stf)
                        net.n_iter[k, n] = wu["n_iter"]
                        net.resid[k, n] = max(wu["resid_geom"], wu["resid_comp"])
                        res_max["geom"] = max(res_max["geom"], wu["resid_geom"])
                        res_max["comp"] = max(res_max["comp"], wu["resid_comp"])
                        if stf != ST_OK:
                            end = {"reason": "wall_failed", "front": k, "detail": STATE_NAMES[stf]}
                            k_cap = k
                        elif k >= k_cap_hard:
                            end = {"reason": "front_cap", "front": k, "detail": f"前線の上限 {k_cap_hard}"}
                            k_cap = k
                    else:
                        if "x" in wu:
                            net.a["x"][k, n] = wu["x"]
                        net.status[k, n] = st
                        net.n_iter[k, n] = wu["n_iter"]
                        end = {"reason": "wall_failed", "front": k, "detail": wu.get("why", STATE_NAMES.get(st))}
                        k_cap = k
                        wall_info[k]["xA"], wall_info[k]["rA"] = float(A[0]), float(A[1])
                # 上限が変わったら同じ段の他の前線 (k > k_cap) は落とす
                keep = ks <= k_cap
                ks, js, mw = ks[keep], js[keep], mw[keep]
            # --- 内点 (1 ≤ j ≤ n−1) ---
            mi = (js >= 1) & (js <= n - 1)
            if np.any(mi):
                ki, ji = ks[mi], js[mi]
                Ast = net.status[ki - 1, ji - 1]
                Bst = net.status[ki, ji + 1]
                good = (Ast == ST_OK) & (Bst == ST_OK)
                beyond = (Ast == ST_BEYOND_EXIT) | (Bst == ST_BEYOND_EXIT)
                net.status[ki[~good & beyond], ji[~good & beyond]] = ST_BEYOND_EXIT
                net.status[ki[~good & ~beyond], ji[~good & ~beyond]] = ST_MISSING
                if np.any(good):
                    kg, jg = ki[good], ji[good]
                    A = [net.a[f][kg - 1, jg - 1] for f in _Net.FIELDS]
                    B = [net.a[f][kg, jg + 1] for f in _Net.FIELDS]
                    thrA = None
                    if analytic:
                        thrA = np.where(jg - 1 == 0, net.thr_axis[kg - 1], np.nan)
                    stt: dict = {}
                    xP, rP, thP, nuP, ok = interior_vec(*A, *B, g, delta, 2, axis_limit=axis_limit,
                                                        thrA=thrA, thrB=None, corrector="converge",
                                                        tol=tol, max_corr=max_corr, stats=stt)
                    st = np.array(stt["status"], dtype=np.int8)
                    conv = st == PAIR_CONVERGED
                    st[conv & ~(rP >= -1e-12)] = ST_BELOW_AXIS
                    fin = np.isfinite(xP) & np.isfinite(rP)
                    out_x = fin & (xP > L)
                    rw = np.full(len(xP), np.nan)
                    inside_x = fin & ~out_x
                    if np.any(inside_x):
                        rw[inside_x] = np.asarray(W.w.r(xP[inside_x], 0), dtype=float)
                    ok2 = st == PAIR_CONVERGED
                    st[ok2 & (nuP < 0.0)] = ST_SUBSONIC
                    st[ok2 & (nuP > nu_max)] = ST_ABOVE_TABLE
                    st[(st == PAIR_CONVERGED) & inside_x & (rP > rw + OUTSIDE_WALL_TOL)] = ST_OUTSIDE_WALL
                    st[(st == PAIR_CONVERGED) & out_x] = ST_BEYOND_EXIT
                    MP = _mach(np.where(np.isfinite(nuP), nuP, 0.0), g)
                    for f, v in zip(_Net.FIELDS, (xP, np.maximum(rP, 0.0), thP, nuP, MP)):
                        net.a[f][kg, jg] = v
                    net.status[kg, jg] = st
                    net.n_iter[kg, jg] = stt["n_iter"]
                    rg, rc = stt["resid_geom"], stt["resid_comp"]
                    net.resid[kg, jg] = np.fmax(rg, rc)
                    inside_conv = (st == ST_OK)
                    if np.any(inside_conv):
                        res_max["geom"] = max(res_max["geom"], float(np.nanmax(rg[inside_conv])))
                        res_max["comp"] = max(res_max["comp"], float(np.nanmax(rc[inside_conv])))
            # --- 軸点 (j = 0) ---
            ma = js == 0
            if np.any(ma):
                k = int(ks[ma][0])
                B_st = net.status[k, 1]
                if B_st != ST_OK:
                    net.status[k, 0] = ST_BEYOND_EXIT if B_st == ST_BEYOND_EXIT else ST_MISSING
                else:
                    B = tuple(net.a[f][k, 1] for f in _Net.FIELDS)
                    srcP = axis_thr_neighbor(B[1], B[2])
                    au = _axis_unit(B, g, delta, srcP, tol, max_corr)
                    axis_info[k] = {kk: vv for kk, vv in au.items() if kk in ("state", "why", "n_iter")}
                    st = au["state"]
                    if st == ST_OK:
                        nuA = au["nu"]
                        st = ST_SUBSONIC if nuA < 0.0 else (ST_ABOVE_TABLE if nuA > nu_max else ST_OK)
                        if st == ST_OK and au["x"] > L:
                            st = ST_BEYOND_EXIT
                        MA = float(_mach(np.array([max(nuA, 0.0)]), g)[0])
                        net.put(k, 0, au["x"], 0.0, 0.0, nuA, MA, st)
                        net.n_iter[k, 0] = au["n_iter"]
                        net.resid[k, 0] = au["resid_comp"]
                        if st == ST_OK:
                            net.thr_axis[k] = au["thr"]
                            res_max["comp"] = max(res_max["comp"], au["resid_comp"])
                    else:
                        net.status[k, 0] = st
                        net.n_iter[k, 0] = au["n_iter"]
                        if np.isfinite(au.get("x", np.nan)):
                            net.a["x"][k, 0] = au["x"]
        lev += 1
        if lev > 2 * k_cap + 1:
            break
    K = int(min(k_cap, net.k_alloc - 1))
    # 最後に計算した前線まで切り詰める
    done = np.any(net.status != ST_NONE, axis=1)
    K = int(np.max(np.flatnonzero(done))) if np.any(done) else 0
    net.trim(K)
    if end["reason"] is None:
        end = {"reason": "levels_exhausted", "front": K, "detail": "段を使い切った"}
    t_march = time.time() - t0
    out = {**base, "n": n, "n_fronts": K, "end": end, "net": net, "wall_info": wall_info, "axis_info": axis_info,
           "resid_max": res_max, "mdot_star": float(_line_flux(x0, r0, th0, M0, g))}
    out["gate"] = _gate(out, W, g, n_cross_sample)
    out["timing"] = {"march_s": t_march, "total_s": time.time() - t0}
    return out


# ------------------------------------------------------------------------------------------------
# 診断とゲート
def _line_flux(x, r, th, M, g) -> float:
    """折れ線 (軸→壁) を横切る質量流量 (moc_inverse._flux_along と同じ式、配列版)。"""
    Mm = 0.5 * (M[1:] + M[:-1])
    thm = 0.5 * (th[1:] + th[:-1])
    rm = 0.5 * (r[1:] + r[:-1])
    d = (_mass_flux_density(Mm, g) * (np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x)) * 2.0 * np.pi * rm)
    return float(np.sum(d))


def _signed_len(xP, xQ, ang):
    """特性線の向き (角 ang) に測った Q→P の符号付き長さ (P−Q は向き (cos, sin) と平行なので (xP−xQ)/cos)。"""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (xP - xQ) / np.cos(ang)


def _evidence_j_summary(ev: list) -> dict:
    """根拠の点の j の要約 (j = 0 が軸、None は位置が分からない根拠)。"""
    js = [e["j"] for e in ev if e.get("j") is not None]
    return {"n": len(ev), "n_unlocated": sum(1 for e in ev if e.get("j") is None),
            "j_min": min(js) if js else None, "j_max": max(js) if js else None,
            "n_j_le_near_axis": sum(1 for j in js if j <= NEAR_AXIS_J)}


def _gate(res: dict, W: _Wall, g, n_cross_sample: int) -> dict:
    net: _Net = res["net"]
    L = res["L"]
    n = res["n"]
    X, R, TH, NU, MM = (net.a[f] for f in _Net.FIELDS)
    S = net.status
    K = S.shape[0] - 1
    codes, reasons = [], []
    evidence = []                                # 不受理の根拠の点 (種類・k・j・x・r; j は軸からの添字)

    def _ev(kind, k, j, x=None, r=None):
        evidence.append({"kind": kind, "k": None if k is None else int(k), "j": None if j is None else int(j),
                         "x": None if x is None or not np.isfinite(x) else float(x),
                         "r": None if r is None or not np.isfinite(r) else float(r)})
    MU = mu_vec(np.where(np.isfinite(MM), MM, 1.0))
    ok = S == ST_OK
    inside = ok & (X <= L)
    # --- 1. 単位過程の失敗 (位置: 計算できた点はその x、できなかった点は担体の x の最大で判定) ---
    fail_in = {}
    fail_loc = []
    for k in range(1, K + 1):
        for j in np.flatnonzero(np.isin(S[k], FAIL_STATES)):
            st = int(S[k, j])
            if j == n:
                carriers = [(k - 1, n - 1)]
            elif j == 0:
                carriers = [(k, 1)]
            else:
                carriers = [(k - 1, j - 1), (k, j + 1)]
            xc = float(np.nanmax([X[a, b] for a, b in carriers]))
            xp = float(X[k, j]) if np.isfinite(X[k, j]) else np.nan
            where_x = xp if np.isfinite(xp) else xc
            if where_x <= L:
                name = STATE_NAMES[st]
                fail_in[name] = fail_in.get(name, 0) + 1
                rr = float(R[k, j]) if np.isfinite(R[k, j]) else float(R[carriers[0]])
                _ev(f"unit_{name}", k, j, where_x, rr)
                if len(fail_loc) < 20:
                    fail_loc.append({"k": k, "j": int(j), "state": name, "x": where_x, "r": rr})
    for name, cnt in sorted(fail_in.items()):
        codes.append(f"unit_{name}")
        label = {"geom_parallel": "幾何的棄却 (特性線が平行)", "iter_nonfinite": "反復の失敗 (非有限)",
                 "iter_maxiter": "反復の失敗 (上限到達)", "below_axis": "新点が軸より下",
                 "subsonic": "亜音速 (ν < 0)", "above_table": "表の外 (ν > ν_max)",
                 "wall_miss": "C⁺ が壁に届かない", "outside_wall": "新点が壁の外"}.get(name, name)
        reasons.append(f"{label} {cnt} 点 (x ≤ L)")
    # --- 2. 最終残差 ---
    rmax = max(res["resid_max"]["geom"], res["resid_max"]["comp"])
    if not rmax <= RESID_TOL:
        codes.append("resid")
        reasons.append(f"最終残差 {rmax:.3e} > {RESID_TOL:g}")
        rz = np.where(inside & np.isfinite(net.resid), net.resid, -np.inf)
        if np.any(np.isfinite(rz)):
            kr, jr = np.unravel_index(int(np.argmax(rz)), rz.shape)
            _ev("resid", kr, jr, X[kr, jr], R[kr, jr])
        else:
            _ev("resid", None, None)
    # --- 3. 網の反転 (担体より特性線の向きで手前) ---
    # C⁺: P = (k, j) ← A = (k−1, j−1)  (j = 1..n)
    Px, Pth, Pmu = X[1:, 1:], TH[1:, 1:], MU[1:, 1:]
    Ax, Ath, Amu = X[:-1, :-1], TH[:-1, :-1], MU[:-1, :-1]
    both_p = ok[1:, 1:] & ok[:-1, :-1] & (X[1:, 1:] <= L)
    ang_p = 0.5 * (Ath + Pth) + 0.5 * (Amu + Pmu)
    d_p = np.where(both_p, _signed_len(Px, Ax, ang_p), np.nan)
    # C⁻: P = (k, j) ← B = (k, j+1)  (j = 0..n−1)
    Bx, Bth, Bmu = X[:, 1:], TH[:, 1:], MU[:, 1:]
    Qx, Qth, Qmu = X[:, :-1], TH[:, :-1], MU[:, :-1]
    both_m = ok[:, 1:] & ok[:, :-1] & (X[:, :-1] <= L)
    ang_m = 0.5 * (Bth + Qth) - 0.5 * (Bmu + Qmu)
    e_m = np.where(both_m, _signed_len(Qx, Bx, ang_m), np.nan)
    rev_p = both_p & ~(d_p > 0.0)
    rev_m = both_m & ~(e_m > 0.0)
    n_rev_p = int(np.sum(rev_p))
    n_rev_m = int(np.sum(rev_m))
    init_rev = []                                # 前線 0→1 の反転 (初期線不整合の判定)
    for (kk, jj) in zip(*np.nonzero(rev_p)):
        _ev("reversal_cplus_step", kk + 1, jj + 1, X[kk + 1, jj + 1], R[kk + 1, jj + 1])
        if kk + 1 <= 1:
            init_rev.append(("reversal_cplus_step", int(kk + 1), int(jj + 1)))
    for (kk, jj) in zip(*np.nonzero(rev_m)):
        _ev("reversal_cminus_step", kk, jj, X[kk, jj], R[kk, jj])
        if kk <= 1:
            init_rev.append(("reversal_cminus_step", int(kk), int(jj)))
    onset = {}
    for name, rv, XX, RR, dk in (("cplus_step", rev_p, Px, R[1:, 1:], 1), ("cminus_step", rev_m, Qx, R[:, :-1], 0)):
        if np.any(rv):
            i = int(np.argmin(np.where(rv, XX, np.inf)))
            ki, ji = np.unravel_index(i, XX.shape)
            onset[name] = {"x": float(XX.flat[i]), "r": float(RR.flat[i]), "k": int(ki + dk), "j": int(ji + dk)}
    if n_rev_p:
        codes.append("reversal_cplus_step")
        reasons.append(f"網の反転: C⁺ 沿いに新点が担体の手前 {n_rev_p} 点 (C⁻ 族の前線どうしの交差)")
    if n_rev_m:
        codes.append("reversal_cminus_step")
        reasons.append(f"網の反転: C⁻ 沿いに新点が担体の手前 {n_rev_m} 点 (C⁺ 族どうしの交差)")
    kw = np.flatnonzero(ok[:, n])
    ka = np.flatnonzero(ok[:, 0] & (X[:, 0] <= L))
    wx = X[kw, n]
    ax = X[ka, 0]
    bad_w = np.flatnonzero(~(np.diff(wx) > 0.0)) if len(wx) > 1 else np.zeros(0, dtype=int)
    bad_a = np.flatnonzero(~(np.diff(ax) > 0.0)) if len(ax) > 1 else np.zeros(0, dtype=int)
    n_wall_nonmono = int(len(bad_w))
    n_axis_nonmono = int(len(bad_a))
    for i in bad_w:
        _ev("reversal_wall", kw[i + 1], n, X[kw[i + 1], n], R[kw[i + 1], n])
        if kw[i + 1] <= 1:
            init_rev.append(("reversal_wall", int(kw[i + 1]), int(n)))
    for i in bad_a:
        _ev("reversal_axis", ka[i + 1], 0, X[ka[i + 1], 0], 0.0)
        if ka[i + 1] <= 1:
            init_rev.append(("reversal_axis", int(ka[i + 1]), 0))
    if n_wall_nonmono:
        codes.append("reversal_wall")
        reasons.append(f"網の反転: 壁点の x が単調でない {n_wall_nonmono} 箇所")
    if n_axis_nonmono:
        codes.append("reversal_axis")
        reasons.append(f"網の反転: 軸点の x が単調でない {n_axis_nonmono} 箇所 (前線が軸の近くで交差)")
    # --- 4. 非隣接の同族特性線の交差 (間引き、moc_diagnostics の折れ線交差) ---
    n_cross_p, n_cross_m, cross_loc = _sampled_crossings(X, R, inside, n, n_cross_sample)
    for c in cross_loc:
        # 交差した線分の端点の j のうち最大を根拠の j にする (両端とも j ≤ 2 のときだけ軸近傍に数える)
        _ev(f"crossing_{c['family']}", c["k"], c["j_max"], c["x"], c["r"])
        evidence[-1]["j_min"] = c["j_min"]
    if n_cross_m:
        codes.append("crossing_cminus")
        reasons.append(f"同族特性線の交差: C⁻ 族 (前線) の非隣接の組 {n_cross_m} (間引き {n_cross_sample} 本)")
    if n_cross_p:
        codes.append("crossing_cplus")
        reasons.append(f"同族特性線の交差: C⁺ 族の非隣接の組 {n_cross_p} (間引き {n_cross_sample} 本)")
    # --- 5. 出口まで届いたか ---
    end = res["end"]
    if end["reason"] != "exit_reached":
        codes.append("exit_not_reached" if end["reason"] != "wall_failed" else "wall_not_reached")
        last_w = float(wx[-1]) if len(wx) else float("nan")
        reasons.append(f"出口まで届かない ({end['reason']}: {end.get('detail')}; 最後の壁点 x = {last_w:.4g} / L = {L:.4g})")
        fk = end.get("front")
        if end["reason"] == "wall_failed":
            _ev("exit_not_reached", fk, n, last_w)
        elif end["reason"] == "carrier_lost":
            _ev("exit_not_reached", fk, n - 1, X[fk, n - 1] if fk is not None else None)
        else:                                    # 前線数の上限・段の使い切り: 位置が分からない
            _ev("exit_not_reached", fk, None)
    # --- 指標 (連続量) ---
    ind = _indicators(X, R, TH, MU, ok, d_p, e_m, n, L)
    mdot = _front_mdot(res, X, R, TH, MM, ok, n, L, g)
    counts = {STATE_NAMES[int(s)]: int(np.sum(S == s)) for s in np.unique(S)}
    # --- 判定 (3 値) ---
    ev_j = _evidence_j_summary(evidence)
    undecidable = None
    if init_rev:
        verdict = "undecidable"
        kinds = sorted({t[0] for t in init_rev})
        undecidable = {"kind": "initial_line",
                       "detail": f"前線 0→1 に網の反転 {len(init_rev)} 点 ({', '.join(kinds)}) — 初期線と前進の適合式の不整合"
                                 f" (壁の不受理として数えない)", "points": init_rev[:20]}
        reasons.append("判定不能 (初期線不整合): " + undecidable["detail"])
    elif not codes:
        verdict = "pass"
    elif ev_j["n"] > 0 and ev_j["n_unlocated"] == 0 and ev_j["j_max"] <= NEAR_AXIS_J:
        verdict = "undecidable"
        undecidable = {"kind": "near_axis",
                       "detail": f"不受理の根拠 {ev_j['n']} 点がすべて軸から j ≤ {NEAR_AXIS_J} (j の最大 {ev_j['j_max']}) — "
                                 f"前進の軸点の θ_r は 1 次の代用なので壁の不受理として数えない"}
        reasons.append("判定不能 (軸近傍): " + undecidable["detail"])
    else:
        verdict = "reject"
    gate = {"pass": verdict == "pass", "verdict": verdict, "undecidable": undecidable,
            "codes": codes, "reasons": reasons, "stage": "march",
            "indicators": {**ind, "n_reversal_cplus_step": n_rev_p, "n_reversal_cminus_step": n_rev_m,
                           "n_wall_nonmonotone": n_wall_nonmono, "n_axis_nonmonotone": n_axis_nonmono,
                           "n_crossing_cplus_sampled": n_cross_p, "n_crossing_cminus_sampled": n_cross_m,
                           "resid_max": rmax, "reversal_onset": onset, **mdot},
            "evidence": evidence[:200], "evidence_j": ev_j,
            "failures_inside": fail_in, "failure_locations": fail_loc, "counts": counts,
            "last_wall_x": float(wx[-1]) if len(wx) else None, "end": end}
    return gate


def _indicators(X, R, TH, MU, ok, d_p, e_m, n, L) -> dict:
    r"""同族特性線の間隔の比 (内部衝撃波の指標)。
    C⁻ 族: $d_k[j]/d_k[n]$ ($d$ = 前線 $k-1\to k$ の C⁺ 沿いの間隔、分母は壁の C⁺ 区間)。
    C⁺ 族: $e_k[j]/e_{k_0}[j_0]$ ($e$ = 隣り合う C⁺ の C⁻ 沿いの間隔、分母はその 2 本が最初にそろった前線での値)。"""
    out = {}
    # C⁻ 族 (d_p[k−1, j−1] が P=(k, j) の値; 壁 j = n は d_p[k−1, n−1])
    dw = d_p[:, n - 1:n]                         # 前線 k の壁の区間
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_m = d_p / dw
    valid = np.isfinite(ratio_m)
    valid[:, n - 1] = False                      # 壁そのもの (=1) は除く
    if np.any(valid):
        i = np.argmin(np.where(valid, ratio_m, np.inf))
        k1, j1 = np.unravel_index(i, ratio_m.shape)
        out.update(cminus_sep_ratio_min=float(ratio_m[k1, j1]),
                   cminus_sep_ratio_at={"front": int(k1 + 1), "j": int(j1 + 1), "x": float(X[k1 + 1, j1 + 1]),
                                        "r": float(R[k1 + 1, j1 + 1])},
                   cminus_sep_ratio_p01=float(np.percentile(ratio_m[valid], 1.0)))
    else:
        out.update(cminus_sep_ratio_min=None, cminus_sep_ratio_at=None, cminus_sep_ratio_p01=None)
    # C⁺ 族: e_m[k, j] は (k, j+1)→(k, j) の C⁻ 沿いの間隔。C⁺ の組 (a = k−j, a−1)
    Kp1, nm = e_m.shape
    kk, jj = np.meshgrid(np.arange(Kp1), np.arange(nm), indexing="ij")
    a = kk - jj
    k0 = np.maximum(a, 0)
    j0 = k0 - a
    ref = np.full(e_m.shape, np.nan)
    inb = (k0 < Kp1) & (j0 < nm)
    ref[inb] = e_m[k0[inb], j0[inb]]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_p = e_m / ref
    valid = np.isfinite(ratio_p) & ~((kk == k0) & (jj == j0))
    if np.any(valid):
        i = np.argmin(np.where(valid, ratio_p, np.inf))
        k1, j1 = np.unravel_index(i, ratio_p.shape)
        out.update(cplus_sep_ratio_min=float(ratio_p[k1, j1]),
                   cplus_sep_ratio_at={"front": int(k1), "j": int(j1), "x": float(X[k1, j1]), "r": float(R[k1, j1])})
    else:
        out.update(cplus_sep_ratio_min=None, cplus_sep_ratio_at=None)
    return out


def _sampled_crossings(X, R, inside, n, n_sample) -> tuple:
    """C⁺ 族 (対角線 k−j 一定) と C⁻ 族 (前線) から n_sample 本ずつ間引き、非隣接 (2 本以上離れた) 組の折れ線交差を数える。
    折れ線は x ≤ L の収束点だけで作る (欠けた点で切る)。交差判定は moc_diagnostics の `_polylines_cross_aabb`。
    戻り (C⁺ の交差組数, C⁻ の交差組数, 交差の位置のリスト)。位置は交差した組ごとに 1 件: 交差した線分 (全組) の端点の
    j の最小・最大と、x 最小の交差線分の始点 (k, j, x, r)。組数の判定は位置の探索と独立 (`_polylines_cross_aabb` のまま)。"""
    from .moc_diagnostics import _polylines_cross_aabb, _segments_cross_matrix
    K1 = X.shape[0]

    def runs(ks, js):
        # 連続した有効区間のうち最長のもの (欠けで切れた線は最長の区間だけを見る)。戻り (点 (x, r), k, j)
        best, cur = None, []
        for kv, jv in zip(ks, js):
            if inside[kv, jv]:
                cur.append((kv, jv))
            else:
                if len(cur) >= 2 and (best is None or len(cur) > len(best)):
                    best = cur
                cur = []
        if len(cur) >= 2 and (best is None or len(cur) > len(best)):
            best = cur
        if best is None:
            return None
        kb = np.array([t[0] for t in best])
        jb = np.array([t[1] for t in best])
        return np.column_stack([X[kb, jb], R[kb, jb]]), kb, jb

    # C⁻ 族: 前線 k (j = n → 0 の順)
    fronts = []
    jrev = np.arange(n, -1, -1)
    for k in np.unique(np.linspace(0, K1 - 1, min(n_sample, K1)).astype(int)):
        ln = runs(np.full(n + 1, k), jrev)
        if ln is not None:
            fronts.append((int(k), ln))
    # C⁺ 族: 線 a = k − j (a ∈ [−n, K]); 点 (k, k − a)
    a_all = np.arange(-n, K1)
    lines_p = []
    for a in np.unique(np.linspace(a_all[0], a_all[-1], min(n_sample, len(a_all))).astype(int)):
        k = np.arange(max(a, 0), min(K1, a + n + 1))
        j = k - a
        m = (j >= 0) & (j <= n)
        k, j = k[m], j[m]
        if len(k) < 2:
            continue
        ln = runs(k, j)
        if ln is not None:
            lines_p.append((int(a), ln))

    def locate(p, q):
        # 交差した線分の組 (全組) → 端点の j の最小・最大と、x 最小の交差線分の始点
        (P, kp, jp), (Q, kq, jq) = p, q
        M = _segments_cross_matrix(P[:-1], P[1:], Q[:-1], Q[1:])
        ip, iq = np.nonzero(M)
        if len(ip) == 0:                          # `_polylines_cross_aabb` と同じ定義なので起きない想定
            return None
        jj = np.concatenate([jp[ip], jp[ip + 1], jq[iq], jq[iq + 1]])
        i0 = int(np.argmin(P[ip, 0]))
        return {"j_min": int(jj.min()), "j_max": int(jj.max()), "k": int(kp[ip[i0]]), "j": int(jp[ip[i0]]),
                "x": float(P[ip[i0], 0]), "r": float(P[ip[i0], 1])}

    def count(lines, min_gap, family, locs):
        c = 0
        boxes = [(l[0][:, 0].min(), l[0][:, 0].max(), l[0][:, 1].min(), l[0][:, 1].max()) for _, l in lines]
        for p in range(len(lines)):
            for q in range(p + 1, len(lines)):
                if abs(lines[q][0] - lines[p][0]) < min_gap:
                    continue                      # 隣接 (共有点を持ちうる) は除く
                b1, b2 = boxes[p], boxes[q]
                if b1[1] < b2[0] or b2[1] < b1[0] or b1[3] < b2[2] or b2[3] < b1[2]:
                    continue
                if _polylines_cross_aabb(lines[p][1][0], lines[q][1][0]):
                    c += 1
                    loc = locate(lines[p][1], lines[q][1])
                    if loc is None:
                        loc = {"j_min": None, "j_max": None, "k": None, "j": None, "x": None, "r": None}
                    locs.append({"family": family, "lines": (lines[p][0], lines[q][0]), **loc})
        return c

    locs: list = []
    n_p = count(lines_p, 2, "cplus", locs)
    n_m = count(fronts, 2, "cminus", locs)
    return n_p, n_m, locs


def _front_mdot(res, X, R, TH, MM, ok, n, L, g) -> dict:
    """完全な前線 (軸→壁の全点が x ≤ L で収束) を横切る質量流量の初期線比。前進 MOC は流量を構成的に保存しないので、
    非循環の精度の指標 (逆 MOC の流線整合残差に当たる)。"""
    ms = res["mdot_star"]
    full = np.all(ok & (X <= L), axis=1)
    vals = []
    for k in np.flatnonzero(full):
        vals.append(_line_flux(X[k], R[k], TH[k], MM[k], g) / ms - 1.0)
    if not vals:
        return {"mdot_err_max": None, "mdot_err_last": None, "n_full_fronts": 0}
    v = np.asarray(vals)
    return {"mdot_err_max": float(np.max(np.abs(v))), "mdot_err_last": float(v[-1]), "n_full_fronts": int(len(v))}


# ------------------------------------------------------------------------------------------------
# 取り出し
def wall_points(res) -> np.ndarray:
    """壁点 [x, r, θ, ν, M] (収束・x ≤ L)。"""
    net = res["net"]
    n = res["n"]
    m = (net.status[:, n] == ST_OK)
    return np.column_stack([net.a[f][m, n] for f in _Net.FIELDS])


def axis_points(res) -> np.ndarray:
    net = res["net"]
    m = (net.status[:, 0] == ST_OK)
    return np.column_stack([net.a[f][m, 0] for f in _Net.FIELDS])


def gate_record(res) -> dict:
    """CSV 1 行ぶんの要約 (ゲートの合否・理由・指標)。"""
    gt = res["gate"]
    ind = gt.get("indicators", {})
    at = ind.get("cminus_sep_ratio_at") or {}
    atp = ind.get("cplus_sep_ratio_at") or {}
    ons = ind.get("reversal_onset") or {}
    on0 = min(ons.values(), key=lambda v: v["x"]) if ons else {}
    evj = gt.get("evidence_j") or {}
    return {"gate_pass": bool(gt["pass"]), "verdict": gt.get("verdict"),
            "undecidable_kind": (gt.get("undecidable") or {}).get("kind"),
            "codes": ";".join(gt["codes"]), "reasons": " / ".join(gt["reasons"]),
            "end_reason": (res.get("end") or {}).get("reason"), "n_fronts": res.get("n_fronts"),
            "last_wall_x": gt.get("last_wall_x"),
            "cminus_sep_ratio_min": ind.get("cminus_sep_ratio_min"), "cminus_sep_x": at.get("x"),
            "cminus_sep_r": at.get("r"), "cminus_sep_j": at.get("j"), "cminus_sep_ratio_p01": ind.get("cminus_sep_ratio_p01"),
            "cplus_sep_ratio_min": ind.get("cplus_sep_ratio_min"), "cplus_sep_x": atp.get("x"), "cplus_sep_j": atp.get("j"),
            "reversal_onset_x": on0.get("x"), "reversal_onset_r": on0.get("r"), "reversal_onset_j": on0.get("j"),
            "n_evidence": evj.get("n"), "evidence_j_min": evj.get("j_min"), "evidence_j_max": evj.get("j_max"),
            "n_evidence_j_le_near_axis": evj.get("n_j_le_near_axis"), "n_evidence_unlocated": evj.get("n_unlocated"),
            "n_reversal_cplus_step": ind.get("n_reversal_cplus_step"),
            "n_reversal_cminus_step": ind.get("n_reversal_cminus_step"),
            "n_crossing_cminus_sampled": ind.get("n_crossing_cminus_sampled"),
            "n_crossing_cplus_sampled": ind.get("n_crossing_cplus_sampled"),
            "resid_max": ind.get("resid_max"), "mdot_err_max": ind.get("mdot_err_max"),
            "time_s": (res.get("timing") or {}).get("total_s")}
