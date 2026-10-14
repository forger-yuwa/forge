r"""逆 MOC の単位過程 (`moc_kernel.interior_vec` を 1 段の対に通す部分) の numba 版。

計画: plans/active/tooling-nozzle-mknot-fast.md §4.4・§4.5 (受入は §6 V5)。

`moc_inverse.InverseMOC._pairs` の 1 段分 — 予測子・修正子 (`_interior_pass`)、最終残差 (`_interior_residual`)、源項の端点評価
(`_src_vec`・`_sin_over_r_vec`)、`mu_vec`、ν→M (`pm_mach_vec`)、源項の分岐 (`source_branches`)、対の分類 — を `pairs_nb` にした。
段の中は対ごとのスカラーのループ。式と演算の順序は numpy の版と同じだが、**ビット同一は求めない** (§4.5、ユーザ決定 2026-10-15):

- ν→M の CPG の Newton は対ごとに止め (|ΔM| < 1e-13、上限 60 回)、前の反復の M から始める (numpy の版は配列全体の収束で止め、
  毎回 1 + max(ν, 1e-6) から始める)。semi-perfect の表は二分探索と線形の式 (範囲外は端の値、NaN は NaN)。
- 段のループ (前線の対・軸より下の棄却・次の前線) と集計 (`MocFillDiag.add`) は Python 側 (`moc_inverse`) にあり、両方の backend が
  同じものを使う。ここは数値と、`add` に渡す対ごとの生の状態・反復回数・最終残差・分岐、次の前線の M だけを返す。

この模組は numba を import する。呼び出し側は numba の版を使うときだけ import する (numba が無ければ ImportError)。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from .moc_kernel import PAIR_CONVERGED, PAIR_DONE, PAIR_GEOM, PAIR_MAXITER, PAIR_MISSING, PAIR_NONFINITE

NEWTON_TOL = 1e-13           # pm_mach_vec の既定 tol
NEWTON_ITERS = 60            # pm_mach_vec の既定 iters
GAS_CPG = 0                  # ν→M = CPG の Newton (γ は float、GasCPG・CPG の Isentrope は gamma_ref)
GAS_TABLE = 1                # ν→M = 表の線形補間 (GasSemiPerfect の ν・M の表)
MODE_LEGACY_FIXED2 = 0       # legacy + fixed2: 分類は moc_inverse._posthoc_stats と同じ
MODE_FIXED2 = 1              # analytic + fixed2: 分類は interior_vec の stats と同じ
MODE_CONVERGE = 2            # converge (legacy/analytic)
BRANCH_KEYS = ("axis_analytic", "axis_borrow", "axis_zero", "frac_A", "frac_B", "frac_PA", "frac_PB", "P_on_axis")


@njit(cache=True)
def _nmax(a, b):
    """`np.maximum` (どちらかが NaN なら NaN)。"""
    if a >= b or a != a:
        return a
    return b


@njit(cache=True)
def _mach(nu, M0, gas_kind, g, tx, ty):
    """ν → M (1 点)。CPG: `pm_mach_vec` と同じ式の Newton を M0 から (M0 が 1 以下・非有限なら 1 + max(ν, 1e-6) から)。
    表: tx (ν) の上の線形補間。NaN は NaN。"""
    if nu != nu:
        return nu
    if gas_kind == GAS_TABLE:
        n = tx.size
        if nu <= tx[0]:
            return ty[0]
        if nu >= tx[n - 1]:
            return ty[n - 1]
        lo = 0
        hi = n - 1                           # tx[lo] ≤ ν < tx[hi]
        while hi - lo > 1:
            mid = (lo + hi) >> 1
            if tx[mid] <= nu:
                lo = mid
            else:
                hi = mid
        return ty[lo] + (ty[hi] - ty[lo]) * ((nu - tx[lo]) / (tx[hi] - tx[lo]))
    a = math.sqrt((g + 1.0) / (g - 1.0))
    c = 0.5 * (g - 1.0)
    M = M0 if (M0 > 1.0 and math.isfinite(M0)) else 1.0 + _nmax(nu, 1e-6)
    for _ in range(NEWTON_ITERS):
        b = math.sqrt(_nmax(M * M - 1.0, 0.0))
        f = (a * math.atan(b / a) - math.atan(b)) - nu
        df = b / (1.0 + c * M * M) / M
        dM = f / _nmax(df, 1e-14)
        M = M - dM
        if M <= 1.0:
            M = 1.0 + 1e-9
        if abs(dM) < NEWTON_TOL:
            break
    return M


@njit(cache=True)
def _mu(M):
    """`mu_vec`。"""
    return math.asin(1.0 / _nmax(M, 1.0 + 1e-12))


@njit(cache=True)
def _src(r_p, th_p, r_o, th_o, thr_p, use_thr, frac, eps):
    """`_src_vec` (`_sin_over_r_vec` に、analytic の真の軸端点の θ_r を足したもの)。"""
    if use_thr and r_p <= eps and math.isfinite(thr_p):
        return thr_p
    if (r_p <= 1e-9) or (r_p < frac * r_o):     # 軸上・AXIS_LIMIT_FRAC の分岐は相手の値
        r_p = r_o
        th_p = th_o
    if r_p > 1e-9:
        return math.sin(th_p) / r_p
    return 0.0


@njit(cache=True)
def _src_P(rPc, sP, Ar, bA, frac):
    """新点 P の端点の sinθ/r (`_src_vec(rPc, thP, Ar, Ath, None, …)`): sP = sin θ_P、bA = 相手の値 (軸上・AXIS_LIMIT_FRAC の分岐)。"""
    if (rPc <= 1e-9) or (rPc < frac * Ar):
        return bA
    if rPc > 1e-9:
        return sP / rPc
    return 0.0                                   # rPc が NaN


@njit(cache=True)
def _pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, thP, nuP, MP, delta, frac):
    """`_interior_pass` の 1 対分 (MP = M(nuP))。反復で変わらない端点の項は外で作る: fA0 = sin μ_A·(sinθ/r)_A、
    fB0 = sin μ_B·(sinθ/r)_B、bA・bB = P が軸上のときに借りる A・B の sinθ/r。戻り: xP, rP, thP, nuP, good, parallel。"""
    muP = _mu(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    mp = math.tan(ang_p)
    mm = math.tan(ang_m)
    den = mm - mp
    good = abs(den) >= 1e-12
    parallel = math.isfinite(den) and not good
    if not good:
        den = 1.0
    xP = (Ar - Br + mm * Bx - mp * Ax) / den
    rP = Ar + mp * (xP - Ax)
    rPc = _nmax(rP, 0.0)
    sP = math.sin(thP)
    smP = math.sin(muP)
    fA = 0.5 * (fA0 + smP * _src_P(rPc, sP, Ar, bA, frac))
    fB = 0.5 * (fB0 + smP * _src_P(rPc, sP, Br, bB, frac))
    Sp = -delta * fA / math.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / math.cos(ang_m) * (xP - Bx)
    Jp = (Ath - Anu) + Sp
    Jm = (Bth + Bnu) + Sm
    return xP, rP, 0.5 * (Jp + Jm), 0.5 * (Jm - Jp), good, parallel


@njit(cache=True)
def _resid(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, xP, rP, thP, nuP, MP, delta, frac):
    """`_interior_residual` の 1 対分 (端点の項は `_pass` と同じ)。戻り: (幾何の残差, 適合式の残差)。"""
    muP = _mu(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    mp = math.tan(ang_p)
    mm = math.tan(ang_m)
    rg = _nmax(abs(rP - (Ar + mp * (xP - Ax))), abs(rP - (Br + mm * (xP - Bx))))
    rPc = _nmax(rP, 0.0)
    sP = math.sin(thP)
    smP = math.sin(muP)
    fA = 0.5 * (fA0 + smP * _src_P(rPc, sP, Ar, bA, frac))
    fB = 0.5 * (fB0 + smP * _src_P(rPc, sP, Br, bB, frac))
    Sp = -delta * fA / math.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / math.cos(ang_m) * (xP - Bx)
    rc = _nmax(abs((thP - nuP) - ((Ath - Anu) + Sp)), abs((thP + nuP) - ((Bth + Bnu) + Sm)))
    return rg, rc


@njit(cache=True)
def pairs_nb(x, r, th, nu, M, thrA, thrB, analytic, mode, n_corr, tol, max_corr, delta, gas_kind, g, tx, ty, frac, eps):
    r"""前線 (x, r, th, nu, M; 長さ m + 1) の隣接対 (A = 前線[j + 1] が C⁺ 担体、B = 前線[j] が C⁻ 担体) を単位過程に通す
    (`moc_inverse.InverseMOC._pairs` の 1 段分)。thrA・thrB は対ごとの θ_r (analytic のときだけ使う、無ければ NaN)。

    戻り: (n_lack, xP, rP, thP, nuP, ok, status, n_iter, resid_geom, resid_comp, branch, MP)。
    n_lack > 0 は analytic で θ_r の無い軸上の既知点の対の数 (呼び出し側が interior_vec と同じ ValueError を出す; 他の欄は未計算)。
    ok・status・n_iter・resid_*・branch は interior_vec (と legacy + fixed2 の `_posthoc_stats`) と同じ意味 (branch の行は BRANCH_KEYS)。
    MP = M(nuP) (ok の対だけ、他は NaN; 次の前線の M)。"""
    m = x.size - 1
    xP = np.full(m, np.nan)
    rP = np.full(m, np.nan)
    thP = np.empty(m)
    nuP = np.empty(m)
    MP = np.full(m, np.nan)
    ok = np.zeros(m, np.bool_)
    status = np.empty(m, np.int8)
    n_it = np.zeros(m, np.int16)
    rg = np.full(m, np.nan)
    rc = np.full(m, np.nan)
    br = np.zeros((8, m), np.int8)
    n_lack = 0
    for j in range(m):
        fin = (math.isfinite(x[j + 1]) and math.isfinite(r[j + 1]) and math.isfinite(th[j + 1]) and math.isfinite(nu[j + 1])
               and math.isfinite(M[j + 1]) and math.isfinite(x[j]) and math.isfinite(r[j]) and math.isfinite(th[j])
               and math.isfinite(nu[j]) and math.isfinite(M[j]))
        if analytic and fin and ((r[j + 1] <= eps and not math.isfinite(thrA[j])) or (r[j] <= eps and not math.isfinite(thrB[j]))):
            n_lack += 1
    if n_lack > 0:
        return n_lack, xP, rP, thP, nuP, ok, status, n_it, rg, rc, br, MP
    for j in range(m):
        Ax, Ar, Ath, Anu, AM = x[j + 1], r[j + 1], th[j + 1], nu[j + 1], M[j + 1]
        Bx, Br, Bth, Bnu, BM = x[j], r[j], th[j], nu[j], M[j]
        tA, tB = thrA[j], thrB[j]
        fin_in = (math.isfinite(Ax) and math.isfinite(Ar) and math.isfinite(Ath) and math.isfinite(Anu) and math.isfinite(AM)
                  and math.isfinite(Bx) and math.isfinite(Br) and math.isfinite(Bth) and math.isfinite(Bnu)
                  and math.isfinite(BM))
        tP = 0.5 * (Ath + Bth)
        nP = 0.5 * (Anu + Bnu)
        muA = _mu(AM)
        muB = _mu(BM)
        fA0 = math.sin(muA) * _src(Ar, Ath, Br, Bth, tA, analytic, frac, eps)    # 反復で変わらない端点の項
        fB0 = math.sin(muB) * _src(Br, Bth, Ar, Ath, tB, analytic, frac, eps)
        bA = _src(Ar, Ath, Ar, Ath, np.nan, False, frac, eps)                     # P が軸上のとき借りる A・B の sinθ/r
        bB = _src(Br, Bth, Br, Bth, np.nan, False, frac, eps)
        Mg = 0.5 * (AM + BM)                     # Newton の初期値 (以後は前の反復の M)
        xp = np.nan
        rp = np.nan
        if mode != MODE_CONVERGE:
            # 予測 1 + 修正 n_corr 回
            good_all = True
            par = False
            for _ in range(1 + n_corr):
                Mg = _mach(nP, Mg, gas_kind, g, tx, ty)
                xp, rp, tP, nP, good, para = _pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, tP, nP, Mg,
                                                   delta, frac)
                good_all = good_all and good
                par = par or para
            fin_out = math.isfinite(xp) and math.isfinite(rp) and math.isfinite(tP) and math.isfinite(nP)
            okj = good_all and fin_out
            st = PAIR_DONE
            if mode == MODE_LEGACY_FIXED2:       # moc_inverse._posthoc_stats (出力の有限性で分ける)
                if fin_in and not okj:
                    st = PAIR_GEOM if fin_out else PAIR_NONFINITE
            else:                                # interior_vec の stats (fixed2)
                if fin_in and ((not good_all) or (not fin_out)):
                    st = PAIR_NONFINITE
                if fin_in and par:
                    st = PAIR_GEOM
            if not fin_in:
                st = PAIR_MISSING
            n_it[j] = n_corr
        else:
            # converge: 修正 1 回の θ・ν の更新量が tol 以下になるまで (上限 max_corr 回)。入力が非有限の対は回さない
            st = PAIR_MISSING
            if fin_in:
                st = -1
                for kk in range(max_corr + 1):   # kk = 0 予測子、kk ≥ 1 修正子
                    Mg = _mach(nP, Mg, gas_kind, g, tx, ty)
                    a0, a1, a2, a3, good, para = _pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, tP, nP,
                                                       Mg, delta, frac)
                    fin = math.isfinite(a0) and math.isfinite(a1) and math.isfinite(a2) and math.isfinite(a3)
                    upd = _nmax(abs(a2 - tP), abs(a3 - nP))
                    xp, rp, tP, nP = a0, a1, a2, a3
                    n_it[j] = kk
                    if para:
                        st = PAIR_GEOM
                    elif (not good) or (not fin):
                        st = PAIR_NONFINITE
                    elif kk >= 1 and upd <= tol:
                        st = PAIR_CONVERGED
                    if st >= 0:
                        break
                if st < 0:
                    st = PAIR_MAXITER
            okj = st == PAIR_CONVERGED
        xP[j] = xp
        rP[j] = rp
        thP[j] = tP
        nuP[j] = nP
        ok[j] = okj
        status[j] = st
        if okj:
            Mg = _mach(nP, Mg, gas_kind, g, tx, ty)
            MP[j] = Mg
            if mode == MODE_CONVERGE:
                rg[j], rc[j] = _resid(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, xp, rp, tP, nP, Mg,
                                      delta, frac)
        # 源項の分岐 (source_branches、最終状態の幾何で)
        rPc = _nmax(rp, 0.0)
        nA = Ar <= eps
        nB = Br <= eps
        hA = nA and analytic and math.isfinite(tA)
        hB = nB and analytic and math.isfinite(tB)
        br[0, j] = np.int8(hA) + np.int8(hB)
        br[1, j] = np.int8(nA and not hA and Br > eps) + np.int8(nB and not hB and Ar > eps)
        br[2, j] = np.int8(nA and not hA and not (Br > eps)) + np.int8(nB and not hB and not (Ar > eps))
        br[3, j] = (Ar > eps) and (Ar < frac * Br)
        br[4, j] = (Br > eps) and (Br < frac * Ar)
        br[5, j] = (rPc > eps) and (rPc < frac * Ar)
        br[6, j] = (rPc > eps) and (rPc < frac * Br)
        br[7, j] = rPc <= eps
    return n_lack, xP, rP, thP, nuP, ok, status, n_it, rg, rc, br, MP
