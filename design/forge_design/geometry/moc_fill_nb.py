r"""逆 MOC の網の充填 (`moc_inverse.InverseMOC.fill_levels`) の段のループの numba 版。

計画: plans/active/tooling-nozzle-mknot-fast.md §4.4・§5.1 #4 (受入は §6 V5)。

`fill_levels` の段のループ全体 — 単位過程 (`moc_kernel.interior_vec` の予測子・修正子 `_interior_pass`、最終残差
`_interior_residual`、源項の端点評価 `_src_vec`・`_sin_over_r_vec`、`mu_vec`、ν→M の `pm_mach_vec`)、源項の分岐
(`source_branches`)、対の集計 (`moc_inverse.MocFillDiag.add`) — を 1 つの `@njit` 関数 `fill_levels_nb` にした。
段の中は対ごとのスカラーのループで、numpy 版と**同じ式・同じ演算の順序** (float64、`fastmath=False`)。
この模組は numba を import する。呼び出し側 (`moc_inverse`) は numba の版を使うときだけこの模組を import する
(numba が無ければ ImportError — 黙って numpy の版に落ちない)。

numpy 版とビット同一にするための約束 (どちらかを変えるときは両方を変える):

- ν→M の CPG の Newton (`pm_mach_vec`) は、**渡された配列全体**で |ΔM| < 1e-13 になるまで回す (最大 60 回)。
  要素ごとに止めると、先に収束した要素の追加の反復 (1 ulp 動きうる) が numpy 版と変わる。そこで Newton だけは numpy 版と
  同じ集合 (fixed2 は段の全部の対、converge は反復中の対、最終残差は収束した対、段の終わりは段の全部の点 [ok でない点は
  ν = 0]) について配列で回し、その後の単位過程は対ごとに回す。NaN の要素はその後もずっと NaN なので計算を飛ばす
  (配列全体の停止が起きないのも numpy 版と同じ: 60 回まで回る)。
- semi-perfect (`GasSemiPerfect`) の ν→M は表の線形補間 `np.interp` を、numpy の C 実装 (`binary_search_with_guess`、
  傾きの式、NaN のときの逆向きの式) のとおりに写した (要素ごとに独立)。
- `np.maximum` は NaN を伝える (`_nmax`)。`np.where` の選択は if 文。
- 超越関数 (sin・cos・tan・arcsin・arctan・sqrt) は numba が libm を呼ぶ。numpy の umath が libm と同じ値を返す環境
  (手元の X86_V3 で確認: 2026-10-15) ではビット同一。numpy が SIMD の別実装 (例: AVX512 の SVML) を使う環境では
  ulp の違いが出うる (受入は plan §6 V5 の許容差)。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from .moc_kernel import (PAIR_BELOW_AXIS, PAIR_CONVERGED, PAIR_DONE, PAIR_GEOM, PAIR_MAXITER, PAIR_MISSING,
                         PAIR_NONFINITE)

NEWTON_TOL = 1e-13           # pm_mach_vec の既定 tol (変えるなら moc_kernel.pm_mach_vec と両方)
NEWTON_ITERS = 60            # pm_mach_vec の既定 iters
GAS_CPG = 0                  # ν→M = CPG の Newton (γ は float、GasCPG・cpg の Isentrope は gamma_ref)
GAS_TABLE = 1                # ν→M = 表の線形補間 (GasSemiPerfect の np.interp(ν, _nu, _M))
MODE_LEGACY_FIXED2 = 0       # 既定 (legacy + fixed2): 分類は moc_inverse._posthoc_stats
MODE_FIXED2 = 1              # analytic + fixed2: 分類は interior_vec の stats
MODE_CONVERGE = 2            # converge (legacy/analytic)
ERR_LACK_THR = 1             # analytic で θ_r の無い軸上の既知点 (interior_vec の ValueError)
_LIKELY_IN_CACHE_SIZE = 8    # numpy の binary_search_with_guess の定数
# 集計の並び (moc_inverse.MocFillDiag と同じ名前・同じ順)
COUNT_KEYS = ("input", "missing", "geom_parallel", "geom_below_axis", "iter_nonfinite", "iter_maxiter", "converged",
              "done_fixed")
BRANCH_KEYS = ("axis_analytic", "axis_borrow", "axis_zero", "frac_A", "frac_B", "frac_PA", "frac_PB", "P_on_axis")
FRAC_KEYS = ("frac_A", "frac_B", "frac_PA", "frac_PB")
GEOM_KINDS = ("parallel", "below_axis")
FAIL_KINDS = ("nonfinite", "maxiter")


def no_trace():
    """`fill_levels_nb` の末尾の引数 (trace を使わないとき): trace=False と大きさ 1 の配列。"""
    return (False, np.zeros((1, 1), np.int8), np.zeros((1, 1), np.int64), np.zeros((1, 1)), np.zeros((1, 1)),
            np.zeros((1, 8, 1), np.int8))


def trace_arrays(n: int):
    """`fill_levels_nb` の末尾の引数 (試験用の trace、n = 初期前線の点の数)。書かれない欄は状態・回数 −1、残差 NaN、分岐 0。"""
    return (True, np.full((n, n), -1, np.int8), np.full((n, n), -1, np.int64), np.full((n, n), np.nan),
            np.full((n, n), np.nan), np.zeros((n, 8, n), np.int8))


@njit(cache=True)
def _nmax(a, b):
    """`np.maximum` (どちらかが NaN なら NaN)。"""
    if a >= b or a != a:
        return a
    return b


@njit(cache=True)
def _mach_newton(nus, cnt, g, Mo):
    """`pm_mach_vec` (CPG) の Newton。対象は nus[:cnt] (numpy 版の配列と同じ集合)。配列全体で |dM| < tol になったら止める。"""
    a = np.sqrt((g + 1.0) / (g - 1.0))
    c = 0.5 * (g - 1.0)
    for t in range(cnt):
        Mo[t] = 1.0 + _nmax(nus[t], 1e-6)
    for _ in range(NEWTON_ITERS):
        conv = True
        for t in range(cnt):
            M = Mo[t]
            if M != M:                       # NaN は以後も NaN (numpy 版も同じ値; 停止しないことだけを写す)
                conv = False
                continue
            b = np.sqrt(_nmax(M * M - 1.0, 0.0))
            f = (a * np.arctan(b / a) - np.arctan(b)) - nus[t]
            df = np.sqrt(_nmax(M * M - 1.0, 0.0)) / (1.0 + c * M * M) / M
            dM = f / _nmax(df, 1e-14)
            M = M - dM
            if M <= 1.0:
                M = 1.0 + 1e-9
            Mo[t] = M
            if not (abs(dM) < NEWTON_TOL):
                conv = False
        if conv:
            break


@njit(cache=True)
def _bsearch_guess(key, arr, n, guess):
    """numpy の `binary_search_with_guess` (compiled_base.c) の写し。arr[i] <= key < arr[i+1] の i (範囲外は -1・n)。"""
    imin = 0
    imax = n
    if key > arr[n - 1]:
        return n
    elif key < arr[0]:
        return -1
    if n <= 4:
        i = 1
        while i < n and key >= arr[i]:
            i += 1
        return i - 1
    if guess > n - 3:
        guess = n - 3
    if guess < 1:
        guess = 1
    if key < arr[guess]:
        if key < arr[guess - 1]:
            imax = guess - 1
            if guess > _LIKELY_IN_CACHE_SIZE and key >= arr[guess - _LIKELY_IN_CACHE_SIZE]:
                imin = guess - _LIKELY_IN_CACHE_SIZE
        else:
            return guess - 1
    else:
        if key < arr[guess + 1]:
            return guess
        else:
            if key < arr[guess + 2]:
                return guess + 1
            else:
                imin = guess + 2
                if guess < n - _LIKELY_IN_CACHE_SIZE - 1 and key < arr[guess + _LIKELY_IN_CACHE_SIZE]:
                    imax = guess + _LIKELY_IN_CACHE_SIZE
    while imin < imax:
        imid = imin + ((imax - imin) >> 1)
        if key >= arr[imid]:
            imin = imid + 1
        else:
            imax = imid
    return imin - 1


@njit(cache=True)
def _interp_np(xs, cnt, xp, fp, out):
    """`np.interp(xs[:cnt], xp, fp)` (left・right は既定) の写し (numpy の arr_interp と同じ分岐・同じ式)。"""
    lenxp = xp.size
    lval = fp[0]
    rval = fp[lenxp - 1]
    if lenxp == 1:
        for i in range(cnt):
            x = xs[i]
            out[i] = lval if x < xp[0] else (rval if x > xp[0] else fp[0])
        return
    j = 0
    for i in range(cnt):
        x = xs[i]
        if x != x:
            out[i] = x
            continue
        j = _bsearch_guess(x, xp, lenxp, j)
        if j == -1:
            out[i] = lval
        elif j == lenxp:
            out[i] = rval
        elif j == lenxp - 1:
            out[i] = fp[j]
        elif xp[j] == x:
            out[i] = fp[j]
        else:
            slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j])
            v = slope * (x - xp[j]) + fp[j]
            if v != v:
                v = slope * (x - xp[j + 1]) + fp[j + 1]
                if v != v and fp[j] == fp[j + 1]:
                    v = fp[j]
            out[i] = v


@njit(cache=True)
def _mach(nus, cnt, gas_kind, g, tx, ty, Mo):
    """`pm_mach_vec(nus[:cnt], g)` → Mo[:cnt]。CPG の Newton と表の補間は別の経路。"""
    if gas_kind == GAS_CPG:
        _mach_newton(nus, cnt, g, Mo)
    else:
        _interp_np(nus, cnt, tx, ty, Mo)


@njit(cache=True)
def _mu(M):
    """`mu_vec`。"""
    return math.asin(1.0 / _nmax(M, 1.0 + 1e-12))


@njit(cache=True)
def _sor(r_p, th_p, r_o, th_o, frac):
    """`_sin_over_r_vec` (点自身で評価、軸上と AXIS_LIMIT_FRAC の分岐だけ相手の値)。"""
    if (r_p <= 1e-9) or (r_p < frac * r_o):
        ra = r_o
        tha = th_o
    else:
        ra = r_p
        tha = th_p
    if ra > 1e-9:
        return math.sin(tha) / ra
    return 0.0


@njit(cache=True)
def _src(r_p, th_p, r_o, th_o, thr_p, use_thr, frac, eps):
    """`_src_vec` (use_thr = analytic かつ θ_r を渡す端点)。"""
    v = _sor(r_p, th_p, r_o, th_o, frac)
    if use_thr and r_p <= eps and math.isfinite(thr_p):
        return thr_p
    return v


@njit(cache=True)
def _pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, thrA, thrB, thP, nuP, MP, delta, analytic, frac, eps):
    """`_interior_pass` の 1 対分 (MP = pm_mach(nuP) は配列の Newton で先に求めたもの)。
    戻り: xP, rP, thP, nuP, good, parallel。"""
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
    fA = 0.5 * (math.sin(muA) * _src(Ar, Ath, Br, Bth, thrA, analytic, frac, eps)
                + math.sin(muP) * _src(rPc, thP, Ar, Ath, np.nan, False, frac, eps))
    fB = 0.5 * (math.sin(muB) * _src(Br, Bth, Ar, Ath, thrB, analytic, frac, eps)
                + math.sin(muP) * _src(rPc, thP, Br, Bth, np.nan, False, frac, eps))
    Sp = -delta * fA / math.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / math.cos(ang_m) * (xP - Bx)
    Jp = (Ath - Anu) + Sp
    Jm = (Bth + Bnu) + Sm
    return xP, rP, 0.5 * (Jp + Jm), 0.5 * (Jm - Jp), good, parallel


@njit(cache=True)
def _resid(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, thrA, thrB, xP, rP, thP, nuP, MP, delta, analytic, frac, eps):
    """`_interior_residual` の 1 対分。戻り: (幾何の残差, 適合式の残差)。"""
    muP = _mu(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    mp = math.tan(ang_p)
    mm = math.tan(ang_m)
    rg = _nmax(abs(rP - (Ar + mp * (xP - Ax))), abs(rP - (Br + mm * (xP - Bx))))
    rPc = _nmax(rP, 0.0)
    fA = 0.5 * (math.sin(muA) * _src(Ar, Ath, Br, Bth, thrA, analytic, frac, eps)
                + math.sin(muP) * _src(rPc, thP, Ar, Ath, np.nan, False, frac, eps))
    fB = 0.5 * (math.sin(muB) * _src(Br, Bth, Ar, Ath, thrB, analytic, frac, eps)
                + math.sin(muP) * _src(rPc, thP, Br, Bth, np.nan, False, frac, eps))
    Sp = -delta * fA / math.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / math.cos(ang_m) * (xP - Bx)
    rc = _nmax(abs((thP - nuP) - ((Ath - Anu) + Sp)), abs((thP + nuP) - ((Bth + Bnu) + Sm)))
    return rg, rc


@njit(cache=True)
def _branch(Ar, Br, rP, thrA, thrB, analytic, frac, eps, br, j):
    """`source_branches` の 1 対分を br[:, j] に (並びは BRANCH_KEYS)。"""
    rPc = _nmax(rP, 0.0)
    nA = Ar <= eps
    nB = Br <= eps
    hA = nA and analytic and math.isfinite(thrA)
    hB = nB and analytic and math.isfinite(thrB)
    br[0, j] = np.int8(hA) + np.int8(hB)
    br[1, j] = np.int8(nA and not hA and Br > eps) + np.int8(nB and not hB and Ar > eps)
    br[2, j] = np.int8(nA and not hA and not (Br > eps)) + np.int8(nB and not hB and not (Ar > eps))
    br[3, j] = (Ar > eps) and (Ar < frac * Br)
    br[4, j] = (Br > eps) and (Br < frac * Ar)
    br[5, j] = (rPc > eps) and (rPc < frac * Ar)
    br[6, j] = (rPc > eps) and (rPc < frac * Br)
    br[7, j] = rPc <= eps


@njit(cache=True)
def fill_levels_nb(out, x0, r0, th0, nu0, M0, thr0, analytic, mode, n_corr, tol, max_corr, delta,
                   gas_kind, g, tx, ty, frac, eps, cap,
                   trace, tr_status, tr_iter, tr_rg, tr_rc, tr_branch):
    r"""`InverseMOC.fill_levels` の段のループ (段 k = 1, 2, …)。`out` (n, n, 5) の段 0 は呼び出し側が埋めておく。

    入力: 初期前線 x0・r0・th0・nu0・M0、θ_r `thr0` (analytic のときの初期前線の値、それ以外は使わない)、
    `mode` (MODE_*)、`n_corr`、converge の `tol`・`max_corr`、`delta`、ν→M (`gas_kind`・γ `g`・表 `tx`=ν・`ty`=M)、
    `frac` = AXIS_LIMIT_FRAC、`eps` = AXIS_R_EPS、`cap` = MocFillDiag.CAP。
    `trace` が真なら段 k・対 j の状態 (軸より下の上書きの後)・反復回数・最終残差・分岐を tr_* [k, j] に書く (試験用)。

    戻り: (err, err_n, counts, hist, hist_order, resid_max, resid_at_i, resid_at_f, resid_has, branch, frac_rng,
    fp_i, fp_f, g_i, g_f, fl_i, fl_f, n_levels) — MocFillDiag の各欄の材料 (moc_inverse._diag_from_nb が組み立てる)。
    err = ERR_LACK_THR のとき err_n はその段の θ_r の無い対の数 (呼び出し側が interior_vec と同じ ValueError を出す)。"""
    n = x0.size
    nb = max(n, 1)
    x = x0.copy()
    r = r0.copy()
    th = th0.copy()
    nu = nu0.copy()
    M = M0.copy()
    live = np.ones(nb, np.bool_)
    pair = np.zeros(nb, np.bool_)
    fin_in = np.zeros(nb, np.bool_)
    xP = np.empty(nb)
    rP = np.empty(nb)
    thP = np.empty(nb)
    nuP = np.empty(nb)
    muA = np.empty(nb)
    muB = np.empty(nb)
    thrA = np.empty(nb)
    thrB = np.empty(nb)
    okk = np.zeros(nb, np.bool_)
    par = np.zeros(nb, np.bool_)
    status = np.zeros(nb, np.int8)
    n_it = np.zeros(nb, np.int64)
    rg = np.empty(nb)
    rc = np.empty(nb)
    br = np.zeros((8, nb), np.int8)
    act = np.zeros(nb, np.int64)
    nus = np.empty(nb)
    Ms = np.empty(nb)
    # --- 集計 (MocFillDiag の欄) ---
    counts = np.zeros(8, np.int64)
    H = max(n_corr, max_corr, 0) + 1
    hist = np.zeros(H, np.int64)
    hist_lvl = np.zeros(H, np.int64)
    hist_order = np.zeros(H, np.int64)
    n_order = 0
    resid_max = np.zeros(2)
    resid_at_i = np.full((2, 2), -1, np.int64)
    resid_at_f = np.full((2, 2), np.nan)
    resid_has = np.zeros(2, np.bool_)
    branch = np.zeros(8, np.int64)
    frac_rng = np.array([np.inf, -np.inf, np.inf, -np.inf])
    fp_i = np.zeros((max(cap, 1), 3), np.int64)
    fp_f = np.zeros((max(cap, 1), 3))
    n_fp = 0
    g_i = np.zeros((64, 3), np.int64)
    g_f = np.zeros((64, 4))
    n_g = 0
    fl_i = np.zeros((max(cap, 1), 4), np.int64)
    fl_f = np.zeros((max(cap, 1), 4))
    n_fl = 0
    n_levels = 0
    err = 0
    err_n = 0
    for k in range(1, n):
        m = n - k                            # 段 k の対の数 (= 前の段の点の数 − 1)
        anyp = False
        for j in range(m):
            pair[j] = live[j] and live[j + 1]
            if pair[j]:
                anyp = True
        if not anyp:
            break
        # --- 入力の有限性と θ_r (interior_vec の冒頭) ---
        lack = 0
        for j in range(m):
            # A = front[1:] (C⁺ 担体)、B = front[:-1] (C⁻ 担体)
            fin_in[j] = (math.isfinite(x[j + 1]) and math.isfinite(r[j + 1]) and math.isfinite(th[j + 1])
                         and math.isfinite(nu[j + 1]) and math.isfinite(M[j + 1]) and math.isfinite(x[j])
                         and math.isfinite(r[j]) and math.isfinite(th[j]) and math.isfinite(nu[j])
                         and math.isfinite(M[j]))
            if k == 1:
                thrA[j] = thr0[j + 1]
                thrB[j] = thr0[j]
            else:                            # 生成した点は軸外 (θ_r なし = NaN)
                thrA[j] = np.nan
                thrB[j] = np.nan
            if analytic and fin_in[j] and ((r[j + 1] <= eps and not math.isfinite(thrA[j]))
                                           or (r[j] <= eps and not math.isfinite(thrB[j]))):
                lack += 1
        if analytic and lack > 0:
            err = ERR_LACK_THR
            err_n = lack
            break
        for j in range(m):
            thP[j] = 0.5 * (th[j + 1] + th[j])
            nuP[j] = 0.5 * (nu[j + 1] + nu[j])
            muA[j] = _mu(M[j + 1])
            muB[j] = _mu(M[j])
        if mode != MODE_CONVERGE:
            # --- fixed2: 予測 1 + 修正 n_corr 回 (全部の対) ---
            for j in range(m):
                xP[j] = 0.0
                rP[j] = 0.0
                okk[j] = True
                par[j] = False
            for _ in range(1 + n_corr):
                for j in range(m):
                    nus[j] = nuP[j]
                _mach(nus, m, gas_kind, g, tx, ty, Ms)
                for j in range(m):
                    a0, a1, a2, a3, good, para = _pass(x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j],
                                                       nu[j], muA[j], muB[j], thrA[j], thrB[j], thP[j], nuP[j], Ms[j],
                                                       delta, analytic, frac, eps)
                    xP[j] = a0
                    rP[j] = a1
                    thP[j] = a2
                    nuP[j] = a3
                    okk[j] = okk[j] and good
                    par[j] = par[j] or para
            for j in range(m):
                fin_out = (math.isfinite(xP[j]) and math.isfinite(rP[j]) and math.isfinite(thP[j])
                           and math.isfinite(nuP[j]))
                ok_ret = okk[j] and fin_out
                fi = fin_in[j]
                st = PAIR_DONE
                if mode == MODE_LEGACY_FIXED2:
                    # moc_inverse._posthoc_stats (出力の有限性で幾何的棄却と反復中の非有限を分ける)
                    if fi and not ok_ret and fin_out:
                        st = PAIR_GEOM
                    if fi and not ok_ret and not fin_out:
                        st = PAIR_NONFINITE
                else:
                    # interior_vec の stats (fixed2)
                    if fi and ((not okk[j]) or (not fin_out)):
                        st = PAIR_NONFINITE
                    if fi and par[j]:
                        st = PAIR_GEOM
                if not fi:
                    st = PAIR_MISSING
                status[j] = st
                n_it[j] = n_corr
                rg[j] = np.nan
                rc[j] = np.nan
                okk[j] = ok_ret
        else:
            # --- converge: 対ごとに収束まで (入力が有限の対だけで判定、収束した対はその時点で凍結) ---
            na = 0
            for j in range(m):
                xP[j] = np.nan
                rP[j] = np.nan
                n_it[j] = 0
                if fin_in[j]:
                    status[j] = -1
                    act[na] = j
                    na += 1
                else:
                    status[j] = PAIR_MISSING
            for kk in range(max_corr + 1):          # kk = 0 予測子、kk ≥ 1 修正子 kk 回目
                if na == 0:
                    break
                for t in range(na):
                    nus[t] = nuP[act[t]]
                _mach(nus, na, gas_kind, g, tx, ty, Ms)
                nn = 0
                for t in range(na):
                    j = act[t]
                    a0, a1, a2, a3, good, para = _pass(x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j],
                                                       nu[j], muA[j], muB[j], thrA[j], thrB[j], thP[j], nuP[j], Ms[t],
                                                       delta, analytic, frac, eps)
                    fin = math.isfinite(a0) and math.isfinite(a1) and math.isfinite(a2) and math.isfinite(a3)
                    upd = _nmax(abs(a2 - thP[j]), abs(a3 - nuP[j]))
                    xP[j] = a0
                    rP[j] = a1
                    thP[j] = a2
                    nuP[j] = a3
                    n_it[j] = kk
                    st = -1
                    if (not para) and ((not good) or (not fin)):
                        st = PAIR_NONFINITE
                    if para:
                        st = PAIR_GEOM
                    if kk >= 1 and good and fin and upd <= tol:
                        st = PAIR_CONVERGED
                    status[j] = st
                    if st < 0:
                        act[nn] = j
                        nn += 1
                na = nn
            for t in range(na):
                status[act[t]] = PAIR_MAXITER
            cnt = 0
            for j in range(m):
                rg[j] = np.nan
                rc[j] = np.nan
                okk[j] = status[j] == PAIR_CONVERGED
                if okk[j]:
                    act[cnt] = j
                    nus[cnt] = nuP[j]
                    cnt += 1
            if cnt > 0:                              # 最終残差 (収束した対の配列で Newton)
                _mach(nus, cnt, gas_kind, g, tx, ty, Ms)
                for t in range(cnt):
                    j = act[t]
                    rg[j], rc[j] = _resid(x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j], nu[j], muA[j],
                                          muB[j], thrA[j], thrB[j], xP[j], rP[j], thP[j], nuP[j], Ms[t], delta,
                                          analytic, frac, eps)
        for j in range(m):
            _branch(r[j + 1], r[j], rP[j], thrA[j], thrB[j], analytic, frac, eps, br, j)
        # --- 軸より下 (fill_levels) と状態の上書き (MocFillDiag.add) ---
        anyok = False
        for j in range(m):
            rpj = rP[j]
            below = pair[j] and math.isfinite(rpj) and not (rpj >= -1e-12)
            st = status[j]
            if below and (st == PAIR_CONVERGED or st == PAIR_DONE):
                status[j] = PAIR_BELOW_AXIS
            okk[j] = okk[j] and pair[j] and (rpj >= -1e-12)
            if okk[j]:
                anyok = True
        # --- 集計 (MocFillDiag.add と同じ順) ---
        counts[0] += m
        for v in range(H):
            hist_lvl[v] = 0
        for j in range(m):
            st = status[j]
            if st == PAIR_MISSING:
                counts[1] += 1
            elif st == PAIR_GEOM:
                counts[2] += 1
            elif st == PAIR_BELOW_AXIS:
                counts[3] += 1
            elif st == PAIR_NONFINITE:
                counts[4] += 1
            elif st == PAIR_MAXITER:
                counts[5] += 1
            elif st == PAIR_CONVERGED:
                counts[6] += 1
            elif st == PAIR_DONE:
                counts[7] += 1
            if st == PAIR_CONVERGED or st == PAIR_DONE:
                hist_lvl[n_it[j]] += 1
        for v in range(H):                           # 辞書の挿入順も numpy 版 (段ごとに np.unique の昇順) と同じ
            if hist_lvl[v] > 0:
                if hist[v] == 0:
                    hist_order[n_order] = v
                    n_order += 1
                hist[v] += hist_lvl[v]
        if m > 0:
            for q in range(2):
                lmax = 0.0
                jb = -1
                for j in range(m):
                    v = rg[j] if q == 0 else rc[j]
                    if not math.isfinite(v):
                        v = -1.0
                    if jb < 0 or v > lmax:           # 最初の最大 (np.argmax)
                        lmax = v
                        jb = j
                if lmax > resid_max[q]:
                    resid_max[q] = lmax
                    resid_has[q] = True
                    resid_at_i[q, 0] = k
                    resid_at_i[q, 1] = jb
                    resid_at_f[q, 0] = xP[jb]
                    resid_at_f[q, 1] = rP[jb]
        for j in range(m):
            if status[j] != PAIR_MISSING:
                branch[0] += br[0, j]
                branch[1] += br[1, j]
                branch[2] += br[2, j]
        for q in range(4):                           # frac_A, frac_B, frac_PA, frac_PB の順、各々 j の昇順
            for j in range(m):
                if br[3 + q, j] and status[j] != PAIR_MISSING:
                    branch[3 + q] += 1
                    if q == 0:
                        xp, rp, ro = x[j + 1], r[j + 1], r[j]
                    elif q == 1:
                        xp, rp, ro = x[j], r[j], r[j + 1]
                    elif q == 2:
                        xp, rp, ro = xP[j], _nmax(rP[j], 0.0), r[j + 1]
                    else:
                        xp, rp, ro = xP[j], _nmax(rP[j], 0.0), r[j]
                    if xp < frac_rng[0]:
                        frac_rng[0] = xp
                    if xp > frac_rng[1]:
                        frac_rng[1] = xp
                    if rp < frac_rng[2]:
                        frac_rng[2] = rp
                    if rp > frac_rng[3]:
                        frac_rng[3] = rp
                    if n_fp < cap:
                        fp_i[n_fp, 0] = k
                        fp_i[n_fp, 1] = j
                        fp_i[n_fp, 2] = q
                        fp_f[n_fp, 0] = xp
                        fp_f[n_fp, 1] = rp
                        fp_f[n_fp, 2] = ro
                        n_fp += 1
        for j in range(m):
            if br[7, j] and (status[j] == PAIR_CONVERGED or status[j] == PAIR_DONE):
                branch[7] += 1
        for q in range(2):                           # 幾何的棄却 (全数): 平行、軸より下 の順
            code = PAIR_GEOM if q == 0 else PAIR_BELOW_AXIS
            for j in range(m):
                if status[j] == code:
                    if n_g == g_i.shape[0]:
                        gi2 = np.zeros((2 * n_g, 3), np.int64)
                        gf2 = np.zeros((2 * n_g, 4))
                        gi2[:n_g] = g_i
                        gf2[:n_g] = g_f
                        g_i = gi2
                        g_f = gf2
                    g_i[n_g, 0] = k
                    g_i[n_g, 1] = j
                    g_i[n_g, 2] = q
                    g_f[n_g, 0] = x[j + 1]
                    g_f[n_g, 1] = r[j + 1]
                    g_f[n_g, 2] = x[j]
                    g_f[n_g, 3] = r[j]
                    n_g += 1
        for q in range(2):                           # 反復の失敗 (先頭 cap 件): 非有限、上限 の順
            code = PAIR_NONFINITE if q == 0 else PAIR_MAXITER
            for j in range(m):
                if status[j] == code and n_fl < cap:
                    fl_i[n_fl, 0] = k
                    fl_i[n_fl, 1] = j
                    fl_i[n_fl, 2] = q
                    fl_i[n_fl, 3] = n_it[j]
                    fl_f[n_fl, 0] = x[j + 1]
                    fl_f[n_fl, 1] = r[j + 1]
                    fl_f[n_fl, 2] = x[j]
                    fl_f[n_fl, 3] = r[j]
                    n_fl += 1
        if k > n_levels:
            n_levels = k
        if trace:
            for j in range(m):
                tr_status[k, j] = status[j]
                tr_iter[k, j] = n_it[j]
                tr_rg[k, j] = rg[j]
                tr_rc[k, j] = rc[j]
                for q in range(8):
                    tr_branch[k, q, j] = br[q, j]
        if not anyok:
            break
        # --- 次の段の前線 (ok でない点は NaN; M は段の全部の点で ν→M、ok でない点は ν = 0) ---
        for j in range(m):
            nus[j] = nuP[j] if okk[j] else 0.0
        _mach(nus, m, gas_kind, g, tx, ty, Ms)
        for j in range(m):
            if okk[j]:
                x[j] = xP[j]
                r[j] = rP[j]
                th[j] = thP[j]
                nu[j] = nuP[j]
                M[j] = Ms[j]
            else:
                x[j] = np.nan
                r[j] = np.nan
                th[j] = np.nan
                nu[j] = np.nan
                M[j] = np.nan
            out[k, j, 0] = x[j]
            out[k, j, 1] = r[j]
            out[k, j, 2] = th[j]
            out[k, j, 3] = nu[j]
            out[k, j, 4] = M[j]
            live[j] = okk[j]
    return (err, err_n, counts, hist, hist_order[:n_order], resid_max, resid_at_i, resid_at_f, resid_has, branch,
            frac_rng, fp_i[:n_fp], fp_f[:n_fp], g_i[:n_g], g_f[:n_g], fl_i[:n_fl], fl_f[:n_fl], n_levels)
