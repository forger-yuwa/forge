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

NEWTON_TOL = 1e-13           # pm_mach_vec の既定 tol (緩めると M → 1 の近くで |ΔM| が 4e-9 に達する、2026-10-15 試験)
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
def _sin_mu(M):
    """sin μ = 1/M (μ = arcsin(1/M) なので sin(arcsin) を省く; numpy の版は sin(μ) を計算する)。"""
    return 1.0 / _nmax(M, 1.0 + 1e-12)


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
    cp = math.cos(ang_p)                         # tan は sin/cos から (cos は源項でも使う)
    cm = math.cos(ang_m)
    mp = math.sin(ang_p) / cp
    mm = math.sin(ang_m) / cm
    den = mm - mp
    good = abs(den) >= 1e-12
    parallel = math.isfinite(den) and not good
    if not good:
        den = 1.0
    xP = (Ar - Br + mm * Bx - mp * Ax) / den
    rP = Ar + mp * (xP - Ax)
    rPc = _nmax(rP, 0.0)
    sP = math.sin(thP)
    smP = _sin_mu(MP)
    fA = 0.5 * (fA0 + smP * _src_P(rPc, sP, Ar, bA, frac))
    fB = 0.5 * (fB0 + smP * _src_P(rPc, sP, Br, bB, frac))
    Sp = -delta * fA / cp * (xP - Ax)
    Sm = +delta * fB / cm * (xP - Bx)
    Jp = (Ath - Anu) + Sp
    Jm = (Bth + Bnu) + Sm
    return xP, rP, 0.5 * (Jp + Jm), 0.5 * (Jm - Jp), good, parallel


@njit(cache=True)
def _resid(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, fA0, fB0, bA, bB, xP, rP, thP, nuP, MP, delta, frac):
    """`_interior_residual` の 1 対分 (端点の項は `_pass` と同じ)。戻り: (幾何の残差, 適合式の残差)。"""
    muP = _mu(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    cp = math.cos(ang_p)                         # tan は sin/cos から (cos は源項でも使う)
    cm = math.cos(ang_m)
    mp = math.sin(ang_p) / cp
    mm = math.sin(ang_m) / cm
    rg = _nmax(abs(rP - (Ar + mp * (xP - Ax))), abs(rP - (Br + mm * (xP - Bx))))
    rPc = _nmax(rP, 0.0)
    sP = math.sin(thP)
    smP = _sin_mu(MP)
    fA = 0.5 * (fA0 + smP * _src_P(rPc, sP, Ar, bA, frac))
    fB = 0.5 * (fB0 + smP * _src_P(rPc, sP, Br, bB, frac))
    Sp = -delta * fA / cp * (xP - Ax)
    Sm = +delta * fB / cm * (xP - Bx)
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
    # 対ごとの前処理 (反復で変わらない量)。反復は「回」のループを外側、対のループを内側に回す (対どうしは独立なので、
    # 1 対の中の直列の依存の待ち時間を CPU が別の対の計算で埋められる)
    fin_in = np.empty(m, np.bool_)
    muA = np.empty(m)
    muB = np.empty(m)
    fA0 = np.empty(m)
    fB0 = np.empty(m)
    bA = np.empty(m)
    bB = np.empty(m)
    Mg = np.empty(m)                             # Newton の初期値 (以後は前の反復の M)
    good_all = np.ones(m, np.bool_)
    par = np.zeros(m, np.bool_)
    act = np.empty(m, np.int64)
    na = 0
    for j in range(m):
        Ar, Ath, AM = r[j + 1], th[j + 1], M[j + 1]
        Br, Bth, BM = r[j], th[j], M[j]
        fin_in[j] = (math.isfinite(x[j + 1]) and math.isfinite(Ar) and math.isfinite(Ath) and math.isfinite(nu[j + 1])
                     and math.isfinite(AM) and math.isfinite(x[j]) and math.isfinite(Br) and math.isfinite(Bth)
                     and math.isfinite(nu[j]) and math.isfinite(BM))
        thP[j] = 0.5 * (Ath + Bth)
        nuP[j] = 0.5 * (nu[j + 1] + nu[j])
        muA[j] = _mu(AM)
        muB[j] = _mu(BM)
        fA0[j] = _sin_mu(AM) * _src(Ar, Ath, Br, Bth, thrA[j], analytic, frac, eps)    # 端点の項
        fB0[j] = _sin_mu(BM) * _src(Br, Bth, Ar, Ath, thrB[j], analytic, frac, eps)
        bA[j] = _src(Ar, Ath, Ar, Ath, np.nan, False, frac, eps)                      # P が軸上のとき借りる sinθ/r
        bB[j] = _src(Br, Bth, Br, Bth, np.nan, False, frac, eps)
        Mg[j] = 0.5 * (AM + BM)
        if fin_in[j]:                            # converge で回す対 (入力が非有限の対は回さない; fixed2 は全部の対を回す)
            act[na] = j
            na += 1
        status[j] = PAIR_MISSING if mode == MODE_CONVERGE else PAIR_DONE
    if mode != MODE_CONVERGE:
        # 予測 1 + 修正 n_corr 回 (全部の対)
        for _ in range(1 + n_corr):
            for j in range(m):
                Mg[j] = _mach(nuP[j], Mg[j], gas_kind, g, tx, ty)
                xP[j], rP[j], thP[j], nuP[j], good, para = _pass(
                    x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j], nu[j], muA[j], muB[j],
                    fA0[j], fB0[j], bA[j], bB[j], thP[j], nuP[j], Mg[j], delta, frac)
                good_all[j] = good_all[j] and good
                par[j] = par[j] or para
        for j in range(m):
            fin_out = math.isfinite(xP[j]) and math.isfinite(rP[j]) and math.isfinite(thP[j]) and math.isfinite(nuP[j])
            okj = good_all[j] and fin_out
            st = PAIR_DONE
            if mode == MODE_LEGACY_FIXED2:       # moc_inverse._posthoc_stats (出力の有限性で分ける)
                if fin_in[j] and not okj:
                    st = PAIR_GEOM if fin_out else PAIR_NONFINITE
            else:                                # interior_vec の stats (fixed2)
                if fin_in[j] and ((not good_all[j]) or (not fin_out)):
                    st = PAIR_NONFINITE
                if fin_in[j] and par[j]:
                    st = PAIR_GEOM
            if not fin_in[j]:
                st = PAIR_MISSING
            status[j] = st
            ok[j] = okj
            n_it[j] = n_corr
    else:
        # converge: 修正 1 回の θ・ν の更新量が tol 以下になるまで (上限 max_corr 回)。止まった対は凍結
        for kk in range(max_corr + 1):           # kk = 0 予測子、kk ≥ 1 修正子
            if na == 0:
                break
            nn = 0
            for t in range(na):
                j = act[t]
                Mg[j] = _mach(nuP[j], Mg[j], gas_kind, g, tx, ty)
                a0, a1, a2, a3, good, para = _pass(x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j], nu[j],
                                                   muA[j], muB[j], fA0[j], fB0[j], bA[j], bB[j], thP[j], nuP[j], Mg[j],
                                                   delta, frac)
                fin = math.isfinite(a0) and math.isfinite(a1) and math.isfinite(a2) and math.isfinite(a3)
                upd = _nmax(abs(a2 - thP[j]), abs(a3 - nuP[j]))
                xP[j], rP[j], thP[j], nuP[j] = a0, a1, a2, a3
                n_it[j] = kk
                st = -1
                if para:
                    st = PAIR_GEOM
                elif (not good) or (not fin):
                    st = PAIR_NONFINITE
                elif kk >= 1 and upd <= tol:
                    st = PAIR_CONVERGED
                status[j] = st
                if st < 0:
                    act[nn] = j
                    nn += 1
            na = nn
        for t in range(na):
            status[act[t]] = PAIR_MAXITER
        for j in range(m):
            ok[j] = status[j] == PAIR_CONVERGED
    for j in range(m):
        if ok[j]:                                # 次の前線の M (と converge の最終残差)
            Mg[j] = _mach(nuP[j], Mg[j], gas_kind, g, tx, ty)
            MP[j] = Mg[j]
            if mode == MODE_CONVERGE:
                rg[j], rc[j] = _resid(x[j + 1], r[j + 1], th[j + 1], nu[j + 1], x[j], r[j], th[j], nu[j], muA[j], muB[j],
                                      fA0[j], fB0[j], bA[j], bB[j], xP[j], rP[j], thP[j], nuP[j], Mg[j], delta, frac)
        # 源項の分岐 (source_branches、最終状態の幾何で)
        Ar, Br = r[j + 1], r[j]
        rPc = _nmax(rP[j], 0.0)
        nA = Ar <= eps
        nB = Br <= eps
        hA = nA and analytic and math.isfinite(thrA[j])
        hB = nB and analytic and math.isfinite(thrB[j])
        br[0, j] = np.int8(hA) + np.int8(hB)
        br[1, j] = np.int8(nA and not hA and Br > eps) + np.int8(nB and not hB and Ar > eps)
        br[2, j] = np.int8(nA and not hA and not (Br > eps)) + np.int8(nB and not hB and not (Ar > eps))
        br[3, j] = (Ar > eps) and (Ar < frac * Br)
        br[4, j] = (Br > eps) and (Br < frac * Ar)
        br[5, j] = (rPc > eps) and (rPc < frac * Ar)
        br[6, j] = (rPc > eps) and (rPc < frac * Br)
        br[7, j] = rPc <= eps
    return n_lack, xP, rP, thP, nuP, ok, status, n_it, rg, rc, br, MP


# --- 網の検査の段ごとの走査 (`moc_diagnostics.characteristic_topology_diagnostics`、2026-10-15 plan §4.6 (2)) -------------------
@njit(cache=True)
def level_scan_nb(xs, rs, n_ax, seam_i_pad, seam_k_max, seam_x_max, has_sx, cap):
    r"""網の各段 (`levels[k, :, 0:2]`、NaN = 欠損) について、(1) 隣接する有限の点の間隔の最小とその位置、(2) 点の堆積 = 同じ段の
    **非隣接** (添字の差 > 2) の点の対で、距離 ≤ 0.2 × (その段の隣接間隔の中央値) (下限 1e-9) のものを数える。
    以前の書き方 (段ごとの cKDTree.query_pairs) と同じ判定 (距離² = Δx² + Δr² ≤ 閾値²)。対は x で並べた掃き出しで探す。
    seam の扱いも同じ: 対の片方の添字が n_ax ± seam_i_pad の中で段 ≤ seam_k_max、または両方が x < seam_x_max (has_sx のとき)。
    戻り: (最小間隔 (無ければ inf), その段, その添字, 堆積の対の数, seam を除いた数, seam を除いた対の位置 (先頭 cap 件、段ごとに
    (a, b) の昇順) [段, x_a, r_a])。"""
    n_lev, n_pt = xs.shape
    min_gap = np.inf
    gk = -1
    gi = -1
    n_pile = 0
    n_pile_int = 0
    locs = np.empty((max(cap, 1), 3))
    n_loc = 0
    idx = np.empty(n_pt, np.int64)
    px = np.empty(n_pt)
    pr = np.empty(n_pt)
    d = np.empty(max(n_pt - 1, 1))
    pa = np.empty(n_pt, np.int64)
    pb = np.empty(n_pt, np.int64)
    for k in range(n_lev):
        m = 0
        for i in range(n_pt):
            if math.isfinite(xs[k, i]) and math.isfinite(rs[k, i]):
                idx[m] = i
                px[m] = xs[k, i]
                pr[m] = rs[k, i]
                m += 1
        if m < 3:
            continue
        for q in range(m - 1):
            d[q] = math.hypot(px[q + 1] - px[q], pr[q + 1] - pr[q])
        jb = 0
        for q in range(1, m - 1):
            if d[q] < d[jb]:
                jb = q
        if d[jb] < min_gap:
            min_gap = d[jb]
            gk = k
            gi = idx[jb]
        if m < 5:
            continue
        thr = max(np.median(d[:m - 1]) * 0.2, 1e-9)
        thr2 = thr * thr
        inc = True                               # x が単調な段 (ほとんど) は並べ替えずに掃く
        dec = True
        for q in range(m - 1):
            if px[q + 1] < px[q]:
                inc = False
            if px[q + 1] > px[q]:
                dec = False
        if inc:
            order = np.arange(m)
        elif dec:
            order = np.arange(m - 1, -1, -1)
        else:
            order = np.argsort(px[:m])
        npair = 0
        for s in range(m):
            a = order[s]
            for u in range(s + 1, m):
                b = order[u]
                if px[b] - px[a] > thr:
                    break
                dx = px[a] - px[b]
                dr = pr[a] - pr[b]
                if dx * dx + dr * dr <= thr2 and abs(idx[a] - idx[b]) > 2:
                    if npair == pa.size:             # 対の置き場を広げる
                        pa2 = np.empty(2 * npair, np.int64)
                        pb2 = np.empty(2 * npair, np.int64)
                        pa2[:npair] = pa[:npair]
                        pb2[:npair] = pb[:npair]
                        pa = pa2
                        pb = pb2
                    pa[npair] = min(a, b)
                    pb[npair] = max(a, b)
                    npair += 1
        if npair == 0:
            continue
        key = np.empty(npair, np.int64)
        for q in range(npair):
            key[q] = pa[q] * m + pb[q]
        for q in np.argsort(key):
            a = pa[q]
            b = pb[q]
            n_pile += 1
            in_seam = ((abs(idx[a] - n_ax) <= seam_i_pad or abs(idx[b] - n_ax) <= seam_i_pad) and k <= seam_k_max)
            if has_sx and px[a] < seam_x_max and px[b] < seam_x_max:
                in_seam = True
            if not in_seam:
                n_pile_int += 1
                if n_loc < cap:
                    locs[n_loc, 0] = k
                    locs[n_loc, 1] = px[a]
                    locs[n_loc, 2] = pr[a]
                    n_loc += 1
    return min_gap, gk, gi, n_pile, n_pile_int, locs[:n_loc]


# --- C⁺ 流束閉包の壁 (`moc_inverse.cplus_flux_wall`、2026-10-15 plan §4.6 (3)) -------------------------------------------------
@njit(cache=True)
def _interp1(x, xp, fp):
    """`np.interp(x, xp, fp)` の 1 点 (xp は狭義単調増加; 範囲外は端の値、NaN は NaN、傾きの式は numpy の C 実装と同じ)。"""
    if x != x:
        return x
    n = xp.size
    if x < xp[0]:
        return fp[0]
    if x > xp[n - 1]:
        return fp[n - 1]
    if x == xp[n - 1]:
        return fp[n - 1]
    lo = 0
    hi = n - 1                                   # xp[lo] ≤ x < xp[hi]
    while hi - lo > 1:
        mid = (lo + hi) >> 1
        if xp[mid] <= x:
            lo = mid
        else:
            hi = mid
    if xp[lo] == x:
        return fp[lo]
    slope = (fp[lo + 1] - fp[lo]) / (xp[lo + 1] - xp[lo])
    v = slope * (x - xp[lo]) + fp[lo]
    if v != v:
        v = slope * (x - xp[lo + 1]) + fp[lo + 1]
        if v != v and fp[lo] == fp[lo + 1]:
            v = fp[lo]
    return v


@njit(cache=True)
def _flux_density(M, gas_kind, g, tx, ty):
    """`moc_inverse._mass_flux_density` の 1 点 (CPG の式、または表の線形補間)。"""
    if gas_kind == GAS_TABLE:
        return _interp1(M, tx, ty)
    t = 1.0 + 0.5 * (g - 1.0) * M * M
    return M * t ** (-0.5 * (g + 1.0) / (g - 1.0))


@njit(cache=True)
def cplus_flux_wall_nb(lx, lr, lth, lM, init_cum, mdot_star, gas_kind, g, tx, ty):
    r"""`cplus_flux_wall` の C⁺ 線ごとのループ (線 = $L_k[m-k]$、$m$ の降順、各線は最初の非有限の点の手前まで)。式と累積の順序は
    numpy の書き方と同じ。戻り: (n, 4) [x, r, θ, M] (線の順; x での並べ替えは呼び出し側)。"""
    n_lev, n = lx.shape
    out = np.empty((n, 4))
    no = 0
    for m in range(n - 1, -1, -1):
        L = 0                                    # 線の長さ (cplus_lines)
        while L <= m and L < n_lev and math.isfinite(lx[L, m - L]) and math.isfinite(lr[L, m - L]):
            L += 1
        if L < 2:
            continue
        c0 = init_cum[m]
        if c0 + 0.0 >= mdot_star:                # 起点が既に壁 (初期値線の壁足)
            out[no, 0] = lx[0, m]
            out[no, 1] = lr[0, m]
            out[no, 2] = lth[0, m]
            out[no, 3] = lM[0, m]
            no += 1
            continue
        cs = 0.0
        prev = c0 + 0.0
        for i in range(L - 1):
            a = m - i
            b = m - i - 1
            Mm = 0.5 * (lM[i + 1, b] + lM[i, a])
            thm = 0.5 * (lth[i + 1, b] + lth[i, a])
            rm = 0.5 * (lr[i + 1, b] + lr[i, a])
            d = (_flux_density(Mm, gas_kind, g, tx, ty)
                 * (math.cos(thm) * (lr[i + 1, b] - lr[i, a]) - math.sin(thm) * (lx[i + 1, b] - lx[i, a]))
                 * 2.0 * math.pi * rm)
            cs += d
            cur = c0 + cs
            if cur >= mdot_star:
                s = (mdot_star - prev) / max(cur - prev, 1e-30)
                out[no, 0] = lx[i, a] + s * (lx[i + 1, b] - lx[i, a])
                out[no, 1] = lr[i, a] + s * (lr[i + 1, b] - lr[i, a])
                out[no, 2] = lth[i, a] + s * (lth[i + 1, b] - lth[i, a])
                out[no, 3] = lM[i, a] + s * (lM[i + 1, b] - lM[i, a])
                no += 1
                break
            prev = cur
    return out[:no]
