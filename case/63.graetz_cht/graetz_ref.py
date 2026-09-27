#!/usr/bin/env python3
r"""Graetz 問題 (管内層流・十分発達した速度・壁温一様) の基準解。forge とは独立の実装。

無次元化: ξ = r/R, x⁺ = x/(D·Pe) (Pe = Re·Pr = U_m D/α), θ = (T − T_w)/(T_in − T_w)。
速度 u/U_m = 2(1 − ξ²)。局所 Nu_x = q_w D / (k (T_w − T_b))、T_b は混合平均 ∫u T r dr / ∫u r dr。

    --mode march   : 軸方向伝導なし (Pe → ∞)。x⁺ 方向 Crank–Nicolson マーチング。古典 Graetz 解。
    --mode ellip   : 軸方向伝導あり (有限 Pe)。上流の断熱区間 (u 発達済み) と下流の断熱区間を含む楕円型。
                     x は加熱開始点に集中させた不等間隔、対流は 2 次風上、伝導は 2 次中心差分。

自己検査 (--selftest): march で Nu_∞ = 3.6568 (第一固有値 λ₀² /2)、Shah の相関式 (±3 %) との差、格子収束。
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

NU_INF = 3.656793  # Shah & London (1978) の UWT 十分発達値


# Graetz 級数 (UWT)。係数は文献値を写さず、march と独立に**固有値問題を直接解いて**作る
# (Chebyshev でなく細かい FV の一般化固有値問題 −4∇²φ = λ̃ u φ、nr=4000)。Nu_x = Σ G_n e^{−λ̃_n x⁺} / Σ (G_n/λ̃_n) e^{−λ̃_n x⁺}
# の形 (G_n は入口 θ=1 の展開係数 × 壁勾配)。第一固有値から Nu_∞ = λ̃₀/4·... を文献 3.656793 と照合する。
def graetz_series(xp, nr=4000, nterm=60):
    xi, h, vol, L = radial_ops(nr)
    u = 2.0 * (1.0 - xi ** 2)
    m = nr
    A = (-4.0 * L[:m, :m]).tocsc()
    B = sp.diags(u[:m] * vol[:m]).tocsc()
    w, V = spla.eigsh(A, k=nterm, M=B, sigma=0.0, which="LM")
    o = np.argsort(w); w, V = w[o], V[:, o]
    # 正規化 V^T B V = I。入口 θ=1 の係数 c_n = V_n^T B 1
    c = V.T @ (B @ np.ones(m))
    full = np.vstack([V, np.zeros((1, nterm))])
    g = (3 * full[-1] - 4 * full[-2] + full[-3]) / (2 * h)       # 各モードの ∂φ/∂ξ|_w
    ub = (u * vol).sum()
    tbn = (full.T * u * vol).sum(axis=1) / ub                    # 各モードの混合平均
    xp = np.atleast_1d(np.asarray(xp, float))
    E = np.exp(-np.outer(xp, w))
    return -2.0 * (E @ (c * g)) / (E @ (c * tbn)), w[0]


def shah_nu_x(xp):
    """Shah (1975) の局所 Nu 相関 (UWT、±3 % 程度)。自己検査の粗い照合用。"""
    xp = np.asarray(xp, float)
    a = 1.077 * xp ** (-1.0 / 3.0) - 0.7
    b = 3.657 + 6.874 * (1e3 * xp) ** (-0.488) * np.exp(-57.2 * xp)
    return np.where(xp <= 0.01, a, b)


def radial_ops(nr):
    """ξ ∈ [0,1] の節点 (nr+1 点、一様)。有限体積 (軸は半 CV) で ∇²θ = (1/ξ)∂(ξ ∂θ/∂ξ) の行列と CV 体積 ∫ξ dξ。"""
    xi = np.linspace(0.0, 1.0, nr + 1)
    h = xi[1] - xi[0]
    n = nr + 1
    # 面 ξ_{i±1/2}
    vol = np.zeros(n)
    L = sp.lil_matrix((n, n))
    for i in range(n):
        lo = max(xi[i] - h / 2, 0.0)
        hi = min(xi[i] + h / 2, 1.0)
        vol[i] = 0.5 * (hi ** 2 - lo ** 2)
        if i > 0:
            L[i, i - 1] += lo / h
            L[i, i] -= lo / h
        if i < n - 1:
            L[i, i + 1] += hi / h
            L[i, i] -= hi / h
    return xi, h, vol, L.tocsr()


def mixing_cup(theta, u, vol):
    return (theta * u * vol).sum(axis=-1) / (u * vol).sum()


def wall_grad(theta, h):
    """∂θ/∂ξ at ξ=1 (壁 = 最後の節点、θ_w=0)、2 次片側差分。"""
    return (3 * theta[..., -1] - 4 * theta[..., -2] + theta[..., -3]) / (2 * h)


def solve_march(nr, xmax, nx, return_tb=False):
    """Pe→∞。無次元式 u (∂θ/∂x⁺) = 4 ∇²_ξ θ  (x⁺ = x/(D Pe), ξ = r/R: α ∂²/∂r² → 4/(D²) ...)。"""
    xi, h, vol, L = radial_ops(nr)
    u = 2.0 * (1.0 - xi ** 2)
    # 壁 (最後の節点) は Dirichlet θ=0 なので自由度から外す
    m = nr
    Lf = L[:m, :m]
    M = sp.diags(u[:m] * vol[:m])
    # x⁺ は入口近傍を細かく: 対数刻み
    xs = np.concatenate([[0.0], np.geomspace(1e-7, xmax, nx)])
    th = np.ones(m)
    out_nu, out_x, out_tb = [], [], []
    for k in range(1, len(xs)):
        dx = xs[k] - xs[k - 1]
        A = (M - 0.5 * dx * 4.0 * Lf).tocsc()
        B = (M + 0.5 * dx * 4.0 * Lf)
        # 初手は後退 Euler (不連続の初期条件での CN の振動を避ける)
        if k <= 20:
            A = (M - dx * 4.0 * Lf).tocsc()
            th = spla.spsolve(A, M @ th)
        else:
            th = spla.spsolve(A, B @ th)
        full = np.append(th, 0.0)
        tb = mixing_cup(full, u, vol)
        g = wall_grad(full, h)          # ∂θ/∂ξ
        nu = -2.0 * g / tb              # Nu = q D/(k ΔT_b) = −(∂θ/∂ξ)·(D/R)/θ_b
        out_x.append(xs[k]); out_nu.append(nu); out_tb.append(tb)
    if return_tb:
        return np.array(out_x), np.array(out_nu), np.array(out_tb)
    return np.array(out_x), np.array(out_nu)


def solve_ellip(nr, pe, xu, xh, xd, nx_h, nx_u, nx_d):
    """有限 Pe。x⁺ ∈ [−xu, xh+xd]、加熱区間 [0, xh] は θ_w=0、それ以外の壁は断熱。入口 θ=1、出口は ∂θ/∂x=0。

    無次元式: u ∂θ/∂x⁺ = 4 ∇²_ξ θ + (1/Pe²) ∂²θ/∂x⁺²   (x⁺ = x/(D Pe), ξ = r/R)。
    """
    xi, h, vol, L = radial_ops(nr)
    u = 2.0 * (1.0 - xi ** 2)
    n_r = nr + 1

    # 点の配置: 加熱開始点 x⁺=0 から両側へ対数刻み (march と同じ考え方。等比を固定すると n を増やしても
    # 細分化にならない — 2026-09-27 の初版はこれで nx 1200 が壊れた)。
    x0 = 1e-6
    xh_pts = np.concatenate([[0.0], np.geomspace(x0, xh, nx_h)])
    xu_pts = -np.geomspace(x0, xu, nx_u)[::-1]
    xd_pts = xh + np.geomspace(x0 * 10, xd, nx_d)
    x = np.concatenate([xu_pts, xh_pts, xd_pts])
    n_x = len(x)
    N = n_x * n_r
    idx = lambda i, j: i * n_r + j
    rows, cols, vals = [], [], []
    rhs = np.zeros(N)
    eps = 1.0 / pe ** 2
    heated = (x >= -1e-15) & (x <= xh + 1e-15)
    for i in range(n_x):
        for j in range(n_r):
            p = idx(i, j)
            if i == 0:                                    # 入口 θ=1
                rows.append(p); cols.append(p); vals.append(1.0); rhs[p] = 1.0; continue
            if j == nr and heated[i]:                     # 加熱壁 θ=0
                rows.append(p); cols.append(p); vals.append(1.0); continue
            # 半径方向 (FV、断熱壁は L の境界で流束 0)
            a = {}
            for jj, v in zip(L.indices[L.indptr[j]:L.indptr[j + 1]], L.data[L.indptr[j]:L.indptr[j + 1]]):
                a[idx(i, jj)] = a.get(idx(i, jj), 0.0) - 4.0 * v / vol[j]
            # 軸方向
            if i == n_x - 1:                              # 出口: θ_i = θ_{i-1} (零勾配)
                rows += [p, p]; cols += [p, idx(i - 1, j)]; vals += [1.0, -1.0]; continue
            dxm, dxp = x[i] - x[i - 1], x[i + 1] - x[i]
            # 2 次中心の 2 階微分 (非等間隔)
            c_m = 2.0 / (dxm * (dxm + dxp)); c_p = 2.0 / (dxp * (dxm + dxp))
            a[idx(i - 1, j)] = a.get(idx(i - 1, j), 0.0) - eps * c_m
            a[idx(i + 1, j)] = a.get(idx(i + 1, j), 0.0) - eps * c_p
            a[p] = a.get(p, 0.0) + eps * (c_m + c_p)
            # 対流 u ∂θ/∂x: 2 次風上 (非等間隔 3 点後退)、i==1 は 1 次
            if i >= 2:
                d1, d2 = x[i] - x[i - 1], x[i] - x[i - 2]
                w0 = (d1 + d2) / (d1 * d2); w1 = -d2 / (d1 * (d2 - d1)); w2 = d1 / (d2 * (d2 - d1))
                for ii, w in ((i, w0), (i - 1, w1), (i - 2, w2)):
                    a[idx(ii, j)] = a.get(idx(ii, j), 0.0) + u[j] * w
            else:
                a[p] = a.get(p, 0.0) + u[j] / dxm
                a[idx(i - 1, j)] = a.get(idx(i - 1, j), 0.0) - u[j] / dxm
            for c, v in a.items():
                rows.append(p); cols.append(c); vals.append(v)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    th = spla.spsolve(A.tocsc(), rhs).reshape(n_x, n_r)
    tb = mixing_cup(th, u, vol)
    g = wall_grad(th, h)
    nu = np.where(heated, -2.0 * g / tb, np.nan)
    return x, nu, tb


def selftest():
    ok = True
    print("== march (Pe→∞) の格子収束と Nu_∞")
    res = {}
    for nr in (100, 200, 400):
        xs, nu = solve_march(nr, 0.3, 3000)
        res[nr] = (xs, nu)
        print(f"  nr={nr}: Nu(x⁺=0.3) = {nu[-1]:.5f}")
    nu_inf = res[400][1][-1]
    e = abs(nu_inf - NU_INF) / NU_INF
    print(f"  |Nu_∞ − 3.656793| / 3.656793 = {e:.2e}  (許容 1e-3)  {'PASS' if e < 1e-3 else 'FAIL'}")
    ok &= e < 1e-3
    xq = np.array([1e-3, 3e-3, 1e-2, 3e-2, 0.1])
    n2 = np.interp(np.log(xq), np.log(res[200][0]), res[200][1])
    n4 = np.interp(np.log(xq), np.log(res[400][0]), res[400][1])
    gs, lam0 = graetz_series(xq)
    print(f"  固有値問題: 第一固有値から Nu_∞ = λ̃₀/2·(…) — 級数の x⁺=10 の値 {graetz_series([10.0])[0][0]:.6f}")
    print("  x⁺        Nu(march nr400)  nr200→400 [%]  級数 (固有値 60 項, x⁺≥3e-3) との差 [%]")
    for a, b, c, d in zip(xq, n4, n2, gs):
        tol_ok = abs(b - c) / b < 1e-3 and (a < 3e-3 or abs(b - d) / d < 1e-3)
        print(f"  {a:8.1e}  {b:12.5f}  {abs(b - c) / b * 100:12.4f}      {(b - d) / d * 100:+9.4f}  {'ok' if tol_ok else 'NG'}")
        ok &= tol_ok
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--mode", choices=["march", "ellip"], default="march")
    ap.add_argument("--nr", type=int, default=200)
    ap.add_argument("--pe", type=float, default=720.0)
    ap.add_argument("--xu", type=float, default=0.01, help="上流断熱区間 [x⁺]")
    ap.add_argument("--xh", type=float, default=0.12, help="加熱区間 [x⁺]")
    ap.add_argument("--xd", type=float, default=0.01, help="下流断熱区間 [x⁺]")
    ap.add_argument("--nx", type=int, default=600)
    ap.add_argument("--out", default=None, help="CSV (x⁺, Nu_x, θ_b)")
    a = ap.parse_args()
    if a.selftest:
        ok = selftest()
        print("VERDICT:", "PASS" if ok else "FAIL")
        return 0 if ok else 1
    if a.mode == "march":
        xs, nu = solve_march(a.nr, a.xh, a.nx)
        tb = np.full_like(xs, np.nan)
    else:
        xs, nu, tb = solve_ellip(a.nr, a.pe, a.xu, a.xh, a.xd, a.nx, a.nx // 6, a.nx // 6)
    if a.out:
        np.savetxt(a.out, np.c_[xs, nu, tb], delimiter=",", header="xplus,Nu_x,theta_b", comments="")
    for xq in (1e-3, 3e-3, 1e-2, 3e-2, 0.1):
        m = np.isfinite(nu) & (xs > 0)
        print(f"x⁺={xq:7.1e}  Nu={np.interp(xq, xs[m], nu[m]):.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
