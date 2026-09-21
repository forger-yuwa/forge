#!/usr/bin/env python3
r"""Langtry–Menter 2009 ($\gamma$–$Re_{\theta t}$) のソース項の**参照実装** (numpy、倍精度)。

SU2 8.x の `trans_sources.hpp` / `trans_correlations.hpp` (`MENTER_LANGTRY`) / `CTransLMSolver.cpp` を式の順序・下限ごと写したもの。
forge のカーネル `transitionSource_d.cu` の単体検査に使う (methods/turbulence/theory.md §11、plan turbulence-transition-lm2009 §6)。
ここを直すときは SU2 の該当行を必ず見ること。

usage (自己検査): python3 solver_density_cuda/tools/lm2009_reference.py
"""
import numpy as np

C_E1, C_A1, C_E2, C_A2, C_THETA = 1.0, 2.0, 50.0, 0.06, 0.03


def re_theta_c(re_t):
    re_t = np.asarray(re_t, float)
    lo = (-396.035e-2 + 10120.656e-4 * re_t - 868.230e-6 * re_t ** 2 + 696.506e-9 * re_t ** 3 - 174.105e-12 * re_t ** 4)
    hi = re_t - (593.11 + 0.482 * (re_t - 1870.0))
    return np.where(re_t <= 1870.0, lo, hi)


def f_length1(re_t):
    re_t = np.asarray(re_t, float)
    return np.select([re_t < 400, re_t < 596, re_t < 1200],
                     [39.8189 - 119.270e-4 * re_t - 132.567e-6 * re_t ** 2,
                      263.404 - 123.939e-2 * re_t + 194.548e-5 * re_t ** 2 - 101.695e-8 * re_t ** 3,
                      0.5 - 3.0e-4 * (re_t - 596.0)], 0.3188)


def re_theta_t_freestream(tu):
    """入口値: $\\lambda_\\theta=0$ の相関。tu は [%]。"""
    tu = np.asarray(tu, float)
    return np.where(tu <= 1.3, 1173.51 - 589.428 * tu + 0.2196 / np.maximum(tu, 0.027) ** 2, 331.5 * np.maximum(tu - 0.5658, 1e-12) ** -0.671)


def re_theta_t_corr(tu, ro, mu, U, duds, iters=100):
    """$Re_{\\theta t}(Tu,\\lambda_\\theta)$。$\\theta$ が $Re_{\\theta t}$ に依るので不動点反復 (SU2 と同じ初期値 20)。"""
    ret = np.full(np.shape(tu), 20.0)
    for _ in range(iters):
        th = ret * mu / ro / U
        lam = np.clip(ro * th * th / mu * duds, -0.1, 0.1)
        fl = np.where(lam <= 0.0,
                      1.0 - (-12.986 * lam - 123.66 * lam ** 2 - 405.689 * lam ** 3) * np.exp(-(tu / 1.5) ** 1.5),
                      1.0 + 0.275 * (1.0 - np.exp(-35.0 * lam)) * np.exp(-tu / 0.5))
        new = np.where(tu <= 1.3, fl * (1173.51 - 589.428 * tu + 0.2196 / tu ** 2), 331.5 * fl * np.maximum(tu - 0.5658, 1e-12) ** -0.671)
        new = np.maximum(new, 20.0)
        if np.all(np.abs(new - ret) <= 1e-7 * np.abs(ret)):
            ret = new; break
        ret = new
    return ret


def sources(ro, U, G, k, om, mu, dist, gam, ret):
    """体積あたりのソースと診断を返す。U[...,3], G[...,i,j]=dUi/dxj。戻り値は dict。"""
    ro, k, om, mu, dist, gam, ret = (np.asarray(a, float) for a in (ro, k, om, mu, dist, gam, ret))
    Umag = np.maximum(np.linalg.norm(U, axis=-1), 1e-30)
    Sij = 0.5 * (G + np.swapaxes(G, -1, -2))
    S = np.sqrt(2.0 * np.sum(Sij * Sij, axis=(-1, -2)))
    W = np.stack([G[..., 2, 1] - G[..., 1, 2], G[..., 0, 2] - G[..., 2, 0], G[..., 1, 0] - G[..., 0, 1]], -1)
    Om = np.linalg.norm(W, axis=-1)
    tu = np.maximum(100.0 * np.sqrt(2.0 * k / 3.0) / Umag, 0.027)
    rec = re_theta_c(ret)
    r_omega = ro * dist ** 2 * om / mu
    f_sub = np.exp(-(r_omega / 200.0) ** 2)
    f_len = f_length1(ret) * (1 - f_sub) + 40.0 * f_sub
    r_t = ro * k / mu / om
    re_v = ro * dist ** 2 * S / mu
    f1 = re_v / (2.193 * rec)
    f_onset = np.maximum(np.minimum(np.maximum(f1, f1 ** 4), 2.0) - np.maximum(1.0 - (r_t / 2.5) ** 3, 0.0), 0.0)
    gradU = np.einsum("...i,...ij->...j", U, G) / Umag[..., None]           # d|U|/dx_j
    duds = np.sum(U * gradU, -1) / Umag
    tscale = 500.0 * mu / ro / Umag ** 2
    delta = 50.0 * Om * dist / Umag * (7.5 * ret * mu / ro / Umag) + 1e-20
    f_wake = np.exp(-(r_omega / 1.0e5) ** 2)
    var1 = (gam - 1.0 / C_E2) / (1.0 - 1.0 / C_E2)
    f_theta = np.minimum(np.maximum(f_wake * np.exp(-(dist / delta) ** 4), 1.0 - var1 ** 2), 1.0)
    f_turb = np.exp(-(r_t / 4.0) ** 4)
    ret_corr = re_theta_t_corr(tu, ro, mu, Umag, duds)
    Pg = f_len * C_A1 * ro * S * np.sqrt(np.maximum(f_onset * gam, 0.0)) * (1.0 - C_E1 * gam)
    Dg = C_A2 * ro * Om * gam * f_turb * (C_E2 * gam - 1.0)
    Pt = C_THETA * ro / tscale * (ret_corr - ret) * (1.0 - f_theta)
    gsep = np.minimum(2.0 * np.maximum(re_v / (3.235 * rec) - 1.0, 0.0) * np.exp(-(r_t / 20.0) ** 4), 2.0) * f_theta
    gsep = np.clip(gsep, 0.0, 2.0)
    # forge の陰的対角 (項ごとの負の部分; methods/turbulence/implementation.md)。SU2 はソース微分そのもの (su2_jac_gamma) を足す。
    jac_g = 1.5 * C_E1 * f_len * C_A1 * S * np.sqrt(np.maximum(f_onset * gam, 0.0)) + np.maximum(C_A2 * Om * f_turb * (2.0 * C_E2 * gam - 1.0), 0.0)
    jac_t = C_THETA / tscale * (1.0 - f_theta)
    with np.errstate(divide='ignore', invalid='ignore'):
        su2_jg = f_len * C_A1 * S * np.sqrt(f_onset) * (0.5 * gam ** -0.5 - 1.5 * C_E1 * gam ** 0.5) - C_A2 * Om * f_turb * (2.0 * C_E2 * gam - 1.0)
    wall = dist <= 1e-10
    z = lambda a: np.where(wall, 0.0, a)
    return dict(src_gamma=z(Pg - Dg), src_reth=z(Pt), gamma_eff=np.maximum(gam, gsep), gamma_sep=gsep, Tu=tu, f_onset=f_onset,
                jac_gamma=z(jac_g), jac_reth=z(jac_t), su2_jac_gamma=z(su2_jg), P_gamma=z(Pg), E_gamma=z(Dg),
                f_length=f_len, f_theta=f_theta, re_theta_c=rec, re_theta_corr=ret_corr, re_v=re_v, r_t=r_t)


def _selftest():
    ok = True
    # 1) 相関の連続性 (区分の継ぎ目)
    for x in (400.0, 596.0, 1200.0):
        a, b = f_length1(x - 1e-6), f_length1(x + 1e-6)
        print(f"  F_length1 at Re_t={x:6.0f}: {a:.5f} / {b:.5f}"); ok &= abs(a - b) < 2e-2 * max(abs(a), 1e-3) + 1e-3
    a, b = re_theta_c(1870.0 - 1e-6), re_theta_c(1870.0 + 1e-6); print(f"  Re_theta_c at 1870: {a:.3f} / {b:.3f}"); ok &= abs(a - b) < 2.0      # 文献の相関自体に 1.2 の段差がある
    # 2) 文献の代表値: Tu=3.3 % (T3A), 6.5 % (T3B), 0.9 % (T3A-) の自由流 Re_theta_t
    for tu, ref in ((3.3, 168.0), (6.5, 100.0), (0.9, 643.0)):
        v = float(re_theta_t_freestream(tu)); print(f"  Re_theta_t(Tu={tu} %) = {v:.1f}  (目安 {ref:.0f})"); ok &= abs(v - ref) / ref < 0.08
    # 3) 平板の層流境界層の 1 点 (Blasius 的な値): 自由流では F_onset=0, 層の中で Re_v が Re_theta_c の 2.193 倍を超えると立ち上がる
    U = np.array([[50.0, 0.0, 0.0]]); G = np.zeros((1, 3, 3)); G[0, 0, 1] = 2.0e4
    r = sources(1.2, U, G, np.array([0.5]), np.array([2000.0]), 1.8e-5, np.array([5e-4]), np.array([0.02]), np.array([160.0]))
    print("  one-point:", {k_: float(np.ravel(v)[0]) for k_, v in r.items()})
    ok &= r["src_gamma"][0] > 0 and r["f_onset"][0] > 0
    print("SELFTEST:", "PASS" if ok else "FAIL"); return ok


if __name__ == "__main__":
    import sys; sys.exit(0 if _selftest() else 1)
