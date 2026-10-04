#!/usr/bin/env python3
"""平板境界層の共通抽出・壁法則ルーチン。

`check_quasisteady.py` (準定常判定) と `case/26.flat_plate_sst/tools/cf_retheta_analysis.py`
(正式後処理) が**同じ式**を使うための共通モジュール。

かつて両者は式を複製していて、fit 窓 (±0.08 vs 実質 ±0.12)・$U_e,\\rho_e,\\mu_e$ の取り方・
音速の評価が食い違い、**同じ run に対して別の値**を報告していた。数値を突き合わせるとき
「どちらのツールで出したか」を気にしなければならない状態は監査上あってはならないので、
実装をここへ一本化する。

含むもの:
  - `stations` / `column` / `profile` … 壁点を必ず含むカラムと $\\theta,\\delta^*$
  - `cf_karman_schoenherr` … 外部相関
  - `uplus` / `solve_utau` … 壁法則 (reichardt / spalding, $\\kappa$=0.41, $B$=5.0)
  - `cf_momentum` … 運動量積分 $C_f$ (壁出力非依存の収支診断)
"""
import numpy as np


def _trapz(y, x):
    """台形積分。**numpy 2 で `np.trapz` が削除された**ので互換に包む (AWS の numpy で落ちた)。"""
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return f(y, x)

KAPPA = 0.41
B_LOG = 5.0
GAMMA = 1.4


# ---------------------------------------------------------------- 壁法則
def _reichardt_uplus(yp):
    return np.log(1.0 + KAPPA * yp) / KAPPA + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


def _spalding_uplus(yp):
    """Spalding は陰形式 $y^+=H(u^+)$ なので $u^+$ を Newton で求める。"""
    up = np.log(max(yp, 1e-12)) / KAPPA + B_LOG
    for _ in range(200):
        e = KAPPA * up
        f = up + np.exp(-KAPPA * B_LOG) * (np.exp(e) - 1 - e - e**2 / 2 - e**3 / 6) - yp
        df = 1 + np.exp(-KAPPA * B_LOG) * KAPPA * (np.exp(e) - 1 - e - e**2 / 2)
        s = f / df
        up -= s
        if abs(s) < 1e-13:
            break
    return up


def uplus(yp, law="reichardt"):
    yp = max(float(yp), 1e-12)
    if law == "spalding":
        return _spalding_uplus(yp)
    if law == "log":
        return np.log(yp) / KAPPA + B_LOG
    if law == "none":
        return float("nan")
    if law == "reichardt":
        return _reichardt_uplus(yp)
    raise ValueError(f"未知の壁法則: {law}")


def solve_utau(Ut, y, nu, law="reichardt"):
    """$U_t/u_\\tau = u^+(u_\\tau y/\\nu)$ の Newton 逆解き。

    **ソルバが使った法則と揃えること** — 揃えないと A/B の符号すら逆に見える。
    """
    if law == "none" or not np.isfinite(Ut) or Ut <= 0:
        return float("nan")
    ut = np.sqrt(max(nu * Ut / max(y, 1e-30), 1e-30))
    for _ in range(200):
        ut = max(ut, 1e-12)
        f = Ut / ut - uplus(ut * y / nu, law)
        d = 1e-6 * ut
        df = ((Ut / (ut + d) - uplus((ut + d) * y / nu, law)) - f) / d
        if abs(df) < 1e-30:
            break
        s = f / df
        ut -= s
        if abs(s) < 1e-14 * max(ut, 1e-12):
            break
    return max(ut, 0.0)


def cf_karman_schoenherr(re_theta):
    """NASA TMR flat-plate validation と同じ Karman-Schoenherr 相関。"""
    if not np.isfinite(re_theta) or re_theta <= 1.0:
        return float("nan")
    l = np.log10(re_theta)
    return 1.0 / (17.08 * l * l + 25.11 * l + 6.012)


# ---------------------------------------------------------------- 抽出
def stations(C, xmin, xmax):
    x = C[:, 0]
    xa = np.unique(np.round(x[x > 1e-6], 6))
    return xa[(xa >= xmin) & (xa <= xmax)]


def column(C, xc):
    x, y = C[:, 0], C[:, 1]
    i = np.sort(np.where(np.abs(x - xc) < 1e-4)[0])
    return i[np.argsort(y[i])]


def profile(C, V, xc, ytop=None, edge="local", inf=None):
    """壁点を必ず含むカラム。cell は $(y,u)=(0,0)$ を補い、node は既にある壁点を使う。

    `V` は少なくとも `u`, `ro`, `mu` を持つ dict。`P` があれば端の音速も返す。
    `edge="global"` のときだけ `inf`=(U_INF, RHO_INF, MU) を使う。
    """
    col = column(C, xc)
    if len(col) < 5:
        return None
    y = C[col, 1].astype(float)
    u = np.asarray(V["u"])[col].astype(float)
    ro = np.asarray(V["ro"])[col].astype(float)
    mu = np.asarray(V["mu"])[col].astype(float)
    P = np.asarray(V["P"])[col].astype(float) if V.get("P") is not None else None
    if y[0] > 1e-12:                               # cell: 壁点を補う
        y = np.r_[0.0, y]; u = np.r_[0.0, u]; ro = np.r_[ro[0], ro]; mu = np.r_[mu[0], mu]
        if P is not None:
            P = np.r_[P[0], P]
    m = np.ones(len(y), bool) if ytop is None else (y <= ytop)
    j = int(np.argmax(u[m]))
    if edge == "global" and inf is not None:
        ue, roe, mue = inf
    else:
        ue, roe, mue = u[m][j], ro[m][j], mu[m][j]
    theta = float(_trapz((ro[m] / roe) * (u[m] / ue) * (1.0 - u[m] / ue), y[m]))
    dstar = float(_trapz(1.0 - (ro[m] * u[m]) / (roe * ue), y[m]))
    ae = float(np.sqrt(GAMMA * P[m][j] / ro[m][j])) if P is not None else float("nan")
    return dict(y=y, u=u, ro=ro, mu=mu, ue=ue, roe=roe, mue=mue, ae=ae,
                theta=theta, dstar=dstar, col=col, y1=y[1], u1=u[1], nu1=mu[1] / ro[1])


def cf_momentum(C, V, xs, xmin, xmax, ytop=None, edge="local", inf=None,
                fit_window=0.08, fit_order=2, ts=None):
    """運動量積分 $C_f$ (壁出力非依存の収支診断) と $Re_\\theta$ を返す。

      $C_f/2 = d\\theta/dx + \\theta (H + 2 - M_e^2)/U_e \\cdot dU_e/dx$

    **「定義非依存の壁摩擦」ではない** — 定常性・境界層近似・$U_e/\\rho_e$ の取り方・
    積分上端・微分 fit に依存する。`ts` を渡すと $M_e$ をその静温から評価し、
    渡さなければ場の $\\sqrt{\\gamma P/\\rho}$ を使う。
    """
    xa = stations(C, xmin, xmax)
    if len(xa) < fit_order + 2:
        return float("nan"), float("nan")
    prof = {}
    for xx in xa:
        p = profile(C, V, xx, ytop, edge, inf)
        if p is not None and np.isfinite(p["theta"]):
            prof[xx] = p
    xa = np.array([xx for xx in xa if xx in prof])
    if len(xa) < fit_order + 2:
        return float("nan"), float("nan")
    th = np.array([prof[xx]["theta"] for xx in xa])
    ue = np.array([prof[xx]["ue"] for xx in xa])
    xc = xa[int(np.argmin(np.abs(xa - xs)))]
    w = np.abs(xa - xc) <= fit_window
    if w.sum() < fit_order + 2:
        w = np.abs(xa - xc) <= 2 * fit_window
    if w.sum() < fit_order + 2:
        return float("nan"), float("nan")
    dthdx = np.polyval(np.polyder(np.polyfit(xa[w], th[w], fit_order)), xc)
    duedx = np.polyval(np.polyder(np.polyfit(xa[w], ue[w], fit_order)), xc)
    p = prof[xc]
    ae = np.sqrt(GAMMA * 287.0 * ts) if ts is not None else p["ae"]
    Me = p["ue"] / ae if np.isfinite(ae) and ae > 0 else 0.0
    H = p["dstar"] / p["theta"]
    pg = p["theta"] * (H + 2.0 - Me * Me) / p["ue"] * duedx
    ret = p["roe"] * p["ue"] * p["theta"] / p["mue"]
    return float(2.0 * (dthdx + pg)), float(ret)
