"""V4 壁 (plan verification-m6-axis-wave-mesh-su2 §5.1 #16, diagnostician 2026-10-05): 位置+壁角の同時当てはめで、始点は r(0)=1 だけ固定し
r′(0)・r″(0) は MOC 点と λ に任せる (x=0 の θ は Hall 足の入力なので当てはめに使わない)。縮流部の quintic Hermite を (r=1, r′(0), r″(0)) に接続して壁全体を C² にする。
生産コードへの組み込み前の case 内実装 (design_chain の壁を差し替える)。
"""
import copy
import numpy as np
from scipy.interpolate import BSpline
from forge_design.geometry.wall import _hermite_quintic, _poly_eval

H0, H1, XG, SR, ST = 0.0125, 0.5, 6.0, 1e-6, 1e-4


def fit_v4(tb, lam=1e-9, k=5):
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]; x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (XG - x0), 1.0); xs.append(xs[-1] + H0 + (H1 - H0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * H1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1; E = np.eye(nc)
    D = lambda xq, d: np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1); w = np.gradient(x); w = w / w.mean()
    wt = w.copy(); wt[0] = 0.0                       # x=0 の θ (Hall 足) は使わない
    xq = np.linspace(x0, xe, 8000); B3 = D(xq, 3)
    A = (B0.T * (w / SR ** 2)) @ B0 + (B1.T * (wt / ST ** 2)) @ B1 + lam * (B3.T * np.gradient(xq)) @ B3 / SR ** 2
    b = B0.T @ (w * r / SR ** 2) + B1.T @ (wt * np.tan(th) / ST ** 2)
    Cm = np.array([D(np.r_[x0], 0)[0], D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]]); dv = np.array([r[0], r[-1], np.tan(th[-1])])
    c = np.linalg.solve(np.block([[A, Cm.T], [Cm, np.zeros((3, 3))]]), np.r_[b, dv])[:nc]
    return BSpline(t, c, k), nc


class UpstreamSlope:
    """縮流部 quintic Hermite: (r_U, 0, 0) @ x=−L_U → (1, s0, k0) @ x=0。UpstreamThroatPoly の端傾きを 0 以外にした版。"""

    def __init__(self, r_U, L_U, s0, k0):
        self.r_U, self.L_U, self.s0, self.k0 = float(r_U), float(L_U), float(s0), float(k0)
        self._c = _hermite_quintic(-self.L_U, 0.0, self.r_U, 0.0, 0.0, 1.0, self.s0, self.k0)

    @property
    def mu(self):
        return self.L_U ** 2 * self.k0 / (self.r_U - 1.0)

    def r(self, x, deriv=0):
        return _poly_eval(self._c, -self.L_U, 0.0, np.asarray(x, dtype=float), deriv)

    def validate(self, n=4001):
        return []


def build_v4(d, lam=1e-9):
    """design_chain の結果 d から V4 の設計壁 (AxisMachCFDWall 互換) を作る。"""
    s, nc = fit_v4(d["wall_inv"], lam)
    w = copy.copy(d["wall"]); w._spl = s
    up0 = d["wall"].up
    w.up = UpstreamSlope(up0.r_U, up0.L_U, float(s(np.r_[0.0], 1)[0]), float(s(np.r_[0.0], 2)[0]))
    w.fit_info = {"kind": "joint_fit_v4", "lam": lam, "h0": H0, "h1": H1, "n_cp": nc,
                  "r1_0": float(s(np.r_[0.0], 1)[0]), "r2_0": float(s(np.r_[0.0], 2)[0])}
    return w
