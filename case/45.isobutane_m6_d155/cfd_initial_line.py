"""CFD ピン (plan tooling-nozzle-cfd-pinned-initial-line §4.1) の試作: node Euler 場からスロート特性線 (MOC の初期線) を抽出する provider。
`HallThroat` 互換の throat_characteristic(n) / mach(x, r) / axis_anchor(x)。生産コード (design/forge_design/feedback/cfd_initial_line.py) へ移す前の case 内実装。

- 場: 構造格子 (ni 断面 × nj、断面は x 一定) の M, θ を (x, η = r/r_w(x)) 平面の 3 次スプライン (RectBivariateSpline) で補間。末尾 n_tail スナップショット平均。
- C⁻: 壁の x=0 の点から dr/dx = tan(θ−μ) を RK4 で軸まで。r 等間隔 n 点に再標本化 (軸→壁)。壁足の θ は 0 (設計スロートの壁接線)。
- 軸: η=0 の断面値をそのまま使う (スプラインは η 方向に軸の値を持つ)。axis_anchor は軸の M(x) の平滑化スプラインの微分。
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import RectBivariateSpline, make_smoothing_spline, CubicSpline


class CFDThroatCharacteristic:
    def __init__(self, X, R, M, TH, x_window=(-2.0, 2.5)):
        xs = X[:, 0]; eta = R[0] / R[0, -1]
        if np.abs(X - xs[:, None]).max() > 1e-12:
            raise ValueError("断面が x 一定でない")
        E = R / R[:, -1:]
        if np.abs(E - eta[None, :]).max() > 1e-6:
            raise ValueError("η 分布が断面ごとに違う (RectBivariateSpline が使えない)")
        m = (xs >= x_window[0]) & (xs <= x_window[1])
        self.xs, self.eta = xs[m], eta
        self._rw = CubicSpline(xs, R[:, -1])          # 壁は断面の壁節点を通る 3 次補間 (直線補間はスロートで弦誤差 κΔx²/8 ≈ 7e-5)
        self._M = RectBivariateSpline(self.xs, eta, M[m], kx=3, ky=3)
        self._T = RectBivariateSpline(self.xs, eta, TH[m], kx=3, ky=3)
        self._axis = make_smoothing_spline(self.xs, M[m, 0], lam=1e-6)

    @classmethod
    def from_run(cls, run_dir, n_tail=3, res=None):
        rd = Path(run_dir); info = json.loads((rd / "prepare_info.json").read_text()); S = float(info["scale_m"]); ni = int(info["mesh"]["ni"])
        with h5py.File(rd / "nozzle.h5") as f:
            nc = f["/MESH/COORD"][:].reshape(-1, 3)
        X = (nc[:, 0] / S).reshape(ni, -1); R = (nc[:, 1] / S).reshape(ni, -1)
        files = [rd / r for r in res] if res else sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-n_tail:]
        Ms, Ts = [], []
        for fp in files:
            with h5py.File(fp) as f:
                Ux, Uy, son = (f["/VALUE/" + k][:] for k in ("Ux", "Uy", "sonic"))
            Ms.append(np.hypot(Ux, Uy) / son); Ts.append(np.arctan2(Uy, Ux))
        obj = cls(X, R, np.mean(Ms, 0).reshape(ni, -1), np.mean(Ts, 0).reshape(ni, -1))
        obj.source = {"run": str(rd), "res": [f.name for f in files]}
        return obj

    def r_wall(self, x):
        return float(self._rw(float(x)))

    def field(self, x, r):
        x = float(x); e = min(max(float(r) / self.r_wall(x), 0.0), 1.0)
        return float(np.ravel(self._M.ev(x, e))[0]), float(np.ravel(self._T.ev(x, e))[0])

    def mach(self, x, r):
        x = np.atleast_1d(np.asarray(x, float)); r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self.field(a, b)[0] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])

    def theta(self, x, r):
        x = np.atleast_1d(np.asarray(x, float)); r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self.field(a, b)[1] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])

    def axis_anchor(self, x):
        return float(self._axis(x)), float(self._axis(x, 1)), float(self._axis(x, 2))

    def throat_characteristic(self, n: int = 61, ds: float = 2e-4, x_start: float = 0.0, **_):
        def slope(p):
            Mv, tv = self.field(p[0], p[1])
            a = tv - np.arcsin(1.0 / max(Mv, 1.0 + 1e-9))
            return np.r_[np.cos(a), np.sin(a)]
        p = np.r_[x_start, self.r_wall(x_start)]; P = [p.copy()]
        while p[1] > 0 and len(P) < 400000:
            k1 = slope(p); k2 = slope(p + 0.5 * ds * k1); k3 = slope(p + 0.5 * ds * k2); k4 = slope(p + ds * k3)
            q = p + ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            if q[1] <= 0:          # 軸を跨いだ: 線形に r=0 へ
                q = p + (q - p) * (p[1] / (p[1] - q[1])); q[1] = 0.0
            P.append(q.copy()); p = q
        P = np.array(P)[::-1]                       # 軸→壁
        rq = np.linspace(0.0, P[-1, 1], n)
        xq = np.interp(rq, P[:, 1], P[:, 0])
        Mq = np.array([self.field(a, b)[0] for a, b in zip(xq, rq)])
        tq = np.array([self.field(a, b)[1] for a, b in zip(xq, rq)])
        tq[0] = 0.0; tq[-1] = 0.0                   # 軸と壁足 (設計スロートの壁接線) の厳密値
        return xq, rq, Mq, tq
