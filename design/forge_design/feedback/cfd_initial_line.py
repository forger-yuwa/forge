r"""CFD ピン: node Euler 場からスロート特性線 (MOC の初期線) と軸アンカーを取る provider。

計画: plans/active/tooling-nozzle-cfd-pinned-initial-line.md §4.1・§4.2 (diagnostician 2026-10-05 採用版)。
仕様: methods/design/overview.md「初期線の出所: Hall / CFD ピン」。case/45 の試作
(`case/45.isobutane_m6_d155/cfd_initial_line.py`) を生産コードへ移したもので、挙動・入力の拒否は試作と同じ。

`CFDPinnedThroat` は `HallThroat` を継承し `throat_characteristic` / `axis_anchor` / `mach` / `theta` だけを
上書きする (`cd_series` 等は Hall のまま)。`design_chain` は `geometry.initial_line: cfd` でこれを使う。

- 場: 構造格子 (断面は x 一定) の M, θ を (x, η = r/r_w(x)) 平面の 3 次スプラインで補間。壁 r_w(x) は断面の壁節点の 3 次補間。
- 線: 壁足は (0, 1, θ=0) に厳密。追跡は壁足 1e-6 内側から dr/dx = tan(θ−μ) を RK4 (ds 2e-4) で軸を跨ぐまで。
  壁足の M は内側 10 点の線形外挿、軸端 (r=0) の x・M・θ は追跡点の r ≤ 0.05 の 2 次多項式 (偶関数当てはめは使わない —
  C⁻ は軸を斜めに横切る)。
- 軸アンカー: M_A = 線の軸端の M (線と同じ出所)、M′_A = 軸の evenfit (各断面で r>0 の 4 点に M=a₀+a₂r²) の 4 次窓フィット
  (窓 ±0.25)、M″_A は Hall の式を x₀ で評価 (CFD の局所 M″ は law を実測軸から離すので使わない)。
- 凍結源: 同じ形・ガスで Hall 初期線の V0 型壁を Euler で解いた場 (ピン壁自身の Euler は循環定義になるので使わない)。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

from ..geometry.transonic import HallThroat


def load_run_field(run_dir, res=None, n_tail: int = 1):
    """run の node 場 (r_t 単位の格子 X, R と M, θ) を読む。res: ファイル名 (str) かその列 (平均する)。
    None なら末尾 n_tail 枚。対応範囲 (codex plan M5): 構造化 node・軸対称・Euler (slip) の run だけ。不適合は拒否する。"""
    import h5py
    rd = Path(run_dir)
    info = json.loads((rd / "prepare_info.json").read_text())
    S = float(info["scale_m"])
    ni = int(info["mesh"]["ni"])
    cfg = (rd / "solverConfig.yaml").read_text()
    if 'discretization: "node"' not in cfg.replace("'", '"') or "isAxisymmetric: 1" not in cfg:
        raise ValueError("CFD ピン: node・軸対称の run でない")
    if info.get("viscous") not in (None, False, "none", "euler"):
        raise ValueError(f"CFD ピン: Euler run でない (viscous={info.get('viscous')})")
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    X = (nc[:, 0] / S).reshape(ni, -1)
    R = (nc[:, 1] / S).reshape(ni, -1)
    if isinstance(res, (str, Path)):
        res = [res]
    files = ([rd / r for r in res] if res else
             sorted(rd.glob("res_[0-9]*.h5"),
                    key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-n_tail:])
    if not files:
        raise ValueError(f"CFD ピン: {rd} に res_*.h5 が無い")
    Ms, Ts = [], []
    for fp in files:
        with h5py.File(fp) as f:
            Ux, Uy, son = (f["/VALUE/" + k][:] for k in ("Ux", "Uy", "sonic"))
        Ms.append(np.hypot(Ux, Uy) / son)
        Ts.append(np.arctan2(Uy, Ux))
    return (X, R, np.mean(Ms, 0).reshape(ni, -1), np.mean(Ts, 0).reshape(ni, -1),
            {"run": str(rd), "res": [f.name for f in files]})


class CFDPinnedThroat(HallThroat):
    """`HallThroat` 互換の CFD ピン provider (初期線・軸アンカー・場の M/θ を CFD 場から)。"""

    def __init__(self, R, gamma, X, Rr, M, TH, source=None, ds=2e-4, x_window=(-2.0, 2.5),
                 anchor_halfwidth=0.25):
        super().__init__(R=R, gamma=gamma)
        xs = X[:, 0]
        eta = Rr[0] / Rr[0, -1]
        if np.abs(X - xs[:, None]).max() > 1e-12 or np.abs(Rr / Rr[:, -1:] - eta[None, :]).max() > 1e-6:
            raise ValueError("断面が x 一定でない / η 分布が断面ごとに違う")
        m = (xs >= x_window[0]) & (xs <= x_window[1])
        self._xs, self._eta = xs[m], eta
        self._rw = CubicSpline(xs, Rr[:, -1])
        self._Ms = RectBivariateSpline(self._xs, eta, M[m], kx=3, ky=3)
        self._Ts = RectBivariateSpline(self._xs, eta, TH[m], kx=3, ky=3)
        # 軸の evenfit: 各断面で r>0 の 4 点に M = a0 + a2 r²
        rr = Rr[m, 1:5]
        mm = M[m, 1:5]
        self._ax_x = self._xs
        self._ax_M = np.array([np.polyfit(rr[i] ** 2, mm[i], 1)[1] for i in range(len(self._xs))])
        self.source = source or {}
        self._build_line(ds)
        self._anchor_hw = anchor_halfwidth

    # -- 場 --
    def _field(self, x, r):
        x = float(x)
        e = min(max(float(r) / float(self._rw(x)), 0.0), 1.0)
        return float(np.ravel(self._Ms.ev(x, e))[0]), float(np.ravel(self._Ts.ev(x, e))[0])

    def _build_line(self, ds):
        def slope(p):
            if not (self._xs[0] <= p[0] <= self._xs[-1]):
                raise ValueError(f"CFD ピン: 追跡点 x={p[0]:.4f} が場の窓 "
                                 f"[{self._xs[0]:.3f}, {self._xs[-1]:.3f}] の外")
            Mv, tv = self._field(p[0], p[1])
            if not (np.isfinite(Mv) and np.isfinite(tv)):
                raise ValueError(f"CFD ピン: 非有限の場 (x={p[0]:.4f}, r={p[1]:.4f})")
            if Mv <= 1.0:
                raise ValueError(f"CFD ピン: 亜音速 M={Mv:.5f} (x={p[0]:.4f}, r={p[1]:.4f}) — "
                                 "スロート特性線が引けない")
            a = tv - np.arcsin(1.0 / Mv)
            return np.r_[np.cos(a), np.sin(a)]
        p = np.r_[0.0, 1.0 - 1e-6]
        P = [p.copy()]
        while p[1] > 0 and len(P) < 400000:
            k1 = slope(p)
            k2 = slope(p + 0.5 * ds * k1)
            k3 = slope(p + 0.5 * ds * k2)
            k4 = slope(p + ds * k3)
            p = p + ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            P.append(p.copy())
        if p[1] > 0:
            raise ValueError("CFD ピン: 軸に到達しない")
        P = np.array(P[:-1])
        P = P[P[:, 1] > 0]                                                 # 軸を跨いだ最後の点は捨てる
        F = np.array([self._field(a, b) for a, b in P])
        nr = P[:, 1] <= 0.05
        ax = [np.polyval(np.polyfit(P[nr, 1], v, 2), 0.0) for v in (P[nr, 0], F[nr, 0], F[nr, 1])]
        top = np.argsort(-P[:, 1])[:10]                                    # 壁足の M: 内側 10 点の線形外挿
        M_foot = float(np.polyval(np.polyfit(P[top, 1], F[top, 0], 1), 1.0))
        o = np.argsort(P[:, 1])
        self._line_r = np.r_[0.0, P[o, 1], 1.0]
        self._line_x = np.r_[ax[0], P[o, 0], 0.0]
        self._line_M = np.r_[ax[1], F[o, 0], M_foot]
        self._line_T = np.r_[0.0, F[o, 1], 0.0]
        self.x0_cfd, self.M_axis_line = float(ax[0]), float(ax[1])

    # -- HallThroat 互換 API の上書き --
    def throat_characteristic(self, n: int = 61, **_):
        r = np.linspace(0.0, 1.0, n)
        return (np.interp(r, self._line_r, self._line_x), r,
                np.interp(r, self._line_r, self._line_M), np.interp(r, self._line_r, self._line_T))

    def axis_anchor(self, x: float) -> tuple:
        x = float(x)
        if abs(x - self.x0_cfd) > 1e-9:
            raise ValueError(f"CFD ピンの軸アンカーは線の軸着地 x0={self.x0_cfd:.6f} でだけ定義 (要求 x={x:.6f})")
        w = np.abs(self._ax_x - x) <= self._anchor_hw
        c = np.polyfit(self._ax_x[w] - x, self._ax_M[w], 4)
        Mpp_hall = super().axis_anchor(x)[2]
        return self.M_axis_line, float(c[-2]), float(Mpp_hall)

    def mach(self, x, r):
        x = np.atleast_1d(np.asarray(x, float))
        r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self._field(a, b)[0] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])

    def theta(self, x, r):
        x = np.atleast_1d(np.asarray(x, float))
        r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self._field(a, b)[1] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])


def pinned_factory(run_dir, res=None):
    """`HallThroat(R=, gamma=)` と同じ呼び方で CFD ピンの throat を返す factory (場は 1 回だけ読む)。"""
    X, Rr, M, TH, src = load_run_field(run_dir, res)

    def make(R, gamma):
        return CFDPinnedThroat(R, gamma, X, Rr, M, TH, source=src)
    return make
