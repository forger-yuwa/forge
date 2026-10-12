"""区分構成ノズル壁 (無次元, r* = 1)。

区間: 入口直管 → 収縮 (両端 C2 の 5 次 Hermite; Bell–Mehta の趣旨を保ちつつ
円弧側は曲率まで一致させる — サーベイ A2.5(E) 「Bell–Mehta の罠」の解) →
上流円弧 R_u → スロート → 下流円弧 R_d (終端角 theta_a) → ベル初期壁
(C1 の 3 次 Hermite; 逆設計の出発点であり最終形状は帰還エンジンが決める) → 出口。

座標: x = 軸 (スロート = 0), r = 半径。全長は r* 単位。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _hermite_quintic(x0, x1, f0, d0, s0, f1, d1, s1):
    """両端の (値, 1階, 2階) を満たす 5 次多項式係数 (ξ=(x-x0)/(x1-x0) 基底)。"""
    h = x1 - x0
    A = np.zeros((6, 6))
    b = np.array([f0, d0 * h, s0 * h * h, f1, d1 * h, s1 * h * h], dtype=float)
    for j in range(6):
        A[0, j] = 1.0 if j == 0 else 0.0
        A[1, j] = 1.0 if j == 1 else 0.0
        A[2, j] = 2.0 if j == 2 else 0.0
        A[3, j] = 1.0
        A[4, j] = j
        A[5, j] = j * (j - 1)
    return np.linalg.solve(A, b)  # c[j] ξ^j


def _poly_eval(c, x0, x1, x, order=0):
    h = x1 - x0
    xi = (np.asarray(x, dtype=float) - x0) / h
    j = np.arange(len(c))
    if order == 0:
        return np.polynomial.polynomial.polyval(xi, c)
    cc = c.copy()
    for _ in range(order):
        cc = cc[1:] * np.arange(1, len(cc))
    return np.polynomial.polynomial.polyval(xi, cc) / h ** order


@dataclass
class WallParams:
    r_inlet: float = 2.5        # 入口半径 / r*
    L_pipe: float = 0.5         # 入口直管長
    L_contract: float = 3.0     # 収縮長 (直管端 → 上流円弧始端)
    Ru: float = 1.5             # 上流円弧曲率半径
    Rd: float = 0.4             # 下流円弧曲率半径
    phi_u: float = np.deg2rad(25.0)   # 上流円弧の入口側角
    theta_a: float = np.deg2rad(20.0) # 下流円弧終端角 (= 初期膨張角 θn)
    L_div: float = 7.0          # スロート → 出口の軸長
    r_exit: float = 3.0         # 出口半径 (= sqrt(面積比))
    theta_e: float = np.deg2rad(8.0)  # 出口壁角
    bell_type: str = "hermite"  # 下流壁の族: "hermite" (3次) | "top" (Rao TOP 放物線)


class NozzleWall:
    """r_w(x) と導関数を返す区分構成壁。"""

    def __init__(self, p: WallParams) -> None:
        self.p = p
        # 上流円弧始端 (角 phi_u)
        self.x_au = -p.Ru * np.sin(p.phi_u)
        self.r_au = 1.0 + p.Ru * (1.0 - np.cos(p.phi_u))
        # 収縮区間
        self.x_cs = self.x_au - p.L_contract
        self.x_in = self.x_cs - p.L_pipe
        curv_au = 1.0 / (p.Ru * np.cos(p.phi_u) ** 3)  # 円弧の d2r/dx2
        self._c_con = _hermite_quintic(
            self.x_cs, self.x_au,
            p.r_inlet, 0.0, 0.0,
            self.r_au, -np.tan(p.phi_u), curv_au,
        )
        # 下流円弧終端
        self.x_a = p.Rd * np.sin(p.theta_a)
        self.r_a = 1.0 + p.Rd * (1.0 - np.cos(p.theta_a))
        # ベル壁: TOP 放物線 (2 次 Bézier) または 3 次 Hermite
        h = p.L_div - self.x_a
        if h <= 0:
            raise ValueError("L_div が下流円弧より短い")
        self.x_e = p.L_div
        if p.bell_type == "top":
            # Rao TOP: P_a で接線 θn(=θa)・出口で接線 θe の放物線 (任意の 2 次 Bézier
            # は放物線)。制御点 P1 = 両接線の交点。整合には θn > 弦勾配 > θe が必要
            s0, s2 = np.tan(p.theta_a), np.tan(p.theta_e)
            sc = (p.r_exit - self.r_a) / h
            if not (s0 > sc > s2):
                raise ValueError(
                    f"TOP 整合 NG: tanθn({s0:.4f}) > 弦勾配({sc:.4f}) > tanθe({s2:.4f}) が必要")
            x1 = (p.r_exit - self.r_a + s0 * self.x_a - s2 * self.x_e) / (s0 - s2)
            r1 = self.r_a + s0 * (x1 - self.x_a)
            # x(t) = x_a + Bx·t + Ax·t², r(t) 同型 (t∈[0,1] の Bernstein 2 次)
            self._top_Bx = 2.0 * (x1 - self.x_a)
            self._top_Ax = self.x_a - 2.0 * x1 + self.x_e
            self._top_Br = 2.0 * (r1 - self.r_a)
            self._top_Ar = self.r_a - 2.0 * r1 + p.r_exit
        elif p.bell_type == "hermite":
            A = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1], [0, 1, 2, 3]], dtype=float)
            b = np.array([self.r_a, np.tan(p.theta_a) * h, p.r_exit, np.tan(p.theta_e) * h])
            self._c_bell = np.linalg.solve(A, b)
        else:
            raise ValueError(f"bell_type '{p.bell_type}' は未実装 (hermite | top)")

    def _top_t(self, x):
        """TOP 放物線の t(x) (Citardauq 型 — A→0 でも安定, x(t) 単調前提)。"""
        dx = np.asarray(x, dtype=float) - self.x_a
        B, A = self._top_Bx, self._top_Ax
        disc = np.sqrt(np.maximum(B * B + 4.0 * A * dx, 0.0))
        return 2.0 * dx / (B + disc)

    # -- 評価 ---------------------------------------------------------------
    def r(self, x):
        x = np.asarray(x, dtype=float)
        p = self.p
        out = np.empty_like(x)
        m_pipe = x < self.x_cs
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_ru = (x >= self.x_au) & (x < 0.0)
        m_rd = (x >= 0.0) & (x < self.x_a)
        m_bell = x >= self.x_a
        out[m_pipe] = p.r_inlet
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con])
        out[m_ru] = 1.0 + p.Ru - np.sqrt(np.maximum(p.Ru ** 2 - x[m_ru] ** 2, 0.0))
        out[m_rd] = 1.0 + p.Rd - np.sqrt(np.maximum(p.Rd ** 2 - x[m_rd] ** 2, 0.0))
        if p.bell_type == "top":
            t = self._top_t(x[m_bell])
            out[m_bell] = self.r_a + self._top_Br * t + self._top_Ar * t * t
        else:
            out[m_bell] = _poly_eval(self._c_bell, self.x_a, self.x_e, x[m_bell])
        return out

    def drdx(self, x):
        x = np.asarray(x, dtype=float)
        p = self.p
        out = np.zeros_like(x)
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_ru = (x >= self.x_au) & (x < 0.0)
        m_rd = (x >= 0.0) & (x < self.x_a)
        m_bell = x >= self.x_a
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con], order=1)
        out[m_ru] = x[m_ru] / np.sqrt(np.maximum(p.Ru ** 2 - x[m_ru] ** 2, 1e-30))
        out[m_rd] = x[m_rd] / np.sqrt(np.maximum(p.Rd ** 2 - x[m_rd] ** 2, 1e-30))
        if p.bell_type == "top":
            t = self._top_t(x[m_bell])
            out[m_bell] = ((self._top_Br + 2.0 * self._top_Ar * t)
                           / (self._top_Bx + 2.0 * self._top_Ax * t))
        else:
            out[m_bell] = _poly_eval(self._c_bell, self.x_a, self.x_e, x[m_bell], order=1)
        return out

    # -- 健全性 (事前フィルタの一部) -----------------------------------------
    def validate(self) -> list:
        """幾何の破綻を検査してメッセージ一覧を返す (空なら OK)。"""
        msgs = []
        xs = np.linspace(self.x_in, self.x_e, 2000)
        r = self.r(xs)
        if np.any(r <= 0.0):
            msgs.append("壁半径が非正")
        rmin = float(r.min())
        if abs(rmin - 1.0) > 5e-3:
            msgs.append(f"最小半径 {rmin:.4f} != 1 (スロートが最小でない)")
        d = self.drdx(xs)
        m_bell = xs >= self.x_a
        if np.any(d[m_bell] < -1e-9):
            msgs.append("ベル区間で半径が非単調 (収縮している)")
        m_con = (xs >= self.x_cs) & (xs < self.x_au)
        if np.any(d[m_con] > 1e-9):
            msgs.append("収縮区間で半径が増加")
        return msgs


class DirectBellWall:
    r"""③ベル (D) 直接モードの壁 (`wall_bell_c2.BellC2Wall`) を runner の壁の口 (`x_in`・`x_e`・`r(x)`・`r(x, 1)`・`drdx`) に
    合わせる薄い包み (plan tooling-nozzle-bell-two-mode §5.1 #3f)。値は `BellC2Wall.r` そのまま (C2 スロート・U→T Hermite の収縮部・
    拡大部の 5 次 B-spline)。スロートは $(0, 1)$ (IC の `_throat_of` の既定と同じ)。

    `BellC2Wall.r` は $[x_{in}, L]$ の外で `ValueError` にする。格子座標を書き戻した節点・双対重心は丸めで端をわずかに出るので、
    端から `END_TOL` ($r_t$ 単位) 以内だけ端へ寄せて評価する (それより外は従来どおり例外: 黙って外挿しない)。"""

    END_TOL = 1.0e-9

    def __init__(self, c2wall) -> None:
        self.c2 = c2wall
        self.x_in = float(c2wall.x_in)
        self.x_e = float(c2wall.L)
        self.x_throat, self.r_throat = 0.0, 1.0

    def _clip(self, x):
        x = np.asarray(x, dtype=float)
        tol = self.END_TOL * max(1.0, abs(self.x_in), abs(self.x_e))
        if np.any(x < self.x_in - tol) or np.any(x > self.x_e + tol):
            raise ValueError(f"DirectBellWall: x が [{self.x_in}, {self.x_e}] の外 (端の許容 {tol:.1e})")
        return np.clip(x, self.x_in, self.x_e)

    def r(self, x, deriv: int = 0):
        return np.asarray(self.c2.r(self._clip(x), deriv), dtype=float)

    def drdx(self, x):
        return self.r(x, 1)


class RaoBellWall:
    r"""③ベル (M) 特性曲線モードの壁 (plan tooling-nozzle-bell-two-mode §4.2・§4.4・§5.1 #6 (W1))。runner の壁の口
    (`x_in`・`x_e`・`r(x)`・`r(x, 1)`・`drdx`・スロート `(x_throat, r_throat)` = (0, 1)) に合わせる。

    - 入口直管 ($x<-L_U$、半径 $r_U=\sqrt{CR}$) → 収縮部 ($-L_U\le x<0$、U→T 5 次 Hermite `wall_walldriven.UpstreamThroatPoly`;
      (D) の `BellC2Wall` と同じ部品・同じ寸法 spec `BellC2Spec`) → 拡大部 ($0\le x\le x_e$、`rao_axisym.RaoWall`: 円弧 kernel
      [曲率半径 $R_t$、閉じた式] + BE [生成器の壁点列の罰則付き 5 次 B-spline])。
    - **スロートで C2**: 収縮部の終端 $(1, 0, 1/R_t)$ と円弧の始端 $(1, 0, 1/R_t)$ が一致する (構成上。`throat_jump` で検査)。
    - **B で C1**: 円弧の終端と BE の始端で $r$・$r'$ が一致 (BE の当てはめの拘束)。曲率は跳ぶ ($\Delta\kappa_B$、`B_jump`) —
      (M) は古典どおり円弧 kernel (plan §4.4「(M) は円弧 kernel・B で C1」)。
    - 出口 $x_e$ = BE の終点 (生成器の E)。$x_e$ より先 (`RaoWall` の接線の延長) は使わない。

    `r(x, deriv)` は deriv 0..2 ($r''$ は記録と試験用: 円弧は閉じた式、BE は B-spline の 2 階微分)。$[x_{in},x_e]$ の外は
    `ValueError` (端から `END_TOL` 以内の丸めだけ端へ寄せる — `DirectBellWall` と同じ)。
    `offset` (`DisplacementOffset`、既定なし): δ* 補正の物理壁 = 非粘性壁 + δ(x) (固定 x の半径方向、拡大部だけ。
    δ はスロートで値・1・2 階が 0 なので、スロートの C2 と物理スロートの r・r′・r″ は変わらない — plan §5.1 #7 (7b-4))。"""

    END_TOL = 1.0e-9

    def __init__(self, rao_wall, spec, offset=None) -> None:
        from .wall_walldriven import UpstreamThroatPoly
        spec.check()
        if abs(float(rao_wall.R_t) - float(spec.R_t)) > 1e-9 * float(spec.R_t):
            raise ValueError(f"RaoBellWall: 壁の R_t {rao_wall.R_t} と寸法 spec の R_t {spec.R_t} が違う (スロートで C2 にならない)")
        self.rao, self.spec, self.offset = rao_wall, spec, offset
        if offset is not None and abs(float(offset.x_end) - float(rao_wall.x_end)) > 1e-6 * max(1.0, float(rao_wall.x_end)):
            raise ValueError(f"RaoBellWall: δ* のずれの終点 {offset.x_end} と非粘性壁の終点 {rao_wall.x_end} が違う")
        self.up = UpstreamThroatPoly(r_U=spec.r_U, R_t=spec.R_t, L_U=spec.L_U, r_t=1.0)
        self.x_in = -float(spec.L_U) - float(spec.L_pipe)
        self.x_e = float(rao_wall.x_end)
        self.x_B = float(rao_wall.x_B)
        self.x_throat, self.r_throat = 0.0, 1.0
        self._d2spl = rao_wall.spl.derivative(2)

    def _clip(self, x):
        x = np.asarray(x, dtype=float)
        tol = self.END_TOL * max(1.0, abs(self.x_in), abs(self.x_e))
        if np.any(x < self.x_in - tol) or np.any(x > self.x_e + tol) or not np.all(np.isfinite(x)):
            raise ValueError(f"RaoBellWall: x が [{self.x_in}, {self.x_e}] の外 (端の許容 {tol:.1e}) または非有限")
        return np.clip(x, self.x_in, self.x_e)

    def _expansion(self, x, deriv: int):
        """拡大部 $[0, x_e]$ の $r^{(deriv)}$ (deriv 0..2)。δ* のずれ (`offset`) があれば足す (物理壁 = 非粘性壁 + δ)。"""
        out = self._expansion_inviscid(x, deriv)
        if self.offset is not None:
            out = out + self.offset(x, deriv)
        return out

    def _expansion_inviscid(self, x, deriv: int):
        if deriv in (0, 1):
            return np.asarray(self.rao.r(x, deriv), dtype=float)
        R = float(self.rao.R_t)
        out = np.empty_like(x)
        k = x <= self.x_B
        out[k] = R * R / np.maximum(R * R - x[k] * x[k], 1e-300) ** 1.5
        out[~k] = self._d2spl(x[~k])
        return out

    def r(self, x, deriv: int = 0):
        if deriv not in (0, 1, 2):
            raise ValueError("RaoBellWall.r: deriv は 0..2")
        x = np.asarray(x, dtype=float)
        xa = np.atleast_1d(self._clip(x))
        out = np.empty_like(xa)
        m_pipe = xa < -self.spec.L_U
        m_con = (xa >= -self.spec.L_U) & (xa < 0.0)
        m_exp = xa >= 0.0
        out[m_pipe] = self.spec.r_U if deriv == 0 else 0.0
        if np.any(m_con):
            out[m_con] = self.up.r(xa[m_con], deriv)
        if np.any(m_exp):
            out[m_exp] = self._expansion(xa[m_exp], deriv)
        return out.reshape(np.shape(x))

    def drdx(self, x):
        return self.r(x, 1)

    def throat_jump(self) -> dict:
        """スロート $x=0$ での収縮部 (左) と円弧 (右) の $r$・$r'$・$r''$ の差 (C2 なら丸めの範囲で 0)。"""
        z = np.array([0.0])
        return {k: float(self._expansion(z, k)[0] - self.up.r(z, k)[0]) for k in (0, 1, 2)}

    def B_jump(self) -> dict:
        r"""B での BE (右) と円弧 (左) の差: $r$・$r'$ (C1 なら丸めの範囲で 0) と $r''$・曲率 $\kappa$ の跳び (円弧 kernel では非零)。
        非粘性壁の値 (δ* のずれは C2 なので跳びを変えない)。"""
        R, xB = float(self.rao.R_t), self.x_B
        left = (1.0 + R - np.sqrt(R * R - xB * xB), xB / np.sqrt(R * R - xB * xB), R * R / (R * R - xB * xB) ** 1.5)
        right = (float(self.rao.spl(xB)), float(self.rao.dspl(xB)), float(self._d2spl(xB)))
        kap = lambda rp, rpp: rpp / (1.0 + rp * rp) ** 1.5  # noqa: E731
        return {0: right[0] - left[0], 1: right[1] - left[1], 2: right[2] - left[2],
                "kappa_left": kap(left[1], left[2]), "kappa_right": kap(right[1], right[2]),
                "dkappa": kap(right[1], right[2]) - kap(left[1], left[2])}

    def validate(self, n: int = 4001) -> list:
        """壁の破綻の検査 (空なら OK): 収縮部の単調性 (`UpstreamThroatPoly.validate`)、拡大部で $r'\\ge0$ (最小半径がスロート)、
        半径が正。"""
        msgs = [f"収縮部: {m}" for m in self.up.validate()]
        xs = np.linspace(0.0, self.x_e, int(n))
        rp = self.r(xs, 1)
        if float(rp.min()) < -1e-12:
            msgs.append(f"拡大部で r' < 0 (最小 {float(rp.min()):.3e} @ x = {float(xs[int(np.argmin(rp))]):.4g})")
        rr = self.r(np.r_[np.linspace(self.x_in, self.x_e, int(n)), 0.0])
        if float(rr.min()) < 1.0 - 1e-12:
            msgs.append(f"最小半径 {float(rr.min()):.12g} がスロートの 1 を下回る")
        return msgs


class DisplacementOffset:
    r"""δ* 補正の半径方向のずれ $\delta(x)$ (plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-4))。物理壁 = 非粘性壁 + $\delta(x)$。

    $\delta(x)=f(x)\,S(t)$、$f$ = 較正した理論式 $k_\delta C(x)$ の 3 次の平滑化スプライン (`make_smoothing_spline`、λ は GCV、
    $x\in[0,x_{end}]$)、$S$ = 5 次の滑らかな立ち上げ $10t^3-15t^4+6t^5$ ($t=(x-a)/(b-a)$ を [0, 1] に切る)。$S$ の 0・1・2 階は
    $t=0$ で 0 なので、$a\ge0$ なら $x\le a$ で $\delta=\delta'=\delta''=0$ (スロートの $r$・$r'$・$r''$ を保つ)。$t\ge1$ では $\delta=f$。
    `__call__(x, deriv)` は deriv 0..2 (積の微分)。保存は npz (スプラインの t・c・k、立ち上げ、終点、来歴の JSON)。"""

    def __init__(self, spl, ramp, x_end: float, meta: dict | None = None, fit_info: dict | None = None) -> None:
        a, b = float(ramp[0]), float(ramp[1])
        if not (0.0 <= a < b <= float(x_end)):
            raise ValueError(f"DisplacementOffset: 立ち上げ [{a}, {b}] は 0 ≤ a < b ≤ x_end ({x_end})")
        self.spl, self.ramp, self.x_end = spl, (a, b), float(x_end)
        self.meta = dict(meta or {})
        self.fit_info = dict(fit_info or {})
        self._d = [spl, spl.derivative(1), spl.derivative(2)]

    END_WEIGHT = 1.0e8   # 終点の重み: 出口の δ を較正値に (ほぼ) 一致させる — ε_inv = (√ε − δ_E)² と物理壁の出口半径 √ε の整合

    @classmethod
    def fit(cls, x, d, ramp, meta: dict | None = None, lam: float | None = None) -> "DisplacementOffset":
        """`x` は 0 から始まる単調増加 (先頭は x = 0 の点)。λ は既定 GCV。終点は重み `END_WEIGHT` で値に寄せる。"""
        from scipy.interpolate import make_smoothing_spline
        x, d = np.asarray(x, float), np.asarray(d, float)
        if np.any(np.diff(x) <= 0.0) or abs(x[0]) > 1e-12 or not np.all(np.isfinite(d)):
            raise ValueError("DisplacementOffset.fit: x は 0 から単調増加・d は有限")
        w = np.ones_like(x); w[-1] = cls.END_WEIGHT
        spl = make_smoothing_spline(x, d, w=w, lam=lam)
        res = spl(x) - d
        info = {"n": int(len(x)), "lam": lam, "lam_rule": "GCV" if lam is None else "given", "end_weight": cls.END_WEIGHT,
                "fit_rms": float(np.sqrt(np.mean(res ** 2))), "fit_max": float(np.max(np.abs(res))),
                "fit_end_resid": float(res[-1]),
                "fit_rel_max_x_ge_0p1": float(np.max(np.abs(res[x >= 0.1] / d[x >= 0.1])))}
        return cls(spl, ramp, float(x[-1]), meta, info)

    @staticmethod
    def _S(t, deriv):
        if deriv == 0:
            return t * t * t * (10.0 - 15.0 * t + 6.0 * t * t)
        if deriv == 1:
            return 30.0 * t * t * (1.0 - t) ** 2
        return 60.0 * t * (1.0 - t) * (1.0 - 2.0 * t)

    def __call__(self, x, deriv: int = 0):
        if deriv not in (0, 1, 2):
            raise ValueError("DisplacementOffset: deriv は 0..2")
        x = np.asarray(x, float)
        a, b = self.ramp
        h = b - a
        xc = np.clip(x, 0.0, self.x_end)
        t = np.clip((xc - a) / h, 0.0, 1.0)
        f = [self._d[k](xc) for k in range(deriv + 1)]
        S = [self._S(t, k) / h ** k for k in range(deriv + 1)]
        if deriv == 0:
            out = f[0] * S[0]
        elif deriv == 1:
            out = f[1] * S[0] + f[0] * S[1]
        else:
            out = f[2] * S[0] + 2.0 * f[1] * S[1] + f[0] * S[2]
        return np.where(x >= 0.0, out, 0.0)   # x > x_end は終点の値で頭打ち (壁の終点との丸めの差で 0 に落とさない)

    def save(self, path) -> None:
        import json
        np.savez(path, t=self.spl.t, c=self.spl.c, k=np.array(self.spl.k), ramp=np.array(self.ramp), x_end=np.array(self.x_end),
                 meta=np.array(json.dumps(self.meta, ensure_ascii=False)), fit_info=np.array(json.dumps(self.fit_info)))

    @classmethod
    def load(cls, path) -> "DisplacementOffset":
        import json
        from scipy.interpolate import BSpline
        with np.load(path) as z:
            spl = BSpline(z["t"], z["c"], int(z["k"]))
            return cls(spl, tuple(z["ramp"]), float(z["x_end"]), json.loads(str(z["meta"])), json.loads(str(z["fit_info"])))

    def record(self) -> dict:
        xs = np.linspace(0.0, self.x_end, 20001)
        d1, d2 = self(xs, 1), self(xs, 2)
        z = np.array([0.0])
        return {**self.fit_info, "ramp": list(self.ramp), "x_end": self.x_end, "meta": self.meta,
                "d0": float(self(z)[0]), "d1_0": float(self(z, 1)[0]), "d2_0": float(self(z, 2)[0]),
                "delta_end": float(self(np.array([self.x_end]))[0]),
                "d1_min": float(d1.min()), "d1_max": float(d1.max()), "x_d1_max": float(xs[int(np.argmax(d1))]),
                "d2_absmax": float(np.max(np.abs(d2))), "x_d2_absmax": float(xs[int(np.argmax(np.abs(d2)))])}
