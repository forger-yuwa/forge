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
