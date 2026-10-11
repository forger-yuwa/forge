r"""③ベル (D) 直接モードの C2 壁: 拡大部の 5 次 B-spline (無次元, $r_t=1$)。

計画: plans/active/tooling-nozzle-bell-two-mode.md §4.2・§4.3・§5.1 #3。

**規約** (§5.1 #3 で固定、2026-10-11):

- 座標は $r_t$ で無次元、$x$ はスロート基準 (幾何スロート = 0)。
- 拡大部 $x\in[0,L]$ を**等間隔 16 区間の clamped 5 次 B-spline** (端ノット 6 重、内部ノット 15 個、係数 21 個) で表す。
- 係数は 7 個の線形拘束の下で $J=\int_0^L (r''')^2dx$ を最小にする問題を KKT で解いて一意に決める:
  - スロート $x=0$: $r=1$, $r'=0$, $r''=1/R_t$
  - 最大壁角の点 $x_m$: $r'=\tan\theta_{max}$, $r''=0$
  - 出口 $x=L$: $r=\sqrt\varepsilon$, $r'=\tan\theta_e$

  $J$ の零空間は 2 次多項式 ($r'''\equiv0$) で、スロートの 3 拘束がそれを除くので、拘束行列が満行ランク (7) なら
  解は一意 (拘束の零空間の上で $J$ のヘッセ行列が正定値)。ランクと正定値性は構築時に検査する。
- 設計変数 $(L,\theta_e,\theta_{max},x_m/L)$。範囲は $\theta_{max}$ 20〜45°、$\theta_e$ 2〜20°、$x_m/L$ 0.05〜0.40
  ($L$ は比較の長さ $L_k$)。範囲外・非有限は構築時に `ValueError` (入力の誤り)。
- **成立条件** (`feasibility`): 2000 点の等間隔標本で (i) $r'\ge0$、(ii) $\theta'$ の符号変化がちょうど 1 回
  ($\theta'=r''/(1+r'^2)$ なので $r''$ の符号で数える)、(iii) $|\kappa|\le1/R_t$
  ($\kappa=r''/(1+r'^2)^{3/2}$)。破れたら INFEASIBLE (違反を文で返す)。
- **収縮部** (§4.2、両モード共通): 入口直管 (長さ $L_{pipe}$、半径 $r_U=\sqrt{CR}$) → U→T 5 次 Hermite
  (`wall_walldriven.UpstreamThroatPoly`、$(r_U,0,0)\to(1,0,1/R_t)$、長さ $L_U$)。スロートで拡大部と
  $(r,r',r'')$ が一致するので構成上 C2。単調収縮 ⇔ $L_U^2/[R_t(r_U-1)]\le20$。

この壁がそのまま物理壁 (δ\* 補正なし、§4.3)。CFD の runner への組み込みは別項目 (#2b と衝突するため後で行う)。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline

from .wall_walldriven import UpstreamThroatPoly

DEGREE = 5
N_SPAN = 16                       # 拡大部の等間隔区間数
N_COEF = N_SPAN + DEGREE          # 21
N_CONSTRAINT = 7
N_SAMPLE = 2000                   # 成立条件の標本点数 (§5.1 #3)
THETA_MAX_RANGE_DEG = (20.0, 45.0)
THETA_E_RANGE_DEG = (2.0, 20.0)
XM_FRAC_RANGE = (0.05, 0.40)
RESID_TOL = 1e-10                 # 拘束の残差の合格上限 (§5.1 #3)
_GAUSS_N = 5                      # Gram 行列の求積点 (r''' は 2 次なので積は 4 次。3 点で厳密、余裕をみて 5 点)


@dataclass(frozen=True)
class BellC2Spec:
    """寸法 spec (§4.2・§8-1、キャンペーン内で固定)。長さは $r_t$ 単位。"""
    eps: float = 50.0                 # 面積比 (出口半径 = √ε)
    R_t: float = 1.0                  # スロートの曲率半径 / r_t (上下流とも; 曲率は 1 つの値 1/R_t)
    contraction_ratio: float = 3.0    # 収縮比 (面積比 A_U/A_t)。r_U = √CR
    L_U: float = 3.0                  # 収縮部 (U→T Hermite) の長さ
    L_pipe: float = 1.0               # 入口直管の長さ

    @property
    def r_exit(self) -> float:
        return float(np.sqrt(self.eps))

    @property
    def r_U(self) -> float:
        return float(np.sqrt(self.contraction_ratio))

    def check(self) -> None:
        for k in ("eps", "R_t", "contraction_ratio", "L_U", "L_pipe"):
            v = getattr(self, k)
            if not np.isfinite(v):
                raise ValueError(f"BellC2Spec: {k} = {v!r} が有限でない")
        if not (self.eps > 1.0 and self.contraction_ratio > 1.0):
            raise ValueError(f"BellC2Spec: eps {self.eps} と contraction_ratio {self.contraction_ratio} は 1 より大")
        if not (self.R_t > 0.0 and self.L_U > 0.0 and self.L_pipe >= 0.0):
            raise ValueError(f"BellC2Spec: R_t {self.R_t}・L_U {self.L_U} は正、L_pipe {self.L_pipe} は非負")


def cone15_length(eps: float, half_angle_deg: float = 15.0) -> float:
    r"""等面積比の円錐ノズルの長さ $(\sqrt\varepsilon-1)/\tan\alpha$ ($r_t$ 単位、スロートの丸めを含めない)。
    plan §8-1 の比較の長さ $L_k$ の基準 (ε 50 で 22.657 $r_t$ = 3.172 m @ $r_t$ 0.14 m)。"""
    return float((np.sqrt(eps) - 1.0) / np.tan(np.deg2rad(half_angle_deg)))


def _sign_changes(v: np.ndarray, rel_zero: float = 1e-10) -> tuple:
    """符号変化の回数と位置 (添字: 変化の直後の標本)。|v| ≤ rel_zero·max|v| の標本は 0 とみなして飛ばす
    (拘束点 x_m で r'' = 0 になる標本を二重に数えないため)。"""
    v = np.asarray(v, dtype=float)
    scale = float(np.max(np.abs(v))) if v.size else 0.0
    idx = np.flatnonzero(np.abs(v) > rel_zero * max(scale, 1e-300))
    if idx.size < 2:
        return 0, []
    s = np.sign(v[idx])
    ch = np.flatnonzero(s[1:] != s[:-1])
    return int(ch.size), [int(idx[c + 1]) for c in ch]


class BellC2Wall:
    r"""(D) の壁。拡大部 $[0,L]$ は拘束付き最小 $\int(r''')^2$ の 5 次 B-spline、上流は直管 + U→T Hermite。

    角度は度で受ける (`theta_e_deg`・`theta_max_deg`)。`r(x, deriv)` は $x\in[-L_{pipe}-L_U,\ L]$ で定義
    (範囲外は `ValueError`)。構築は決定的 (同じ入力で係数はビット同一)。"""

    def __init__(self, L: float, theta_e_deg: float, theta_max_deg: float, xm_frac: float,
                 spec: BellC2Spec | None = None):
        spec = spec or BellC2Spec()
        spec.check()
        for k, v in (("L", L), ("theta_e_deg", theta_e_deg), ("theta_max_deg", theta_max_deg), ("xm_frac", xm_frac)):
            if isinstance(v, bool) or not np.isfinite(v):
                raise ValueError(f"BellC2Wall: {k} = {v!r} が有限でない")
        L, te, tm, xf = float(L), float(theta_e_deg), float(theta_max_deg), float(xm_frac)
        if not L > 0.0:
            raise ValueError(f"BellC2Wall: L = {L} ≤ 0")
        for name, v, (lo, hi) in (("theta_max_deg", tm, THETA_MAX_RANGE_DEG), ("theta_e_deg", te, THETA_E_RANGE_DEG),
                                  ("xm_frac", xf, XM_FRAC_RANGE)):
            if not (lo <= v <= hi):
                raise ValueError(f"BellC2Wall: {name} = {v} が設計変数の範囲 [{lo}, {hi}] の外 (plan §5.1 #3)")
        self.spec = spec
        self.L, self.theta_e_deg, self.theta_max_deg, self.xm_frac = L, te, tm, xf
        self.x_m = xf * L
        self.up = UpstreamThroatPoly(r_U=spec.r_U, R_t=spec.R_t, L_U=spec.L_U, r_t=1.0)
        self.x_in = -spec.L_U - spec.L_pipe
        # --- clamped 5 次 B-spline (等間隔 16 区間) ---
        brk = np.linspace(0.0, L, N_SPAN + 1)
        self.knots = np.concatenate([np.zeros(DEGREE), brk, np.full(DEGREE, L)])
        self._basis = BSpline(self.knots, np.eye(N_COEF), DEGREE, extrapolate=False)
        # --- 拘束 C c = d ---
        rows = [(0.0, 0), (0.0, 1), (0.0, 2), (self.x_m, 1), (self.x_m, 2), (L, 0), (L, 1)]
        self.constraint_def = rows
        self.d = np.array([1.0, 0.0, 1.0 / spec.R_t, np.tan(np.deg2rad(tm)), 0.0, spec.r_exit, np.tan(np.deg2rad(te))])
        self.C = np.vstack([self._basis_row(x, k) for (x, k) in rows])
        # --- 目的 J = c^T Q c (Q = ∫ B_i''' B_j''' dx、区間ごとの Gauss–Legendre) ---
        gx, gw = np.polynomial.legendre.leggauss(_GAUSS_N)
        xs, ws = [], []
        for a, b in zip(brk[:-1], brk[1:]):
            xs.append(0.5 * (a + b) + 0.5 * (b - a) * gx)
            ws.append(0.5 * (b - a) * gw)
        xs, ws = np.concatenate(xs), np.concatenate(ws)
        B3 = self._basis_matrix(xs, 3)
        self.Q = (B3 * ws[:, None]).T @ B3
        # --- 階数と一意性の検査 ---
        sv = np.linalg.svd(self.C, compute_uv=False)
        self.constraint_rank = int(np.sum(sv > sv[0] * 1e-12))
        if self.constraint_rank != N_CONSTRAINT:
            raise ValueError(f"BellC2Wall: 拘束行列の階数 {self.constraint_rank} ≠ {N_CONSTRAINT} (特異値 {sv})")
        _, _, Vt = np.linalg.svd(self.C)
        N = Vt[N_CONSTRAINT:].T                          # 拘束の零空間の基底 (21 × 14)
        ev = np.linalg.eigvalsh(N.T @ self.Q @ N)
        self.reduced_hessian_eig = (float(ev[0]), float(ev[-1]))
        if not ev[0] > 1e-12 * ev[-1]:
            raise ValueError(f"BellC2Wall: 拘束の零空間で ∫(r''')² が正定値でない (固有値 {ev[0]:.3e}, {ev[-1]:.3e})")
        # --- KKT ---
        n, m = N_COEF, N_CONSTRAINT
        K = np.zeros((n + m, n + m))
        K[:n, :n] = 2.0 * self.Q
        K[:n, n:] = self.C.T
        K[n:, :n] = self.C
        rhs = np.concatenate([np.zeros(n), self.d])
        sol = np.linalg.solve(K, rhs)
        self.coef = sol[:n].copy()
        self.multipliers = sol[n:].copy()
        self.kkt_cond = float(np.linalg.cond(K))
        self.expansion = BSpline(self.knots, self.coef, DEGREE, extrapolate=False)
        self.J = float(self.coef @ self.Q @ self.coef)

    # ------------------------------------------------------------------ 基底
    def _basis_matrix(self, x, deriv: int) -> np.ndarray:
        x = np.atleast_1d(np.asarray(x, dtype=float))
        out = np.asarray(self._basis(x, nu=deriv) if deriv else self._basis(x), dtype=float)
        if not np.all(np.isfinite(out)):
            raise ValueError("BellC2Wall: 基底の評価が非有限 (x が [0, L] の外)")
        return out

    def _basis_row(self, x: float, deriv: int) -> np.ndarray:
        return self._basis_matrix(np.array([x]), deriv)[0]

    # ------------------------------------------------------------------ 評価
    def r_expansion(self, x, deriv: int = 0):
        """拡大部 [0, L] の r^{(deriv)} (deriv 0..5)。"""
        x = np.asarray(x, dtype=float)
        if np.any((x < 0.0) | (x > self.L)) or not np.all(np.isfinite(x)):
            raise ValueError(f"BellC2Wall.r_expansion: x は [0, {self.L}] (範囲外または非有限)")
        return np.asarray(self.expansion(x, nu=deriv) if deriv else self.expansion(x), dtype=float)

    def r(self, x, deriv: int = 0):
        """全体の壁 r^{(deriv)}(x)、x ∈ [x_in, L] (入口直管 → U→T Hermite → 拡大部)。deriv 0..3。"""
        if deriv not in (0, 1, 2, 3):
            raise ValueError("deriv は 0..3")
        x = np.asarray(x, dtype=float)
        if np.any((x < self.x_in) | (x > self.L)) or not np.all(np.isfinite(x)):
            raise ValueError(f"BellC2Wall.r: x は [{self.x_in}, {self.L}] (範囲外または非有限)")
        xa = np.atleast_1d(x)
        out = np.empty_like(xa)
        m_pipe = xa < -self.spec.L_U
        m_con = (xa >= -self.spec.L_U) & (xa < 0.0)
        m_exp = xa >= 0.0
        out[m_pipe] = self.spec.r_U if deriv == 0 else 0.0
        if np.any(m_con):
            out[m_con] = self.up.r(xa[m_con], deriv)
        if np.any(m_exp):
            out[m_exp] = self.r_expansion(xa[m_exp], deriv)
        return out.reshape(x.shape) if x.ndim else float(out[0])

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def kappa(self, x):
        rp, rpp = self.r(x, 1), self.r(x, 2)
        return rpp / (1.0 + rp * rp) ** 1.5

    # ------------------------------------------------------------------ 検査
    def constraint_residuals(self) -> np.ndarray:
        """7 拘束の残差 (評価値 − 指定値)。順序は `constraint_def`。"""
        got = np.array([float(self.r_expansion(np.array([x]), k)[0]) for (x, k) in self.constraint_def])
        return got - self.d

    def throat_jump(self) -> dict:
        """スロート x=0 での収縮部 (左) と拡大部 (右) の r・r'・r''・r''' の差 (C2 なら 0..2 階が丸めの範囲で 0)。"""
        out = {}
        for k in (0, 1, 2, 3):
            left = float(self.up.r(np.array([0.0]), k)[0])
            right = float(self.r_expansion(np.array([0.0]), k)[0])
            out[k] = right - left
        return out

    def feasibility(self, n: int = N_SAMPLE) -> dict:
        r"""成立条件 (§5.1 #3): 等間隔 n 点の標本で $r'\ge0$、$\theta'$ の符号変化ちょうど 1 回、$|\kappa|\le1/R_t$。

        許容差は拘束の残差の上限 (`RESID_TOL`) に揃える: $r'\ge-10^{-10}$、$|\kappa|R_t\le1+10^{-10}$
        (スロートでは $r'=0$・$\kappa=1/R_t$ が拘束そのものなので丸めの分だけ緩める)。"""
        x = np.linspace(0.0, self.L, int(n))
        rp = self.r_expansion(x, 1)
        rpp = self.r_expansion(x, 2)
        kap = rpp / (1.0 + rp * rp) ** 1.5
        nchg, at = _sign_changes(rpp)
        kmax = float(np.max(np.abs(kap)) * self.spec.R_t)
        rp_min = float(np.min(rp))
        th = np.degrees(np.arctan(rp))
        v = []
        if rp_min < -RESID_TOL:
            v.append(f"r' < 0 (最小 {rp_min:.4g} @ x = {float(x[int(np.argmin(rp))]):.4g})")
        if nchg != 1:
            v.append(f"θ' の符号変化が {nchg} 回 (1 回であること)")
        if kmax > 1.0 + RESID_TOL:
            v.append(f"|κ| R_t の最大 {kmax:.6g} > 1 (@ x = {float(x[int(np.argmax(np.abs(kap)))]):.4g})")
        return {"ok": not v, "violations": v, "n_sample": int(n), "rp_min": rp_min,
                "theta_sign_changes": nchg, "x_sign_change": [float(x[i]) for i in at],
                "kappa_Rt_max": kmax, "theta_max_sampled_deg": float(th.max()),
                "x_theta_max_sampled": float(x[int(np.argmax(th))]), "r_max": float(np.max(self.r_expansion(x)))}

    def summary(self) -> dict:
        res = self.constraint_residuals()
        tj = self.throat_jump()
        return {"L": self.L, "theta_e_deg": self.theta_e_deg, "theta_max_deg": self.theta_max_deg,
                "xm_frac": self.xm_frac, "x_m": self.x_m, "spec": {**self.spec.__dict__, "r_U": self.spec.r_U,
                                                                    "r_exit": self.spec.r_exit},
                "contraction_mu": float(self.up.mu), "J_rppp2": self.J, "constraint_rank": self.constraint_rank,
                "constraint_resid_max": float(np.max(np.abs(res))), "kkt_cond": self.kkt_cond,
                "reduced_hessian_eig": list(self.reduced_hessian_eig),
                "throat_jump_r_rp_rpp": [tj[0], tj[1], tj[2]], "throat_jump_rppp": tj[3]}
