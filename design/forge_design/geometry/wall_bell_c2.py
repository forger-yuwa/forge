r"""③ベル (D) 直接モードの C2 壁: 拡大部の 5 次 B-spline (無次元, $r_t=1$)。

計画: plans/active/tooling-nozzle-bell-two-mode.md §4.2・§4.3・§5.1 #3・#3c・#3d。

**規約** (§5.1 #3 で固定、2026-10-11):

- 座標は $r_t$ で無次元、$x$ はスロート基準 (幾何スロート = 0)。
- 拡大部 $x\in[0,L]$ を**等間隔 16 区間の clamped 5 次 B-spline** (端ノット 6 重、内部ノット 15 個、係数 21 個) で表す。
- 係数は 7 個の線形拘束の下で $J=\int_0^L (r''')^2dx$ を最小にする問題を KKT で解いて一意に決める:
  - スロート $x=0$: $r=1$, $r'=0$, $r''=1/R_t$
  - 最大壁角の点 $x_m$: $r'=\tan\theta_{max}$, $r''=0$
  - 出口 $x=L$: $r=\sqrt\varepsilon$, $r'=\tan\theta_e$

  $J$ の零空間は 2 次多項式 ($r'''\equiv0$) で、スロートの 3 拘束がそれを除くので、拘束行列が満行ランク (7) なら
  解は一意 (拘束の零空間の上で $J$ のヘッセ行列が正定値)。ランクと正定値性は構築時に検査する。
- **設計変数は $(L,\theta_e,\theta_{max})$ の 3 個** (§4.3、diagnostician 2026-10-11)。範囲は $\theta_{max}$ 20〜45°、
  $\theta_e$ 2〜20° ($L$ は比較の長さ $L_k$)。範囲外・非有限は `ValueError` (入力の誤り)。
  **$x_m$ は設計変数にせず `design_bell_c2` が内部で決める** (下の「$x_m$ の内部決定」)。成立する $x_m$ は $L$ に比例せず
  絶対値 1.1〜2.1 $r_t$ の帯 ($R_t$ に縛られる) にあり、$x_m/L$ で振ると長い $L$ で帯を外す (V0-a で成立 8/135 の原因)。
- `BellC2Wall` は $x_m$ を与えて 1 本の壁を作る下位の部品。$x_m$ は絶対値 `x_m=` ($0<x_m<L$) か、試験・照合用に残した
  旧来の比 `xm_frac` (0.05〜0.40) のどちらか一方で与える。
- **成立条件** (`feasibility`): 2000 点の等間隔標本で (i) $r'\ge0$、(ii) $r''\ge-tol$ on $[0,x_m]$ かつ
  $r''\le tol$ on $[x_m,L]$ ($tol$ = `RESID_TOL`。$\theta'=r''/(1+r'^2)$ なので $\theta$ が $x_m$ まで非減少・その後は非増加、
  壁角の最大が $x_m$)、(iii) $|\kappa|\le1/R_t$ ($\kappa=r''/(1+r'^2)^{3/2}$)。破れたら INFEASIBLE (違反を文で返す)。
  (ii) は旧条件「$\theta'$ の符号変化がちょうど 1 回」の置き換え (#3d): 旧条件は $x_m$ で $r''$ が符号を変えずに 0 に接する壁
  (例: $L$ 60 %・$\theta_{max}$ 20°・$\theta_e$ 10°・$x_m$ 1.2 は $x_m$ 以降に $r''$ が +0.056 まで戻り、最大壁角 30.18°) を通していた。
- **収縮部** (§4.2、両モード共通): 入口直管 (長さ $L_{pipe}$、半径 $r_U=\sqrt{CR}$) → U→T 5 次 Hermite
  (`wall_walldriven.UpstreamThroatPoly`、$(r_U,0,0)\to(1,0,1/R_t)$、長さ $L_U$)。スロートで拡大部と
  $(r,r',r'')$ が一致するので構成上 C2。単調収縮 ⇔ $L_U^2/[R_t(r_U-1)]\le20$。

**$x_m$ の内部決定** (`design_bell_c2`、§4.3・§5.1 #3・#3c・#3d):

1. **走査**: $x_m$ を 0.5〜4.0 $r_t$ の等間隔 141 点 (0.025 $r_t$ 刻み、$x_m<L$ の点だけ) で振り、各点で `BellC2Wall` を作って
   成立条件 (`feasibility`、2000 点の標本) と $J$ を求める。
2. **成立区間**: 成立した走査点の連続した並び (隣の走査点も成立なら同じ区間)。`intervals` の端は成立した走査点そのもの
   (区間の選び方はこの走査点で決める)。
3. **区間の選び方** (複数に分かれる場合): **成立した走査点のうち $J$ が最小の点を含む区間** (= $J$ 最小の区間)。
   $J$ が同じなら $x_m$ の小さい側。
4. **区間の端を詰める** (#3d): 選んだ区間の両端それぞれで、端の走査点 (成立) と隣の走査点 (不成立) の間を**二分法**
   (幅 1e-6 $r_t$ まで) で詰め、成立側の点を区間の端とする (`interval_refined`、記録は `bisect`)。隣が走査範囲の外
   (走査の端) なら詰めない。これで $x_m$ が走査の 0.025 $r_t$ 刻みに量子化されなくなる ($J$ は不成立の側へ下がり続けることが
   多く、旧手順では $x_m$ が区間の端の走査点に留まっていた)。
5. **確定**: その最小の走査点 $x_i$ と、区間内の両隣の走査点 $x_{i\pm1}$ で挟んだ範囲 ($x_i$ が区間の端なら、その側は
   4. で詰めた端まで) を黄金分割 (幅 1e-6 $r_t$ まで、$J$ が同じなら左側を残す) で詰め、評価した点 (範囲内の走査点・詰めた端を
   含む) のうち $J$ 最小の点を採る (**$J$ 最小の規則は #3d の前と同じ**)。採った点が成立を確かめていない黄金分割の点で
   不成立なら (走査点の間の穴)、$J$ の小さい順に次の点へ移る (`golden["fallback"]`・`n_rejected`)。
6. **成立区間が空**なら例外にせず `status = "INFEASIBLE"` と理由 (走査範囲と不成立の内訳) を返す。走査点の間にしか無い
   幅 0.025 $r_t$ 未満の成立区間は拾わない (V0-b の 0.001 刻みの確認で見つかった細い区間 10 点は、すべて (ii) の旧条件の
   穴による接触型で、#3d の条件では不成立)。

同じ入力ならビット同一 (乱数・並列なし、走査点・二分法・黄金分割の手順が固定)。成立区間が走査の端点 (0.5 か 4.0 $r_t$、または
$L$ の手前の最後の走査点) を含むときは `edge_touch` を立てる (帯が走査範囲の外へ続いているかもしれない — 範囲の見直しの合図)。

この壁がそのまま物理壁 (δ\* 補正なし、§4.3)。CFD の runner への組み込みは別項目 (#2b と衝突するため後で行う)。
"""
from __future__ import annotations

from dataclasses import dataclass, field

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
XM_FRAC_RANGE = (0.05, 0.40)     # 旧来の x_m/L 指定 (試験・照合用) の範囲
XM_SCAN_RANGE = (0.5, 4.0)       # x_m の内部決定の走査範囲 [r_t] (§4.3)
XM_SCAN_N = 141                  # 走査点数 (0.025 r_t 刻み)
XM_GOLDEN_TOL = 1e-6             # 黄金分割の打ち切り幅 [r_t]
XM_BISECT_TOL = 1e-6             # 成立区間の端の二分法の打ち切り幅 [r_t] (#3d)
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


def _check_design_vars(L, theta_e_deg, theta_max_deg, who: str) -> tuple:
    """設計変数 (L, θ_e, θ_max) の検査 (非有限・範囲外は ValueError = 入力の誤り)。戻り (L, θ_e, θ_max) の float。"""
    for k, v in (("L", L), ("theta_e_deg", theta_e_deg), ("theta_max_deg", theta_max_deg)):
        if isinstance(v, bool) or not np.isfinite(v):
            raise ValueError(f"{who}: {k} = {v!r} が有限でない")
    L, te, tm = float(L), float(theta_e_deg), float(theta_max_deg)
    if not L > 0.0:
        raise ValueError(f"{who}: L = {L} ≤ 0")
    for name, v, (lo, hi) in (("theta_max_deg", tm, THETA_MAX_RANGE_DEG), ("theta_e_deg", te, THETA_E_RANGE_DEG)):
        if not (lo <= v <= hi):
            raise ValueError(f"{who}: {name} = {v} が設計変数の範囲 [{lo}, {hi}] の外 (plan §5.1 #3)")
    return L, te, tm


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

    def __init__(self, L: float, theta_e_deg: float, theta_max_deg: float, xm_frac: float | None = None,
                 spec: BellC2Spec | None = None, *, x_m: float | None = None):
        if (xm_frac is None) == (x_m is None):
            raise ValueError("BellC2Wall: x_m (絶対値) か xm_frac (試験・照合用の比) のどちらか一方を与える "
                             "(設計の入口は design_bell_c2 — x_m を内部で決める)")
        spec = spec or BellC2Spec()
        spec.check()
        L, te, tm = _check_design_vars(L, theta_e_deg, theta_max_deg, "BellC2Wall")
        if xm_frac is not None:
            if isinstance(xm_frac, bool) or not np.isfinite(xm_frac):
                raise ValueError(f"BellC2Wall: xm_frac = {xm_frac!r} が有限でない")
            xf = float(xm_frac)
            lo, hi = XM_FRAC_RANGE
            if not (lo <= xf <= hi):
                raise ValueError(f"BellC2Wall: xm_frac = {xf} が範囲 [{lo}, {hi}] の外 (旧来の指定、plan §5.1 #3)")
            xm = xf * L
            self.xm_mode = "frac"
        else:
            if isinstance(x_m, bool) or not np.isfinite(x_m):
                raise ValueError(f"BellC2Wall: x_m = {x_m!r} が有限でない")
            xm = float(x_m)
            if not (0.0 < xm < L):
                raise ValueError(f"BellC2Wall: x_m = {xm} が (0, L = {L}) の外")
            xf = xm / L
            self.xm_mode = "abs"
        self.spec = spec
        self.L, self.theta_e_deg, self.theta_max_deg, self.xm_frac = L, te, tm, xf
        self.x_m = xm
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
        r"""成立条件 (§5.1 #3・#3d): 等間隔 n 点の標本で (i) $r'\ge0$、(ii) $r''\ge0$ on $[0,x_m]$ かつ
        $r''\le0$ on $[x_m,L]$ ($\theta$ が $x_m$ まで非減少・その後は非増加、壁角の最大が $x_m$)、(iii) $|\kappa|\le1/R_t$。

        許容差は拘束の残差の上限 (`RESID_TOL`) に揃える: $r'\ge-10^{-10}$、$r''\ge-10^{-10}$ / $\le10^{-10}$、
        $|\kappa|R_t\le1+10^{-10}$ (スロートでは $r'=0$・$\kappa=1/R_t$、$x_m$ では $r''=0$ が拘束そのものなので丸めの分だけ緩める)。
        (ii) は旧条件「$\theta'$ の符号変化がちょうど 1 回」を置き換えたもの (#3d、diagnostician 2026-10-11): 旧条件は
        $x_m$ で $r''$ が符号を変えずに 0 に接し、下流で一度だけ負になる壁 (壁角の最大が $x_m$ の下流で $\theta_{max}$ を超える)
        を通していた。符号変化の回数と位置 (`theta_sign_changes`・`x_sign_change`) は記録として残し、判定には使わない。"""
        x = np.linspace(0.0, self.L, int(n))
        rp = self.r_expansion(x, 1)
        rpp = self.r_expansion(x, 2)
        kap = rpp / (1.0 + rp * rp) ** 1.5
        nchg, at = _sign_changes(rpp)
        kmax = float(np.max(np.abs(kap)) * self.spec.R_t)
        rp_min = float(np.min(rp))
        th = np.degrees(np.arctan(rp))
        # (ii) x_m の前後で r'' の符号 (標本 x = x_m は両側に入れる。0 < x_m < L なので両側とも空でない)
        i_pre = np.flatnonzero(x <= self.x_m)
        i_post = np.flatnonzero(x >= self.x_m)
        k_pre = int(i_pre[int(np.argmin(rpp[i_pre]))])
        k_post = int(i_post[int(np.argmax(rpp[i_post]))])
        rpp_pre_min, rpp_post_max = float(rpp[k_pre]), float(rpp[k_post])
        v = []
        if rp_min < -RESID_TOL:
            v.append(f"r' < 0 (最小 {rp_min:.4g} @ x = {float(x[int(np.argmin(rp))]):.4g})")
        if rpp_pre_min < -RESID_TOL:
            v.append(f"x_m の前で r'' < 0 (最小 {rpp_pre_min:.4g} @ x = {float(x[k_pre]):.4g}、θ が x_m の前で減る)")
        if rpp_post_max > RESID_TOL:
            v.append(f"x_m の後で r'' > 0 (最大 {rpp_post_max:.4g} @ x = {float(x[k_post]):.4g}、θ が x_m の後で増える)")
        if kmax > 1.0 + RESID_TOL:
            v.append(f"|κ| R_t の最大 {kmax:.6g} > 1 (@ x = {float(x[int(np.argmax(np.abs(kap)))]):.4g})")
        return {"ok": not v, "violations": v, "n_sample": int(n), "rp_min": rp_min,
                "rpp_min_before_xm": rpp_pre_min, "rpp_max_after_xm": rpp_post_max,
                "theta_sign_changes": nchg, "x_sign_change": [float(x[i]) for i in at],
                "kappa_Rt_max": kmax, "theta_max_sampled_deg": float(th.max()),
                "x_theta_max_sampled": float(x[int(np.argmax(th))]), "r_max": float(np.max(self.r_expansion(x)))}

    def summary(self) -> dict:
        res = self.constraint_residuals()
        tj = self.throat_jump()
        return {"L": self.L, "theta_e_deg": self.theta_e_deg, "theta_max_deg": self.theta_max_deg,
                "xm_frac": self.xm_frac, "x_m": self.x_m, "xm_mode": self.xm_mode, "spec": {**self.spec.__dict__, "r_U": self.spec.r_U,
                                                                    "r_exit": self.spec.r_exit},
                "contraction_mu": float(self.up.mu), "J_rppp2": self.J, "constraint_rank": self.constraint_rank,
                "constraint_resid_max": float(np.max(np.abs(res))), "kkt_cond": self.kkt_cond,
                "reduced_hessian_eig": list(self.reduced_hessian_eig),
                "throat_jump_r_rp_rpp": [tj[0], tj[1], tj[2]], "throat_jump_rppp": tj[3]}


# ---------------------------------------------------------------------------------------------- x_m の内部決定
def _violation_kind(v: str) -> str:
    """`feasibility` の違反文を種類に (集計用)。"""
    if v.startswith("r' < 0"):
        return "r' < 0"
    if v.startswith("x_m の前で r'' < 0"):
        return "x_m の前で r'' < 0"
    if v.startswith("x_m の後で r'' > 0"):
        return "x_m の後で r'' > 0"
    if "|κ|" in v:
        return "|κ| > 1/R_t"
    return v


@dataclass
class BellC2Design:
    r"""`design_bell_c2` の戻り。`status` は `"OK"` か `"INFEASIBLE"` (成立区間が空、例外にしない)。

    `wall` は確定した壁 (`INFEASIBLE` なら None)、`reason` は不成立の理由 (`OK` なら空)。走査の記録
    (`scan_x`・`scan_ok`・`scan_J`・`scan_violations`)、成立区間 (`intervals`: (始点 x, 終点 x, 走査点数, 区間内の J の最小))、
    選んだ区間の番号 (`chosen_interval`)、選んだ区間の端を二分法で詰めた範囲 (`interval_refined`: (始点 x, 終点 x)) と
    その記録 (`bisect`: 端 `lo`・`hi` ごとに詰めた成立側 `x_ok`・不成立側 `x_bad`・評価回数・不成立側の違反の種類 `binding`、
    走査の端で隣が無ければ None)、黄金分割の記録 (`golden`)、走査の端に触れたか (`edge_touch`) を持つ。"""
    status: str
    reason: str
    L: float
    theta_e_deg: float
    theta_max_deg: float
    spec: BellC2Spec
    x_m: float | None = None
    wall: BellC2Wall | None = None
    scan_x: np.ndarray = field(default_factory=lambda: np.zeros(0))
    scan_ok: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    scan_J: np.ndarray = field(default_factory=lambda: np.zeros(0))
    scan_violations: dict = field(default_factory=dict)
    intervals: list = field(default_factory=list)
    chosen_interval: int | None = None
    edge_touch: bool = False
    interval_refined: tuple | None = None
    bisect: dict = field(default_factory=dict)
    golden: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "OK"

    def record(self) -> dict:
        """CSV 1 行ぶんの要約。"""
        f = self.wall.feasibility() if self.wall is not None else None
        return {"design_status": self.status, "design_reason": self.reason,
                "x_m": self.x_m, "xm_over_L": (self.x_m / self.L) if self.x_m is not None else None,
                "J_rppp2": self.wall.J if self.wall is not None else None,
                "n_scan": int(len(self.scan_x)), "n_scan_ok": int(np.sum(self.scan_ok)),
                "intervals": ";".join(f"{a:.3f}-{b:.3f}" for (a, b, _, _) in self.intervals),
                "n_intervals": len(self.intervals), "chosen_interval": self.chosen_interval,
                "edge_touch": self.edge_touch,
                "interval_refined": (f"{self.interval_refined[0]:.6f}-{self.interval_refined[1]:.6f}"
                                     if self.interval_refined is not None else None),
                "bisect_n_eval": (sum(v["n_eval"] for v in self.bisect.values() if v) if self.interval_refined else None),
                "binding_lo": (self.bisect.get("lo") or {}).get("binding"),
                "binding_hi": (self.bisect.get("hi") or {}).get("binding"),
                "golden_n_eval": self.golden.get("n_eval"),
                "golden_fallback": self.golden.get("fallback"), "golden_at_scan_point": self.golden.get("at_scan_point"),
                "golden_at_refined_end": self.golden.get("at_refined_end"),
                "rp_min": f["rp_min"] if f else None, "kappa_Rt_max": f["kappa_Rt_max"] if f else None,
                "x_theta_max_sampled": f["x_theta_max_sampled"] if f else None,
                "theta_max_sampled_deg": f["theta_max_sampled_deg"] if f else None}


def _golden_min(fun, a: float, b: float, tol: float) -> list:
    r"""$[a,b]$ の黄金分割 (幅 ≤ tol まで、$f(c)\le f(d)$ なら左を残す)。評価した (x, f) の列を返す (決定的)。"""
    invphi = 0.5 * (np.sqrt(5.0) - 1.0)
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, fd = fun(c), fun(d)
    ev = [(c, fc), (d, fd)]
    while (b - a) > tol:
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = fun(c)
            ev.append((c, fc))
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = fun(d)
            ev.append((d, fd))
    return ev


def design_bell_c2(L: float, theta_e_deg: float, theta_max_deg: float,
                   spec: BellC2Spec | None = None) -> BellC2Design:
    r"""(D) の設計の入口: 設計変数 $(L,\theta_e,\theta_{max})$ から $x_m$ を内部で決めて壁を作る (§4.3、手順はモジュールの
    docstring「$x_m$ の内部決定」)。設計変数の非有限・範囲外は `ValueError` (入力の誤り)。成立区間が空なら例外にせず
    `status = "INFEASIBLE"` と理由を返す。"""
    spec = spec or BellC2Spec()
    spec.check()
    L, te, tm = _check_design_vars(L, theta_e_deg, theta_max_deg, "design_bell_c2")
    out = BellC2Design(status="INFEASIBLE", reason="", L=L, theta_e_deg=te, theta_max_deg=tm, spec=spec)
    xs_all = np.linspace(XM_SCAN_RANGE[0], XM_SCAN_RANGE[1], XM_SCAN_N)
    xs = xs_all[xs_all < L]
    oks = np.zeros(len(xs), dtype=bool)
    Js = np.full(len(xs), np.nan)
    viol: dict = {}
    for i, x in enumerate(xs):
        try:
            w = BellC2Wall(L, te, tm, spec=spec, x_m=float(x))
        except ValueError as e:                      # 拘束の階数・正定値性の破れ (入力の誤りは上で検査済み)
            viol["構築の失敗"] = viol.get("構築の失敗", 0) + 1
            viol.setdefault("_first_error", str(e))
            continue
        f = w.feasibility()
        Js[i] = w.J
        oks[i] = f["ok"]
        for v in f["violations"]:
            k = _violation_kind(v)
            viol[k] = viol.get(k, 0) + 1
    out.scan_x, out.scan_ok, out.scan_J, out.scan_violations = xs, oks, Js, viol
    # --- 成立区間 (成立した走査点の連続した並び) ---
    runs = []
    i = 0
    while i < len(xs):
        if oks[i]:
            j = i
            while j + 1 < len(xs) and oks[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    out.intervals = [(float(xs[a]), float(xs[b]), int(b - a + 1), float(np.min(Js[a:b + 1]))) for (a, b) in runs]
    if not runs:
        lo, hi = XM_SCAN_RANGE
        detail = "・".join(f"{k} {v} 点" for k, v in viol.items() if not k.startswith("_"))
        out.reason = (f"x_m {lo:g}〜{hi:g} r_t の走査 {len(xs)} 点 (x_m < L) に成立点なし"
                      + (f" — 不成立の内訳 (重複あり): {detail}" if detail else ""))
        return out
    # --- 区間の選択: 成立した走査点で J 最小 (同値は x の小さい側) を含む区間 ---
    Jok = np.where(oks, Js, np.inf)
    i_best = int(np.argmin(Jok))
    ci = next(c for c, (a, b) in enumerate(runs) if a <= i_best <= b)
    a, b = runs[ci]
    out.chosen_interval = ci
    out.edge_touch = any(a_ == 0 or b_ == len(xs) - 1 for (a_, b_) in runs)

    def probe(xm):
        """x_m で壁を作り (成立, J, 違反の種類) を返す。構築の失敗は不成立。"""
        try:
            w = BellC2Wall(L, te, tm, spec=spec, x_m=float(xm))
        except ValueError:
            return False, float("inf"), ["構築の失敗"]
        f = w.feasibility()
        return bool(f["ok"]), float(w.J), sorted({_violation_kind(v) for v in f["violations"]})

    # --- 選んだ区間の端を二分法で詰める (#3d): 端の走査点 (成立) と隣の走査点 (不成立) の間 ---
    verified = {float(xs[k]): float(Js[k]) for k in range(a, b + 1)}     # 成立を確かめた x_m → J
    bis = {}
    lo_x, hi_x = float(xs[a]), float(xs[b])
    for side, k_ok, k_bad in (("lo", a, a - 1), ("hi", b, b + 1)):
        if not (0 <= k_bad < len(xs)):                                   # 走査の端 (隣が無い) は詰めない
            bis[side] = None
            continue
        x_ok, J_ok, x_bad, n_b = float(xs[k_ok]), float(Js[k_ok]), float(xs[k_bad]), 0
        kinds_bad = None
        while abs(x_bad - x_ok) > XM_BISECT_TOL:
            xm = 0.5 * (x_ok + x_bad)
            ok_m, J_m, kinds_m = probe(xm)
            n_b += 1
            if ok_m:
                x_ok, J_ok = xm, J_m
            else:
                x_bad, kinds_bad = xm, kinds_m
        if kinds_bad is None:                                            # 一度も詰めて不成立にならなかった
            kinds_bad = probe(x_bad)[2]
        verified[x_ok] = J_ok
        bis[side] = {"x_scan_ok": float(xs[k_ok]), "x_scan_bad": float(xs[k_bad]), "x_ok": x_ok, "x_bad": x_bad,
                     "J_ok": J_ok, "n_eval": n_b, "binding": ";".join(kinds_bad)}
        if side == "lo":
            lo_x = x_ok
        else:
            hi_x = x_ok
    out.bisect, out.interval_refined = bis, (lo_x, hi_x)
    # --- 黄金分割: 最小の走査点と区間内の両隣で挟む。区間の端なら詰めた端まで広げる ---
    lo_b = lo_x if i_best == a else float(xs[i_best - 1])
    hi_b = hi_x if i_best == b else float(xs[i_best + 1])
    cand = [(x_, J_) for x_, J_ in verified.items() if lo_b <= x_ <= hi_b]
    n_eval = 0
    if hi_b > lo_b:
        def Jfun(xm):
            try:
                return float(BellC2Wall(L, te, tm, spec=spec, x_m=float(xm)).J)
            except ValueError:
                return float("inf")
        ev = _golden_min(Jfun, lo_b, hi_b, XM_GOLDEN_TOL)
        n_eval = len(ev)
        cand += [(float(x_), float(J_)) for x_, J_ in ev]
    cand.sort(key=lambda t: (t[1], t[0]))
    # J の小さい順に、成立を確かめた点か、確かめて成立した点を採る (黄金分割の点が区間の穴に落ちた場合の保険)
    x_star, n_rej = None, 0
    for x_, _ in cand:
        if x_ in verified or probe(x_)[0]:
            x_star = x_
            break
        n_rej += 1
    if x_star is None:                                                   # 起こらない (走査点は成立を確認済み)
        x_star = float(xs[i_best])
    wall = BellC2Wall(L, te, tm, spec=spec, x_m=x_star)
    at_scan = bool(np.any(xs[a:b + 1] == x_star))
    at_end = any(s is not None and s["x_ok"] == x_star and x_star != s["x_scan_ok"] for s in bis.values())
    out.golden = {"bracket": (lo_b, hi_b), "x_scan_best": float(xs[i_best]),
                  "J_scan_best": float(Js[i_best]), "n_eval": n_eval, "fallback": n_rej > 0, "n_rejected": n_rej,
                  "at_scan_point": at_scan, "at_refined_end": at_end}
    out.status, out.x_m, out.wall = "OK", float(x_star), wall
    return out
