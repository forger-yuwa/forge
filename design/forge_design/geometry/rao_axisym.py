r"""③ベル (M) モード: 凍結 TP (と CPG) の軸対称 Rao 輪郭の生成器 (S1)。

計画: plans/active/tooling-nozzle-bell-two-mode.md §4.4 (S1 の導出)・§5.1 #5・§6 V1 (T1–T6)。導出と構成手順は
diagnostician の記録 notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md。座標は $r_t$ で無次元 ($x$ はスロート基準)、
状態は音速単位 (`gas/isentrope.py`: $V/a^*$、$\rho/\rho^*$、$p/(\rho^*a^{*2})$)。推力・流量は $\rho^*a^{*2}r_t^2$・$\rho^*a^*r_t^2$ の単位
($2\pi$ を含む)。

**条件** (凍結 TP でも古典と同形、比熱比は陽に現れない): 制御面 DE は C⁺ で、DE 上で
$\Phi_2=V\cos(\theta-\mu)/\cos\mu$ と $\Phi_3=r\rho V^2\sin^2\theta\tan\mu$ が一定。固定 $(L,\varepsilon)$ では縁条件は課さない
(端点項が消える。`edge_function` $=r_E[(p_E-p_a)-\tfrac12\rho_EV_E^2\tan\mu_E\sin2\theta_E]$ は $\partial F^*/\partial r_E|_L/2\pi$
として報告)。自由 $\varepsilon$ (T1c) では縁条件 `edge_function` $=0$ を課す。

**構成** (`design_rao`):

1. kernel (`kernel_wall`): 曲率 $\kappa(s)=1/R_t+c\,s$ ($c$ = 曲率の勾配、$c=0$ が円弧 = 古典で (M) の既定)。記録の族
   $\kappa(s)=1/R_t+(\kappa_B-1/R_t)\,s/s_B$ と同じ集合で、$(\kappa_B,\theta_B)\leftrightarrow(c,s_B)$ は
   $\kappa_B=1/R_t+cs_B$、$\theta_B=s_B/R_t+cs_B^2/2$ の 1 対 1 ($c$ でパラメータ化するのは、同じ $c$ の中で B を動かせるので
   kernel 場が 1 本で済むため。$c\ne0$ は記録用に明示して与える — 外の根は撤回、下の 5)。場を作るために B の先は**同じ曲線の延長** (曲率が連続、BD の下流側も滑らか — BD の上流側の
   片側極限は延長に依らない) を `s_lin` まで続け、その先は C2 の 5 次式で直線 $\theta_{st}$ へつないでどこまでも直線
   (`KernelWall` の docstring: 前進 MOC の壁の交点の検索と「壁の外」の検査に合わせた形。$\theta_{st}$ は初期線の壁点の
   マッハ角 − 5° 以下)。
2. kernel 場 (`kernel_field`): `moc_forward.forward_moc` (`interior_vec` analytic・converge)。前線 $k$ (C⁻) × 点 $j$ (0 = 軸、$n$ = 壁) の
   網を、列 $j$ ごとの $k$ 方向の 3 次スプラインで連続化する (`KernelNet`)。網の論理座標 $(\xi,\eta)=(k,j)$ は特性座標で、
   $\xi$ 一定が C⁻、$(\xi,\eta)$ が同じだけ増える方向が C⁺。
3. BD (`trace_bd`): B ($s_B$) から下流へ出る C⁻ = $\xi=\xi_B$ の線 ($\xi_B$ は壁の節点の弧長 $s(k)$ のスプラインの逆)。
   値は列のスプラインを $\xi_B$ で評価し、$\eta$ 方向にもう一度 3 次スプライン (`BDLine`)。BD に沿った流量・推力は
   区間ごとの Gauss–Legendre。
4. DE (`rao_line`、2026-10-11 diagnostician の訂正で **$r$ をパラメータに**): D の状態から $\Phi_2,\Phi_3$ を決め、各 $r$ で
   2 式 ($\Phi_2$・$\Phi_3$ 一定) を $(M,\theta)$ について Newton で解いて $r$ の増える向き (C⁺ の下流) に引き、流量が
   $\dot m_{BD}(D\to B)$ に達した点を E。$\theta=\mu$ では $M$ が極値をとる ($M$ のパラメータ化の折り返し) が、$r$ では正則で
   $d\theta/dr=-\tan\theta/(2r)$ (旧版の「$|c|>1$・$dr/dM$ の符号反転は解なし」は撤回 — $c=\cos(\theta-\mu)\le1$ は恒等的)。
   解なしは $\theta\le0$・$\theta+\mu\ge90°$・$M$ が表の外 (TP は `M_max` 16)・`M_cap` まで流量に達しない、の 4 条件。
   陰関数の行列式 $\det$ の符号は D で決まり、DE の族を分ける: $\det_D>0$ は $\theta-\mu$ が減る向き ($\theta_D>\mu_D$ なら
   θ = μ を跨ぐ)、$\det_D<0$ ($\theta_D$ が $\mu_D$ より十分大きい D) は $M$ と $\theta-\mu$ がともに増える向き。$\det_D=0$ を
   跨いで D を動かすと E が跳ぶ (内の根の検算で跳びを根と取り違えない)。$\det_D\approx0$ の D では線の出だしが $r$ に対して急
   なので、区間幅を半分ずつ (`SUB_MIN` まで) 縮めて引く。積分は $r$ 格子 (`dr`、軸の近くでは $r$ に比例) の区間ごとの
   Gauss–Legendre (T2 の「格子を倍に」の格子)。
5. 入れ子の根: 内 `solve_D` ($\eta_D$: $r_E=\sqrt\varepsilon$、自由 $\varepsilon$ では縁条件)、中 `solve_thetaB` ($s_B$: $x_E=L$)。
   すべて走査で符号変化を全部拾い、brentq で詰め、根で残差を検算する (値の跳びは根にしない)。**中の根の枝は
   $\eta_D$ の連続性で追い** (添字では結ばない)、**受理条件は $\det_D>0$** (DE が D で流れを軸向きに曲げ始める族 = 古典の枝の
   延長。2026-10-11 diagnostician、plan §4.4「主枝の受理規則」)。$\det_D\le0$ の候補は理由つきで捨てて記録する。旧条件
   $\theta_E\le\theta_B$ と $de=(d\theta/d\sigma)_{DE,D}<0$ は候補ごとのフラグ (記録) で、受理には使わない。$\det_D>0$ なのに
   $\theta_E>\theta_B$ の候補は `policy_holes` に出す (方針の穴として諮る)。候補が複数なら MOC 推力 (制御面の会計) が最大のもの。
   **外の根 (C2) は撤回** (2026-10-11 diagnostician、plan §4.4「(M) は円弧 kernel・B で C1」): $\kappa_B$ を振る族には
   $\tilde g=(d\theta/d\sigma)_{DE,D}-(d\theta/d\sigma)_{kernel\,C^+,D}=0$ の根が無い。`solve_kappaB` は $c$ を明示パラメータとして
   受け取って中の根を解くだけ (走査も brentq もしない)。既定は円弧 ($c=0$)、$\tilde g$ と曲率の跳び
   $\Delta\kappa_B=\kappa_W(B^+)-\kappa_B$ (円弧では $\kappa_W(B^+)-1/R_t$) は記録。
6. 領域 BDE を Goursat 問題で埋め (`goursat_fill`、`interior_vec` を符号付きで — C⁻ を逆向きに積分)、壁 BE は各 C⁺ 線上で
   流量が $\dot m_{BD}(b_i\to B)$ に達する点 (`wall_from_massflow`、`moc_inverse.cplus_flux_wall` と同じ閉包)。
7. 推力の 2 会計 (`thrust_accounts`): (a) 制御面 = (軸 → D) + DE、(b) 初期線 + kernel 壁 + BE 壁。(a) の軸 → D は BD の代わりに
   D を通る kernel 場の C⁺ を遡って閉じる (`cplus_closure`。θ_B > μ_B の BD は出口よりはるかに下流 [θ_B 36° で x ≈ 72 r_t、
   軸の M ≈ 10] で軸に着くので網が粗く、BD と C⁺ の間には壁が無いので運動量は等しい)。BD は B から、kernel 場の列が
   使える深さ (`KernelNet.j_low`) までだけを使う。

**名前**: 記録が「最大推力」と呼ぶのは S1 を通ってから。それまでは「指定した kernel と拘束の下の停留輪郭」。
"""
from __future__ import annotations

import time
import warnings

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from ..gas.isentrope import Isentrope
from .moc_kernel import PAIR_CONVERGED, interior_vec, mu_vec

GL_N = 6
_GLX, _GLW = np.polynomial.legendre.leggauss(GL_N)
TWO_PI = 2.0 * np.pi

# 探索域の既定 (S1 の定義: plan §4.4・§5.1 #4)
THETA_B_RANGE_DEG = (12.0, 50.0)      # B の壁角の探索域
KAPPA_B_MIN = -0.9                    # κ_B ≥ KAPPA_B_MIN / R_t (κ_B → −1/R_t で s_B → ∞)
N_SCAN_SB = 14                       # 中の走査点数
N_SCAN_D = 36                         # 内の走査点数
XTOL = 1e-12
THETA_LIN_MAX_DEG = 75.0              # kernel 場の線形曲率の延長の終わりの壁角の上限
MARGIN_K = 4.0                        # BD の列のスプラインに要る ξ_B より先の前線数
DR_DE = 0.02                          # DE の r 格子の既定 [r_t] (区間ごとの Gauss–Legendre)。区間幅は dr·min(1, r/R_CELL)
R_CELL = 1.0                          # これより軸寄りでは区間幅を r に比例させる (D が軸の近くでも相対幅 dr) [r_t]
NEWTON_TOL = 1e-13                    # DE の各 r での (M, θ) の Newton の収束 (相対 M・絶対 θ)
NEWTON_FTOL = 2e-15                   # 同じく残差 F₁・F₂ (対数の量) での収束 (det が小さい点用)
SUB_MIN = 2.0 ** -20                  # DE の区間幅の縮めの下限 (D が det ≈ 0 の近くで線の出だしが急なとき)
ETA_LINK_FRAC = 0.15                  # 中の根の枝を結ぶ η_D の跳びの上限 (網の列数 n に対する比)
ROOT_RESID_TOL = 1e-8                 # 内・中の根の検算 (根で r_E − √ε・x_E − L が尺度のこの倍以下、超えたら値の跳び)


def _gl(a, b):
    a = np.asarray(a, dtype=float)[..., None]
    b = np.asarray(b, dtype=float)[..., None]
    h = 0.5 * (b - a)
    return 0.5 * (a + b) + h * _GLX, h * _GLW


def as_isentrope(gas) -> Isentrope:
    """ガスを `Isentrope` に (TP の表は M 16 まで — BD の軸寄りは出口より下流で高マッハになる)。"""
    return gas if isinstance(gas, Isentrope) else Isentrope(gas, M_max=16.0)


def default_L_march(L: float) -> float:
    """kernel 場の前進 MOC の打ち切り (既定 2.5 L)。D の探索域 (BD の下端) を決める。"""
    return 2.5 * float(L)


# ================================================================================================
# kernel 壁
class KernelWall:
    r"""kernel 壁 $\theta(s)$ (区分多項式、$s$ は壁の弧長):

    1. $[0,s_{lin}]$: $\theta=s/R_t+c\,s^2/2$ (曲率 $1/R_t+cs$) — kernel 本体と、BD の下流側を滑らかにするための延長
    2. $[s_{lin},s_{lin}+\ell_1]$: $(\theta,\theta',\theta'')$ を両端で合わせた 5 次式で直線 $\theta_{st}$ へ
       ($\theta_{st}=\min(\theta_{lin}+\kappa_{lin}\ell_1/2,\ \theta_{cap})$、$\theta_{cap}$ は初期線の壁点のマッハ角 − 5°)
    3. その先は直線 $\theta_{st}$ (どこまでも)

    2・3 は場を作るためだけのもので、使う前線 ($s\le s_{lin}$ の壁点から出る C⁻) には影響しない。ただし前進 MOC の
    幾何の検査 (新点が壁の外) と壁の交点の検索 (`moc_forward._Wall.intersect` は打ち切り $L$ で壁が C⁺ の下にあることを前提に
    する) に効く: 直線が C⁺ の最小の傾き (スロートの壁点の $\tan\mu^*$) より緩ければ後者を満たし、θ_B > μ_B で上向きに出る
    BD の C⁻ (傾き $\theta-\mu$) より急なら前者を満たす。
    $(x,r)(s)$ は区間ごとの Gauss–Legendre の累積と 3 次 Hermite (微分 $\cos\theta,\sin\theta$ は厳密)。
    `r(x, deriv)` (deriv 0・1) は `moc_forward.forward_moc` の壁の約束。"""

    def __init__(self, R_t: float, c: float = 0.0, s_lin: float = 1.0, l1: float | None = None,
                 theta_cap: float = np.deg2rad(45.0), l_tab: float = 6.0, h_tab: float = 5e-4):
        self.R_t, self.c, self.s_lin = float(R_t), float(c), float(s_lin)
        if not (self.R_t > 0 and self.s_lin > 0):
            raise ValueError("KernelWall: R_t・s_lin は正")
        s0 = self.s_lin
        self.th_lin = s0 / self.R_t + 0.5 * self.c * s0 * s0
        self.k_lin = 1.0 / self.R_t + self.c * s0
        if not (abs(self.th_lin) < np.deg2rad(80.0)):
            raise ValueError(f"KernelWall: θ(s_lin) = {np.rad2deg(self.th_lin):.1f}° (80° 以内)")
        l1 = float(l1) if l1 is not None else 0.25 * self.R_t
        th_st = min(self.th_lin + 0.5 * self.k_lin * l1, float(theta_cap))
        segs = [(0.0, [0.0, 1.0 / self.R_t, 0.5 * self.c, 0.0, 0.0, 0.0]),
                (s0, self._quintic(self.th_lin, self.k_lin, self.c, th_st, 0.0, 0.0, l1)),
                (s0 + l1, [th_st, 0.0, 0.0, 0.0, 0.0, 0.0])]
        self.segs = segs
        self._bp = np.array([sg[0] for sg in segs])
        self._co = np.array([sg[1] for sg in segs])
        self.th_st = float(th_st)
        self.theta_inf = float(th_st)
        s_end = s0 + l1 + float(l_tab)
        self.s_end = s_end
        tt = self.theta(np.linspace(0.0, s_end, 20001))
        if np.max(tt) >= np.deg2rad(85.0):
            raise ValueError(f"KernelWall: 延長で壁角が {np.rad2deg(np.max(tt)):.1f}° になる")
        # 表 (区分の境に節点を置く)。表の先 (s > s_end) は直線 θ_st を解析的に
        nodes = [np.array([0.0])]
        for a, b in zip(self._bp, np.append(self._bp[1:], s_end)):
            nseg = max(int(np.ceil((b - a) / h_tab)), 8)
            nodes.append(np.linspace(a, b, nseg + 1)[1:])
        s = np.concatenate(nodes)
        sg, wg = _gl(s[:-1], s[1:])
        thg = self.theta(sg)
        dx = np.sum(np.cos(thg) * wg, axis=-1)
        dr = np.sum(np.sin(thg) * wg, axis=-1)
        self._s = s
        self._x = np.concatenate([[0.0], np.cumsum(dx)])
        self._r = 1.0 + np.concatenate([[0.0], np.cumsum(dr)])
        th = self.theta(s)
        self._cx, self._cr = np.cos(th), np.sin(th)
        self.x_end, self.r_end = float(self._x[-1]), float(self._r[-1])
        if np.any(np.diff(self._x) <= 0.0):
            raise ValueError("KernelWall: x(s) が単調でない (壁角が ±90° を超える)")

    @staticmethod
    def _quintic(t0, d0, dd0, t1, d1, dd1, ell):
        """Δ ∈ [0, ℓ] で (θ, θ', θ'') が両端 (t0, d0, dd0)・(t1, d1, dd1) の 5 次式の係数 [a0..a5]。"""
        A = np.array([[ell ** 3, ell ** 4, ell ** 5], [3 * ell ** 2, 4 * ell ** 3, 5 * ell ** 4],
                      [6 * ell, 12 * ell ** 2, 20 * ell ** 3]])
        rhs = np.array([t1 - (t0 + d0 * ell + 0.5 * dd0 * ell * ell), d1 - (d0 + dd0 * ell), dd1 - dd0])
        a3, a4, a5 = np.linalg.solve(A, rhs)
        return [t0, d0, 0.5 * dd0, a3, a4, a5]

    def _seg(self, s):
        s = np.asarray(s, dtype=float)
        i = np.clip(np.searchsorted(self._bp, s, side="right") - 1, 0, len(self._bp) - 1)
        return i, s - self._bp[i]

    # --- 壁角と曲率 ---
    def theta(self, s):
        i, d = self._seg(s)
        a = self._co[i]
        return a[..., 0] + d * (a[..., 1] + d * (a[..., 2] + d * (a[..., 3] + d * (a[..., 4] + d * a[..., 5]))))

    def kappa(self, s):
        i, d = self._seg(s)
        a = self._co[i]
        return a[..., 1] + d * (2 * a[..., 2] + d * (3 * a[..., 3] + d * (4 * a[..., 4] + d * 5 * a[..., 5])))

    # --- 位置 ---
    def xr(self, s):
        s = np.atleast_1d(np.asarray(s, dtype=float))
        x = np.empty_like(s)
        r = np.empty_like(s)
        inside = (s >= 0.0) & (s <= self.s_end)
        if np.any(inside):
            si = s[inside]
            i = np.clip(np.searchsorted(self._s, si, side="right") - 1, 0, len(self._s) - 2)
            h = self._s[i + 1] - self._s[i]
            t = (si - self._s[i]) / h
            h00, h10, h01, h11 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t, -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
            x[inside] = h00 * self._x[i] + h10 * h * self._cx[i] + h01 * self._x[i + 1] + h11 * h * self._cx[i + 1]
            r[inside] = h00 * self._r[i] + h10 * h * self._cr[i] + h01 * self._r[i + 1] + h11 * h * self._cr[i + 1]
        beyond = s > self.s_end
        if np.any(beyond):
            ds = s[beyond] - self.s_end
            x[beyond] = self.x_end + ds * np.cos(self.theta_inf)
            r[beyond] = self.r_end + ds * np.sin(self.theta_inf)
        before = s < 0.0                        # スロートの上流 (使わない想定。曲率 1/R_t の円で延長)
        if np.any(before):
            a = s[before] / self.R_t
            x[before] = self.R_t * np.sin(a)
            r[before] = 1.0 + self.R_t * (1.0 - np.cos(a))
        return x, r

    def s_of_x(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        s = np.empty_like(x)
        inside = (x >= 0.0) & (x <= self.x_end)
        if np.any(inside):
            xi = x[inside]
            i = np.clip(np.searchsorted(self._x, xi, side="right") - 1, 0, len(self._s) - 2)
            h = self._s[i + 1] - self._s[i]
            x0, x1, d0, d1 = self._x[i], self._x[i + 1], self._cx[i], self._cx[i + 1]
            t = (xi - x0) / np.maximum(x1 - x0, 1e-300)
            for _ in range(30):
                h00, h10, h01, h11 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t, -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
                f = h00 * x0 + h10 * h * d0 + h01 * x1 + h11 * h * d1 - xi
                df = (6 * t * t - 6 * t) * x0 + (3 * t * t - 4 * t + 1) * h * d0 + (-6 * t * t + 6 * t) * x1 + (3 * t * t - 2 * t) * h * d1
                dt = f / df
                t = t - dt
                if np.all(np.abs(dt) < 1e-15):
                    break
            s[inside] = self._s[i] + t * h
        beyond = x > self.x_end
        if np.any(beyond):
            s[beyond] = self.s_end + (x[beyond] - self.x_end) / np.cos(self.theta_inf)
        before = x < 0.0
        if np.any(before):
            s[before] = self.R_t * np.arcsin(np.clip(x[before] / self.R_t, -1.0, 1.0))
        return s

    def r(self, x, deriv: int = 0):
        s = self.s_of_x(x)
        if deriv == 0:
            return self.xr(s)[1].reshape(np.shape(x))
        if deriv == 1:
            th = np.where(s < 0.0, s / self.R_t, self.theta(np.maximum(s, 0.0)))
            return np.tan(th).reshape(np.shape(x))
        raise ValueError("KernelWall.r: deriv は 0 か 1")


def kernel_wall(R_t: float, *, c: float | None = None, kappa_B: float | None = None, theta_B: float | None = None,
                kind: str = "linear", s_lin: float | None = None) -> KernelWall:
    r"""kernel 壁。`kind="arc"` は円弧 ($c=0$)。`kind="linear"` は $c$ を直接、または記録の $(\kappa_B,\theta_B)$ で与える
    ($s_B=2\theta_B/(1/R_t+\kappa_B)$、$c=(\kappa_B-1/R_t)/s_B$)。`s_lin` は線形曲率の区間の終わり (既定 $s_B$ の 1.5 倍か 1)。"""
    if kind == "arc":
        c = 0.0
    elif kind == "linear":
        if c is None:
            if kappa_B is None or theta_B is None:
                raise ValueError("kernel_wall: linear は c か (kappa_B, theta_B)")
            if not (1.0 / R_t + kappa_B > 0.0):
                raise ValueError(f"kernel_wall: 1/R_t + κ_B = {1.0 / R_t + kappa_B} ≤ 0")
            sB = 2.0 * theta_B / (1.0 / R_t + kappa_B)
            c = (kappa_B - 1.0 / R_t) / sB
            s_lin = s_lin if s_lin is not None else 1.5 * sB
    else:
        raise ValueError("kernel_wall: kind は 'arc' か 'linear'")
    return KernelWall(R_t, c=float(c), s_lin=float(s_lin if s_lin is not None else 1.0))


def _s_of_theta(R_t: float, c: float, theta: float) -> float:
    """θ(s) = s/R_t + c s²/2 が最初に theta に達する s (達しなければ極大の s)。"""
    s = np.linspace(0.0, 20.0 * R_t, 200001)
    th = s / R_t + 0.5 * c * s * s
    i = int(np.argmax(th >= theta)) if np.any(th >= theta) else int(np.argmax(th))
    return float(s[i])


def sB_range(R_t: float, c: float, theta_range_deg=THETA_B_RANGE_DEG, kappa_min=KAPPA_B_MIN) -> tuple:
    r"""$c$ を固定したとき B の弧長 $s_B$ の探索域: $\theta_B\in$ `theta_range`、$\kappa_B\ge$ `kappa_min`$/R_t$、$\theta$ が
    最初に `theta_range` の下限を越えてから。戻り (s_lo, s_hi) (空なら None)。"""
    lo, hi = np.deg2rad(theta_range_deg[0]), np.deg2rad(theta_range_deg[1])
    s = np.linspace(0.0, 20.0 * R_t, 40001)
    th = s / R_t + 0.5 * c * s * s
    kap = 1.0 / R_t + c * s
    ok = (th >= lo) & (th <= hi) & (kap >= kappa_min / R_t)
    if not np.any(ok):
        return None
    i0 = int(np.argmax(ok))
    i1 = i0
    while i1 + 1 < len(s) and ok[i1 + 1]:
        i1 += 1
    return float(s[i0]), float(s[i1])


# ================================================================================================
# kernel 場
class KernelNet:
    r"""前進 MOC の網 (前線 $k$ × 点 $j$) の連続化。列 $j$ ごとに、その列で収束した前線 $k\in[0,k_{max}(j)]$ の範囲で
    $k$ 方向の 3 次スプライン ($x,r,\theta,\nu$ を一括)。$x>$ `L_march` の点は「出口の外」で計算されないので、軸に近い列ほど
    $k_{max}$ が小さい (前線は下流ほど軸の着地点が遠い)。"""

    def __init__(self, res: dict, wall: KernelWall, ig: Isentrope, line: dict):
        from .moc_forward import ST_OK
        net = res["net"]
        n = res["n"]
        S = net.status
        ok = S == ST_OK
        Kc = S.shape[0] - 1
        kmax = np.empty(n + 1, dtype=int)
        for j in range(n + 1):
            bad = np.flatnonzero(~ok[:, j])
            kmax[j] = (bad[0] - 1) if len(bad) else Kc
        self.kmax = kmax
        self.K = int(kmax[n])                      # 壁の節点が続く前線
        if self.K < 8:
            raise RuntimeError(f"KernelNet: 壁の節点が {self.K + 1} 個しかない ({res['end']})")
        self.n = n
        self.A = np.stack([net.a["x"], net.a["r"], net.a["th"], net.a["nu"]], axis=-1)
        self.wall, self.ig, self.line, self.res = wall, ig, line, res
        self.s_wall = wall.s_of_x(net.a["x"][: self.K + 1, n])
        if np.any(np.diff(self.s_wall) <= 0.0):
            raise RuntimeError("KernelNet: 壁の節点の弧長が単調でない")
        self.cols = [None] * (n + 1)
        for j in range(n + 1):
            if kmax[j] >= 4:
                self.cols[j] = CubicSpline(np.arange(kmax[j] + 1, dtype=float), self.A[: kmax[j] + 1, j, :], axis=0)
        self.dcols = [None if c is None else c.derivative() for c in self.cols]
        kk = np.arange(self.K + 1, dtype=float)
        self.s_spl = CubicSpline(kk, self.s_wall)
        self.nu_wall_s = CubicSpline(self.s_wall, net.a["nu"][: self.K + 1, n])

    def xi_of_s(self, s_B: float) -> float:
        k = int(np.searchsorted(self.s_wall, s_B) - 1)
        if k < 0 or k >= self.K:
            raise ValueError(f"xi_of_s: s_B {s_B} が網の壁の範囲 [{self.s_wall[0]}, {self.s_wall[-1]}] の外")
        return float(brentq(lambda xi: float(self.s_spl(xi)) - s_B, k, k + 1, xtol=1e-14, rtol=4 * np.finfo(float).eps))

    def j_low(self, xi: float, margin: float = MARGIN_K) -> int:
        """ξ で列のスプラインが使える最も軸寄りの j (それより上の列はすべて k_max ≥ ξ + margin)。"""
        good = self.kmax >= xi + margin
        if not good[self.n]:
            raise ValueError(f"j_low: ξ {xi} で壁の列も使えない")
        bad = np.flatnonzero(~good)
        return int(bad[-1] + 1) if len(bad) else 0

    def eval_cols(self, xi: float, js) -> tuple:
        """列 js を ξ で評価: (値 (len, 4), ∂_ξ (len, 4))。"""
        js = list(js)
        v = np.array([self.cols[j](xi) for j in js])
        d = np.array([self.dcols[j](xi) for j in js])
        return v, d

    def eval_point(self, xi, eta) -> tuple:
        r"""任意の $(\xi,\eta)$ での値と $\partial_\xi$・$\partial_\eta$ ($\eta$ は 4 点の 3 次 Lagrange、$\xi$ は列のスプライン)。"""
        j0 = int(np.clip(np.floor(eta) - 1, 0, self.n - 3))
        js = np.arange(j0, j0 + 4)
        v, d = self.eval_cols(float(xi), js)
        w = np.empty(4)
        dw = np.empty(4)
        for a in range(4):
            others = [b for b in range(4) if b != a]
            den = np.prod([js[a] - js[b] for b in others])
            w[a] = np.prod([float(eta) - js[b] for b in others]) / den
            dw[a] = sum(np.prod([float(eta) - js[b] for b in others if b != c]) for c in others) / den
        return w @ v, w @ d, dw @ v



def kernel_field(ig, R_t: float, c: float, line: dict, s_need: float, *, L_march: float, n_after: int = 8,
                 K_guess: int | None = None) -> KernelNet:
    r"""$c$ の kernel 場。$s\le$ `s_need` の B を扱えるよう、壁の節点が `s_need` を越えてから `n_after` 本の前線まで
    線形曲率の区間 (`s_lin`) に入るよう、`s_lin` と前線数を必要に応じて広げて `forward_moc` を回し直す。
    `L_march`: 前進 MOC の打ち切り ($x>L$ の点は出口の外として計算しない)。設計の $L$ より十分長くとる (D の探索域を BD の
    下の方まで確保するため)。"""
    from .moc_forward import forward_moc
    ig = as_isentrope(ig)
    n = len(line["x"]) - 1
    th_cap = float(np.arcsin(1.0 / float(line["M"][-1]))) - np.deg2rad(5.0)
    s_cap = _s_of_theta(R_t, c, np.deg2rad(THETA_LIN_MAX_DEG))       # s_lin の上限 (θ_lin ≤ 75°)
    s_need = min(float(s_need), s_cap - 0.02 * R_t)
    s_lin = min(s_need + 0.06 * R_t, s_cap)
    K = K_guess if K_guess is not None else n
    last = None
    for _ in range(12):
        wall = KernelWall(R_t, c=c, s_lin=s_lin, theta_cap=th_cap)
        res = forward_moc(wall, L_march, ig.moc, line, max_fronts=K)
        try:
            kn = KernelNet(res, wall, ig, line)
        except RuntimeError as e:
            last = str(e)
            K = int(K * 1.5)
            continue
        sw = kn.s_wall
        if sw[-1] <= s_need:
            if res["end"]["reason"] not in ("front_cap",):
                last = f"壁の節点が s_need に届く前に前進 MOC が止まった ({res['end']})"
                break
            K = int(K * 1.6) + 10
            continue
        k_need = int(np.searchsorted(sw, s_need))
        if k_need + n_after > kn.K:
            if res["end"]["reason"] not in ("front_cap",):
                last = f"前線 {k_need + n_after} まで届かない ({res['end']})"
                break
            K = int(max(K, k_need + n_after + 5) * 1.2)
            continue
        if sw[k_need + n_after] > s_lin:
            if s_lin >= s_cap - 1e-12:
                # 延長を伸ばせない: B の探索域の上端を下げる (s_need = n_after 本手前の壁点)
                k_need = int(np.searchsorted(sw, s_lin)) - n_after - 1
                if k_need < 4:
                    last = "θ_lin の上限で B の探索域が空"
                    break
                s_need = float(sw[k_need])
            else:
                s_lin = min(float(sw[k_need + n_after]) + 0.02 * R_t, s_cap)
                continue
        kn.k_need = k_need
        kn.s_need = s_need
        kn.c = c
        kn.L_march = L_march
        return kn
    raise RuntimeError(f"kernel_field: c = {c} で s_need = {s_need} までの場が作れない ({last})")


# ================================================================================================
# BD と DE
class BDLine:
    r"""C⁻ 上の状態 ($\eta\in[\eta_0,N]$、N が B。kernel 場由来では $\eta=j$、$\eta_0$ は列が使える最も軸寄りの j) の
    3 次スプライン表現と、$\eta$ に沿った流量・推力 (下 → 上の向き)。
    `vals`: (N−η0+1, 4) [x, r, θ, ν]、`dxi`: 同じ形の $\partial_\xi$ (kernel 場由来のとき。C2 の関数に使う)。"""

    def __init__(self, vals, ig, dxi=None, xi_B=None, s_B=None, eta0: int = 0):
        self.ig = as_isentrope(ig)
        v = np.asarray(vals, dtype=float)
        self.eta0 = int(eta0)
        self.N = self.eta0 + v.shape[0] - 1
        self.eta = np.arange(self.eta0, self.N + 1, dtype=float)
        self.spl = CubicSpline(self.eta, v, axis=0)
        self.dspl = self.spl.derivative()
        self.dxi = None if dxi is None else CubicSpline(self.eta, np.asarray(dxi, dtype=float), axis=0)
        self.xi_B, self.s_B = xi_B, s_B
        self._cum = {}

    def at(self, eta):
        """η での [x, r, θ, ν]。"""
        return self.spl(eta)

    def point(self, eta) -> dict:
        x, r, th, nu = (float(v) for v in self.spl(float(eta)))
        M = float(np.atleast_1d(self.ig.mach_of_nu(np.array([nu])))[0])
        return {"x": x, "r": r, "th": th, "nu": nu, "M": M, "eta": float(eta)}

    def _flux_density(self, eta, p_a):
        v = self.spl(eta)
        d = self.dspl(eta)
        return flux_density(self.ig, v, d, p_a)

    def _cumulative(self, p_a):
        key = float(p_a)
        if key not in self._cum:
            eg, wg = _gl(self.eta[:-1], self.eta[1:])
            dm, dT = self._flux_density(eg, p_a)
            self._cum[key] = (np.concatenate([[0.0], np.cumsum(np.sum(dm * wg, axis=-1))]),
                              np.concatenate([[0.0], np.cumsum(np.sum(dT * wg, axis=-1))]))
        return self._cum[key]

    def integrals_from_bottom(self, eta, p_a: float = 0.0) -> tuple:
        """η0 から η までの (流量, 推力)。区間内は部分区間の Gauss–Legendre (η について滑らか)。"""
        cm, cT = self._cumulative(p_a)
        eta = float(np.clip(eta, self.eta0, self.N))
        i = min(int(np.floor(eta)) - self.eta0, self.N - self.eta0 - 1)
        j = self.eta0 + i
        if eta == j:
            return float(cm[i]), float(cT[i])
        eg, wg = _gl(float(j), eta)
        dm, dT = self._flux_density(eg, p_a)
        return float(cm[i] + np.sum(dm * wg)), float(cT[i] + np.sum(dT * wg))

    def integrals_to_B(self, eta, p_a: float = 0.0) -> tuple:
        """η から B (N) までの (流量, 推力)。"""
        cm, cT = self._cumulative(p_a)
        m0, T0 = self.integrals_from_bottom(eta, p_a)
        return float(cm[-1] - m0), float(cT[-1] - T0)


def flux_density(ig, v, d, p_a):
    r"""曲線 (パラメータ $t$、下 → 上) の $d\dot m/dt$ と $dT/dt$: 値 `v` [x, r, θ, ν] と $t$ 微分 `d` から
    $2\pi r\rho V(\cos\theta\,r_t-\sin\theta\,x_t)$ と $2\pi r[\rho V^2\cos\theta(\cos\theta\,r_t-\sin\theta\,x_t)+(p-p_a)r_t]$。"""
    x, r, th, nu = v[..., 0], v[..., 1], v[..., 2], v[..., 3]
    xe, re = d[..., 0], d[..., 1]
    M = ig.mach_of_nu(nu)
    st = ig.state(np.maximum(np.ravel(M), 1.0))
    shp = np.shape(M)
    rhoV = (st["rho"] * st["V"]).reshape(shp)
    V = st["V"].reshape(shp)
    p = st["p"].reshape(shp)
    nflux = np.cos(th) * re - np.sin(th) * xe
    return TWO_PI * r * rhoV * nflux, TWO_PI * r * (rhoV * V * np.cos(th) * nflux + (p - p_a) * re)


def trace_bd(net: KernelNet, s_B: float) -> BDLine:
    r"""B ($s_B$) から下流へ出る C⁻ (BD) を kernel 場の補間で取る: $\xi_B$ = 壁の弧長のスプラインの逆、
    各列 $j\ge j_{low}$ を $\xi_B$ で評価 (値と $\partial_\xi$)。θ_B > μ_B の作動域では網の前線は BD に使えない (逆向き) ので、
    前線そのものではなく網を通る C⁻ の族の 1 本として取る (記録の構成手順)。"""
    xi = net.xi_of_s(s_B)
    jl = net.j_low(xi)
    v, d = net.eval_cols(xi, range(jl, net.n + 1))
    return BDLine(v, net.ig, dxi=d, xi_B=xi, s_B=s_B, eta0=jl)


def cplus_closure(net: KernelNet, xi_B: float, eta_D: float, p_a: float = 0.0, n_per: int = 1) -> dict:
    r"""D から軸 (または初期線) へ、D を通る kernel 場の C⁺ を遡る曲線 $(\xi_B-t,\eta_D-t)$ に沿った流量と推力
    (軸 → D の向き)。C⁺ が初期線 ($\xi=0$) に先に当たるときは、初期線の軸 → 当たり点を足す。BD (軸 → D) の代わりの閉じ方
    (BD と C⁺ の間に壁は無いので運動量は等しい)。"""
    ig = net.ig
    t_end = min(eta_D, xi_B)
    nseg = max(int(np.ceil(t_end)) * n_per, 1)
    tt = np.linspace(0.0, t_end, nseg + 1)
    tg, wg = _gl(tt[:-1], tt[1:])
    vals = np.empty(tg.shape + (4,))
    ders = np.empty(tg.shape + (4,))
    for idx in np.ndindex(tg.shape):
        t = float(tg[idx])
        v, dxi, deta = net.eval_point(xi_B - t, eta_D - t)
        vals[idx] = v
        ders[idx] = dxi + deta             # 軸 → D の向き (τ = t_end − t) の微分 d/dτ = ∂_ξ + ∂_η
    dm, dT = flux_density(ig, vals, ders, p_a)
    m = float(np.sum(dm * wg))
    T = float(np.sum(dT * wg))
    out = {"mdot": m, "T": T, "t_end": t_end, "hits": "axis" if eta_D <= xi_B else "initial_line"}
    if eta_D > xi_B:
        j_hit = eta_D - xi_B
        ln = net.line
        v0 = np.column_stack([ln["x"], ln["r"], ln["theta"], ig.nu(ln["M"])])
        sp = CubicSpline(np.arange(len(v0), dtype=float), v0, axis=0)
        jj = np.linspace(0.0, j_hit, max(int(np.ceil(j_hit)), 1) + 1)
        jg, wj = _gl(jj[:-1], jj[1:])
        dm2, dT2 = flux_density(ig, sp(jg), sp.derivative()(jg), p_a)
        out["mdot"] += float(np.sum(dm2 * wj))
        out["T"] += float(np.sum(dT2 * wj))
    return out


def edge_function(ig, M, th, r, p_a: float = 0.0) -> float:
    r"""縁の関数 $r[(p-p_a)-\tfrac12\rho V^2\tan\mu\sin2\theta]$ ($=\partial F^*/\partial r_E|_L/2\pi$、0 が自由 $\varepsilon$ の
    縁条件 = Cain 2010 式 4.3)。平面では $r$ を外した同じ式 (根は `rao_planar.theta_lip` と同じ)。"""
    ig = as_isentrope(ig)
    st = ig.state(np.array([float(M)]))
    p, rho, V, mu = (float(st[k][0]) for k in ("p", "rho", "V", "mu"))
    return float(r) * ((p - p_a) - 0.5 * rho * V * V * np.tan(mu) * np.sin(2.0 * float(th)))


class RaoLine:
    r"""制御面 DE (C⁺): D の状態から $\Phi_2,\Phi_3$ 一定で **$r$ をパラメータに**引く (2026-10-11 diagnostician の訂正、
    plan §4.4)。各 $r$ で
    $F_1=\ln[V\cos(\theta-\mu)/\cos\mu]-\ln\Phi_2=0$、$F_2=\ln[r\rho V^2\sin^2\theta\tan\mu]-\ln\Phi_3=0$
    を $(M,\theta)$ について Newton で解く。陰関数の微分 ($a=\partial F_1/\partial M$、$c=\partial F_2/\partial M$):
    $dM/dr=-\tan(\theta-\mu)/(r\det)$、$d\theta/dr=-a/(r\det)$、$\det=2a\cot\theta+c\tan(\theta-\mu)$。
    $\theta=\mu$ では $dM/dr=0$ ($M$ の極値 = $M$ のパラメータ化の折り返し) で $d\theta/dr=-\tan\theta/(2r)$ (C⁺ の適合式と
    $\Phi_3$ 一定が同じ値を与える) — 線は正則に続き、$\arccos$ の枝は $r$ の連続性で選ばれる。

    $\det=-\sin(\theta+\mu)\,(\nu'\tan(\theta-\mu)-a)/(\sin\mu\sin\theta)$ ($\nu'=d\nu/dM$) の符号は D で決まり、
    D から $r$ が増える向き (C⁺ の下流、$dr/d\sigma=\sin(\theta+\mu)>0$) が DE の向き。$\det_D<0$ ($\theta_D$ が $\mu_D$ より
    十分大きい D) では $M$ と $\theta-\mu$ がともに増える向きに進む。

    **解なし** (記録の 4 条件): $\theta\le0$、$\theta+\mu\ge90°$、$M$ が表の外 ($M<1$、TP は表の上限 `M_max` 16)、`M_cap` まで
    流量に達しない。以下は起きないはずの数値上の終わりで、起きたら理由として返す (方針の穴の検出用): $\det$ の符号反転
    ($r$ の折り返し)、Newton の不収束、区間内の追跡の不連続。$\det_D\approx0$ の D では出だしの $dM/dr\propto1/\det_D$ が大きく、
    区間幅を `SUB_MIN` 近くまで縮めて引く (2026-10-11 の基準ケース θ_B 36° で $|\det_D|$ 0.0013 まで、独立な擬似弧長の積分と
    E が 4 桁一致)。

    積分は $r$ 格子の区間ごとの Gauss–Legendre (区間幅 `dr`$\cdot\min(1,r/$`R_CELL`$)$ — 軸の近くの D では相対幅 `dr`。
    Newton が予測子から収束しない区間は区間幅を半分ずつ [`SUB_MIN` まで] にしてやり直す。被積分 $d\dot m/dr=2\pi r\rho V\sin\mu/\sin(\theta+\mu)$、
    $dx/dr=1/\tan(\theta+\mu)$、$dT/dr=2\pi r[\rho V^2\cos\theta\sin\mu/\sin(\theta+\mu)+(p-p_a)]$)。T2 の「格子を倍に」は `dr`。"""

    def __init__(self, ig, D: dict, dr: float = DR_DE, p_a: float = 0.0):
        self.ig = as_isentrope(ig)
        self.D = dict(D)
        self.dr = float(dr)
        self.p_a = float(p_a)
        MD, thD, rD = float(D["M"]), float(D["th"]), float(D["r"])
        st = self.ig.state(np.array([MD]))
        V, rho, mu = float(st["V"][0]), float(st["rho"][0]), float(st["mu"][0])
        self.mu_D = mu
        self.Phi2 = V * np.cos(thD - mu) / np.cos(mu)
        self.Phi3 = rD * rho * V * V * np.sin(thD) ** 2 * np.tan(mu)
        self.why_D = None
        if not thD > 0.0:
            self.why_D = "θ ≤ 0 (D が軸上か壁角が負)"
        elif not thD + mu < 0.5 * np.pi:
            self.why_D = "θ+μ ≥ 90° (D)"
        self.ok_D = self.why_D is None
        self.M_tab = self.ig.M_max if not self.ig.cpg else np.inf
        if self.ok_D:
            self._lnP2, self._lnP3 = np.log(self.Phi2), np.log(self.Phi3)
            fD = self._eval(np.array([MD]), np.array([thD]), np.array([rD]))
            self.det_D = float(fD["det"][0])
            self.sgn = 1.0 if self.det_D > 0.0 else -1.0
            self._fD = fD
        self._bnd = None

    # ------------------------------------------------------------------ 点の式
    def _eval(self, M, th, r) -> dict:
        """(M, θ, r) での残差 F₁・F₂、Jacobian、r 微分と被積分 (配列)。"""
        M, th, r = (np.asarray(v, dtype=float) for v in (M, th, r))
        st = self.ig.state(M.ravel())
        V = st["V"].reshape(M.shape)
        rho = st["rho"].reshape(M.shape)
        p = st["p"].reshape(M.shape)
        mu = st["mu"].reshape(M.shape)
        dlnV = st["dlnV_dM"].reshape(M.shape)
        tmu = np.tan(mu)
        psi = th - mu
        tpsi = np.tan(psi)
        a = dlnV - (tpsi + tmu) * tmu / M                       # ∂F1/∂M (dμ/dM = −tanμ/M)
        b = -tpsi                                               # ∂F1/∂θ
        c = (2.0 - M * M) * dlnV - 1.0 / (M * np.cos(mu) ** 2)  # ∂F2/∂M (d ln ρ/dM = −M² d ln V/dM)
        d = 2.0 / np.tan(th)                                    # ∂F2/∂θ
        det = a * d - b * c
        with np.errstate(divide="ignore", invalid="ignore"):
            F1 = np.log(V) + np.log(np.cos(psi)) - np.log(np.cos(mu)) - self._lnP2
            F2 = np.log(r) + np.log(rho) + 2.0 * np.log(V) + 2.0 * np.log(np.sin(th)) + np.log(tmu) - self._lnP3
            dM_dr = -tpsi / (r * det)
            dth_dr = -a / (r * det)
        spm = np.sin(th + mu)
        return {"M": M, "th": th, "r": r, "mu": mu, "V": V, "rho": rho, "p": p, "dlnV": dlnV,
                "a": a, "b": b, "c": c, "d": d, "det": det, "F1": F1, "F2": F2, "dM_dr": dM_dr, "dth_dr": dth_dr,
                "dx": 1.0 / np.tan(th + mu), "dm": TWO_PI * r * rho * V * np.sin(mu) / spm,
                "dT": TWO_PI * r * (rho * V * V * np.cos(th) * np.sin(mu) / spm + (p - self.p_a))}

    def _solve(self, r, M0, th0, maxit: int = 40) -> dict:
        """各 r で (M, θ) を Newton で解く (予測子 M0・θ0 から、歩幅を制限)。戻りは収束点の `_eval` に `conv` を足したもの。"""
        r = np.asarray(r, dtype=float)
        M = np.array(M0, dtype=float, copy=True)
        th = np.array(th0, dtype=float, copy=True)
        M = np.clip(M, 1.0, self.M_tab)
        conv = np.zeros(r.shape, dtype=bool)
        for _ in range(maxit):
            f = self._eval(M, th, r)
            with np.errstate(divide="ignore", invalid="ignore"):
                dM = (f["d"] * f["F1"] - f["b"] * f["F2"]) / f["det"]
                dth = (f["a"] * f["F2"] - f["c"] * f["F1"]) / f["det"]
                lam = np.minimum(1.0, np.minimum(0.2 * np.maximum(M - 1.0, 0.05) / np.abs(dM), 0.05 / np.abs(dth)))
            fin = np.isfinite(dM) & np.isfinite(dth) & np.isfinite(lam)
            lam = np.where(fin, lam, 0.0)
            M = np.clip(M - lam * np.where(fin, dM, 0.0), 1.0, self.M_tab)
            mu = np.arcsin(1.0 / np.maximum(M, 1.0 + 1e-15))
            th = np.clip(th - lam * np.where(fin, dth, 0.0), 1e-12, np.minimum(np.pi - 1e-9, mu + 0.5 * np.pi - 1e-9))
            # 収束: 全歩幅で歩幅が NEWTON_TOL 以下、または残差が丸めの床 (det が小さい点では歩幅の精度の床が eps/det なので
            # 歩幅だけでは判定できない)
            small_F = (np.abs(f["F1"]) <= NEWTON_FTOL) & (np.abs(f["F2"]) <= NEWTON_FTOL)
            conv = fin & (lam == 1.0) & (((np.abs(dM) <= NEWTON_TOL * M) & (np.abs(dth) <= NEWTON_TOL)) | small_F)
            if np.all(conv | ~fin):
                break
        f = self._eval(M, th, r)
        f["conv"] = conv & (np.abs(f["F1"]) <= 1e-12) & (np.abs(f["F2"]) <= 1e-12)
        return f

    def _classify(self, f, k) -> tuple:
        """点 k の分類: ("ok", None)・("stop", 理由)・("retry", 理由 or None)。不収束の点は予測子・区間幅をやり直す (retry)。
        不収束で表の端・θ の下限に張り付いた点 (Newton は M を [1, 表の上限] に、θ を (0, …) に切り詰める) は、予測子の外れと
        区別できないので retry とし、区間幅を下限まで縮めても張り付いたままなら 4 条件の「解なし」として返す (理由を添える)。"""
        M, th, mu = float(f["M"].flat[k]), float(f["th"].flat[k]), float(f["mu"].flat[k])
        rr = float(f["r"].flat[k])
        bound = None
        if M >= self.M_cap * (1.0 - 1e-12) or M >= self.M_tab * (1.0 - 1e-12):
            bound = f"M_cap {min(self.M_cap, self.M_tab):g} まで流量に達しない (M が表の外、r {rr:.5g})"
        elif M <= 1.0 + 1e-12:
            bound = f"M が表の外 (M ≤ 1、r {rr:.5g})"
        elif th <= 1e-12 * (1.0 + 1e-6):
            bound = f"θ ≤ 0 (r {rr:.5g})"
        if not bool(f["conv"].flat[k]):
            return "retry", bound
        if bound is not None:
            return "stop", bound
        if not th > 0.0:
            return "stop", f"θ ≤ 0 (r {rr:.5g})"
        if not th + mu < 0.5 * np.pi:
            return "stop", f"θ+μ ≥ 90° (r {rr:.5g})"
        if float(f["det"].flat[k]) * self.sgn <= 0.0:
            return "stop", f"r の折り返し (det の符号反転、M {M:.5g}・θ {np.rad2deg(th):.4f}°・r {rr:.5g}) — 4 条件の外"
        return "ok", None

    # ------------------------------------------------------------------ 区間の境の記録と予測子
    def _hermite(self, i, r):
        """境 i と i+1 の (M, θ) と r 微分の 3 次 Hermite で r の (M, θ) を予測する。"""
        R, Mb, Tb, dM, dT = (np.asarray(self._bnd[k]) for k in ("r", "M", "th", "dM", "dth"))
        i = np.asarray(i)
        r0, r1 = R[i], R[i + 1]
        h = r1 - r0
        t = (r - r0) / np.where(h > 0.0, h, 1.0)
        h00, h10, h01, h11 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t, -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
        M = h00 * Mb[i] + h10 * h * dM[i] + h01 * Mb[i + 1] + h11 * h * dM[i + 1]
        th = h00 * Tb[i] + h10 * h * dT[i] + h01 * Tb[i + 1] + h11 * h * dT[i + 1]
        return M, th

    def _taylor(self, r):
        """最後の境から 2 次の Taylor (2 階微分は直前の 2 境の 1 階微分の差) で予測する。"""
        B = self._bnd
        r0, M0, th0, dM0, dth0 = B["r"][-1], B["M"][-1], B["th"][-1], B["dM"][-1], B["dth"][-1]
        if len(B["r"]) >= 2:
            h = B["r"][-1] - B["r"][-2]
            d2M = (B["dM"][-1] - B["dM"][-2]) / h
            d2th = (B["dth"][-1] - B["dth"][-2]) / h
        else:
            d2M = d2th = 0.0
        s = r - r0
        return M0 + s * dM0 + 0.5 * s * s * d2M, th0 + s * dth0 + 0.5 * s * s * d2th

    def _push(self, r, M, th, dM, dth, cm, cx, cT):
        B = self._bnd
        for k, v in (("r", r), ("M", M), ("th", th), ("dM", dM), ("dth", dth), ("cm", cm), ("cx", cx), ("cT", cT)):
            B[k].append(float(v))

    # ------------------------------------------------------------------ 前進
    def march(self, mdot_target: float, M_cap: float | None = None, max_cells: int = 200000) -> dict:
        r"""流量 $\dot m_{DE}(D\to r)$ が `mdot_target` に達する E を求める。戻り dict: `ok`、`why` (失敗の理由)、
        E の `M`・`x`・`r`・`th`、`T` (DE の推力)、`mdot`、`n_cells`、`dth_dr`・`dM_dr`。"""
        self.M_cap = float(M_cap) if M_cap is not None else (self.ig.M_max if not self.ig.cpg else 30.0)
        if not self.ok_D:
            return {"ok": False, "why": f"DE: {self.why_D}"}
        MD, thD, rD = float(self.D["M"]), float(self.D["th"]), float(self.D["r"])
        if not 1.0 < MD < self.M_cap:
            return {"ok": False, "why": f"DE: M が表の外 (D の M {MD:.5g})"}
        fD = self._fD
        self._bnd = {k: [] for k in ("r", "M", "th", "dM", "dth", "cm", "cx", "cT")}
        self._push(rD, MD, thD, fD["dM_dr"][0], fD["dth_dr"][0], 0.0, 0.0, 0.0)
        if mdot_target <= 0.0:
            return {"ok": True, "M": MD, "x": float(self.D["x"]), "r": rD, "th": thD, "T": 0.0, "mdot": 0.0, "n_cells": 0,
                    "dth_dr": float(fD["dth_dr"][0]), "dM_dr": float(fD["dM_dr"][0])}
        chunk = 8
        sub = 1.0                                   # 区間幅の縮め (不収束のとき半分ずつ)
        n = 0
        while n < max_cells:
            B = self._bnd
            r0 = B["r"][-1]
            rb = np.empty(chunk)
            rr = r0
            for k in range(chunk):
                rr = rr + sub * self.dr * min(1.0, rr / R_CELL)
                rb[k] = rr
            rg, wg = _gl(np.concatenate([[r0], rb[:-1]]), rb)
            pts = np.concatenate([rb, rg.ravel()])
            Mp, thp = self._taylor(pts)
            f = self._solve(pts, Mp, thp)
            # 区間 k の可否 (右端と GL 点): 4 条件・det の符号 (stop)、不収束 (retry)
            bad_k, bad_kind, bad_why = chunk, None, None
            for k in range(chunk):
                for q in [k] + list(chunk + k * GL_N + np.arange(GL_N)):
                    kind, why = self._classify(f, q)
                    if kind != "ok":
                        bad_k, bad_kind, bad_why = k, kind, why
                        break
                if bad_kind is not None:
                    break
            # 受理した区間を積む (境の右端の値・微分、累積)。区間の GL 点が両端の Hermite と 1e-5 以上ずれたら追跡の跳び
            # (別の解へ移った) とみなして retry
            for k in range(bad_k):
                gi = chunk + k * GL_N + np.arange(GL_N)
                dm = float(np.sum(f["dm"][gi] * wg[k]))
                cm = B["cm"][-1]
                if cm + dm >= mdot_target:
                    return self._finish(float(rb[k]), f["M"][k], f["th"][k], f["dM_dr"][k], f["dth_dr"][k], dm, mdot_target)
                self._push(rb[k], f["M"][k], f["th"][k], f["dM_dr"][k], f["dth_dr"][k], cm + dm,
                           B["cx"][-1] + float(np.sum(f["dx"][gi] * wg[k])), B["cT"][-1] + float(np.sum(f["dT"][gi] * wg[k])))
                Mh, thh = self._hermite(len(B["r"]) - 2, rg[k])
                jump = max(float(np.max(np.abs(Mh - f["M"][gi]))), float(np.max(np.abs(thh - f["th"][gi]))))
                if jump > 1e-5:
                    for kk in B:
                        B[kk].pop()
                    bad_k, bad_kind, bad_why = k, "retry", None
                    break
                n += 1
            if bad_k == chunk:
                chunk = min(2 * chunk, 64)
                sub = min(1.0, 2.0 * sub)
                continue
            if bad_kind == "retry":
                if bad_k > 0:
                    chunk = max(bad_k, 1)
                    continue
                if chunk > 1:
                    chunk = max(chunk // 2, 1)
                    continue
                if sub > SUB_MIN:
                    sub *= 0.5
                    continue
                if bad_why is not None:
                    return {"ok": False, "why": f"DE: {bad_why}", "n_cells": n}
                return {"ok": False, "why": (f"DE: r での追跡が収束しない (r {r0:.5g}、区間幅 {SUB_MIN:.1e}·dr まで縮めた、"
                                             f"det_D {self.det_D:+.3g}) — 4 条件の外"), "n_cells": n}
            return {"ok": False, "why": f"DE: {bad_why}", "n_cells": n}
        return {"ok": False, "why": "DE: 流量に達しない (セル数の上限)", "n_cells": n}

    def _finish(self, rb, Mb, thb, dMb, dthb, dm_cell, target) -> dict:
        """最後の区間 [r_a, r_b] の中で累積流量が target になる r_E を Newton (括り付き) で求める。"""
        B = self._bnd
        self._push(rb, Mb, thb, dMb, dthb, np.nan, np.nan, np.nan)          # Hermite の予測子のための仮の右端
        ia = len(B["r"]) - 2
        ra, cm0, cx0, cT0 = B["r"][ia], B["cm"][ia], B["cx"][ia], B["cT"][ia]
        lo, hi = ra, rb
        r = ra + (target - cm0) / dm_cell * (rb - ra)
        f = None
        for _ in range(60):
            rg, wg = _gl(ra, r)
            pts = np.append(rg, r)
            Mp, thp = self._hermite(np.full(pts.shape, ia), pts)
            f = self._solve(pts, Mp, thp)
            m = cm0 + float(np.sum(f["dm"][:GL_N] * wg))
            g = m - target
            if g > 0.0:
                hi = r
            else:
                lo = r
            # 止める: 流量の残差が丸めの床 (目標の 8 ulp)、または r の更新・括りが 1e-14 r 以下
            if abs(g) <= 8.0 * np.finfo(float).eps * abs(target) or hi - lo <= 1e-14 * r:
                break
            r_new = r - g / float(f["dm"][GL_N])
            if not lo < r_new < hi:
                r_new = 0.5 * (lo + hi)
            if abs(r_new - r) <= 1e-14 * r:
                r = r_new
                break
            r = r_new
        rg, wg = _gl(ra, r)
        pts = np.append(rg, r)
        Mp, thp = self._hermite(np.full(pts.shape, ia), pts)
        f = self._solve(pts, Mp, thp)
        dm = float(np.sum(f["dm"][:GL_N] * wg))
        dx = float(np.sum(f["dx"][:GL_N] * wg))
        dT = float(np.sum(f["dT"][:GL_N] * wg))
        for k in B:
            B[k].pop()
        ME, thE = float(f["M"][GL_N]), float(f["th"][GL_N])
        self._push(r, ME, thE, f["dM_dr"][GL_N], f["dth_dr"][GL_N], cm0 + dm, cx0 + dx, cT0 + dT)
        self.M_E, self.r_E = ME, r
        self.n_cells = len(B["r"]) - 1
        return {"ok": True, "M": ME, "x": float(self.D["x"]) + cx0 + dx, "r": r, "th": thE, "T": cT0 + dT,
                "mdot": cm0 + dm, "n_cells": self.n_cells, "dth_dr": float(f["dth_dr"][GL_N]),
                "dM_dr": float(f["dM_dr"][GL_N]), "conv": bool(np.all(f["conv"]))}

    # ------------------------------------------------------------------ 前進の後の問い合わせ
    def state_at(self, rs) -> dict:
        """DE 上の r 列の状態と r 微分 (`_eval` の dict。march で通った範囲 [r_D, r_E])。"""
        rs = np.atleast_1d(np.asarray(rs, dtype=float))
        B = self._bnd
        R = np.asarray(B["r"])
        i = np.clip(np.searchsorted(R, rs, side="right") - 1, 0, len(R) - 2)
        Mp, thp = self._hermite(i, rs)
        return self._solve(rs, Mp, thp)

    def points(self, rs) -> np.ndarray:
        """DE 上の r 列の点 [x, r, θ, ν, M] (march の後。x は境の累積 + 部分区間の Gauss–Legendre)。"""
        rs = np.atleast_1d(np.asarray(rs, dtype=float))
        B = self._bnd
        R = np.asarray(B["r"])
        CX = np.asarray(B["cx"])
        i = np.clip(np.searchsorted(R, rs, side="right") - 1, 0, len(R) - 2)
        rg, wg = _gl(R[i], rs)                               # (n, GL_N)
        allr = np.concatenate([rs, rg.ravel()])
        alli = np.concatenate([i, np.repeat(i, GL_N)])
        Mp, thp = self._hermite(alli, allr)
        f = self._solve(allr, Mp, thp)
        nq = len(rs)
        dx = np.sum(f["dx"][nq:].reshape(nq, GL_N) * wg, axis=1)
        out = np.empty((nq, 5))
        out[:, 0] = float(self.D["x"]) + CX[i] + dx
        out[:, 1] = rs
        out[:, 2] = f["th"][:nq]
        out[:, 4] = f["M"][:nq]
        out[:, 3] = self.ig.nu(out[:, 4])
        return out

    def dtheta_dsigma_D(self) -> float:
        """DE に沿った dθ/dσ (D で、σ は C⁺ の弧長、下流向き) = (dθ/dr)·sin(θ+μ)。"""
        return float(self._fD["dth_dr"][0] * np.sin(float(self.D["th"]) + self.mu_D))


def rao_line(ig, D: dict, mdot_target: float, *, p_a: float = 0.0, dr: float = DR_DE) -> tuple:
    """DE を引いて E を求める。戻り (RaoLine, E の dict)。"""
    rl = RaoLine(ig, D, dr=dr, p_a=p_a)
    return rl, rl.march(mdot_target)


# ================================================================================================
# 入れ子の根
def _de_at(bd: BDLine, eta: float, p_a: float, dr: float) -> dict:
    """D = BD(η) から DE を引き、E と推力を返す (失敗は E["ok"] False と理由)。`F_part` = BD (η0 → D) + DE の推力
    (B を固定した族では全推力と定数差。T2 と同じ B の中の比較に使う)。"""
    D = bd.point(eta)
    m_DB, _ = bd.integrals_to_B(eta, p_a)
    m_bot, T_bot = bd.integrals_from_bottom(eta, p_a)
    rl, E = rao_line(bd.ig, D, m_DB, p_a=p_a, dr=dr)
    out = {"eta": float(eta), "D": D, "E": E, "rl": rl, "mdot_DB": m_DB, "T_BD_bottom_D": T_bot}
    if E.get("ok"):
        out["F_part"] = T_bot + E["T"]
        out["edge"] = edge_function(bd.ig, E["M"], E["th"], E["r"], p_a)
    return out


def _refine_edges(xs, vals, fun, n_bisect: int = 12):
    """走査の隣どうしで片方だけ有限 (成立の境) なら、境へ二分法で詰めて有限側の端の点を走査に足す
    (成立域の端の近くにある符号変化を拾うため)。`fun(x)` は値 (非有限 = 不成立) を返す。戻り (xs, vals) を x 昇順で。"""
    xs, vals = list(xs), list(vals)
    add = []
    for i in range(len(xs) - 1):
        a, b, fa, fb = xs[i], xs[i + 1], vals[i], vals[i + 1]
        if np.isfinite(fa) == np.isfinite(fb):
            continue
        good, bad, fg = (a, b, fa) if np.isfinite(fa) else (b, a, fb)
        for _ in range(n_bisect):
            mid = 0.5 * (good + bad)
            fm = fun(mid)
            if np.isfinite(fm):
                good, fg = mid, fm
            else:
                bad = mid
        add.append((good, fg))
    for x, v in add:
        xs.append(x)
        vals.append(v)
    o = np.argsort(xs)
    return np.asarray(xs)[o], np.asarray(vals)[o]


def solve_D(bd: BDLine, *, eps: float | None = None, p_a: float = 0.0, dr: float = DR_DE, n_scan: int = N_SCAN_D,
            eta_lo: float = 0.25, xtol: float = XTOL) -> dict:
    r"""内の根: D の位置 $\eta_D$ を、固定 $\varepsilon$ なら $r_E=\sqrt\varepsilon$、自由 ($\varepsilon$ = None) なら縁条件
    `edge_function` = 0 で決める。$\eta\in[\eta_{lo}, N)$ を `n_scan` 点で走査して符号変化をすべて候補にし、brentq で詰める。
    戻り dict: `roots` (各根の D・E・推力)、`scan` (走査の記録)、`why` (根が無いときの理由)。"""
    N = bd.N
    etas = np.linspace(bd.eta0 + eta_lo, N - 0.25, n_scan)
    target = None if eps is None else float(np.sqrt(eps))

    def fun(eta):
        d = _de_at(bd, eta, p_a, dr)
        if not d["E"].get("ok"):
            return np.nan, d
        return ((d["E"]["r"] - target) if target is not None else d["edge"]), d

    vals, recs = [], []
    for e in etas:
        v, d = fun(float(e))
        vals.append(v)
        recs.append({"eta": float(e), "val": v, "ok": bool(d["E"].get("ok")), "why": d["E"].get("why"),
                     "r_E": d["E"].get("r"), "x_E": d["E"].get("x")})
    vals = np.asarray(vals)
    etas, vals = _refine_edges(etas, vals, lambda t: fun(float(t))[0])
    roots, jumps = [], []
    for i in range(len(etas) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b <= 0.0 and not (a == 0.0 and b == 0.0):
            try:
                e = brentq(lambda t: fun(t)[0], float(etas[i]), float(etas[i + 1]), xtol=xtol,
                           rtol=4 * np.finfo(float).eps, maxiter=200)
            except ValueError:          # 区間の内側に不成立の点がある (成立域が走査点の間で途切れる)
                continue
            v, d = fun(e)
            # 根の検算: 符号変化が値の跳び (D の det の符号が変わる点で DE の向きが切り替わる — RaoLine の docstring) なら
            # brentq は跳びに収束するので、残差が小さい点だけを根にし、跳びは記録する
            if np.isfinite(v) and abs(v) <= ROOT_RESID_TOL * max(1.0, abs(target) if target is not None else 1.0):
                roots.append(d)
            else:
                jumps.append({"eta": float(e), "val": float(v) if np.isfinite(v) else None,
                              "det_D": float(d["rl"].det_D) if d["rl"].ok_D else None})
    why = None
    if not roots:
        fin = np.isfinite(vals)
        why = ("DE が全走査点で成り立たない: " + str(recs[0]["why"])) if not np.any(fin) else \
            f"r_E − √ε (または縁の関数) の符号が変わらない (有限 {int(np.sum(fin))}/{len(vals)}、範囲 {np.nanmin(vals):.4g}〜{np.nanmax(vals):.4g})"
    return {"roots": roots, "scan": recs, "why": why, "jumps": jumps}


def _middle_eval(net: KernelNet, s_B: float, eps, p_a, dr, n_scan_D):
    bd = trace_bd(net, s_B)
    sd = solve_D(bd, eps=eps, p_a=p_a, dr=dr, n_scan=n_scan_D)
    return bd, sd


def _family(det_D: float) -> int:
    """DE の族 (D での陰関数の行列式の符号、`RaoLine` の docstring)。"""
    return 1 if det_D > 0.0 else -1


def _link_branches(root_lists, tol: float) -> list:
    r"""走査点ごとの D の根 [(η_D, 族)] のリストを、隣の走査点どうし $|\Delta\eta_D|$ の小さい組から貪欲に (`tol` 以下、かつ
    同じ族 [$\det_D$ の符号] だけ) 結んで枝に分ける。戻り: 各走査点の各根の枝番号。添字 (根の並び順) では結ばない — 根が
    1 つ増えると添字がずれて別の枝をつなぐため (plan §4.4 の訂正)。族の違う根を結ばないのは、$\det_D=0$ を跨ぐと E が跳ぶ
    (DE の向きが切り替わる) ので、連続な根の枝は $\det_D$ の符号を変えられないため。"""
    labels, prev, nid = [], [], 0
    for cur in root_lists:
        lab = [None] * len(cur)
        pairs = sorted((abs(e - pe), ip, ic) for ic, (e, fam) in enumerate(cur)
                       for ip, (pe, pfam, _) in enumerate(prev) if fam == pfam)
        used_p, used_c = set(), set()
        for dd, ip, ic in pairs:
            if dd > tol or ip in used_p or ic in used_c:
                continue
            lab[ic] = prev[ip][2]
            used_p.add(ip)
            used_c.add(ic)
        for ic in range(len(cur)):
            if lab[ic] is None:
                lab[ic] = nid
                nid += 1
        labels.append(lab)
        prev = [(cur[ic][0], cur[ic][1], lab[ic]) for ic in range(len(cur))]
    return labels


def _nearest_root(sd: dict, eta_ref: float, tol: float, fam: int):
    """`solve_D` の根のうち族 `fam` で η_D が eta_ref に最も近いもの (|Δη_D| ≤ tol、無ければ None)。"""
    best = None
    for r in sd["roots"]:
        d = abs(r["eta"] - eta_ref)
        if _family(r["rl"].det_D) == fam and d <= tol and (best is None or d < best[0]):
            best = (d, r)
    return None if best is None else best[1]


class _BranchLost(Exception):
    """brentq の途中の s_B で枝の D の根が `eta_tol` 以内に見つからない。"""


def _accept_candidate(cand: dict, theta_B: float, theta_E: float, det_D: float, de: float) -> bool:
    r"""中の根の候補の受理 (plan §4.4「主枝の受理規則」、diagnostician 2026-10-11): **受理条件は $\det_D>0$ だけ**
    ($\det_D>0\iff$ DE が D で流れを軸向きに曲げ始める [$a>0$ の域で $d\theta/dr=-a/(r\det)$] — 古典の枝の延長)。
    $\theta_E\le\theta_B$ と $de=(d\theta/d\sigma)_{DE,D}<0$ は記録のフラグ。$\det_D>0$ なのに $\theta_E>\theta_B$ なら
    `policy_hole` に理由を書く (受理はする — 方針の穴として呼び出し側が止めて諮る)。`cand` に記録を足し、受理なら True。"""
    cand["det_D"] = float(det_D)
    cand["de"] = float(de)
    cand["theta_E_le_theta_B"] = bool(theta_E <= theta_B)
    cand["de_lt_0"] = bool(de < 0.0)
    if det_D > 0.0:
        if not theta_E <= theta_B:
            cand["policy_hole"] = (f"det_D {det_D:+.4g} > 0 なのに θ_E {np.rad2deg(theta_E):.3f}° > θ_B {np.rad2deg(theta_B):.3f}° "
                                   f"(plan §4.4 の受理規則の穴 — 諮る)")
        return True
    cand["why"] = (f"det_D {det_D:+.4g} ≤ 0 (受理条件 det_D > 0 — D で DE が流れを軸から離す向きに曲げる族、"
                   f"θ_E {np.rad2deg(theta_E):.3f}°・θ_B {np.rad2deg(theta_B):.3f}°・de {de:+.4g})")
    return False


def solve_thetaB(net: KernelNet, L: float, *, eps: float | None, p_a: float = 0.0, dr: float = DR_DE,
                 s_range: tuple | None = None, n_scan: int = N_SCAN_SB, n_scan_D: int = N_SCAN_D,
                 xtol: float = XTOL, eta_tol: float | None = None, n_bisect: int = 10) -> dict:
    r"""中の根: B の弧長 $s_B$ (= $\theta_B$) を $x_E=L$ で決める (kernel 場 1 本の中で B を動かす)。

    **枝は $\eta_D$ の連続性で追う** (plan §4.4 の訂正、diagnostician 2026-10-11): 走査点ごとの D の根を、隣の走査点どうし
    $|\Delta\eta_D|\le$ `eta_tol` (既定 `ETA_LINK_FRAC`$\cdot n$) かつ同じ族 ($\det_D$ の符号) で結んで枝を作る (`_link_branches`)。枝の端は二分法で詰め
    (直前の $\eta_D$ に最も近い根が `eta_tol` 以内にある限り同じ枝)、枝ごとに $x_E-L$ の符号変化を brentq で詰める。途中の
    $s_B$ では括りの両端の $\eta_D$ の線形補間に最も近い根を取る (`eta_tol` 以内に無ければその括りは `lost` に記録して捨てる)。
    根の検算は $|x_E-L|\le$ `ROOT_RESID_TOL`$\cdot\max(1,L)$ (超えたら値の跳びとして `lost`)。

    **受理条件は $\det_D>0$** (`_accept_candidate`、plan §4.4「主枝の受理規則」)。$\det_D\le0$ の候補は `rejected` に理由つきで
    残す。$\theta_E\le\theta_B$ と $de<0$ は候補ごとのフラグ (記録)。$\det_D>0$ なのに $\theta_E>\theta_B$ の候補は受理したうえで
    `policy_holes` にも出す (方針の穴として呼び出し側が止めて諮る)。
    戻り dict: `roots` (受理)、`rejected`、`lost`、`policy_holes`、`scan` (走査点ごとの全根)、`branches` (枝ごとの点列)、
    `eta_tol`、`why`。"""
    s_lo, s_hi = s_range if s_range is not None else (net.s_wall[2], net.s_wall[net.k_need])
    tol = float(eta_tol) if eta_tol is not None else ETA_LINK_FRAC * net.n
    ss = np.linspace(s_lo, s_hi, n_scan)
    cache = {}

    def mid(s):
        key = float(s)
        if key not in cache:
            cache[key] = _middle_eval(net, key, eps, p_a, dr, n_scan_D)
        return cache[key]

    def rec_root(r):
        E = r["E"]
        return {"eta_D": float(r["eta"]), "x_E": float(E["x"]), "r_E": float(E["r"]),
                "theta_E_deg": float(np.rad2deg(E["th"])), "det_D": float(r["rl"].det_D)}

    scan = []
    for s in ss:
        bd, sd = mid(s)
        scan.append({"s_B": float(s), "theta_B_deg": float(np.rad2deg(net.wall.theta(s))), "eta0": int(bd.eta0),
                     "N": int(bd.N), "roots": [rec_root(r) for r in sd["roots"]], "jumps": sd.get("jumps", []),
                     "why": sd["why"]})
    labels = _link_branches([[(r["eta_D"], _family(r["det_D"])) for r in sc["roots"]] for sc in scan], tol)
    for sc, lab in zip(scan, labels):
        for r, b in zip(sc["roots"], lab):
            r["branch"] = int(b)
    # 枝ごとの点列 [s_B, η_D, x_E − L] と族
    br_pts, br_fam = {}, {}
    for sc, lab in zip(scan, labels):
        for r, b in zip(sc["roots"], lab):
            br_pts.setdefault(int(b), []).append([sc["s_B"], r["eta_D"], r["x_E"] - L])
            br_fam[int(b)] = _family(r["det_D"])
    # 枝の端を二分法で詰める (隣の走査点にこの枝の根が無い側)
    s_list = [float(s) for s in ss]
    for b, pts in br_pts.items():
        pts.sort()
        add = []
        for end, step in ((pts[0], -1), (pts[-1], +1)):
            i_end = s_list.index(end[0])
            i_nb = i_end + step
            if not 0 <= i_nb < len(s_list):
                continue
            good_s, good_eta, good_v = end
            bad_s = s_list[i_nb]
            for _ in range(n_bisect):
                m = 0.5 * (good_s + bad_s)
                r = _nearest_root(mid(m)[1], good_eta, tol, br_fam[b])
                if r is None:
                    bad_s = m
                else:
                    good_s, good_eta, good_v = m, float(r["eta"]), float(r["E"]["x"]) - L
            if good_s != end[0]:
                add.append([good_s, good_eta, good_v])
        pts.extend(add)
        pts.sort()
    roots, rejected, lost = [], [], []
    for b, pts in sorted(br_pts.items()):
        for pa, pb in zip(pts[:-1], pts[1:]):
            (sa, ea, va), (sb, eb, vb) = pa, pb
            if not (np.isfinite(va) and np.isfinite(vb) and va * vb <= 0.0) or (va == 0.0 and vb == 0.0):
                continue

            def eref(s, sa=sa, sb=sb, ea=ea, eb=eb):
                return ea + (eb - ea) * (s - sa) / (sb - sa)

            def f(s, fam=br_fam[b]):
                r = _nearest_root(mid(s)[1], eref(s), tol, fam)
                if r is None:
                    raise _BranchLost(s)
                return float(r["E"]["x"]) - L
            try:
                s0 = brentq(f, sa, sb, xtol=xtol, rtol=4 * np.finfo(float).eps, maxiter=200)
            except _BranchLost as e:
                lost.append({"branch": b, "theta_B_deg": [float(np.rad2deg(net.wall.theta(sa))), float(np.rad2deg(net.wall.theta(sb)))],
                             "why": f"brentq の途中 (s_B {float(e.args[0]):.6f}) で枝の D の根が η_D の {tol:g} 以内に無い"})
                continue
            bd, sd = mid(s0)
            r = _nearest_root(sd, eref(s0), tol, br_fam[b])
            v = float(r["E"]["x"]) - L
            thB, thE = float(net.wall.theta(s0)), float(r["E"]["th"])
            if not abs(v) <= ROOT_RESID_TOL * max(1.0, abs(L)):
                lost.append({"branch": b, "theta_B_deg": float(np.rad2deg(thB)), "eta_D": float(r["eta"]), "x_E_minus_L": v,
                             "why": "x_E − L の値の跳び (根の検算で残差が大きい)"})
                continue
            if any(q["branch"] == b and abs(q["s_B"] - s0) <= 1e-10 for q in roots + rejected):
                continue
            cand = {"s_B": float(s0), "theta_B": thB, "kappa_B": float(net.wall.kappa(s0)), "bd": bd, "D_root": r,
                    "branch": b, "x_E_minus_L": v}
            if _accept_candidate(cand, thB, thE, float(r["rl"].det_D), r["rl"].dtheta_dsigma_D()):
                roots.append(cand)
            else:
                rejected.append(cand)
    holes = [c for c in roots if c.get("policy_hole")]
    why = None
    if not roots:
        if rejected:
            why = (f"受理条件 det_D > 0 を満たす根なし (捨てた候補 {len(rejected)} 個: "
                   + "、".join(f"θ_B {np.rad2deg(c['theta_B']):.3f}°・det_D {c['det_D']:+.4g}・θ_E {np.rad2deg(c['D_root']['E']['th']):.3f}°"
                              for c in rejected) + ")")
        else:
            why = (f"x_E − L の符号が変わる枝なし (θ_B {scan[0]['theta_B_deg']:.2f}〜{scan[-1]['theta_B_deg']:.2f}°、枝 "
                   f"{len(br_pts)} 本: " + "; ".join(
                       f"枝 {b}: θ_B {np.rad2deg(net.wall.theta(p[0][0])):.2f}〜{np.rad2deg(net.wall.theta(p[-1][0])):.2f}°・"
                       f"x_E − L {min(q[2] for q in p):+.3g}〜{max(q[2] for q in p):+.3g}" for b, p in sorted(br_pts.items()))
                   + (f"; lost {len(lost)}" if lost else "") + ")")
    branches = {b: {"family_det_sign": br_fam[b],
                    "points": [{"s_B": p[0], "theta_B_deg": float(np.rad2deg(net.wall.theta(p[0]))), "eta_D": p[1],
                                "x_E_minus_L": p[2]} for p in pts]} for b, pts in br_pts.items()}
    return {"roots": roots, "rejected": rejected, "lost": lost, "policy_holes": holes, "scan": scan, "branches": branches,
            "eta_tol": tol, "why": why}


def candidate_summary(c: dict) -> dict:
    """中の根の候補 (`solve_thetaB` の `roots`・`rejected` の要素) から表に載せる量だけを抜く (BD・DE のオブジェクトを外す)。"""
    E = c["D_root"]["E"]
    out = {"theta_B_deg": float(np.rad2deg(c["theta_B"])), "s_B": float(c["s_B"]), "kappa_B": float(c["kappa_B"]),
           "eta_D": float(c["D_root"]["eta"]), "x_E": float(E["x"]), "r_E": float(E["r"]),
           "theta_E_deg": float(np.rad2deg(E["th"])), "M_E": float(E["M"]), "branch": c.get("branch"),
           "det_D": float(c["D_root"]["rl"].det_D), "de": float(c["D_root"]["rl"].dtheta_dsigma_D()),
           "theta_E_le_theta_B": bool(E["th"] <= c["theta_B"]), "de_lt_0": bool(c["D_root"]["rl"].dtheta_dsigma_D() < 0.0)}
    for k in ("why", "policy_hole"):
        if c.get(k):
            out[k] = c[k]
    return out


def gtilde(net: KernelNet, mid_root: dict) -> dict:
    r"""C2 の判定量 $\tilde g=(d\theta/d\sigma)_{DE,D}-(d\theta/d\sigma)_{kernel\,C^+,D}$ [rad/$r_t$] (0 なら B で曲率連続と等価。
    2026-10-11 に C2 を撤回したので根は取らず記録する — plan §4.4)。kernel 側は網の論理座標で
    C⁺ が $(d\xi,d\eta)=(1,1)$ の向きなので $(\theta_\xi+\theta_\eta)/|(x_\xi+x_\eta,\,r_\xi+r_\eta)|$ (BD の両側で滑らかな延長場)。"""
    bd = mid_root["bd"]
    rt = mid_root["D_root"]
    eta = rt["eta"]
    d_eta = bd.dspl(eta)
    d_xi = bd.dxi(eta)
    tx = d_xi[0] + d_eta[0]
    tr = d_xi[1] + d_eta[1]
    tth = d_xi[2] + d_eta[2]
    k_side = float(tth / np.hypot(tx, tr))
    de_side = rt["rl"].dtheta_dsigma_D()
    # C⁺ の向きの検算: kernel 側の C⁺ の傾きと tan(θ+μ)
    D = rt["D"]
    mu = float(mu_vec(np.array([D["M"]]))[0])
    slope_err = float(np.arctan2(tr, tx) - (D["th"] + mu))
    return {"g": de_side - k_side, "de": de_side, "kernel": k_side, "cplus_angle_err": slope_err}


# ================================================================================================
# Goursat と壁
def goursat_fill(ig, bd: BDLine, D_root: dict, m: int, q: int, ext: float = 0.08) -> np.ndarray:
    r"""領域 BDE の Goursat 充填。`P[i, j]` ($i=0..m$ が BD 上 D → B、$j=0..q$ が DE 上 D → E [$r$ に等間隔 — DE の $M$ は
    θ = μ で極値をとり単調でない]、その先 E を越える延長) = [x, r, θ, ν, M]。$P_{i,j}$ = interior(A = $P_{i,j-1}$ [C⁺ 担体],
    B = $P_{i-1,j}$ [C⁻ 担体、符号付きで逆向き])、段 $i+j$ ごとに一括。DE は流量で `ext` だけ E の先へ延ばす (流量で切る壁が
    E の近くで網の外に出ないように — 中点則の流量は厳密な目標と $O(h^2)$ ずれる)。延長が作れなければ延長なし。"""
    ig = as_isentrope(ig)
    eta_D = D_root["eta"]
    etas = np.linspace(eta_D, bd.N, m + 1)
    v = bd.at(etas)
    rE = float(D_root["E"]["r"])
    rD = float(D_root["D"]["r"])
    rs = np.linspace(rD, rE, q + 1)
    rl = D_root["rl"]
    q_ext = 0
    if ext > 0.0:
        rl2 = RaoLine(ig, D_root["D"], dr=rl.dr, p_a=rl.p_a)
        E2 = rl2.march(D_root["E"]["mdot"] * (1.0 + ext))
        if E2.get("ok"):
            drs = (rE - rD) / q
            q_ext = max(int(np.floor((E2["r"] - rE) / drs)), 0)
            if q_ext > 0:
                rs = np.concatenate([rs, rE + drs * np.arange(1, q_ext + 1)])
                rl = rl2
    qq = q + q_ext
    P = np.full((m + 1, qq + 1, 5), np.nan)
    P[:, 0, :4] = v
    P[:, 0, 4] = ig.mach_of_nu(v[:, 3])
    P[0, :, :] = rl.points(rs)
    # D は BD と DE の共通点: BD 側の値を使う (DE は同じ状態から作った式なので一致)
    P[0, 0, :4] = v[0]
    P[0, 0, 4] = ig.mach_of_nu(np.array([v[0, 3]]))[0]
    nfail = 0
    for lev in range(2, m + qq + 1):
        i = np.arange(max(1, lev - qq), min(m, lev - 1) + 1)
        j = lev - i
        A = [P[i, j - 1, k] for k in range(5)]
        B = [P[i - 1, j, k] for k in range(5)]
        st: dict = {}
        xP, rP, thP, nuP, ok = interior_vec(*A, *B, ig.moc, 1.0, 2, axis_limit="analytic", thrA=None, thrB=None,
                                            corrector="converge", stats=st)
        good = ok & (st["status"] == PAIR_CONVERGED)
        nfail += int(np.sum(~good & np.isfinite(A[0]) & np.isfinite(B[0])))
        P[i, j, 0] = np.where(good, xP, np.nan)
        P[i, j, 1] = np.where(good, rP, np.nan)
        P[i, j, 2] = np.where(good, thP, np.nan)
        P[i, j, 3] = np.where(good, nuP, np.nan)
        P[i, j, 4] = np.where(good, ig.mach_of_nu(np.where(good, nuP, 0.0)), np.nan)
    goursat_fill.last_nfail = nfail
    goursat_fill.last_q_ext = q_ext
    return P


def wall_from_massflow(ig, bd: BDLine, D_root: dict, P: np.ndarray, p_a: float = 0.0) -> dict:
    r"""壁 BE = 各 C⁺ 線 $i$ (BD 上の $b_i$ から) で、$b_i$ から数えた流量が $\dot m_{BD}(b_i\to B)$ に達する点
    (`moc_inverse.cplus_flux_wall` と同じ: 中点則の流量、跨いだ線分を線形補間)。戻り: `wall` (B → E の順、[x, r, θ, M])、
    `closure` (線 $i=0$ = DE 上の壁点と E の差)、`n_missing` (壁に届かない線)。"""
    ig = as_isentrope(ig)
    m = P.shape[0] - 1
    eta_D = D_root["eta"]
    etas = np.linspace(eta_D, bd.N, m + 1)
    targets = np.array([bd.integrals_to_B(e, p_a)[0] for e in etas])
    targets[-1] = 0.0
    out = []
    miss = 0
    for i in range(m, -1, -1):
        line = P[i]
        ok = np.isfinite(line[:, 0])
        cut = int(np.argmin(ok)) if not np.all(ok) else len(ok)
        ln = line[:cut]
        if targets[i] <= 0.0:
            out.append((ln[0, 0], ln[0, 1], ln[0, 2], ln[0, 4]))
            continue
        if len(ln) < 2:
            miss += 1
            continue
        Mm = 0.5 * (ln[1:, 4] + ln[:-1, 4])
        thm = 0.5 * (ln[1:, 2] + ln[:-1, 2])
        rm = 0.5 * (ln[1:, 1] + ln[:-1, 1])
        d = ig.mass_flux_density(Mm) * (np.cos(thm) * np.diff(ln[:, 1]) - np.sin(thm) * np.diff(ln[:, 0])) * TWO_PI * rm
        cum = np.concatenate([[0.0], np.cumsum(d)])
        if cum[-1] < targets[i]:
            miss += 1
            continue
        jj = int(np.argmax(cum >= targets[i]))
        s = (targets[i] - cum[jj - 1]) / max(cum[jj] - cum[jj - 1], 1e-300)
        a, b = ln[jj - 1], ln[jj]
        out.append(tuple(a[[0, 1, 2, 4]] + s * (b[[0, 1, 2, 4]] - a[[0, 1, 2, 4]])))
    w = np.asarray(out, dtype=float)
    w = w[np.argsort(w[:, 0])]
    E = D_root["E"]
    closure = {"dx": float(w[-1, 0] - E["x"]), "dr": float(w[-1, 1] - E["r"]), "dth": float(w[-1, 2] - E["th"])}
    # 独立の検算: 壁は流線でもあるはず (折れ線の傾きと tanθ)
    slope = np.diff(w[:, 1]) / np.diff(w[:, 0])
    tmid = np.tan(0.5 * (w[1:, 2] + w[:-1, 2]))
    return {"wall": w, "closure": closure, "n_missing": miss,
            "streamline_resid_max": float(np.max(np.abs(np.arctan(slope) - np.arctan(tmid))))}


def wall_curvature_B(wall: np.ndarray, n_fit: int = 6) -> float:
    """BE の出だしの曲率 dθ/ds (B から n_fit 点の 2 次の当てはめ)。"""
    w = wall[:n_fit]
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(w[:, 0]), np.diff(w[:, 1])))])
    c = np.polyfit(s, w[:, 2], 2)
    return float(c[1])


# ================================================================================================
# 推力の会計
def line_integrals(ig, x, r, th, M, p_a: float = 0.0, spline: bool = True) -> tuple:
    """折れ線 (下 → 上) を横切る流量と推力。`spline` は点列の添字の 3 次スプライン + GL (初期線用)、False は中点則。"""
    ig = as_isentrope(ig)
    x, r, th, M = (np.asarray(v, dtype=float) for v in (x, r, th, M))
    if spline:
        t = np.arange(len(x), dtype=float)
        sp = CubicSpline(t, np.stack([x, r, th, ig.nu(M)], axis=-1), axis=0)
        dsp = sp.derivative()
        tg, wg = _gl(t[:-1], t[1:])
        v = sp(tg)
        d = dsp(tg)
        Mg = ig.mach_of_nu(v[..., 3])
        st = ig.state(np.maximum(Mg, 1.0).ravel())
        rhoV = (st["rho"] * st["V"]).reshape(Mg.shape)
        V = st["V"].reshape(Mg.shape)
        p = st["p"].reshape(Mg.shape)
        nflux = np.cos(v[..., 2]) * d[..., 1] - np.sin(v[..., 2]) * d[..., 0]
        dm = TWO_PI * v[..., 1] * rhoV * nflux
        dT = TWO_PI * v[..., 1] * (rhoV * V * np.cos(v[..., 2]) * nflux + (p - p_a) * d[..., 1])
        return float(np.sum(dm * wg)), float(np.sum(dT * wg))
    Mm = 0.5 * (M[1:] + M[:-1])
    st = ig.state(np.maximum(Mm, 1.0))
    thm = 0.5 * (th[1:] + th[:-1])
    rm = 0.5 * (r[1:] + r[:-1])
    nflux = np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x)
    rhoV = st["rho"] * st["V"]
    return (float(np.sum(TWO_PI * rm * rhoV * nflux)),
            float(np.sum(TWO_PI * rm * (rhoV * st["V"] * np.cos(thm) * nflux + (st["p"] - p_a) * np.diff(r)))))


def wall_pressure_force(ig, x, r, M, p_a: float = 0.0) -> float:
    """壁の点列 (x 昇順) の圧力の軸力 ∫(p − p_a) 2πr dr (台形)。"""
    ig = as_isentrope(ig)
    p = ig.state(np.maximum(np.asarray(M, dtype=float), 1.0))["p"]
    r = np.asarray(r, dtype=float)
    return float(np.sum(0.5 * ((p[1:] - p_a) * r[1:] + (p[:-1] - p_a) * r[:-1]) * TWO_PI * np.diff(r)))


def thrust_accounts(net: KernelNet, mid_root: dict, wall: dict, p_a: float = 0.0) -> dict:
    r"""推力の 2 会計: (a) 制御面 = (軸 → D) + DE、(b) 初期線 + kernel 壁 (0 → $s_B$) + BE 壁。(a) の軸 → D は BD の代わりに
    D を通る kernel 場の C⁺ を遡って閉じる (`cplus_closure`。BD は θ_B > μ_B では軸に着くのが出口よりずっと下流で、網が粗い)。
    差は MOC の離散化の運動量の非保存 (n とともに減る)。流量も 2 通り (初期線と、閉じた曲線 + DE) を返す。"""
    ig = net.ig
    line = net.line
    rt = mid_root["D_root"]
    bd = mid_root["bd"]
    cl = cplus_closure(net, bd.xi_B, rt["eta"], p_a)
    F_cs = cl["T"] + rt["E"]["T"]
    m_cs = cl["mdot"] + rt["E"]["mdot"]
    m_line, T_line = line_integrals(ig, line["x"], line["r"], line["theta"], line["M"], p_a, spline=True)
    s_B = mid_root["s_B"]
    sgrid = np.linspace(0.0, s_B, 401)
    sg, wg = _gl(sgrid[:-1], sgrid[1:])
    nu = net.nu_wall_s(sg)
    st = ig.state(np.maximum(ig.mach_of_nu(nu), 1.0).ravel())
    p = st["p"].reshape(sg.shape)
    rr = net.wall.xr(sg.ravel())[1].reshape(sg.shape)
    T_kernel = float(np.sum((p - p_a) * TWO_PI * rr * np.sin(net.wall.theta(sg)) * wg))
    w = wall["wall"]
    T_BE = wall_pressure_force(ig, w[:, 0], w[:, 1], w[:, 3], p_a)
    F_wall = T_line + T_kernel + T_BE
    return {"F_cs": F_cs, "F_wall": F_wall, "rel_diff": (F_cs - F_wall) / F_cs, "T_line": T_line, "T_kernel": T_kernel,
            "T_BE": T_BE, "T_closure": cl["T"], "closure_hits": cl["hits"], "T_DE": rt["E"]["T"], "mdot_line": m_line,
            "mdot_cs": m_cs, "mdot_rel_diff": (m_cs - m_line) / m_line}


# ================================================================================================
# 全体
def _hall_line(R_t, ig, n):
    from .moc_forward import hall_initial_line
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return hall_initial_line(R_t, ig.moc, n)


def solve_kappaB(ig, R_t: float, L: float, *, c: float, eps: float | None, line: dict, p_a: float = 0.0,
                 dr: float = DR_DE, n_scan_sb: int = N_SCAN_SB, n_scan_D: int = N_SCAN_D, L_march: float | None = None,
                 log=None) -> dict:
    r"""与えた $c$ ($\kappa(s)=1/R_t+cs$、$c=0$ が円弧) の kernel 場を作り、中の根を解いて各根の $\tilde g$ を記録する。
    **外の根 (C2) は取らない** (2026-10-11 diagnostician、plan §4.4「(M) は円弧 kernel・B で C1」: $\kappa_B$ を振る族には
    $\tilde g=0$ の根が無い — $\kappa_B$ の影響は BD を伝わる間にほぼ失われ、D での DE の横断微分は D の状態だけで決まる)。
    $c$ は明示パラメータで、走査も brentq もしない (関数名は旧版の外の根の名残)。
    戻り dict: `c`、`ok` (受理した中の根があるか)、`why`、`net`、`mid` (`solve_thetaB` の戻り。kernel 場が作れなければ無し)、
    `roots` (各根の `c`・`net`・`mid_root`・`g`・`gt` [`gtilde` の dict])。"""
    ig = as_isentrope(ig)
    c = float(c)
    L_march = L_march if L_march is not None else default_L_march(L)
    rng = sB_range(R_t, c)
    if rng is None:
        return {"c": c, "ok": False, "why": "s_B の探索域が空", "roots": []}
    try:
        net = kernel_field(ig, R_t, c, line, rng[1], L_march=L_march)
    except RuntimeError as e:
        return {"c": c, "ok": False, "why": f"kernel 場: {e}", "roots": []}
    s_lo = max(rng[0], float(net.s_wall[2]))
    mid = solve_thetaB(net, L, eps=eps, p_a=p_a, dr=dr, s_range=(s_lo, min(rng[1], net.s_need)), n_scan=n_scan_sb,
                       n_scan_D=n_scan_D)
    roots = []
    for r in mid["roots"]:
        gt = gtilde(net, r)
        roots.append({"c": c, "net": net, "mid_root": r, "g": gt["g"], "gt": gt})
    if log:
        log(f"  c·R_t² {c * R_t * R_t:+.6f}: θ_B {[round(np.rad2deg(r['theta_B']), 4) for r in mid['roots']]}° "
            f"κ_B {[round(r['kappa_B'], 5) for r in mid['roots']]} g̃ {[f'{r['g']:+.4e}' for r in roots]}")
    return {"c": c, "ok": bool(roots), "why": None if roots else "中の根なし: " + str(mid["why"]), "net": net, "mid": mid,
            "roots": roots}


def design_rao(gas, *, R_t: float, L: float, eps: float | None, kernel: str = "arc", c: float | None = None,
               n: int = 321, line: dict | None = None, p_a: float = 0.0, dr: float | None = None, m: int | None = None,
               q: int | None = None, L_march: float | None = None, n_scan_sb: int = N_SCAN_SB, n_scan_D: int = N_SCAN_D,
               log=None) -> dict:
    r"""与えた $(L,\varepsilon)$ の停留輪郭 (固定 $\varepsilon$。`eps=None` は自由 $\varepsilon$ で縁条件)。

    `gas`: γ・`GasCPG`・`GasSemiPerfect`・`Isentrope`。`kernel`: **既定 `"arc"`** (古典、曲率半径 $R_t$ の円弧 = $c=0$。
    (M) の標準 — plan §4.4「(M) は円弧 kernel・B で C1」) / `"linear"` (曲率 $1/R_t+cs$ の kernel、`c` を明示して与える。
    C2 の根は取らない [撤回] — 記録用)。`line`: 初期線 (既定 Hall の $n$ 点。$R_t\le1$ は $n\le321$)。
    `m`・`q`: Goursat の BD・DE の点数 (**既定 n** — 2026-10-12 に n/2 から変更、plan §5.1 #5 (g1): T6 で m・q 倍は
    費用 2 s で推力 2 会計の差を約 1/3 にし、θ_max・F_wall・Δκ_B は m に律速されていた。根 [θ_B・θ_E・η_D] は m・q に
    依らない)、`dr`: DE の r 格子 (既定 3.2/n — n 161 で 0.02)。
    中の根の受理条件は $\det_D>0$ (`solve_thetaB`)。B で C1 (接線連続) で、曲率の跳び `dkappa_B` $=\kappa_W(B^+)-\kappa_B$
    (円弧では $\kappa_W(B^+)-1/R_t$) と C2 の判定量 `g_tilde` を記録する。
    戻り dict: `status` (`ok` / `infeasible`)、`why`、`wall` (kernel + BE の [x, r, θ, M])、`theta_B`・`theta_max`・`theta_E`
    [rad]・`kappa_B`・`c`・`s_B`・`eta_D`・`D`・`E`・`thrust`・`closure`・`g_tilde`・`kappaW_B`・`dkappa_B`・`det_D`・`de`・
    `theta_E_le_theta_B`・`de_lt_0`・`candidates`・`timing`、中の根の記録 (`rejected` = 受理条件 det_D > 0 で捨てた候補、
    `policy_holes` = det_D > 0 なのに θ_E > θ_B の受理候補 [方針の穴、空でなければ諮る]、`lost`、`branches`、`scan_mid`)。"""
    if kernel == "arc":
        if c is not None and float(c) != 0.0:
            raise ValueError("design_rao: kernel='arc' は c = 0 (c を与えるなら kernel='linear')")
        c = 0.0
    elif kernel == "linear":
        if c is None:
            raise ValueError("design_rao: kernel='linear' は c (曲率の勾配) を明示する — 外の根 (C2) は撤回 (plan §4.4)")
        c = float(c)
    else:
        raise ValueError("design_rao: kernel は 'arc' か 'linear'")
    t0 = time.time()
    ig = as_isentrope(gas)
    L_march = L_march if L_march is not None else default_L_march(L)
    if line is None:
        line = _hall_line(R_t, ig, n)
    n = len(line["x"])
    dr = dr if dr is not None else 3.2 / n
    m = m if m is not None else max(n, 20)
    q = q if q is not None else max(n, 20)
    base = {"R_t": R_t, "L": L, "eps": eps, "kernel": kernel, "c": c, "n": n, "p_a": p_a, "dr": dr, "m": m, "q": q,
            "gas": ig.kind, "gamma_star": ig.gamma_star, "accept_rule": "det_D > 0"}
    ko = solve_kappaB(ig, R_t, L, c=c, eps=eps, line=line, p_a=p_a, dr=dr, L_march=L_march, n_scan_sb=n_scan_sb,
                      n_scan_D=n_scan_D, log=log)
    mid = ko.get("mid")
    mid_rec = {} if mid is None else {
        "rejected": [candidate_summary(cc) for cc in mid["rejected"]], "lost": mid["lost"],
        "policy_holes": [candidate_summary(cc) for cc in mid["policy_holes"]],
        "branches": mid["branches"], "scan_mid": mid["scan"], "eta_tol": mid["eta_tol"]}
    if not ko["roots"]:
        return {**base, "status": "infeasible", "why": ko["why"], **mid_rec, "timing": {"total_s": time.time() - t0}}
    cands = ko["roots"]
    # 候補ごとに Goursat → 壁 → 推力。複数なら制御面の推力が最大
    built = []
    for cd in cands:
        net, mr = cd["net"], cd["mid_root"]
        P = goursat_fill(ig, mr["bd"], mr["D_root"], m, q)
        nfail = goursat_fill.last_nfail
        wl = wall_from_massflow(ig, mr["bd"], mr["D_root"], P, p_a)
        th_acc = thrust_accounts(net, mr, wl, p_a)
        built.append({"cd": cd, "P": P, "nfail": nfail, "wl": wl, "thrust": th_acc})
    best = max(built, key=lambda b: b["thrust"]["F_cs"])
    cd, net, mr = best["cd"], best["cd"]["net"], best["cd"]["mid_root"]
    s_B = mr["s_B"]
    sk = np.linspace(0.0, s_B, 801)
    xk, rk = net.wall.xr(sk)
    thk = net.wall.theta(sk)
    Mk = net.ig.mach_of_nu(net.nu_wall_s(sk))
    kern = np.column_stack([xk, rk, thk, Mk])
    BE = best["wl"]["wall"]
    full = np.vstack([kern, BE[1:]]) if abs(BE[0, 0] - xk[-1]) < 1e-6 else np.vstack([kern, BE])
    E = mr["D_root"]["E"]
    thmax_k = float(np.max(thk))
    kW = wall_curvature_B(BE)
    out = {**base, "status": "ok", "why": None, "wall": full, "kernel_wall": kern, "BE": BE,
           "theta_B": mr["theta_B"], "theta_max": max(thmax_k, float(np.max(BE[:, 2]))), "theta_E": E["th"],
           "kappa_B": mr["kappa_B"], "c": cd["c"], "s_B": s_B, "x_B": float(xk[-1]), "r_B": float(rk[-1]),
           "D": mr["D_root"]["D"], "E": {k: v for k, v in E.items() if k in ("M", "x", "r", "th", "T", "mdot")},
           "edge": mr["D_root"].get("edge"), "Phi2": mr["D_root"]["rl"].Phi2, "Phi3": mr["D_root"]["rl"].Phi3,
           "thrust": best["thrust"], "closure": best["wl"]["closure"], "streamline_resid_max": best["wl"]["streamline_resid_max"],
           "n_wall_missing": best["wl"]["n_missing"], "goursat_nfail": best["nfail"],
           "g_tilde": cd["gt"], "kappaW_B": kW, "dkappa_B": kW - mr["kappa_B"],
           "eta_D": float(mr["D_root"]["eta"]), "eta0_BD": int(mr["bd"].eta0), "N_BD": int(mr["bd"].N),
           "det_D": mr["det_D"], "de": mr["de"], "theta_E_le_theta_B": mr["theta_E_le_theta_B"], "de_lt_0": mr["de_lt_0"],
           "branch": mr.get("branch"),
           "candidates": [{**candidate_summary(b["cd"]["mid_root"]), "c": b["cd"]["c"], "F_cs": b["thrust"]["F_cs"],
                           "rel_diff": b["thrust"]["rel_diff"], "g_tilde": b["cd"]["g"]} for b in built],
           **mid_rec,
           "net_K": net.K, "line_n": n, "_net": net, "_mid": mr, "_goursat": best["P"],
           "timing": {"total_s": time.time() - t0}}
    F1 = ig.ideal_1d(eps if eps is not None else E["r"] ** 2, p_a)
    out["F_1D"] = F1["F"]
    out["CF"] = best["thrust"]["F_cs"] / (np.pi * ig.p0)
    out["isp_ratio"] = (best["thrust"]["F_cs"] / best["thrust"]["mdot_line"]) / (F1["F"] / F1["mdot"])
    return out


# ================================================================================================
# 与えた壁の推力 (前進 MOC で場を解き直す) — T4 停留性・T1a′ (plan §5.1 #5 (g)・§6 V1)
class ExtendedWall:
    r"""壁 `wall` (`r(x, deriv)`、$x\le L$) を $x>L$ で接線に延ばしたもの。前進 MOC を $L$ より先まで回して出口 $x=L$ を
    跨ぐ壁点を得るために使う ($x\le L$ の場は延長に依らない — 依存は下流へしか伝わらない)。"""

    def __init__(self, wall, L: float):
        self.w, self.L = wall, float(L)
        self.rL = float(np.asarray(wall.r(np.array([self.L]), 0), dtype=float).ravel()[0])
        self.rpL = float(np.asarray(wall.r(np.array([self.L]), 1), dtype=float).ravel()[0])

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        xa = np.atleast_1d(x)
        out = np.empty_like(xa)
        ins = xa <= self.L
        if np.any(ins):
            out[ins] = np.asarray(self.w.r(xa[ins], deriv), dtype=float).ravel()
        if np.any(~ins):
            out[~ins] = (self.rL + self.rpL * (xa[~ins] - self.L)) if deriv == 0 else self.rpL
        return out.reshape(np.shape(x))


def _pspline_fit(x, r, start: tuple, end: tuple, n_int: int, lam: float | None, deg: int = 5):
    r"""罰則付き最小二乗の 5 次 B-spline (`RaoWall` の `"pspline"`): 節点は折れ線 $(x_i,r_i)$ の弧長に等間隔の `n_int` 区間、
    罰則 $\lambda\int(r''')^2dx$、端の拘束 $r(x_0),r'(x_0)$ = `start`・$r(x_N),r'(x_N)$ = `end` (拘束の零空間で解く)。
    `lam` が None なら GCV ($n\,\mathrm{RSS}/(n-\mathrm{edf})^2$ 最小) で $10^{-12}$〜$10^{2}$ の 1/4 桁刻みから選ぶ。
    戻り (BSpline, 記録 dict)。"""
    from scipy.interpolate import BSpline
    from scipy.linalg import null_space
    x, r = np.asarray(x, dtype=float), np.asarray(r, dtype=float)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    xb = np.interp(np.linspace(0.0, s[-1], n_int + 1), s, x)
    xb[0], xb[-1] = x[0], x[-1]
    t = np.concatenate([[x[0]] * deg, xb, [x[-1]] * deg])
    nc = len(t) - deg - 1
    basis = BSpline(t, np.eye(nc), deg, extrapolate=False)
    Bm = basis(x)
    gx, gw = np.polynomial.legendre.leggauss(8)
    xg = (0.5 * (xb[:-1] + xb[1:])[:, None] + 0.5 * (xb[1:] - xb[:-1])[:, None] * gx).ravel()
    wg = (0.5 * (xb[1:] - xb[:-1])[:, None] * gw).ravel()
    D3 = basis.derivative(3)(xg)
    Om = D3.T @ (D3 * wg[:, None])
    C = np.vstack([basis(np.array([x[0]])), basis.derivative(1)(np.array([x[0]])),
                   basis(np.array([x[-1]])), basis.derivative(1)(np.array([x[-1]]))])
    dv = np.array([start[0], start[1], end[0], end[1]], dtype=float)
    c0 = np.linalg.lstsq(C, dv, rcond=None)[0]
    Z = null_space(C)
    BZ = Bm @ Z
    G = BZ.T @ BZ
    W = Z.T @ Om @ Z
    y = r - Bm @ c0
    n = len(r)

    def solve(lm):
        A = G + lm * W
        a = np.linalg.solve(A, BZ.T @ y - lm * (Z.T @ Om @ c0))
        c = c0 + Z @ a
        res = Bm @ c - r
        edf = float(np.trace(BZ @ np.linalg.solve(A, BZ.T))) + C.shape[0]
        return c, res, edf, n * float(np.sum(res * res)) / (n - edf) ** 2
    if lam is None:
        grid = 10.0 ** np.arange(-12.0, 2.0 + 1e-9, 0.25)
        gcv = [solve(lm)[3] for lm in grid]
        lam = float(grid[int(np.argmin(gcv))])
        how = "gcv"
    else:
        how = "given"
    c, res, edf, g = solve(float(lam))
    sp = BSpline(t, c, deg, extrapolate=False)
    return sp, {"fit": "pspline", "deg": deg, "n_int": int(n_int), "n_points": int(n), "lam": float(lam), "lam_how": how,
                "edf": edf, "gcv": g, "resid_max": float(np.max(np.abs(res))), "resid_rms": float(np.sqrt(np.mean(res * res)))}


class RaoWall:
    r"""(M) の壁 $r(x)$: 円弧 kernel (曲率半径 $R_t$、$0\le x\le x_B$、閉じた式) + BE (生成器の壁点列の表現) + 終点の先は接線。

    `BE_xr`: B → E の壁点 (x, r)。先頭は厳密な B ($x_B=R_t\sin\theta_B$、$r_B=1+R_t(1-\cos\theta_B)$) に置き換え、末尾は
    呼び出し側が渡した終点をそのまま使う (生成器の E [DE の流量の根] を渡せば $(L,\sqrt\varepsilon)$ が厳密)。
    `x_end`: 壁の終わり (既定は点列の末尾。負対照のように点列の途中で切るときに与える)。

    BE の表現 (`fit`):

    - `"pspline"` (既定): 5 次 B-spline (節点は点列の折れ線の弧長に等間隔の `n_int` 区間) の罰則付き最小二乗
      $\sum_i(r(x_i)-r_i)^2+\lambda\int(r''')^2dx$ を、端の拘束 [$r(x_B)=r_B$・$r'(x_B)=\tan\theta_B$ (B で C1)、点列の末尾で
      終点の $r$ と勾配 $\tan\theta_E$] の下で解く ($\lambda$ は GCV で $10^{-12}$〜$10^{2}$ の 1/4 桁刻みから)。
      生成器の壁点は C⁺ 線上の流量の線形補間 (`wall_from_massflow`) なので位置に $O(h^2)$ の雑音 (n 121 で rms 1e-5) があり、
      B の近くの密な点 (間隔 0.015) をそのまま補間すると曲率が ±0.1 振動し、前進 MOC が BDE の中で C⁻ 前線の交差を出す
      (2026-10-12 実測。生成器の Goursat 網には反転なし)。当てはめの記録は `fit_info` (λ・有効自由度・残差・曲率の範囲)。
      最終の物理壁の当てはめ (plan §5.1 #7) とは別物 (T4・T1a′ の評価用の滑らかな表現)。
    - `"interp"`: 点列を通る 3 次スプライン (両端の勾配を固定)。比較用。

    摂動 (`PerturbedWall`) のために BE の弧長 $s(x)$ ($x_B\to x_{end}$) を 3 次 Hermite の表で持つ ($ds/dx=\sqrt{1+r'^2}$ は厳密)。"""

    def __init__(self, R_t: float, theta_B: float, BE_xr, theta_E: float, x_end: float | None = None,
                 fit: str = "pspline", n_int: int = 48, lam: float | None = None, h_tab: float = 2e-3):
        from scipy.interpolate import CubicHermiteSpline
        self.R_t, self.theta_B, self.theta_E = float(R_t), float(theta_B), float(theta_E)
        self.x_B = self.R_t * np.sin(self.theta_B)
        self.r_B = 1.0 + self.R_t * (1.0 - np.cos(self.theta_B))
        P = np.array(BE_xr, dtype=float)[:, :2].copy()
        self.B_replaced = (float(P[0, 0] - self.x_B), float(P[0, 1] - self.r_B))
        P[0] = (self.x_B, self.r_B)
        if np.any(np.diff(P[:, 0]) <= 0.0):
            raise ValueError("RaoWall: BE の x が単調増加でない")
        self.fit = fit
        if fit == "interp":
            self.spl = CubicSpline(P[:, 0], P[:, 1], bc_type=((1, np.tan(self.theta_B)), (1, np.tan(self.theta_E))))
            self.fit_info = {"fit": "interp", "n_points": int(len(P))}
        elif fit == "pspline":
            self.spl, self.fit_info = _pspline_fit(P[:, 0], P[:, 1], (self.r_B, np.tan(self.theta_B)),
                                                   (float(P[-1, 1]), np.tan(self.theta_E)), int(n_int), lam)
        else:
            raise ValueError(f"RaoWall: fit は 'pspline' か 'interp' ({fit!r})")
        self.dspl = self.spl.derivative()
        self.x_end = float(x_end) if x_end is not None else float(P[-1, 0])
        if not self.x_B < self.x_end <= P[-1, 0] + 1e-12:
            raise ValueError(f"RaoWall: x_end {self.x_end} が (x_B, 点列の末尾] の外")
        self.r_end = float(self.spl(self.x_end))
        self.rp_end = float(self.dspl(self.x_end))
        self.theta_end = float(np.arctan(self.rp_end))
        # BE の弧長の表
        nt = max(int(np.ceil((self.x_end - self.x_B) / h_tab)), 16)
        xt = np.linspace(self.x_B, self.x_end, nt + 1)
        xg, wg = _gl(xt[:-1], xt[1:])
        ds = np.sum(np.sqrt(1.0 + self.dspl(xg) ** 2) * wg, axis=-1)
        st = np.concatenate([[0.0], np.cumsum(ds)])
        self.S = float(st[-1])
        self._s = CubicHermiteSpline(xt, st, np.sqrt(1.0 + self.dspl(xt) ** 2))
        xs = np.linspace(self.x_B, self.x_end, 20001)
        rp, rpp = self.dspl(xs), self.spl.derivative(2)(xs)
        kap = rpp / (1.0 + rp * rp) ** 1.5
        dk = np.diff(kap)
        self.fit_info.update({"kappa_B_plus": float(kap[0]), "kappa_min": float(kap.min()), "kappa_max": float(kap.max()),
                              "n_dkappa_sign_changes": int(np.sum(np.diff(np.sign(dk[np.abs(dk) > 1e-10])) != 0)),
                              "theta_max_BE_deg": float(np.rad2deg(np.arctan(rp.max()))),
                              "end_r_minus_point": float(self.spl(P[-1, 0]) - P[-1, 1])})

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        xa = np.atleast_1d(x)
        out = np.empty_like(xa)
        R = self.R_t
        k = xa <= self.x_B
        if np.any(k):
            xk = xa[k]
            root = np.sqrt(np.maximum(R * R - xk * xk, 0.0))
            out[k] = (1.0 + R - root) if deriv == 0 else xk / root
        b = (~k) & (xa <= self.x_end)
        if np.any(b):
            out[b] = self.spl(xa[b]) if deriv == 0 else self.dspl(xa[b])
        e = xa > self.x_end
        if np.any(e):
            out[e] = (self.r_end + self.rp_end * (xa[e] - self.x_end)) if deriv == 0 else self.rp_end
        if deriv not in (0, 1):
            raise ValueError("RaoWall.r: deriv は 0 か 1")
        return out.reshape(np.shape(x))

    def t_of_x(self, x):
        """BE の弧長の正規化 t ∈ [0, 1] (B → 終点)。"""
        return np.asarray(self._s(np.clip(np.asarray(x, dtype=float), self.x_B, self.x_end)), dtype=float) / self.S

    def dt_dx(self, x):
        return np.sqrt(1.0 + self.dspl(np.asarray(x, dtype=float)) ** 2) / self.S


def perturbation_family(name: str, k: int):
    r"""摂動の族 $\eta_k(t)$ と $\eta_k'(t)$ ($t\in[0,1]$、両端で 0)。`"sin"`: $\sin(k\pi t)$ (端で勾配が変わる — B に角が
    できる)。`"sin_c1"`: $\sin(\pi t)\sin(k\pi t)$ (端で勾配も 0 — B・E の接線連続を保つ)。"""
    k = int(k)
    if name == "sin":
        return (lambda t: np.sin(k * np.pi * t)), (lambda t: k * np.pi * np.cos(k * np.pi * t))
    if name == "sin_c1":
        return ((lambda t: np.sin(np.pi * t) * np.sin(k * np.pi * t)),
                (lambda t: np.pi * np.cos(np.pi * t) * np.sin(k * np.pi * t) + k * np.pi * np.sin(np.pi * t) * np.cos(k * np.pi * t)))
    raise ValueError(f"perturbation_family: {name!r} は 'sin' か 'sin_c1'")


class PerturbedWall:
    r"""(M) の壁 BE に $\epsilon\,\eta_k(t(x))$ を半径方向 ($x$ 一定) に足した壁。B・終点 ($t=0,1$) で $\eta_k=0$ なので
    B と終点の位置 (したがって $L$・$\varepsilon$) は変わらず、kernel と初期線 (流量) も同じ。終点の先は元の壁の接線。"""

    def __init__(self, base: RaoWall, family: str, k: int, eps: float):
        self.base, self.family, self.k, self.eps = base, family, int(k), float(eps)
        self.eta, self.deta = perturbation_family(family, k)
        self.x_B, self.x_end = base.x_B, base.x_end

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        xa = np.atleast_1d(x)
        out = np.asarray(self.base.r(xa, deriv), dtype=float).copy()
        m = (xa > self.x_B) & (xa < self.x_end)
        if np.any(m) and self.eps != 0.0:
            t = self.base.t_of_x(xa[m])
            if deriv == 0:
                out[m] += self.eps * self.eta(t)
            else:
                out[m] += self.eps * self.deta(t) * self.base.dt_dx(xa[m])
        return out.reshape(np.shape(x))


N_BEYOND = 2                          # forward_thrust: 出口 x = L の先で補間に使う壁点の数


def _wall_dod_check(res: dict, x_last: float) -> dict:
    r"""前進 MOC の網で、推力に使った壁点 ($x\le$ `x_last`) の依存域に網の反転・単位過程の失敗が無いか。
    点 $(k,j)$ は $A=(k-1,j-1)$ (同じ C⁺ の対角線 $a=k-j$) と $B=(k,j+1)$ ($a-1$) に依るので、壁点 $(K,n)$ の依存域は
    添字で $\{k'\le K,\ k'-j'\le K-n\}$。$K$ は使った最後の壁点の前線。反転の定義は `moc_forward` のゲートと同じ
    (C⁺ 沿い・C⁻ 沿いの符号付きの長さ ≤ 0)。軸の近くの内部衝撃波 (DE より下で C⁻ が軸へ集まる — Rao 輪郭の場で x 15〜18 r_t)
    は壁の依存域の外なので、ゲートが不受理でも壁の圧力は使える — その区別をここで記録する。"""
    from .moc_forward import FAIL_STATES, ST_OK, _signed_len
    net = res["net"]
    n = res["n"]
    X, R, TH, MM = (net.a[f] for f in ("x", "r", "th", "M"))
    S = net.status
    ok = S == ST_OK
    MU = mu_vec(np.where(np.isfinite(MM), MM, 1.0))
    kw = np.flatnonzero(ok[:, n] & (X[:, n] <= x_last + 1e-12))
    K = int(kw[-1])
    kk, jj = np.meshgrid(np.arange(S.shape[0]), np.arange(S.shape[1]), indexing="ij")
    D = (kk <= K) & (kk - jj <= K - n)
    both_p = ok[1:, 1:] & ok[:-1, :-1]
    d_p = np.where(both_p, _signed_len(X[1:, 1:], X[:-1, :-1], 0.5 * (TH[:-1, :-1] + TH[1:, 1:]) + 0.5 * (MU[:-1, :-1] + MU[1:, 1:])),
                   np.nan)
    rev_p = np.zeros(S.shape, dtype=bool)
    rev_p[1:, 1:] = both_p & ~(d_p > 0.0)
    both_m = ok[:, 1:] & ok[:, :-1]
    e_m = np.where(both_m, _signed_len(X[:, :-1], X[:, 1:], 0.5 * (TH[:, 1:] + TH[:, :-1]) - 0.5 * (MU[:, 1:] + MU[:, :-1])), np.nan)
    rev_m = np.zeros(S.shape, dtype=bool)
    rev_m[:, :-1] = both_m & ~(e_m > 0.0)
    fail = np.isin(S, FAIL_STATES)
    bad = (rev_p | rev_m | fail) & D
    allbad = rev_p | rev_m | fail
    out = {"clean": not bool(np.any(bad)), "K_front": K, "n_bad_in_dod": int(np.sum(bad)), "n_bad_total": int(np.sum(allbad)),
           "n_missing_in_dod": int(np.sum(D & (S == 1)))}
    if np.any(allbad):
        xb = np.where(allbad, X, np.nan)
        out["bad_x_min"] = float(np.nanmin(xb))
        out["bad_r_max"] = float(np.nanmax(np.where(allbad, R, np.nan)))
    if np.any(bad):
        ib = np.argwhere(bad)[:5]
        out["bad_in_dod"] = [{"k": int(a), "j": int(b), "x": float(X[a, b]), "r": float(R[a, b])} for a, b in ib]
    return out


def forward_thrust(wall, L: float, gas, line: dict, *, p_a: float = 0.0, splits=(), ext_frac: float = 0.1,
                   gate: bool = True) -> dict:
    r"""与えた壁の推力を、壁を与えた前進 MOC (`moc_forward.forward_moc`、analytic・converge) で場を解き直して求める:
    $F$ = 初期線を横切る運動量 + 圧力 (`line_integrals`、生成器の `thrust_accounts` の `T_line` と同じ) +
    壁の圧力の軸力 $\int_{x_0}^L(p-p_a)\,2\pi r\,r'\,dx$ ($x_0$ = 初期線の壁点、スロートなら 0)。

    壁の圧力は前進 MOC の壁点の $p$ (等エントロピー、`Isentrope.state`) を $x$ の 3 次スプラインにし、壁の幾何 $r(x)$・$r'(x)$ は
    与えた壁の式のまま、壁点の間ごとの Gauss–Legendre で積分する。`splits` ($x$ の列) で区間を分け、区間ごとに別の
    スプラインにする (壁の勾配・曲率が跳ぶ点 — (M) の B — を跨いで補間しない。区間の端は各側の点列からの外挿)。
    出口 $x=L$ を跨ぐ壁点を得るため、前進 MOC は壁を接線で延ばして $L(1+$`ext_frac`$)$ まで回す ($x\le L$ の場は延長に
    依らない。出口近くの壁点の間隔は数 $r_t$ になるので、跨ぐ壁点が無ければ延長を倍にして回し直す。補間に使う $L$ の先の
    壁点は `N_BEYOND` 個まで)。`gate` なら同じ壁・初期線で $L$ までの前進 MOC をもう 1 本回してゲートの判定 (`gate_record`)
    を返し、$x\le L$ の壁点が 2 本で一致するか (`nodes_max_abs_diff`) も記録する。

    推力に使った壁点の依存域に網の反転・失敗が無いかを `dod_clean` (`_wall_dod_check`) に返す (ゲートの不受理が壁の依存域の
    外 [軸の近くの内部衝撃波] だけなら、壁の圧力は使える)。
    戻り dict: `ok`・`why`、`F` (= `T_line` + `T_wall`)、`T_line`・`mdot_line`・`T_wall`・`T_wall_seg` (区間ごと)・
    `T_wall_trap` (壁点の台形 [r で] の検算、出口は補間点)、`n_nodes_le_L`・`x_last_node_le_L`・`gate` (dict)・
    `mdot_err_max`・`end_reason`・`time_s`。"""
    from .moc_forward import forward_moc, gate_record, wall_points
    t0 = time.time()
    ig = as_isentrope(gas)
    L = float(L)
    ew = ExtendedWall(wall, L)
    ext = float(ext_frac)
    for _ in range(6):                      # 出口近くの壁点の間隔は数 r_t になるので、L を跨ぐ壁点が取れるまで延長を倍に
        Lx = L * (1.0 + ext)
        res = forward_moc(ew, Lx, ig.moc, line)
        out = {"ok": False, "L": L, "L_eval": Lx, "end_reason": (res.get("end") or {}).get("reason")}
        if res.get("net") is None:
            out["why"] = f"前進 MOC が初期線で止まった ({res['gate'].get('reasons')})"
            return out
        w = wall_points(res)
        w = w[np.argsort(w[:, 0])]
        if len(w) >= 8 and w[-1, 0] > L:
            break
        if out["end_reason"] != "exit_reached":
            break
        ext *= 2.0
    if not (len(w) >= 8 and w[-1, 0] > L):
        out["why"] = (f"壁点が x = L を越えない (壁点 {len(w)}、最後 x {w[-1, 0] if len(w) else np.nan:.6g}、"
                      f"終わり {res.get('end')})")
        return out
    p = ig.state(np.maximum(w[:, 4], 1.0))["p"]
    x0 = float(line["x"][-1])                     # 初期線の壁点 (スロートの壁なら 0)
    edges = [x0] + sorted(float(s) for s in splits if x0 < float(s) < L) + [L]
    # 最後の区間のスプラインに使う L より先の壁点は N_BEYOND 個まで (それより先は延長の場で、x ≤ L の補間に要らない)
    i_beyond = np.flatnonzero(w[:, 0] > L)
    use = np.ones(len(w), dtype=bool)
    use[i_beyond[N_BEYOND:]] = False
    w, p = w[use], p[use]
    segs = []
    T_wall = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        last = b >= L
        sel = (w[:, 0] >= a - 1e-14) & ((w[:, 0] <= b) if not last else np.ones(len(w), dtype=bool))
        xs, ps = w[sel, 0], p[sel]
        if len(xs) < 4:
            out["why"] = f"区間 [{a:.5g}, {b:.5g}] の壁点が {len(xs)} 個 (4 個以上要る)"
            return out
        sp = CubicSpline(xs, ps)
        knots = np.concatenate([[a], xs[(xs > a) & (xs < b)], [b]])
        xg, wg = _gl(knots[:-1], knots[1:])
        rr = np.asarray(wall.r(xg.ravel(), 0), dtype=float).reshape(xg.shape)
        rp = np.asarray(wall.r(xg.ravel(), 1), dtype=float).reshape(xg.shape)
        Ts = float(np.sum((sp(xg) - p_a) * TWO_PI * rr * rp * wg))
        segs.append({"a": a, "b": b, "n_nodes": int(np.sum((xs >= a) & (xs <= b))), "T": Ts,
                     "p_a_end": float(sp(a)), "p_b_end": float(sp(b))})
        T_wall += Ts
    m_line, T_line = line_integrals(ig, line["x"], line["r"], line["theta"], line["M"], p_a, spline=True)
    # 検算: 壁点の台形 (r で、出口は最後の区間のスプラインで補間した点)
    inL = w[:, 0] <= L
    rr_t = np.append(w[inL, 1], float(np.asarray(wall.r(np.array([L]), 0)).ravel()[0]))
    pp = np.append(p[inL], segs[-1]["p_b_end"])
    T_trap = float(np.sum(0.5 * ((pp[1:] - p_a) * rr_t[1:] + (pp[:-1] - p_a) * rr_t[:-1]) * TWO_PI * np.diff(rr_t)))
    dod = _wall_dod_check(res, float(w[-1, 0]))
    out.update({"ok": True, "why": None, "F": T_line + T_wall, "T_line": T_line, "mdot_line": m_line, "T_wall": T_wall,
                "dod_clean": dod["clean"], "dod": dod,
                "T_wall_seg": segs, "T_wall_trap": T_trap, "F_trap": T_line + T_trap, "n_nodes_le_L": int(np.sum(inL)),
                "x_last_node_le_L": float(w[inL, 0][-1]), "n_fronts_eval": res.get("n_fronts"),
                "mdot_err_max_eval": res["gate"]["indicators"].get("mdot_err_max")})
    if gate:
        rg = forward_moc(ew, L, ig.moc, line)
        out["gate"] = gate_record(rg)
        out["mdot_err_max"] = rg["gate"]["indicators"].get("mdot_err_max") if rg.get("net") is not None else None
        if rg.get("net") is not None:
            wg_ = wall_points(rg)
            wg_ = wg_[np.argsort(wg_[:, 0])]
            ne = min(len(wg_), int(np.sum(inL)))
            out["nodes_max_abs_diff"] = float(np.max(np.abs(wg_[:ne] - w[:ne]))) if ne else None
            out["n_nodes_gate_run"] = int(len(wg_))
    out["time_s"] = time.time() - t0
    return out


# ================================================================================================
# 整合な前線の再標本化 — T4 の評価器の初期線 (plan §5.1 #5 (h)、diagnostician 2026-10-12)
def resampled_front_line(net, k0: int, N: int) -> dict:
    r"""前進 MOC の網 `net` (`moc_forward.forward_moc` の `res["net"]`) の前線 $k_0$ (壁点 → 軸の C⁻) を、$r$ 等間隔の
    $N+1$ 点 ($0\le r\le r_w$) へ 3 次補間した初期線 (軸→壁、`forward_moc` の `line` の形)。

    **目的**: T4 (停留性) の推力の評価で、BE の壁節点の密度を初期線の点数 $N$ で上げる。Hall の初期線は級数の不整合
    (~1e-2 rad) のため $R_t\le1$ で $n\le321$ に制限される (`moc_forward.hall_initial_line` のガード)。そこで n 321 の Hall
    初期線から前進 MOC を回し、その**網の前線** (マーチの適合式を満たす整合な C⁻) を新しい初期線として細かく取り直す
    (diagnostician 2026-10-12、notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md)。前線は
    スロートの近く (壁点が円弧 kernel の上、B より上流) を選ぶので、BE の摂動は前線にも初期線の流量にも効かない。

    補間: $x(r)$・$\theta(r)$・$M(r)$ をそれぞれ前線の点 ($r$ は非等間隔) を通る 3 次スプライン (not-a-knot) にし、$r$ 等間隔の
    点で評価する (補間誤差 $O(h^4)$、$h$ = 前線の点の間隔)。両端は前線の値をそのまま使う (軸端 $r=0,\ \theta=0$、壁端は壁点
    そのもの — 壁が閉じた式なら `_check_line` の壁の検査 [`WALL_ON_TOL`] にそのまま乗る)。軸端の $\theta_r$ は
    `net.thr_axis[k0]` (前進 MOC が軸点に持たせた代用値 `axis_thr_neighbor`)。

    前線の点の状態が軸 ($j=0$) から壁 ($j=n$) まで全部 `ST_OK` でないとき、$r$ が狭義単調増加でないとき、軸端が $r=0$ で
    ないとき、`thr_axis` が非有限のときは `ValueError`。戻り dict: `x`・`r`・`M`・`theta`・`thr_axis` と来歴
    (`source` = `"resampled_front"`・`k0`・`N`・`n_front` [前線の区間数]・`front_dr_min`・`front_dr_max`・`x_wall`・`r_wall`・
    `x_axis`)。"""
    from .moc_forward import ST_OK, STATE_NAMES
    k0, N = int(k0), int(N)
    if N < 3:
        raise ValueError(f"resampled_front_line: N = {N} (3 以上)")
    S = net.status
    if not 0 <= k0 < S.shape[0]:
        raise ValueError(f"resampled_front_line: 前線 k0 = {k0} が網 (前線 0..{S.shape[0] - 1}) の外")
    st = S[k0]
    if not np.all(st == ST_OK):
        bad = np.flatnonzero(st != ST_OK)
        names = sorted({STATE_NAMES.get(int(st[j]), int(st[j])) for j in bad})
        raise ValueError(f"resampled_front_line: 前線 {k0} の状態が軸まで OK でない ({len(bad)} 点、j {bad[:5].tolist()}…、{names})")
    x, r, th, M = (np.asarray(net.a[f][k0], dtype=float) for f in ("x", "r", "th", "M"))
    thr = float(net.thr_axis[k0])
    if not np.isfinite(thr):
        raise ValueError(f"resampled_front_line: 前線 {k0} の軸端の θ_r が非有限")
    if r[0] != 0.0:
        raise ValueError(f"resampled_front_line: 前線 {k0} の軸端が r = 0 でない (r = {r[0]:.3e})")
    dr = np.diff(r)
    if not np.all(dr > 0.0):
        raise ValueError(f"resampled_front_line: 前線 {k0} の r が軸→壁で狭義単調増加でない")
    rq = np.linspace(0.0, float(r[-1]), N + 1)
    out = {"r": rq}
    for name, v in (("x", x), ("theta", th), ("M", M)):
        q = np.asarray(CubicSpline(r, v)(rq), dtype=float)
        q[0], q[-1] = v[0], v[-1]
        out[name] = q
    out.update({"thr_axis": thr, "source": "resampled_front", "k0": k0, "N": N, "n_front": int(len(r) - 1),
                "front_dr_min": float(dr.min()), "front_dr_max": float(dr.max()), "x_wall": float(x[-1]), "r_wall": float(r[-1]),
                "x_axis": float(x[0])})
    return out
