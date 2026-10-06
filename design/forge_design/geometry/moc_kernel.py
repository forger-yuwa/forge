"""kernel MOC: 円弧スロート下流の特性曲線法マーチ (軸対称・非回転・等エントロピー)。

目的: スロート円弧 (R_d, theta_a) が支配する kernel 領域の軸上マッハ分布と、
kernel 境界 x_k (円弧終端 P_a からの最終 C- が軸に当たる点) での M, M', M''
(= 目標 Bézier のアンカー) を、設計変数のみの決定的関数として与える
(親計画 §4.7(b))。

特性方程式 (導出は plans §4.7 / 標準教科書 Zucrow-Hoffman):
    C±: dr/dx = tan(θ±μ)
    C+ に沿って d(θ-ν) = -δ (sinμ sinθ / r) dσ+
    C- に沿って d(θ+ν) = +δ (sinμ sinθ / r) dσ-
δ=1 (軸対称)。マーチは「連続する C- フロント」方式: 各フロント = 壁点から軸まで
の 1 本の C-。壁点は前フロントの壁直下点からの C+ で決める。予測子-修正子
(係数平均) の 2 次精度単位プロセス。
"""
from __future__ import annotations

import numpy as np


# --- Prandtl-Meyer -----------------------------------------------------------
def _is_gas(g) -> bool:
    """`g` が γ (float) でなくガスモデル (`forge_design.gas.GasSemiPerfect` 等) か。
    MOC カーネル全体で「γ を渡す場所にガスモデルも渡せる」規約 (semi-perfect 対応、
    2026-08-17)。float なら従来どおり完全気体の閉形式。"""
    return hasattr(g, "nu") and hasattr(g, "mach_of_nu")


def pm_nu(M, g=1.4):
    if _is_gas(g):
        return g.nu(M)
    M = np.asarray(M, dtype=float)
    a = np.sqrt((g + 1.0) / (g - 1.0))
    b = np.sqrt(np.maximum(M * M - 1.0, 0.0))
    return a * np.arctan(b / a) - np.arctan(b)


def pm_mach(nu, g=1.4, tol=1e-12):
    """ν → M (Newton)。"""
    if _is_gas(g):
        return float(g.mach_of_nu(float(nu)))
    nu = float(nu)
    M = 1.0 + max(nu, 1e-6)  # 初期値
    for _ in range(60):
        f = pm_nu(M, g) - nu
        df = np.sqrt(M * M - 1.0) / (1.0 + 0.5 * (g - 1.0) * M * M) / M
        dM = f / max(df, 1e-14)
        M -= dM
        if M <= 1.0:
            M = 1.0 + 1e-9
        if abs(dM) < tol:
            break
    return M


def dnu_dM(M, g=1.4, order: int = 1):
    r"""Prandtl–Meyer 関数の $M$ 微分 ($d\nu/dM$ または $d^2\nu/dM^2$)。

    計画: plans/active/tooling-nozzle-axislaw-smoothness.md (軸則 C: 非負 $d\nu/dx$
    B-spline が Hall アンカー $M'_A, M''_A$ を $q(x_A)=\nu_M M'_A$,
    $q'(x_A)=\nu_{MM}(M'_A)^2+\nu_M M''_A$ に chain rule 変換するのに使う)。

    - CPG (`g` が float、または `GasCPG` — テーブルを持たず定数 γ のみのガスモデル):
      閉形式 $d\nu/dM=\sqrt{M^2-1}/[M(1+\frac{\gamma-1}2M^2)]$ (`pm_mach` の Newton
      反復が使う式と同一、既に検証済み)。$d^2\nu/dM^2$ は $h(M)=M(1+\frac{\gamma-1}2M^2)$,
      $b=\sqrt{M^2-1}$ として $[Mh-(M^2-1)h']/(bh^2)$ (h' = 1+\frac32(\gamma-1)M^2)。
    - semi-perfect (`g` が内部 $(M,\nu)$ テーブル `_M`/`_nu` を持つガスモデル):
      そのテーブルを 3 次スプラインにしてその導関数を返す (数値差分でなく spline 微分。
      テーブルにキャッシュを 1 個だけ持たせる)。
    """
    if hasattr(g, "_M") and hasattr(g, "_nu"):
        spl = getattr(g, "_dnu_dM_spline", None)
        if spl is None:
            from scipy.interpolate import CubicSpline
            spl = CubicSpline(g._M, g._nu)
            g._dnu_dM_spline = spl
        return spl(np.asarray(M, dtype=float), order)
    gamma = float(getattr(g, "gamma_ref", g))    # GasCPG (テーブルなし) or 素の float
    M = np.asarray(M, dtype=float)
    b = np.sqrt(np.maximum(M * M - 1.0, 0.0))
    h = M * (1.0 + 0.5 * (gamma - 1.0) * M * M)
    if order == 1:
        return np.where(b > 1e-12, b / np.maximum(h, 1e-300), 0.0)
    if order == 2:
        hp = 1.0 + 1.5 * (gamma - 1.0) * M * M
        num = M * h - (M * M - 1.0) * hp
        return np.where(b > 1e-9, num / np.maximum(b * h * h, 1e-300), 0.0)
    raise ValueError("dnu_dM: order は 1 か 2")


def axis_theta_r(M, Mp, g=1.4):
    r"""軸上の解析極限 $\theta_r=\lim_{r\to0}\sin\theta/r=\partial\theta/\partial r|_{r=0}$。

    計画: plans/active/discretization-moc-axis-limit-and-corrector.md §3・§4.1。
    軸近傍の質量保存 $2\rho u\,\theta_r=-d(\rho u)/dx$ と、等エントロピー流の
    $d\nu=\sqrt{M^2-1}\,d\ln u$・$d\ln(\rho u)=(1-M^2)\,d\ln u$ から

    $$\theta_r=\tfrac12\sqrt{M^2-1}\,\nu_M(M)\,M'(x)$$

    ($M'$ は軸上の $dM/dx$)。比熱一定なら $\tfrac12(M^2-1)M'/[M(1+\frac{\gamma-1}2M^2)]$。
    semi-perfect は `dnu_dM` (MOC 本体と同じ ν(M) 表の 3 次スプライン微分。本体の ν と逆関数は
    線形補間なので、同じ表でも補間関数は違う — その差は熱力学近似誤差として扱う)。
    $M'$ の代わりに軸上の $d\nu/dx$ を持っているときは `0.5*sqrt(M²-1)*dnu_dx` を直接使い、
    $\nu_M$ を掛けないこと (`moc_inverse.axis_theta_r_init` の代替経路)。"""
    M = np.asarray(M, dtype=float)
    return 0.5 * np.sqrt(np.maximum(M * M - 1.0, 0.0)) * dnu_dM(M, g) * np.asarray(Mp, dtype=float)


def _mu(M):
    return np.arcsin(1.0 / max(M, 1.0 + 1e-12))


class _Pt:
    __slots__ = ("x", "r", "th", "nu", "M")

    def __init__(self, x, r, th, nu, g, M=None):
        self.x = x
        self.r = r
        self.th = th
        self.nu = nu
        # M は ν から一意に決まるが、既に計算済みなら渡して Newton を省ける
        # (ベクトル充填が 100 万点規模の _Pt を作るため — 2026-08-15)
        self.M = pm_mach(nu, g) if M is None else M


# 軸に寄りすぎた点は θ の誤差が 1/r で増幅されるため、相手点から極限を代用する。
# しきい値は「相手点の半径に対する比」(スケール不変)。0 = 常に自身で評価。
AXIS_LIMIT_FRAC = 0.05


def _sin_over_r(p: "_Pt", other: "_Pt"):
    r"""軸対称項 $\sin\theta/r$ を **$p$ 自身**で評価する。

    適合式の源項は特性線分の両端 ($A$ と $P$ など) での値を台形則で平均する。
    したがってここは**その点自身の値**でなければならない。

    **2026-08-16 の修正**: 旧実装は「$p$ と $other$ のうち軸から遠い方」で評価して
    いた (軸上の $0/0$ 回避が目的)。しかし非軸点でも隣接点の値へ置換されるため
    台形則の片端がずれ、**スキーム全体が 1 次精度に落ちていた**。放射源流の厳密解で
    実測 (壁誤差, C⁺ 流束閉包):

    | n_axis | 旧 (遠い方で評価) | 新 (自身で評価) |
    | --- | --- | --- |
    | 140 | 1.131e-3 | 1.06e-4 |
    | 1120 | 1.329e-4 | 3.3e-6 |
    | 観測次数 | 1.02 | 1.6–1.7 |

    軸上 ($r\approx0$) だけは $0/0$ なので $other$ から極限を代用する
    (ここが残る 1 次要因。連続の式から $\partial\theta/\partial r|_{r=0}
    = -\frac12 d\ln F(M_{\rm axis})/dx$ を与えれば除ける)。
    逆 MOC の配列版 (`interior_vec`) には `axis_limit="analytic"` でこの解析極限
    (`axis_theta_r`) を入れた (2026-10-07, plans/active/discretization-moc-axis-limit-and-corrector.md)。
    本スカラー版 (kernel マーチ `KernelMOC` 専用) は従来のまま。
    """
    if p.r > 1e-9 and p.r >= AXIS_LIMIT_FRAC * other.r:
        return np.sin(p.th) / p.r
    if other.r > 1e-9:
        return np.sin(other.th) / other.r
    return 0.0


def trace_cminus(field_M, field_th, x_w: float, r_w: float, n: int, g: float):
    """解析場 (M(x,r), θ(x,r)) の中で C- を壁点から軸まで逆トレースする。

    初期値線を鉛直線でなく C- 特性線にするための前処理。鉛直初期線だと最初の
    C- フロントとの楔領域を長い C+ セグメントで跨ぎ O(Δ²) バイアスが解像度に
    依らず残る (放射流照合で ~1% を確認)。r を独立変数に RK2。
    戻り値: [(x, r, th, nu)] 軸→壁の順。
    """
    pts = [(x_w, r_w)]
    x, r = x_w, r_w
    dr = -r_w / n
    for _ in range(n):
        M1 = field_M(x, r)
        th1 = field_th(x, r)
        s1 = 1.0 / np.tan(th1 - _mu(M1))
        xm = x + s1 * dr * 0.5
        rm = r + dr * 0.5
        M2 = field_M(xm, max(rm, 1e-12))
        th2 = field_th(xm, max(rm, 1e-12))
        x = x + dr / np.tan(th2 - _mu(M2))
        r = r + dr
        pts.append((x, max(r, 0.0)))
    out = []
    for (x, r) in pts[::-1]:
        th = 0.0 if r < 1e-12 else float(field_th(x, r))
        out.append((x, r, th, float(pm_nu(field_M(x, r), g))))
    return out


class KernelMOC:
    """円弧 (または一般壁曲線) 下の kernel を C- フロントでマーチする。"""

    def __init__(self, gamma=1.4, delta=1.0, n_corr=2):
        self.g = gamma
        self.delta = delta
        self.n_corr = n_corr

    # -- 単位プロセス --------------------------------------------------------
    def _interior(self, A: _Pt, B: _Pt) -> _Pt:
        """A (下, C+ 担体) と B (上, C- 担体) から新点。"""
        g, d = self.g, self.delta
        thP, nuP = 0.5 * (A.th + B.th), 0.5 * (A.nu + B.nu)
        xP = rP = None
        for _ in range(1 + self.n_corr):
            MP = pm_mach(nuP, g)
            muP = _mu(MP)
            mp = np.tan(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP))
            mm = np.tan(0.5 * (B.th + thP) - 0.5 * (_mu(B.M) + muP))
            if abs(mm - mp) < 1e-12:
                raise RuntimeError("特性線が平行 (M~1?)")
            xP = (A.r - B.r + mm * B.x - mp * A.x) / (mm - mp)
            rP = A.r + mp * (xP - A.x)
            Pn = _Pt(xP, max(rP, 0.0), thP, max(nuP, 0.0), g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, B) + np.sin(muP) * _sin_over_r(Pn, A))
            fB = 0.5 * (np.sin(_mu(B.M)) * _sin_over_r(B, A) + np.sin(muP) * _sin_over_r(Pn, B))
            Sp = -d * fA / np.cos(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP)) * (xP - A.x)
            Sm = +d * fB / np.cos(0.5 * (B.th + thP) - 0.5 * (_mu(B.M) + muP)) * (xP - B.x)
            Jp = (A.th - A.nu) + Sp   # θP - νP
            Jm = (B.th + B.nu) + Sm   # θP + νP
            thP, nuP = 0.5 * (Jp + Jm), 0.5 * (Jm - Jp)
        return _Pt(xP, rP, thP, nuP, self.g)

    def _axis(self, B: _Pt) -> _Pt:
        """B から C- で軸 (r=0, θ=0) へ。"""
        g, d = self.g, self.delta
        nuP = B.th + B.nu
        xP = B.x
        for _ in range(1 + self.n_corr):
            MP = pm_mach(nuP, g)
            muP = _mu(MP)
            mm = np.tan(0.5 * B.th - 0.5 * (_mu(B.M) + muP))
            xP = B.x - B.r / mm
            Pn = _Pt(xP, 0.0, 0.0, nuP, g)
            fB = 0.5 * (np.sin(_mu(B.M)) * _sin_over_r(B, B) + np.sin(muP) * _sin_over_r(Pn, B))
            Sm = +d * fB / np.cos(0.5 * B.th - 0.5 * (_mu(B.M) + muP)) * (xP - B.x)
            nuP = (B.th + B.nu) + Sm
        return _Pt(xP, 0.0, 0.0, nuP, self.g)

    def _wall_prescribed(self, s: float, wall_xy, wall_th, front) -> _Pt:
        """規定の壁ステーション s の壁点。C+ の足を前フロント折れ線上に補間で求める。

        壁ステップを規定側で制御する (前フロントの点間隔任せだと近壁解像度で
        ステップが決まり kernel 全域で数フロントしか取れず破綻する)。
        """
        g, d = self.g, self.delta
        xw, rw = wall_xy(s)
        thw = wall_th(s)
        nuW = max(front[-1].nu + (thw - front[-1].th), 1e-8)  # 2D 単純波の初期推定
        A, i_foot = front[-2], len(front) - 2
        for _ in range(2 + self.n_corr):
            muW = _mu(pm_mach(nuW, g))
            m = np.tan(0.5 * (A.th + thw) + 0.5 * (_mu(A.M) + muW))
            # W から後方 C+ 直線 r = rw + m (x - xw) と front 折れ線の交点
            A, i_foot = self._foot_on_front(front, xw, rw, m)
            Pn = _Pt(xw, rw, thw, nuW, g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, Pn) + np.sin(muW) * _sin_over_r(Pn, A))
            Sp = -d * fA / np.cos(0.5 * (A.th + thw) + 0.5 * (_mu(A.M) + muW)) * (xw - A.x)
            nuW = max(thw - (A.th - A.nu) - Sp, 1e-8)
        return _Pt(xw, rw, thw, nuW, g), i_foot

    def _foot_on_front(self, front, xw, rw, m):
        """直線 r = rw + m (x - xw) と front 折れ線 (軸→壁) の交点を線形補間。

        戻り値 (点, セグメント下端 index)。interior ペアリングは index 以下の
        front 点のみ使う (足より上の点の C+ は新 C- でなく壁に当たるため)。
        """
        for i in range(len(front) - 2, -1, -1):
            P0, P1 = front[i], front[i + 1]
            dx, drr = P1.x - P0.x, P1.r - P0.r
            den = drr - m * dx
            if abs(den) < 1e-14:
                continue
            t = (rw + m * (P0.x - xw) - P0.r) / den
            if -1e-9 <= t <= 1.0 + 1e-9:
                t = min(max(t, 0.0), 1.0)
                return _Pt(P0.x + t * dx, P0.r + t * drr,
                           P0.th + t * (P1.th - P0.th), P0.nu + t * (P1.nu - P0.nu), self.g), i
        # 交点なし → 壁直下点で代用 (ステップが極端に小さい場合)
        return front[-2], len(front) - 2

    def _wall(self, A: _Pt, wall_xy, wall_th, s0: float) -> tuple:
        """A から C+ で壁曲線 (param s) へ。戻り (点, s)。"""
        g, d = self.g, self.delta
        s = s0
        thP = wall_th(s)
        nuP = A.nu
        for _ in range(2 + self.n_corr):
            MP = pm_mach(max(nuP, 1e-8), g)
            muP = _mu(MP)
            mp = np.tan(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP))
            # C+ 直線と壁曲線の交点 (s を Newton)
            for _ in range(30):
                xw, rw = wall_xy(s)
                f = rw - (A.r + mp * (xw - A.x))
                h = 1e-7
                xw2, rw2 = wall_xy(s + h)
                df = (rw2 - (A.r + mp * (xw2 - A.x)) - f) / h
                if abs(df) < 1e-14:
                    break
                ds = f / df
                s -= ds
                if abs(ds) < 1e-13:
                    break
            xw, rw = wall_xy(s)
            thP = wall_th(s)
            Pn = _Pt(xw, rw, thP, max(nuP, 1e-8), g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, Pn) + np.sin(muP) * _sin_over_r(Pn, A))
            Sp = -d * fA / np.cos(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP)) * (xw - A.x)
            nuP = thP - (A.th - A.nu) - Sp
        return _Pt(xw, rw, thP, nuP, self.g), s

    # -- kernel マーチ (円弧壁) ----------------------------------------------
    def march_arc(self, throat, theta_a: float, n_start: int = 41, n_wall: int = 80):
        """throat: SauerThroat。円弧終端角 theta_a まで規定壁ステーションでマーチ。

        手順: (1) Sauer 有効域内の鉛直 starting line (軸上 M=1.05) を初期値線に、
        **三角形充填** (interior 列 + 軸反射) で初期壁点からの C- (楔境界) を
        細かい折れ線として構築 (Sauer 場を有効域外へ外挿しない)。
        (2) その C- をフロントに壁ステーションマーチ (θa まで)。
        戻り値 dict: axis_x, axis_M, x_k, anchor (M, M', M''), n_fronts。
        """
        g = self.g
        Rd = throat.R

        def wall_xy(s):
            return Rd * np.sin(s), 1.0 + Rd * (1.0 - np.cos(s))

        def wall_th(s):
            return s

        x0, rr, MM, tt = throat.starting_line(M_start=1.05, n=n_start)
        if x0 >= 0.8 * Rd * np.sin(theta_a):
            raise ValueError(
                f"遷音速パッチ (x0={x0:.3f}) が円弧終端に迫る: R={Rd} は Sauer 一次の"
                "適用域外 (Kliegel-Levine 高次が必要 — 既知の制約)")
        s_w0 = float(np.arcsin(min(x0 / Rd, 1.0)))
        front = [_Pt(x0, rr[i], tt[i], float(pm_nu(MM[i], g)), g) for i in range(n_start)]
        front[-1].th = wall_th(s_w0)
        axis_pts = [(front[0].x, front[0].M)]

        # (1) 楔充填: 初期値線の各点 k から出る C- を下から順にフロント化する
        #     (front_k = init[k] から軸までの C-)。最後のフロント = 初期壁点からの
        #     C- = 楔境界で、以後の壁ステーションマーチの初期フロントになる。
        wedge = [front[0]]
        for k in range(1, n_start):
            new = [front[k]]
            for j in range(len(wedge) - 1, -1, -1):
                try:
                    P = self._interior(wedge[j], new[-1])
                except RuntimeError:
                    continue
                if P.r < 0.0 or P.x < min(wedge[j].x, new[-1].x) - 1e-9:
                    continue
                new.append(P)
            ax = self._axis(new[-1])
            new.append(ax)
            wedge = new[::-1]
            axis_pts.append((ax.x, ax.M))
        front = wedge

                # (2) 壁ステーションマーチ
        stations = np.linspace(s_w0, theta_a, n_wall + 1)[1:]
        front = self._march_stations(front, stations, wall_xy, wall_th, axis_pts)
        axis_x = np.array([p[0] for p in axis_pts])
        axis_M = np.array([p[1] for p in axis_pts])
        x_k = float(axis_x[-1])
        n_fit = max(6, int(0.3 * len(axis_x)))
        c = np.polyfit(axis_x[-n_fit:] - x_k, axis_M[-n_fit:], 3)
        anchor = (float(np.polyval(c, 0.0)),
                  float(np.polyval(np.polyder(c), 0.0)),
                  float(np.polyval(np.polyder(c, 2), 0.0)))
        return {"axis_x": axis_x, "axis_M": axis_M, "x_k": x_k, "anchor": anchor,
                "n_fronts": len(axis_x)}

    def _march_stations(self, front, stations, wall_xy, wall_th, axis_pts):
        """規定壁ステーション列で C- フロントを順に構築する共通ループ。"""
        for s_k in stations:
            W, i_foot = self._wall_prescribed(float(s_k), wall_xy, wall_th, front)
            new = [W]
            for j in range(i_foot, -1, -1):
                try:
                    P = self._interior(front[j], new[-1])
                except RuntimeError:
                    continue
                if P.r < 0.0 or P.x < min(front[j].x, new[-1].x) - 1e-9:
                    continue
                new.append(P)
            ax = self._axis(new[-1])
            new.append(ax)
            front = new[::-1]
            axis_pts.append((ax.x, ax.M))
        return front

    # -- 検証用: 一般壁 (円錐等) --------------------------------------------
    def march_wall(self, front_pts, wall_xy, wall_th, stations):
        """既製フロントから一般壁下を規定ステーションでマーチ (テスト用)。"""
        g = self.g
        front = [_Pt(*p, g) for p in front_pts]  # (x, r, th, nu)
        axis_pts = [(front[0].x, front[0].M)]
        self._march_stations(front, np.asarray(stations, dtype=float), wall_xy, wall_th, axis_pts)
        return np.array(axis_pts)


# --- ベクトル化ユーティリティ (fill の O(n²) スカラーループ解消, 2026-08-15) -----
def pm_mach_vec(nu, g=1.4, tol=1e-13, iters=60):
    """ν → M の配列版 Newton (`pm_mach` と同じ反復・同じ初期値)。"""
    if _is_gas(g):
        return np.asarray(g.mach_of_nu(np.asarray(nu, dtype=float)), dtype=float)
    nu = np.asarray(nu, dtype=float)
    M = 1.0 + np.maximum(nu, 1e-6)
    for _ in range(iters):
        f = pm_nu(M, g) - nu
        df = np.sqrt(np.maximum(M * M - 1.0, 0.0)) / (1.0 + 0.5 * (g - 1.0) * M * M) / M
        dM = f / np.maximum(df, 1e-14)
        M = M - dM
        M = np.where(M <= 1.0, 1.0 + 1e-9, M)
        if np.all(np.abs(dM) < tol):
            break
    return M


def mu_vec(M):
    return np.arcsin(1.0 / np.maximum(M, 1.0 + 1e-12))


def _sin_over_r_vec(r_p, th_p, r_o, th_o):
    """`_sin_over_r` の配列版 (点自身で評価、軸上のみ相手から極限を代用)。"""
    on_axis = (r_p <= 1e-9) | (r_p < AXIS_LIMIT_FRAC * r_o)
    ra = np.where(on_axis, r_o, r_p)
    tha = np.where(on_axis, th_o, th_p)
    ok = ra > 1e-9
    return np.where(ok, np.sin(tha) / np.where(ok, ra, 1.0), 0.0)


# --- 軸上の解析極限と予測修正の収束 (2026-10-07) ---------------------------------
# 計画: plans/active/discretization-moc-axis-limit-and-corrector.md §4.1・§4.2。
# 問題 YAML の geometry.moc_axis_limit / geometry.moc_corrector で選ぶ。既定 (legacy + fixed2) は従来とビット同一。
AXIS_R_EPS = 1e-9            # 軸上の判定 r ≤ これ (`_sin_over_r_vec` と同じ値)
AXIS_LIMIT_MODES = ("legacy", "analytic")
CORRECTOR_MODES = ("fixed2", "converge")
CORR_TOL = 1e-12             # converge: 修正 1 回あたりの θ・ν の更新量の上限 [rad]
CORR_MAX = 50                # converge: 修正子の上限回数
RESID_TOL = 1e-10            # converge: 最終状態の残差 (幾何の交点式 [r_t]・適合式 [rad]) の合格上限

# 対の分類 (`interior_vec(..., stats=)` の stats["status"]、§4.2 の 5 分類)。
# 入力時点の対象 = 全対。そのうち:
PAIR_CONVERGED = 0           # 収束 (converge)
PAIR_MISSING = 1             # もともとの欠損 (入力が非有限 — fill_levels の死んだ対)
PAIR_GEOM = 2                # 幾何的棄却: 特性線が平行 (|tan(θ−μ) − tan(θ+μ)| < 1e-12)
PAIR_NONFINITE = 3           # 反復の失敗: 反復中に NaN・Inf
PAIR_MAXITER = 4             # 反復の失敗: 上限到達
PAIR_DONE = 5                # fixed2: 固定回数を回し終えた (収束は判定しない)
PAIR_BELOW_AXIS = 6          # 幾何的棄却: 新点が軸より下 (r < −1e-12、呼び出し側の判定)


def _src_vec(r_p, th_p, r_o, th_o, thr_p, analytic: bool):
    r"""適合式の源項 $\sin\theta/r$ の端点評価。`analytic=False` は `_sin_over_r_vec` そのもの (従来とビット同一)。

    `analytic=True`: **真の軸端点** ($r_p\le$ `AXIS_R_EPS`) で $\theta_r$ (`thr_p`、有限) が与えられていれば
    それを使う (相手の値を借りない)。軸外の点と $\theta_r$ を持たない点は従来どおり
    (`AXIS_LIMIT_FRAC` の分岐を含む。真の軸端点の判定が先)。"""
    v = _sin_over_r_vec(r_p, th_p, r_o, th_o)
    if not analytic or thr_p is None:
        return v
    return np.where((r_p <= AXIS_R_EPS) & np.isfinite(thr_p), thr_p, v)


def _interior_pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, thrA, thrB, thP, nuP,
                   g, delta: float, analytic: bool):
    """予測子・修正子の 1 回分 (P の現在値 thP・nuP → 新しい xP・rP・thP・nuP)。
    `analytic=False` の演算は旧 `interior_vec` のループ本体と同一 (同じ式・同じ順序)。
    戻り: xP, rP, thP, nuP, good (従来の判定 |den| ≥ 1e-12。den が NaN でも False), parallel (den が有限で
    |den| < 1e-12 = 特性線が平行。NaN の den は平行に数えない — 反復中の非有限と区別するため)。"""
    MP = pm_mach_vec(nuP, g)
    muP = mu_vec(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    mp, mm = np.tan(ang_p), np.tan(ang_m)
    den = mm - mp
    good = np.abs(den) >= 1e-12
    parallel = np.isfinite(den) & ~good
    den = np.where(good, den, 1.0)
    xP = (Ar - Br + mm * Bx - mp * Ax) / den
    rP = Ar + mp * (xP - Ax)
    rPc = np.maximum(rP, 0.0)
    fA = 0.5 * (np.sin(muA) * _src_vec(Ar, Ath, Br, Bth, thrA, analytic)
                + np.sin(muP) * _src_vec(rPc, thP, Ar, Ath, None, analytic))
    fB = 0.5 * (np.sin(muB) * _src_vec(Br, Bth, Ar, Ath, thrB, analytic)
                + np.sin(muP) * _src_vec(rPc, thP, Br, Bth, None, analytic))
    Sp = -delta * fA / np.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / np.cos(ang_m) * (xP - Bx)
    Jp = (Ath - Anu) + Sp
    Jm = (Bth + Bnu) + Sm
    return xP, rP, 0.5 * (Jp + Jm), 0.5 * (Jm - Jp), good, parallel


def _interior_residual(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB, thrA, thrB,
                       xP, rP, thP, nuP, g, delta: float, analytic: bool):
    r"""最終状態 $(x_P,r_P,\theta_P,\nu_P)$ で幾何の交点式と適合式を再評価した残差 (§4.2)。

    幾何: $r_P-[r_A+\tan(\bar\theta_+ +\bar\mu_+)(x_P-x_A)]$ と C⁻ 側 (長さ, $r_t$)。
    適合式: $(\theta_P\mp\nu_P)-[(\theta\mp\nu)_{A,B}+S_\pm]$ ($S_\pm$ は最終値の源項, rad)。
    戻り: (幾何の残差の絶対値の大きい方, 適合式の残差の絶対値の大きい方)。"""
    MP = pm_mach_vec(nuP, g)
    muP = mu_vec(MP)
    ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
    ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
    mp, mm = np.tan(ang_p), np.tan(ang_m)
    rg = np.maximum(np.abs(rP - (Ar + mp * (xP - Ax))), np.abs(rP - (Br + mm * (xP - Bx))))
    rPc = np.maximum(rP, 0.0)
    fA = 0.5 * (np.sin(muA) * _src_vec(Ar, Ath, Br, Bth, thrA, analytic)
                + np.sin(muP) * _src_vec(rPc, thP, Ar, Ath, None, analytic))
    fB = 0.5 * (np.sin(muB) * _src_vec(Br, Bth, Ar, Ath, thrB, analytic)
                + np.sin(muP) * _src_vec(rPc, thP, Br, Bth, None, analytic))
    Sp = -delta * fA / np.cos(ang_p) * (xP - Ax)
    Sm = +delta * fB / np.cos(ang_m) * (xP - Bx)
    rc = np.maximum(np.abs((thP - nuP) - ((Ath - Anu) + Sp)), np.abs((thP + nuP) - ((Bth + Bnu) + Sm)))
    return rg, rc


def source_branches(Ar, Br, rP, thrA=None, thrB=None, analytic: bool = False) -> dict:
    r"""源項の端点評価が**どの分岐に入ったか** (対ごと、最終状態の幾何で判定)。記録専用 (値は変えない)。

    - `axis_analytic`: 真の軸端点 (A・B) で $\theta_r$ を使った数 (0〜2)
    - `axis_borrow`: 真の軸端点で相手 (軸外) の値を借りた数 — legacy の代用
    - `axis_zero`: 真の軸端点で相手も軸上なので 0 にした数 — legacy
    - `frac_A`・`frac_B`・`frac_PA`・`frac_PB`: 軸外の点 p が `AXIS_LIMIT_FRAC`·r_相手 より軸に近く、
      相手の値を借りた (A|B、B|A、P|A、P|B の各評価)
    - `P_on_axis`: 新点 P が軸上 ($r_P\le$ `AXIS_R_EPS`。θ_r を持たないので従来どおり相手を借りる)"""
    rPc = np.maximum(rP, 0.0)
    out = {}
    nA = Ar <= AXIS_R_EPS
    nB = Br <= AXIS_R_EPS
    hA = nA & (np.isfinite(thrA) if (analytic and thrA is not None) else False)
    hB = nB & (np.isfinite(thrB) if (analytic and thrB is not None) else False)
    out["axis_analytic"] = hA.astype(np.int8) + hB.astype(np.int8)
    out["axis_borrow"] = ((nA & ~hA & (Br > AXIS_R_EPS)).astype(np.int8)
                          + (nB & ~hB & (Ar > AXIS_R_EPS)).astype(np.int8))
    out["axis_zero"] = ((nA & ~hA & ~(Br > AXIS_R_EPS)).astype(np.int8)
                        + (nB & ~hB & ~(Ar > AXIS_R_EPS)).astype(np.int8))
    out["frac_A"] = (Ar > AXIS_R_EPS) & (Ar < AXIS_LIMIT_FRAC * Br)
    out["frac_B"] = (Br > AXIS_R_EPS) & (Br < AXIS_LIMIT_FRAC * Ar)
    out["frac_PA"] = (rPc > AXIS_R_EPS) & (rPc < AXIS_LIMIT_FRAC * Ar)
    out["frac_PB"] = (rPc > AXIS_R_EPS) & (rPc < AXIS_LIMIT_FRAC * Br)
    out["P_on_axis"] = rPc <= AXIS_R_EPS
    return out


def _check_corr_params(tol, max_corr):
    """converge の許容差 (有限・正) と上限回数 (1 以上の整数) を検査する。"""
    if isinstance(tol, bool) or not isinstance(tol, (int, float, np.floating)) or not (np.isfinite(tol) and tol > 0.0):
        raise ValueError(f"corrector='converge' の tol は有限の正の数 (受け取った値: {tol!r})")
    if isinstance(max_corr, bool) or not isinstance(max_corr, (int, np.integer)) or max_corr < 1:
        raise ValueError(f"corrector='converge' の max_corr は 1 以上の整数 (受け取った値: {max_corr!r})")


def _finite_all(*arrs):
    m = np.isfinite(arrs[0])
    for a in arrs[1:]:
        m = m & np.isfinite(a)
    return m


def interior_vec(Ax, Ar, Ath, Anu, AM, Bx, Br, Bth, Bnu, BM,
                 g: float, delta: float, n_corr: int, *, axis_limit: str = "legacy",
                 thrA=None, thrB=None, corrector: str = "fixed2", tol: float = CORR_TOL,
                 max_corr: int = CORR_MAX, stats: dict | None = None):
    r"""`KernelMOC._interior` の**配列版** (A=C⁺担体, B=C⁻担体 の対を一括処理)。

    スカラー版と同一の予測子-修正子・同一の係数平均・同一の軸対称源項評価。
    戻り: (xP, rP, thP, nuP, ok) — ok=False は特性線が平行/非有限 (converge では未収束も) で棄却すべき対。

    **選択肢** (plans/active/discretization-moc-axis-limit-and-corrector.md §4。既定はビット同一):

    - `axis_limit`: `"legacy"` = 真の軸端点の $\sin\theta/r$ は相手の値で代用 (相手も軸上なら 0) /
      `"analytic"` = その点の解析極限 $\theta_r$ (`thrA`・`thrB`、対ごと。軸外の点は NaN でよい)。
      既知の軸端点の $\theta_r$ は予測・修正を通じて固定 (反復中の θ/r から更新しない)。
      `analytic` で入力が有限の軸端点に $\theta_r$ が無ければ ValueError。
    - `corrector`: `"fixed2"` = 予測 1 回 + 修正 `n_corr` 回 (従来) / `"converge"` = 修正 1 回の
      θ・ν の更新量が全て `tol` 以下になるまで (上限 `max_corr` 回)。対ごとに、収束した時点の値で止める
      (同じ段の他の対の反復回数に結果が依存しない)。停止判定は入力が有限の対だけで行う。
    - `stats` (dict) を渡すと書き込む: `status` (対ごとの分類 `PAIR_*`)、`n_iter` (修正子の回数、
      反復の失敗は失敗時点)、`resid_geom`・`resid_comp` (converge で収束した対の最終残差、他は NaN)、
      `branch` (`source_branches`)。
    """
    if axis_limit not in AXIS_LIMIT_MODES:
        raise ValueError(f"axis_limit は {AXIS_LIMIT_MODES} のどれか (受け取った値: {axis_limit!r})")
    if corrector not in CORRECTOR_MODES:
        raise ValueError(f"corrector は {CORRECTOR_MODES} のどれか (受け取った値: {corrector!r})")
    if corrector == "converge":
        _check_corr_params(tol, max_corr)
    analytic = axis_limit == "analytic"
    fin_in = _finite_all(Ax, Ar, Ath, Anu, AM, Bx, Br, Bth, Bnu, BM)
    if analytic:
        thrA = np.broadcast_to(np.asarray(np.nan if thrA is None else thrA, dtype=float), np.shape(Ath))
        thrB = np.broadcast_to(np.asarray(np.nan if thrB is None else thrB, dtype=float), np.shape(Bth))
        lack = fin_in & (((Ar <= AXIS_R_EPS) & ~np.isfinite(thrA)) | ((Br <= AXIS_R_EPS) & ~np.isfinite(thrB)))
        if np.any(lack):
            raise ValueError(f"axis_limit='analytic': 軸上の既知点 {int(np.sum(lack))} 対に θ_r が無い "
                             "(初期前線の軸端点には axis_thr を渡すこと)")
    else:
        thrA = thrB = None
    thP = 0.5 * (Ath + Bth)
    nuP = 0.5 * (Anu + Bnu)
    muA, muB = mu_vec(AM), mu_vec(BM)
    if corrector == "fixed2":
        # 従来の固定回数 (旧 `interior_vec` と同じ式・同じ順序)
        xP = np.zeros_like(thP)
        rP = np.zeros_like(thP)
        ok = np.ones(thP.shape, dtype=bool)
        par = np.zeros(thP.shape, dtype=bool)
        for _ in range(1 + n_corr):
            xP, rP, thP, nuP, good, parallel = _interior_pass(Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB,
                                                              thrA, thrB, thP, nuP, g, delta, analytic)
            ok &= good
            par |= parallel
        fin_out = np.isfinite(xP) & np.isfinite(rP) & np.isfinite(thP) & np.isfinite(nuP)
        if stats is not None:
            st = np.full(thP.shape, PAIR_DONE, dtype=np.int8)
            st[fin_in & (~ok | ~fin_out)] = PAIR_NONFINITE      # 反復中に非有限 (平行でない対の den の NaN を含む)
            st[fin_in & par] = PAIR_GEOM
            st[~fin_in] = PAIR_MISSING
            stats.update(status=st, n_iter=np.full(thP.shape, int(n_corr), dtype=np.int16),
                         resid_geom=np.full(thP.shape, np.nan), resid_comp=np.full(thP.shape, np.nan),
                         branch=source_branches(Ar, Br, rP, thrA, thrB, analytic))
        ok &= fin_out
        return xP, rP, thP, nuP, ok

    # --- converge: 対ごとに収束まで (有効な対だけで判定、収束した対はその時点で凍結) ---
    shape = np.shape(thP)
    flat = [np.ravel(np.broadcast_to(a, shape)).astype(float) for a in
            (Ax, Ar, Ath, Anu, Bx, Br, Bth, Bnu, muA, muB)]
    tA = None if thrA is None else np.ravel(np.broadcast_to(thrA, shape)).astype(float)
    tB = None if thrB is None else np.ravel(np.broadcast_to(thrB, shape)).astype(float)
    thP, nuP = np.ravel(thP).astype(float), np.ravel(nuP).astype(float)
    n = thP.size
    xP = np.full(n, np.nan)
    rP = np.full(n, np.nan)
    status = np.full(n, -1, dtype=np.int8)
    n_it = np.zeros(n, dtype=np.int16)
    status[~np.ravel(np.broadcast_to(fin_in, shape))] = PAIR_MISSING
    act = np.flatnonzero(status < 0)
    for k in range(int(max_corr) + 1):              # k = 0 予測子、k ≥ 1 修正子 k 回目
        if act.size == 0:
            break
        sub = [a[act] for a in flat]
        xs, rs, ths, nus, good, par = _interior_pass(*sub, None if tA is None else tA[act],
                                                     None if tB is None else tB[act],
                                                     thP[act], nuP[act], g, delta, analytic)
        fin = np.isfinite(xs) & np.isfinite(rs) & np.isfinite(ths) & np.isfinite(nus)
        upd = np.maximum(np.abs(ths - thP[act]), np.abs(nus - nuP[act]))
        xP[act], rP[act], thP[act], nuP[act] = xs, rs, ths, nus
        n_it[act] = k
        st = np.full(act.size, -1, dtype=np.int8)
        st[~par & (~good | ~fin)] = PAIR_NONFINITE      # 反復中に非有限 (den の NaN を含む)
        st[par] = PAIR_GEOM
        if k >= 1:
            st[good & fin & (upd <= tol)] = PAIR_CONVERGED
        status[act] = st
        act = act[st < 0]
    status[act] = PAIR_MAXITER
    ok = status == PAIR_CONVERGED
    if stats is not None:
        rg = np.full(n, np.nan)
        rc = np.full(n, np.nan)
        if np.any(ok):
            i = np.flatnonzero(ok)
            rg[i], rc[i] = _interior_residual(*[a[i] for a in flat], None if tA is None else tA[i],
                                              None if tB is None else tB[i], xP[i], rP[i], thP[i], nuP[i],
                                              g, delta, analytic)
        stats.update(status=status.reshape(shape), n_iter=n_it.reshape(shape),
                     resid_geom=rg.reshape(shape), resid_comp=rc.reshape(shape),
                     branch={k: v.reshape(shape) for k, v in
                             source_branches(flat[1], flat[5], rP, tA, tB, analytic).items()})
    return xP.reshape(shape), rP.reshape(shape), thP.reshape(shape), nuP.reshape(shape), ok.reshape(shape)
