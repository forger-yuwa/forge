r"""平面 2D 非対称 2 壁 (SERN) の MOC と key point 逆設計。

plan: plans/active/tooling-nozzle-sern-chain.md §4.3–4.4。無次元 (p_in = 1, H = 1)、
完全気体 γ 一定 (semi-perfect は圧力比の写像を差し替えれば載る — 未対応)。

## kernel (前進 MOC、逆 [格子] 法)

x 一定の station を進め、各 station の y 格子点 P について C⁺/C⁻ を前 station へ逆トレースし、
足の不変量 J⁺ = θ−ν (C⁺ 上で保存)、K⁻ = θ+ν (C⁻ 上で保存) を線形補間して
θ_P = ½(J⁺+K⁻)、ν_P = ½(K⁻−J⁺)。足の角度は P と足の平均 (2 次精度)。
- ランプ角部扇 / カウル角部扇: J⁺ (K⁻) が入口値のままの点は「その扇の単純波域 or 一様域」
  なので**解析閉包** (θ−μ = atan2(y−1,x) 等) を使い、格子補間の初期スミアを避ける。
- カウル後縁 (TE) 以降の下境界は等圧自由境界 (ν 固定、θ は K⁻ から)。TE の膨張扇は上流の
  K⁻ が一様でないため解析閉包できない → **C⁺ レイ束** (各レイは J⁺ 保存の特性線、状態は局所
  K⁻ との組合せ) として station ごとに前進させ、足の補間データに合流させる (格子だけでは扇が
  最初の station で 1 セルに潰れて伝播しない — 2026-09-04 実測)。p_ext > p_TE (圧縮 = 衝撃)
  は表現できないので警告フラグ。
- 上壁 (ランプ) は kernel では全長 θ_r0 の直線。実際のランプは逆設計で a より下流を置換する。

## key point 逆設計 (平面)

制御面 = c から出る C⁺ (c–e)。平面では制御面上の状態が一様 (rao_planar.py) なので、
c を「入口からの質量流量比 f の流線上で M = M_c となる点」に置き、c から上流へ C⁻ を辿って
ランプ足 a を得る。a–c (C⁻, K⁻ = 一定) と c–e (C⁺, 一様) に挟まれた領域は K⁻ が全域で
等しい = C⁺ 族の単純波 (各 C⁺ は直線で状態一定)。壁 a–e はその中の流線: a–c 上の各点から
C⁺ の直線を張り、隣接レイ間を平均壁角で進めて交点を取る。e = c から出るレイと壁の交点。
dv = (M_c, f) が c を一意に決め、θ_c は従属 (Rao 最適は縁条件で θ_c を選ぶ特殊解)。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .moc_kernel import pm_nu, pm_mach_vec, mu_vec
from .rao_planar import (p_over_pt, mach_from_p_over_pt, massflux_ratio,
                         ideal_gross_thrust, inlet_stream_thrust)
from .sern_geometry import SernDesign

_TOL_INV = 1e-10


@dataclass
class SernKernelSpec:
    M_in: float
    theta_r0: float                 # ランプ初期角 [rad] (>0 で上へ膨張)
    theta_c0: float                 # カウル角 [rad] (>0 で下へ膨張)
    L_cowl: float
    gamma: float = 1.4
    p_ext_over_p_in: float | None = None   # None → カウル後縁圧に一致 (無波)
    x_max: float = 10.0
    nj: int = 301
    dx: float = 2e-3
    dx0: float = 1e-4
    dx_growth: float = 1.06
    n_te_rays: int = 64                     # TE 膨張扇のレイ本数
    te_ray_emit_every: int = 4              # 自由境界からの反射 C⁺ レイ射出間隔 [station]

    def __post_init__(self):
        for k in ("M_in", "L_cowl", "x_max", "dx", "dx0"):
            v = getattr(self, k)
            if not np.isfinite(v) or v <= 0.0:
                raise ValueError(f"SernKernelSpec: {k} = {v!r} は正の有限値でない")
        if self.M_in <= 1.0:
            raise ValueError("SernKernelSpec: M_in は超音速 (>1) が必要 (音速スロート給気は内部ノズルで接続)")
        if not (0.0 <= self.theta_r0 < 0.5 * np.pi) or not (-0.25 * np.pi < self.theta_c0 < 0.5 * np.pi):
            raise ValueError("SernKernelSpec: theta_r0 ∈ [0, π/2), theta_c0 ∈ (−π/4, π/2) が必要")
        if self.L_cowl > self.x_max:
            raise ValueError("SernKernelSpec: L_cowl ≤ x_max が必要")
        if self.nj < 21:
            raise ValueError("SernKernelSpec: nj ≥ 21")


def _solve_nu_minus_mu(target, g):
    """ν(M) − μ(M) = target を M について解く (単調増加、ベクトル二分法)。"""
    target = np.asarray(target, dtype=float)
    lo = np.full(target.shape, 1.0 + 1e-9)
    hi = np.full(target.shape, 400.0)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        v = pm_nu(mid, g) - mu_vec(mid)
        take = v < target
        lo = np.where(take, mid, lo)
        hi = np.where(take, hi, mid)
    return 0.5 * (lo + hi)


class PlanarMOC:
    def __init__(self, spec: SernKernelSpec):
        self.s = spec
        self.g = spec.gamma
        self.warnings: list[str] = []
        self._marched = False

    # ------------------------------------------------------------------ station 配置
    def _stations(self):
        s = self.s
        def graded(x_from, x_to):
            xs = [x_from]
            h = s.dx0
            while xs[-1] + h < x_to - 1e-12:
                xs.append(xs[-1] + h)
                h = min(h * s.dx_growth, s.dx)
            xs.append(x_to)
            return xs
        xs = graded(0.0, s.L_cowl)
        if s.L_cowl < s.x_max:
            xs = xs + graded(s.L_cowl, s.x_max)[1:]
        return np.asarray(xs, dtype=float)

    # ------------------------------------------------------------------ march
    def march(self, stop_at=None):
        """stop_at=(f, M_c): 質量流量比 f の流線を march と同時に積分し、その上で M ≥ M_c に達した station の
        少し先 (余裕 `stop_margin` = 0.3H) で march を打ち切る。設計に要る場は key point c の C⁻ の依存領域
        (c より上流) だけで、c より下流の kernel は使わない (2026-09-05、ユーザ指摘: 解析領域が必要以上に広い)。
        x_max は上限として残る。"""
        s, g = self.s, self.g
        X = self._stations()
        stop_margin = 0.3
        y_s = None; x_stop = None
        if stop_at is not None:
            f_s, M_stop = float(stop_at[0]), float(stop_at[1]); y_s = 1.0 - f_s
        n_st = len(X)
        nj = s.nj
        sj = np.linspace(0.0, 1.0, nj)
        Y = np.zeros((n_st, nj)); TH = np.zeros_like(Y); NU = np.zeros_like(Y); MM = np.zeros_like(Y)
        ylo = np.zeros(n_st); yup = np.zeros(n_st)
        nu_in = float(pm_nu(s.M_in, g))
        Jp_in, Km_in = -nu_in, nu_in
        # 入口 station
        Y[0] = sj; TH[0] = 0.0; NU[0] = nu_in; MM[0] = s.M_in; ylo[0] = 0.0; yup[0] = 1.0
        tan_r, tan_c = np.tan(s.theta_r0), np.tan(s.theta_c0)
        i_te = int(np.argmin(np.abs(X - s.L_cowl)))
        assert abs(X[i_te] - s.L_cowl) < 1e-12
        self.i_te = i_te
        nu_ext = None; th_b_prev = None; y_te = None
        rays_y = rays_J = rays_th = rays_mu = None
        n_since_emit = 0
        self._aug = None
        self.p_te_over_p_in = None
        for n in range(n_st - 1):
            x0, x1 = X[n], X[n + 1]; dx = x1 - x0
            y0 = Y[n]; th0 = TH[n]; nu0 = NU[n]; mu0 = mu_vec(MM[n])
            Jp0 = th0 - nu0; Km0 = th0 + nu0
            yup1 = 1.0 + x1 * tan_r
            lower_wall = x1 <= s.L_cowl + 1e-12
            if lower_wall:
                ylo1 = -x1 * tan_c
            else:
                if nu_ext is None:  # TE 直後: 自由境界の状態を確定
                    y_te = ylo[n]
                    p_te = float(p_over_pt(MM[n, 0], g) / p_over_pt(s.M_in, g))
                    self.p_te_over_p_in = p_te
                    p_ext = p_te if s.p_ext_over_p_in is None else float(s.p_ext_over_p_in)
                    self.p_ext_over_p_in = p_ext
                    if p_ext > p_te * (1.0 + 1e-9):
                        self.warnings.append(
                            f"p_ext/p_in={p_ext:.5f} > p_TE/p_in={p_te:.5f}: カウル後縁で圧縮 (衝撃) — MOC は等エントロピー圧縮で近似")
                    M_ext = float(mach_from_p_over_pt(p_ext * p_over_pt(s.M_in, g), g))
                    nu_ext = float(pm_nu(M_ext, g))
                    # 自由境界は TE から扇通過後の角度で出る (K⁻_TE − ν_ext)。扇前の壁角を
                    # 使うと最初の station で境界が浅く動き、扇の後端レイを境界外として殺す
                    th_b_prev = float(TH[n, 0] + NU[n, 0] - nu_ext)
                ylo1 = ylo[n] + dx * np.tan(th_b_prev)
            # 足の補間データ (TE 扇のレイ束があれば合流済みの配列)
            if self._aug is not None:
                y0, th0, mu0, Jp0, Km0 = self._aug
            th = TH[n].copy(); nu = NU[n].copy(); mu = mu_vec(MM[n])
            for _it in range(4):
                y1 = ylo1 + sj * (yup1 - ylo1)
                # C⁺ 足 (下左へ)
                ang_p = th + mu
                yA = y1 - dx * np.tan(ang_p)
                for _k in range(2):
                    thA = np.interp(yA, y0, th0); muA = np.interp(yA, y0, mu0)
                    yA = y1 - dx * np.tan(0.5 * (ang_p + thA + muA))
                # C⁻ 足 (上左へ)
                ang_m = th - mu
                yB = y1 - dx * np.tan(ang_m)
                for _k in range(2):
                    thB = np.interp(yB, y0, th0); muB = np.interp(yB, y0, mu0)
                    yB = y1 - dx * np.tan(0.5 * (ang_m + thB - muB))
                Jp = np.interp(yA, y0, Jp0)   # 端でクランプ = 壁/境界点の値
                Km = np.interp(yB, y0, Km0)
                # --- 解析閉包 (角部扇の単純波域) ---
                ramp_simple = np.abs(Jp - Jp_in) < _TOL_INV
                cowl_simple = np.abs(Km - Km_in) < _TOL_INV
                # K⁻ が入口値のままでも TE 扇 (C⁺ 源) の下流ではカウル扇の単純波域ではない。
                # TE 扇の先頭レイ (扇前状態の θ+μ) より上流側 (上方) に限定する
                if nu_ext is not None:
                    lead = float(TH[i_te, 0] + mu_vec(MM[i_te, 0]))
                    above_lead = np.arctan2(y1 - y_te, x1 - s.L_cowl) > lead
                    cowl_simple = cowl_simple & above_lead
                both = ramp_simple & cowl_simple
                if s.theta_r0 > 0.0 and ramp_simple.any():
                    phi = np.arctan2(y1 - 1.0, x1)
                    phi_lo = -float(mu_vec(s.M_in))
                    M_r = float(pm_mach_vec(np.array([nu_in + s.theta_r0]), g)[0])
                    phi_hi = s.theta_r0 - float(mu_vec(M_r))
                    Mf = _solve_nu_minus_mu(phi - Jp_in, g)
                    nuf = pm_nu(Mf, g)
                    nu_r = np.where(phi <= phi_lo, nu_in, np.where(phi >= phi_hi, nu_in + s.theta_r0, nuf))
                    th_r = Jp_in + nu_r
                    Km = np.where(ramp_simple, th_r + nu_r, Km)
                if s.theta_c0 > 0.0 and cowl_simple.any():
                    phi = np.arctan2(y1, x1)
                    phi_hi = float(mu_vec(s.M_in))
                    M_c0 = float(pm_mach_vec(np.array([nu_in + s.theta_c0]), g)[0])
                    phi_lo = -s.theta_c0 + float(mu_vec(M_c0))
                    Mf = _solve_nu_minus_mu(Km_in - phi, g)
                    nuf = pm_nu(Mf, g)
                    nu_c = np.where(phi >= phi_hi, nu_in, np.where(phi <= phi_lo, nu_in + s.theta_c0, nuf))
                    th_c = Km_in - nu_c
                    Jp = np.where(cowl_simple, th_c - nu_c, Jp)
                if both.any():  # 一様入口域
                    Jp = np.where(both, Jp_in, Jp); Km = np.where(both, Km_in, Km)
                # --- TE 膨張扇: TE station からの最初の march のみ解析閉包 (以降はレイ束) ---
                if (not lower_wall) and n == i_te:
                    below = yA < y0[0] + 1e-15
                    if below.any():
                        phi = np.arctan2(y1 - y_te, x1 - s.L_cowl)
                        Mf = _solve_nu_minus_mu(Km - phi, g)
                        nuf = pm_nu(Mf, g)
                        th_f = Km - nuf
                        Jp = np.where(below, th_f - nuf, Jp)
                th_new = 0.5 * (Jp + Km); nu_new = 0.5 * (Km - Jp)
                # --- 境界点 ---
                th_new[-1] = s.theta_r0; nu_new[-1] = s.theta_r0 - Jp[-1]
                if lower_wall:
                    th_new[0] = -s.theta_c0; nu_new[0] = Km[0] + s.theta_c0
                else:
                    nu_new[0] = nu_ext; th_new[0] = Km[0] - nu_ext
                    ylo1 = ylo[n] + dx * np.tan(0.5 * (th_b_prev + th_new[0]))
                nu_new = np.maximum(nu_new, 1e-9)
                M_new = pm_mach_vec(nu_new, g)
                th, nu, mu = th_new, nu_new, mu_vec(M_new)
            if not lower_wall:
                th_b_prev = th[0]
            y1 = ylo1 + sj * (yup1 - ylo1)
            Y[n + 1] = y1; TH[n + 1] = th; NU[n + 1] = nu; MM[n + 1] = M_new
            ylo[n + 1] = ylo1; yup[n + 1] = yup1
            # --- key point 到達判定 (f 流線を同時積分) ---
            if y_s is not None and x_stop is None:
                th_s0 = float(np.interp(y_s, y0, th0)); y_pred = y_s + dx * np.tan(th_s0)
                th_s1 = float(np.interp(y_pred, y1, th)); y_s = y_s + dx * np.tan(0.5 * (th_s0 + th_s1))
                M_s = float(np.interp(y_s, y1, M_new))
                if M_s >= M_stop:
                    x_stop = x1 + stop_margin
            if x_stop is not None and x1 >= x_stop:
                n_last = n + 1
                X = X[: n_last + 1]; Y = Y[: n_last + 1]; TH = TH[: n_last + 1]; NU = NU[: n_last + 1]; MM = MM[: n_last + 1]
                ylo = ylo[: n_last + 1]; yup = yup[: n_last + 1]
                break
            # --- TE 扇のレイ束 (C⁺, J⁺ 保存) を station n+1 へ前進し補間データに合流 ---
            if not lower_wall:
                if rays_y is None:  # TE 直後: 扇の状態範囲 [ν_TE, ν_ext] を nr 本に分割
                    nr = s.n_te_rays
                    nu_te = NU[i_te, 0]; K_te = TH[i_te, 0] + nu_te
                    nu_r = np.linspace(nu_te, nu_ext, nr)
                    rays_J = K_te - 2.0 * nu_r
                    th_r = K_te - nu_r; mu_r = mu_vec(pm_mach_vec(nu_r, g))
                    rays_y = y_te + dx * np.tan(th_r + mu_r)
                    rays_th, rays_mu = th_r, mu_r
                else:
                    ang = rays_th + rays_mu
                    yg = rays_y + dx * np.tan(ang)
                    for _k in range(2):
                        Kr = np.interp(yg, y1, th + nu)
                        thr = 0.5 * (rays_J + Kr); nur = np.maximum(0.5 * (Kr - rays_J), 1e-9)
                        mur = mu_vec(pm_mach_vec(nur, g))
                        yg = rays_y + dx * np.tan(0.5 * (ang + thr + mur))
                    rays_y, rays_th, rays_mu = yg, thr, mur
                # 自由境界からの反射 C⁺ (J⁺ = θ_b − ν_ext) も定期的にレイとして射出する。
                # 後端レイと境界に挟まれた楔の内部は、これが無いと格子点が古い値のまま残る
                # (両端にしか正しいデータが無い — 2026-09-04 実測)。
                n_since_emit += 1
                if n_since_emit >= s.te_ray_emit_every:
                    n_since_emit = 0
                    rays_y = np.append(rays_y, ylo1 + 1e-9); rays_J = np.append(rays_J, th[0] - nu_ext)
                    rays_th = np.append(rays_th, th[0]); rays_mu = np.append(rays_mu, mu[0])
                keep = (rays_y < yup1 - 1e-12) & (rays_y > ylo1 + 1e-12)
                rays_y, rays_J, rays_th, rays_mu = rays_y[keep], rays_J[keep], rays_th[keep], rays_mu[keep]
                if len(rays_y):
                    ya = np.concatenate([y1, rays_y]); order = np.argsort(ya, kind="stable")
                    tha = np.concatenate([th, rays_th])[order]; mua = np.concatenate([mu, rays_mu])[order]
                    Jpa = np.concatenate([th - nu, rays_J])[order]
                    Kma = np.concatenate([th + nu, 2.0 * rays_th - rays_J])[order]
                    self._aug = (ya[order], tha, mua, Jpa, Kma)
                else:
                    self._aug = None
        self.X, self.Y, self.TH, self.NU, self.M = X, Y, TH, NU, MM
        self.ylo, self.yup = ylo, yup
        self.stopped_at = x_stop
        if self.p_te_over_p_in is None:  # カウルが x_max まで
            self.p_ext_over_p_in = (float(s.p_ext_over_p_in) if s.p_ext_over_p_in is not None
                                    else float(p_over_pt(MM[-1, 0], g) / p_over_pt(s.M_in, g)))
        self._marched = True
        return self

    # ------------------------------------------------------------------ 場の問合せ
    def _bounds(self, x):
        i = int(np.clip(np.searchsorted(self.X, x), 1, len(self.X) - 1))
        t = (x - self.X[i - 1]) / (self.X[i] - self.X[i - 1])
        return i, float(np.clip(t, 0.0, 1.0))

    def y_up(self, x):
        return 1.0 + x * np.tan(self.s.theta_r0)

    def y_lo(self, x):
        i, t = self._bounds(x)
        return (1 - t) * self.ylo[i - 1] + t * self.ylo[i]

    def state_at(self, x, y):
        """(θ, ν, M) at (x,y)。station 間は x 線形、station 上は y 線形補間。"""
        i, t = self._bounds(x)
        out = []
        for k in (i - 1, i):
            yy = self.Y[k]
            out.append((np.interp(y, yy, self.TH[k]), np.interp(y, yy, self.NU[k])))
        th = (1 - t) * out[0][0] + t * out[1][0]
        nu = (1 - t) * out[0][1] + t * out[1][1]
        M = float(pm_mach_vec(np.array([max(nu, 1e-9)]), self.g)[0])
        return float(th), float(nu), M

    def trace(self, x0, y0, kind: str, backward: bool = False, ds: float = 1e-3, max_steps: int = 200000):
        """kind ∈ {'C+','C-','stream'}。x を独立変数に RK2。境界/範囲で停止し (n,4)[x,y,θ,M] を返す。
        終端は境界と線形に交差させる。"""
        sign = -1.0 if backward else 1.0
        def slope(x, y):
            th, _, M = self.state_at(x, y)
            mu = float(mu_vec(M))
            if kind == "C+":
                return np.tan(th + mu)
            if kind == "C-":
                return np.tan(th - mu)
            return np.tan(th)
        pts = []
        x, y = float(x0), float(y0)
        th, _, M = self.state_at(x, y)
        pts.append((x, y, th, M))
        for _ in range(max_steps):
            h = sign * ds
            k1 = slope(x, y)
            xm, ym = x + 0.5 * h, y + 0.5 * h * k1
            if xm < self.X[0] or xm > self.X[-1]:
                break
            k2 = slope(xm, ym)
            xn, yn = x + h, y + h * k2
            if xn < self.X[0] - 1e-15 or xn > self.X[-1] + 1e-15:
                break
            up, lo = self.y_up(xn), self.y_lo(xn)
            if yn >= up or yn <= lo:
                # 境界との交点 (線形)
                f0 = (y - self.y_up(x)) if yn >= up else (y - self.y_lo(x))
                f1 = (yn - up) if yn >= up else (yn - lo)
                a = f0 / (f0 - f1) if abs(f0 - f1) > 1e-300 else 1.0
                xb, yb = x + a * (xn - x), y + a * (yn - y)
                th, _, M = self.state_at(xb, yb)
                pts.append((xb, yb, th, M))
                return np.asarray(pts)
            x, y = xn, yn
            th, _, M = self.state_at(x, y)
            pts.append((x, y, th, M))
        return np.asarray(pts)

    # ------------------------------------------------------------------ key point
    def key_point(self, M_c: float, f: float):
        """入口から質量流量比 f (ランプ側) の流線上で最初に M ≥ M_c となる点 (x_c,y_c,θ_c)。"""
        if not (0.0 < f < 1.0):
            raise ValueError("f ∈ (0,1)")
        line = self.trace(0.0, 1.0 - f, "stream")
        Ms = line[:, 3]
        idx = np.argmax(Ms >= M_c)
        if Ms[idx] < M_c:
            raise ValueError(f"key point 不成立: 流線上の M は最大 {Ms.max():.4f} < M_c={M_c} (kernel を伸ばすか M_c を下げる)")
        if idx == 0:
            raise ValueError(f"M_c={M_c} ≤ 入口 M")
        a = (M_c - Ms[idx - 1]) / (Ms[idx] - Ms[idx - 1])
        p = line[idx - 1] + a * (line[idx] - line[idx - 1])
        return float(p[0]), float(p[1]), float(p[2]), line

    # ------------------------------------------------------------------ 逆設計
    def design_ramp(self, M_c: float, f: float, ds: float = 1e-3) -> SernDesign:
        if not self._marched:
            self.march()
        s, g = self.s, self.g
        x_c, y_c, th_c, stream = self.key_point(M_c, f)
        th_c, nu_c, M_c_act = self.state_at(x_c, y_c)
        # c から上流へ C⁻ → ランプ足 a
        cm = self.trace(x_c, y_c, "C-", backward=True, ds=ds)
        xa, ya = cm[-1, 0], cm[-1, 1]
        if abs(ya - self.y_up(xa)) > 1e-6 * max(1.0, ya):
            if xa < 1e-3 and abs(ya - 1.0) < 2e-3:
                raise ValueError(
                    f"c=({x_c:.3f},{y_c:.3f}) からの C⁻ がランプ角部 (0,1) に達した: c はランプ角部扇の"
                    " 最終レイより上流 (kernel の直線ランプが実現されない)。M_c を上げるか f を下げる")
            if xa <= self.X[0] + 1e-9:
                raise ValueError(
                    f"c=({x_c:.3f},{y_c:.3f}) からの C⁻ が入口面 (x=0, y={ya:.3f}) に達した: c はランプ角部扇の"
                    " 影響域外 (壁で制御できない)。M_c を上げるか f を下げる")
            raise ValueError(f"c からの C⁻ が上壁に届かず ({xa:.3f},{ya:.3f}) で停止。M_c/f を見直す")
        ac = cm[::-1]  # a → c の順
        # 単純波のレイ: 各点から C⁺ 直線
        xs, ys, ths, Ms = ac[:, 0], ac[:, 1], ac[:, 2], ac[:, 3]
        mus = mu_vec(Ms)
        beta = ths + mus
        W = np.zeros((len(ac), 2)); W[0] = (xa, ya)
        for i in range(len(ac) - 1):
            tb = 0.5 * (ths[i] + ths[i + 1])
            d = np.array([np.cos(tb), np.sin(tb)])
            o = np.array([xs[i + 1], ys[i + 1]]); r = np.array([np.cos(beta[i + 1]), np.sin(beta[i + 1])])
            A = np.array([[d[0], -r[0]], [d[1], -r[1]]])
            rhs = o - W[i]
            det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
            if abs(det) < 1e-14:
                raise ValueError("単純波レイと壁が平行 (退化)")
            t, u = np.linalg.solve(A, rhs)
            if t < 0:
                raise ValueError("壁構築で後退 (レイ順序が不正)")
            W[i + 1] = W[i] + t * d
        # 壁 = 角部 → (kernel 直線ランプ) → a → W → e
        xk = self.X[self.X <= xa]
        ramp_k = np.column_stack([xk, self.y_up(xk)])
        ramp_k_M = self.M[: len(xk), -1]; ramp_k_th = self.TH[: len(xk), -1]
        if len(xk) == 0 or xk[-1] < xa - 1e-12:
            tha, _, Ma = self.state_at(xa, ya)
            ramp_k = np.vstack([ramp_k, [xa, ya]]); ramp_k_M = np.append(ramp_k_M, Ma); ramp_k_th = np.append(ramp_k_th, s.theta_r0)
        ramp_xy = np.vstack([ramp_k, W[1:]])
        ramp_M = np.concatenate([ramp_k_M, Ms[1:]]); ramp_th = np.concatenate([ramp_k_th, ths[1:]])
        # カウル
        xc_ = self.X[self.X <= s.L_cowl + 1e-12]
        cowl_xy = np.column_stack([xc_, -xc_ * np.tan(s.theta_c0)])
        cowl_M = self.M[: len(xc_), 0]; cowl_th = self.TH[: len(xc_), 0]
        # 質量検算: c–e (一様、直線 Mach 線) を横切る流量 / 入口流量
        L_ce = float(np.hypot(W[-1, 0] - x_c, W[-1, 1] - y_c))
        f_chk = float(massflux_ratio(M_c_act, s.M_in, g)) * np.sin(float(mu_vec(M_c_act))) * L_ce
        return SernDesign(ramp_xy=ramp_xy, ramp_M=ramp_M, ramp_theta=ramp_th,
                          cowl_xy=cowl_xy, cowl_M=cowl_M, cowl_theta=cowl_th,
                          key_point=(x_c, y_c, M_c_act, th_c), foot_a=(float(xa), float(ya)),
                          lip_e=(float(W[-1, 0]), float(W[-1, 1])), mass_fraction=f,
                          mass_fraction_check=f_chk, p_ext_over_p_in=float(self.p_ext_over_p_in),
                          info={"warnings": list(self.warnings), "n_rays": int(len(ac)),
                                "theta_e": float(ths[-1]), "K_minus_spread": float(np.ptp(ths + pm_nu(Ms, g)))})


    # ------------------------------------------------------------------ 直線ランプ (検証用・逆設計なし)
    def straight_design(self, L_ramp: float) -> SernDesign:
        """kernel の直線ランプ (θ_r0) を x = L_ramp で切った「設計しない」形状。NASA TM X-71972 の
        平板 SERN のような検証用。壁状態は kernel の境界行から取る。"""
        if not self._marched:
            self.march()
        s = self.s
        if L_ramp > self.X[-1] + 1e-12:
            raise ValueError(f"L_ramp={L_ramp} > x_max={self.X[-1]}: kernel を伸ばす")
        m = self.X <= L_ramp + 1e-12
        xk = self.X[m]
        ramp_xy = np.column_stack([xk, self.y_up(xk)])
        ramp_M = self.M[: len(xk), -1]; ramp_th = self.TH[: len(xk), -1]
        xc_ = self.X[self.X <= s.L_cowl + 1e-12]
        cowl_xy = np.column_stack([xc_, -xc_ * np.tan(s.theta_c0)])
        the, _, Me = self.state_at(float(xk[-1]), float(self.y_up(xk[-1])))
        return SernDesign(ramp_xy=ramp_xy, ramp_M=ramp_M, ramp_theta=ramp_th, cowl_xy=cowl_xy,
                          cowl_M=self.M[: len(xc_), 0], cowl_theta=self.TH[: len(xc_), 0],
                          key_point=(float(xk[-1]), float(self.y_up(xk[-1])), float(Me), float(the)),
                          foot_a=(0.0, 1.0), lip_e=(float(xk[-1]), float(self.y_up(xk[-1]))),
                          mass_fraction=float("nan"), mass_fraction_check=float("nan"),
                          p_ext_over_p_in=float(self.p_ext_over_p_in),
                          info={"warnings": list(self.warnings), "n_rays": 0, "theta_e": float(s.theta_r0),
                                "K_minus_spread": 0.0, "mode": "straight"})


# ---------------------------------------------------------------------- 力積分
def wall_forces(design: SernDesign, M_in: float, g: float = 1.4, pa_over_pin: float | None = None,
                x_ref: float = 0.0, y_ref: float = 0.0) -> dict:
    """壁圧 (p − p_a) の積分から推力/揚力/ピッチモーメント (無次元: p_in=1, H=1)。
    推力 = 壁力の −x 成分 (前向き)。揚力 = +y。モーメントは頭上げ正 (= −z)。
    総推力 = 入口流れ推力 + 壁推力、C_T = 総推力 / 理想総推力 (p_a まで等エントロピー膨張)。"""
    pa = design.p_ext_over_p_in if pa_over_pin is None else float(pa_over_pin)
    p0 = float(p_over_pt(M_in, g))

    def integ(xy, M, upper: bool):
        p = p_over_pt(M, g) / p0 - pa
        d = np.diff(xy, axis=0)
        ln = np.hypot(d[:, 0], d[:, 1])
        n_out = np.column_stack([-d[:, 1], d[:, 0]]) if upper else np.column_stack([d[:, 1], -d[:, 0]])
        pm = 0.5 * (p[1:] + p[:-1])
        F = pm[:, None] * n_out              # 長さ込み (n_out は |dl| 倍)
        xm = 0.5 * (xy[1:] + xy[:-1])
        Mz = (xm[:, 0] - x_ref) * F[:, 1] - (xm[:, 1] - y_ref) * F[:, 0]
        return F.sum(axis=0), -Mz.sum(), ln.sum()

    Fr, Mr, Lr = integ(design.ramp_xy, design.ramp_M, True)
    Fc, Mc, Lc = integ(design.cowl_xy, design.cowl_M, False)
    F = Fr + Fc
    T_wall = -F[0]; L = F[1]; Mn = Mr + Mc
    F_ideal, M_e_ideal = ideal_gross_thrust(M_in, pa, g)
    F_gross = inlet_stream_thrust(M_in, pa, g) + T_wall
    return {"T_wall": float(T_wall), "L": float(L), "M_noseup": float(Mn),
            "F_gross": float(F_gross), "F_ideal": float(F_ideal), "M_e_ideal": float(M_e_ideal),
            "C_T": float(F_gross / F_ideal), "C_T_wall": float(T_wall / F_ideal),
            "C_L": float(L / F_ideal), "C_M": float(Mn / F_ideal),
            "T_ramp": float(-Fr[0]), "T_cowl": float(-Fc[0]), "L_ramp": float(Fr[1]), "L_cowl": float(Fc[1]),
            "pa_over_pin": float(pa)}
