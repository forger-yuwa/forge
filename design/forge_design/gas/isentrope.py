r"""等エントロピー膨張の状態 (音速単位) と、MOC 用の高精度な $\nu(M)$・$M(\nu)$・質量流束密度。

計画: plans/active/tooling-nozzle-bell-two-mode.md §4.4 (S1) ・§5.1 #5。③ベルの (M) モード (TP の Rao、
`geometry/rao_axisym.py`) が、制御面 DE の式 ($\Phi_2,\Phi_3$ が一定) とその解析微分、B の曲率の判定量 $\tilde g$ (記録) に使う。

**入力は既存のガスモデルの公開メソッドだけ** (`gas/semiperfect.py` は他ブランチも触るので本体は変えない):

- `GasSemiPerfect`: `cp_mass(T)`・`h_mass(T)`・`R`・`Tt`。その `nu`・`mach_of_nu` の表 (線形補間) は使わない。
- `GasCPG` と素の γ (float): 閉形式。`moc` は γ (float) を返すので、MOC の単位過程は従来の閉形式の経路を通る。

**音速単位** (`state(M)` の戻り): $T/T^*$、$V/a^*$、$a/a^*$、$\rho/\rho^*$、`p` $=p/(\rho^*a^{*2})$ (運動量流束 $\rho V^2$ と
同じ単位。$p/p^*$ は `p_over_pstar`)、$\mu=\arcsin(1/M)$、$k=d\ln a/d\ln V$。

TP の $k$: $a^2=\gamma(T)RT$、$dT=-V\,dV/c_p$ より $k=-\dfrac{V^2}{2c_p}\left(\dfrac1T+\dfrac{\gamma'}{\gamma}\right)$、
$\gamma'/\gamma=-Rc_p'/[c_p(c_p-R)]$。CPG では $k=-\tfrac{\gamma-1}2M^2$。$c_p'$ は `cp_mass` の中心差分
(NASA-9 の区切り 200・1000・6000 K をまたがない片側差分に切り替える)。

使う関係 (任意の EOS で成立): $d\ln V/dM=1/[M(1-k)]$、$d\nu/dM=\sqrt{M^2-1}/[M(1-k)]$、
$d\ln\rho/dM=-M/(1-k)$、$d\ln(\rho V)/dM=(1-M^2)/[M(1-k)]$。

**TP の表** (`n_tab` 節点、$u=\operatorname{arccosh}M$ に等間隔): 節点ごとに $T(M)$ を Newton で厳密に解き
($2(h_0-h(T))=M^2\gamma(T)RT$)、$\nu$ は $d\nu/du=\sinh^2u/[M(1-k)]$ を区間ごとの Gauss–Legendre で累積 ($u$ では $M=1$ の
特異性が消える)、$\ln(p/p^*)=\int_{T^*}^{T}c_p/(RT)\,dT$ は $T$ の区間ごとの Gauss–Legendre (区切り温度で分割) で累積する。
問い合わせは厳密な微分を使う 3 次 Hermite 補間 (誤差 $\sim h^4$)、$M(\nu)$ はその補間式の区間内 Newton の逆 (補間式と整合)。
"""
from __future__ import annotations

import numpy as np

# NASA-9 の区切り温度 (内蔵種の標準区間。cp の微分がここで不連続になる)
T_BREAKS = (200.0, 1000.0, 6000.0)
GL_N = 8                                    # 区間ごとの Gauss–Legendre の点数
_GLX, _GLW = np.polynomial.legendre.leggauss(GL_N)


def _gl_on(a, b):
    """区間 [a, b] (配列可) の Gauss–Legendre 節点 (…, GL_N) と重み (…, GL_N)。"""
    a = np.asarray(a, dtype=float)[..., None]
    b = np.asarray(b, dtype=float)[..., None]
    h = 0.5 * (b - a)
    return 0.5 * (a + b) + h * _GLX, h * _GLW


def _is_cpg(gas) -> bool:
    if isinstance(gas, (int, float, np.floating)):
        return True
    return getattr(gas, "kind", None) == "cpg" or (hasattr(gas, "gamma_ref") and not hasattr(gas, "h_mass"))


def _hermite(u, h, f0, f1, d0, d1):
    """区間 [0, h] の 3 次 Hermite (u は区間内の位置)。"""
    t = u / h
    t2, t3 = t * t, t * t * t
    return ((2 * t3 - 3 * t2 + 1) * f0 + (t3 - 2 * t2 + t) * h * d0
            + (-2 * t3 + 3 * t2) * f1 + (t3 - t2) * h * d1)


def _hermite_d(u, h, f0, f1, d0, d1):
    t = u / h
    t2 = t * t
    return ((6 * t2 - 6 * t) * f0 / h + (3 * t2 - 4 * t + 1) * d0
            + (-6 * t2 + 6 * t) * f1 / h + (3 * t2 - 2 * t) * d1)


class Isentrope:
    r"""等エントロピー膨張の accessor (音速単位) と MOC 用のガス (`nu`・`mach_of_nu`・`mass_flux_density`・`gamma_throat`)。

    `gas`: γ (float)・`GasCPG`・`GasSemiPerfect`。`moc`: MOC の単位過程に渡すガス (CPG は γ の float、TP は self)。
    `M_max`: 表の上限 (TP。既定 8、表の端の温度 T(M_max) が 60 K を割るなら下げる)。`n_tab`: 表の節点数 (TP)。"""

    def __init__(self, gas, M_max: float = 8.0, n_tab: int = 4001):
        self.source = gas
        self.cpg = _is_cpg(gas)
        if self.cpg:
            self.gamma_ref = float(gas) if isinstance(gas, (int, float, np.floating)) else float(gas.gamma_ref)
            g = self.gamma_ref
            if not (1.0 < g < 3.0):
                raise ValueError(f"Isentrope: γ = {g} が範囲外")
            self.kind = "cpg"
            self.moc = g
            self.gamma_star = g
            self.M_max = float(np.inf)
            # p0/(ρ* a*²) = (p0/p*)/γ
            self.p0 = ((0.5 * (g + 1.0)) ** (g / (g - 1.0))) / g
            return
        self.kind = "semiperfect_isentrope"
        self.moc = self
        self.R = float(gas.R)
        self.Tt = float(gas.Tt)
        self._cp = lambda T: np.asarray(gas.cp_mass(np.atleast_1d(np.asarray(T, dtype=float))), dtype=float)
        self._h = lambda T: np.asarray(gas.h_mass(np.atleast_1d(np.asarray(T, dtype=float))), dtype=float)
        self.h0 = float(self._h(self.Tt)[0])
        from scipy.optimize import brentq
        R = self.R

        def f_star(T):
            cp = float(self._cp(T)[0])
            return 2.0 * (self.h0 - float(self._h(T)[0])) - cp / (cp - R) * R * T
        self.T_star = float(brentq(f_star, 0.3 * self.Tt, self.Tt * (1.0 - 1e-12), xtol=1e-14 * self.Tt,
                                   rtol=4 * np.finfo(float).eps, maxiter=500))
        cps = float(self._cp(self.T_star)[0])
        self.gamma_star = cps / (cps - R)
        self.a_star = float(np.sqrt(self.gamma_star * R * self.T_star))
        # 表の上限: T(M_max) ≥ 60 K
        M_max = float(M_max)
        while M_max > 1.5 and float(self._T_newton(np.array([M_max]))[0]) < 60.0:
            M_max -= 0.25
        self.M_max = M_max
        self._build(int(n_tab))
        # p0/(ρ* a*²): ln(p0/p*) = ∫_{T*}^{Tt} cp/(R T) dT
        self.p0 = float(np.exp(self._lnp_T(self.T_star, self.Tt)[0])) / self.gamma_star

    # ------------------------------------------------------------------ TP の厳密な関数 (T 基準)
    def _cp_prime(self, T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        d = np.maximum(1e-5 * T, 1e-4)
        out = (self._cp(T + d) - self._cp(T - d)) / (2.0 * d)
        for tb in T_BREAKS:
            near = np.abs(T - tb) < d
            if np.any(near):
                Tn, dn = T[near], d[near]
                sgn = np.where(Tn >= tb, 1.0, -1.0)
                f0 = self._cp(Tn)
                f1 = self._cp(Tn + sgn * dn)
                f2 = self._cp(Tn + 2.0 * sgn * dn)
                out[near] = sgn * (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * dn)
        return out

    def _T_newton(self, M):
        """2(h0 − h(T)) = M² γ(T) R T を T について解く (超音速の枝、M ≥ 1)。"""
        M = np.atleast_1d(np.asarray(M, dtype=float))
        R = self.R
        g = self.gamma_star
        T = self.T_star * 0.5 * (g + 1.0) / (1.0 + 0.5 * (g - 1.0) * M * M)
        for _ in range(60):
            cp = self._cp(T)
            gam = cp / (cp - R)
            dgam = -R * self._cp_prime(T) / (cp - R) ** 2
            f = 2.0 * (self.h0 - self._h(T)) - M * M * gam * R * T
            df = -2.0 * cp - M * M * R * (gam + T * dgam)
            dT = f / df
            T = T - dT
            if np.all(np.abs(dT) <= 1e-14 * T):
                break
        return T

    def _lnp_T(self, T0, T1):
        """∫_{T0}^{T1} cp/(R T) dT (配列可、区切り温度で分割した Gauss–Legendre)。"""
        T0 = np.atleast_1d(np.asarray(T0, dtype=float))
        T1 = np.atleast_1d(np.asarray(T1, dtype=float))
        T0, T1 = np.broadcast_arrays(T0, T1)
        out = np.zeros(T0.shape)
        lo, hi = np.minimum(T0, T1), np.maximum(T0, T1)
        edges = [lo]
        for tb in T_BREAKS:
            edges.append(np.clip(np.full(lo.shape, tb), lo, hi))
        edges.append(hi)
        edges = np.sort(np.stack(edges, axis=0), axis=0)
        for a, b in zip(edges[:-1], edges[1:]):
            x, w = _gl_on(a, b)
            f = self._cp(x.ravel()).reshape(x.shape) / (self.R * x)
            out += np.sum(f * w, axis=-1)
        return np.where(T1 >= T0, out, -out)

    def _exact_T(self, T, M):
        """T・M から V/a*, a/a*, k, cp, γ (TP、厳密)。"""
        R = self.R
        cp = self._cp(T)
        gam = cp / (cp - R)
        V = np.sqrt(np.maximum(2.0 * (self.h0 - self._h(T)), 0.0))
        gp_g = -R * self._cp_prime(T) / (cp * (cp - R))
        k = -(V * V) / (2.0 * cp) * (1.0 / T + gp_g)
        return V / self.a_star, V / self.a_star / np.maximum(M, 1e-300), k, cp, gam

    # ------------------------------------------------------------------ TP の表
    def _build(self, n):
        u_max = float(np.arccosh(self.M_max))
        u = np.linspace(0.0, u_max, n)
        self._u = u
        self._hu = u[1] - u[0]
        M = np.cosh(u)
        T = self._T_newton(M)
        T[0] = self.T_star
        V, a, k, cp, gam = self._exact_T(T, M)
        # ν: dν/du = sinh²u / (M (1 − k))
        ug, wg = _gl_on(u[:-1], u[1:])
        Mg = np.cosh(ug)
        Tg = self._T_newton(Mg.ravel()).reshape(Mg.shape)
        _, _, kg, _, _ = self._exact_T(Tg.ravel(), Mg.ravel())
        kg = kg.reshape(Mg.shape)
        dnu_g = np.sinh(ug) ** 2 / (Mg * (1.0 - kg))
        nu = np.concatenate([[0.0], np.cumsum(np.sum(dnu_g * wg, axis=-1))])
        # ln(p/p*) と ln(ρ/ρ*) = ln(p/p*) − ln(T/T*)
        lnp = np.concatenate([[0.0], np.cumsum(self._lnp_T(T[:-1], T[1:]))])
        lnrho = lnp - np.log(T / self.T_star)
        sh = np.sinh(u)
        dlnV_dM = 1.0 / (M * (1.0 - k))
        self._tab = {
            "M": M, "T": T, "V": V, "k": k, "cp": cp, "gam": gam, "nu": nu, "lnp": lnp, "lnrho": lnrho,
            "dnu_du": sh * sh / (M * (1.0 - k)),
            # ln(ρV/ρ*a*) と ln(p/p*) の u 微分
            "lnG": lnrho + np.log(V), "dlnG_du": (1.0 - M * M) * dlnV_dM * sh,
            "dlnp_du": -(V * self.a_star) ** 2 * dlnV_dM * sh / (self.R * T),
            "dT_du": -(V * self.a_star) ** 2 * dlnV_dM * sh / cp,
        }
        self._M = M                                   # moc_kernel.dnu_dM・moc_forward.nu_max_of が読む
        self._nu = nu

    def _cell(self, u):
        u = np.clip(u, 0.0, self._u[-1])
        i = np.minimum((u / self._hu).astype(int), len(self._u) - 2)
        return i, u - self._u[i]

    def _interp(self, key, dkey, u):
        i, du = self._cell(u)
        t = self._tab
        return _hermite(du, self._hu, t[key][i], t[key][i + 1], t[dkey][i], t[dkey][i + 1])

    # ------------------------------------------------------------------ MOC のガス API
    def nu(self, M):
        if self.cpg:
            from ..geometry.moc_kernel import pm_nu
            return pm_nu(M, self.gamma_ref)
        M = np.asarray(M, dtype=float)
        u = np.arccosh(np.maximum(M, 1.0))
        return self._interp("nu", "dnu_du", u)

    def mach_of_nu(self, nu):
        if self.cpg:
            from ..geometry.moc_kernel import pm_mach_vec
            return pm_mach_vec(np.atleast_1d(nu), self.gamma_ref)
        nu = np.asarray(nu, dtype=float)
        shp = nu.shape
        v = np.clip(np.ravel(nu), 0.0, self._nu[-1])
        t = self._tab
        i = np.clip(np.searchsorted(self._nu, v, side="right") - 1, 0, len(self._u) - 2)
        f0, f1, d0, d1 = t["nu"][i], t["nu"][i + 1], t["dnu_du"][i], t["dnu_du"][i + 1]
        h = self._hu
        du = h * (v - f0) / np.maximum(f1 - f0, 1e-300)       # 線形の初期値
        for _ in range(40):
            f = _hermite(du, h, f0, f1, d0, d1) - v
            df = _hermite_d(du, h, f0, f1, d0, d1)
            step = f / np.where(df > 0.0, df, 1e-300)
            du_new = np.clip(du - step, 0.0, h)
            if np.all(np.abs(du_new - du) <= 1e-15 * (1.0 + self._u[i])):
                du = du_new
                break
            du = du_new
        return np.cosh(self._u[i] + du).reshape(shp)

    def mass_flux_density(self, M):
        """ρV/(ρ* a*) (M ≥ 1)。"""
        if self.cpg:
            g = self.gamma_ref
            M = np.asarray(M, dtype=float)
            return M * ((2.0 + (g - 1.0) * M * M) / (g + 1.0)) ** (-0.5 * (g + 1.0) / (g - 1.0))
        u = np.arccosh(np.maximum(np.asarray(M, dtype=float), 1.0))
        return np.exp(self._interp("lnG", "dlnG_du", u))

    def gamma_throat(self, Tt=None):
        return self.gamma_star

    def area_ratio(self, M):
        return 1.0 / self.mass_flux_density(M)

    # ------------------------------------------------------------------ 状態 (音速単位)
    def state(self, M) -> dict:
        r"""等エントロピー膨張の状態 (M ≥ 1、配列可)。戻り dict: `T` (T/T*)・`V` (V/a*)・`a` (a/a*)・`rho` (ρ/ρ*)・
        `p` (p/(ρ* a*²))・`p_over_pstar`・`mu`・`k` (d ln a/d ln V)・`dlnV_dM` (= 1/[M(1−k)])・`gamma`。
        TP は T を Newton で厳密に解き、V・a・k は厳密、ρ と p は表の Hermite 補間 (誤差 ~h⁴)。"""
        M = np.atleast_1d(np.asarray(M, dtype=float))
        if self.cpg:
            g = self.gamma_ref
            Tr = 0.5 * (g + 1.0) / (1.0 + 0.5 * (g - 1.0) * M * M)
            a = np.sqrt(Tr)
            V = M * a
            rho = Tr ** (1.0 / (g - 1.0))
            pps = Tr ** (g / (g - 1.0))
            k = -0.5 * (g - 1.0) * M * M
            gam = np.full(M.shape, g)
        else:
            if np.any(M < 1.0 - 1e-12) or np.any(M > self.M_max * (1.0 + 1e-12)):
                raise ValueError(f"Isentrope.state: M は [1, {self.M_max}] (受け取った範囲 {M.min()}〜{M.max()})")
            T = self._T_newton(M)
            V, a, k, cp, gam = self._exact_T(T, M)
            Tr = T / self.T_star
            u = np.arccosh(np.maximum(M, 1.0))
            lnp = self._interp("lnp", "dlnp_du", u)
            pps = np.exp(lnp)
            rho = pps / Tr
        return {"M": M, "T": Tr, "V": V, "a": a, "rho": rho, "p": pps / self.gamma_star, "p_over_pstar": pps,
                "mu": np.arcsin(1.0 / np.maximum(M, 1.0 + 1e-15)), "k": k, "dlnV_dM": 1.0 / (M * (1.0 - k)),
                "gamma": gam}

    def mach_of_area(self, eps: float) -> float:
        """超音速の枝で A/A* = eps の M (1D)。"""
        from scipy.optimize import brentq
        hi = 2.0
        while float(self.mass_flux_density(hi)) > 1.0 / eps:
            hi *= 1.5
            if hi > self.M_max:
                raise ValueError(f"mach_of_area: ε {eps} は表の上限 M {self.M_max} を超える")
        return float(brentq(lambda m: float(self.mass_flux_density(m)) - 1.0 / eps, 1.0 + 1e-12, hi, xtol=1e-14))

    def ideal_1d(self, eps: float, p_a: float = 0.0) -> dict:
        """1D 理想値 (スロートで一様音速、出口は面積比 eps、音速単位・A* = 1): ṁ = 1、F = u_e + (p_e − p_a) eps。"""
        Me = self.mach_of_area(eps)
        s = self.state(Me)
        F = float(s["V"][0]) + (float(s["p"][0]) - p_a) * eps
        return {"M_e": Me, "u_e": float(s["V"][0]), "p_e": float(s["p"][0]), "F": F, "mdot": 1.0,
                "CF": F / self.p0, "p0": self.p0}


def isentrope(gas, M, **kw) -> dict:
    """関数形の accessor: `Isentrope(gas, **kw).state(M)` (毎回表を作るので、繰り返すなら `Isentrope` を使い回す)。"""
    return Isentrope(gas, **kw).state(M)


def constant_cp_gas(gamma: float, Tt: float = 3000.0, MW: float = 0.0289645, name: str = "CPX"):
    """**定比熱の擬似種**の `GasSemiPerfect` (T1d の回帰用): NASA-9 の a2 = γ/(γ−1) だけを持つ 1 種 (温度区間の制限なし)。
    TP の経路で CPG の閉形式と一致するかを確かめるためのもの。"""
    from .composition import ResolvedSpeciesDB, SpeciesEntry
    from .semiperfect import GasSemiPerfect
    a2 = gamma / (gamma - 1.0)
    coef = [0.0, 0.0, a2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    e = SpeciesEntry(name, float(MW), list(coef), list(coef), Tlo=1.0, Tmid=1000.0, Thi=1.0e5, T_eval_max=None,
                     source="constant_cp_gas (T1d)")
    return GasSemiPerfect({name: 1.0}, Tt=float(Tt), db=ResolvedSpeciesDB({name: e}))
