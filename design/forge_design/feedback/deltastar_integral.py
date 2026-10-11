r"""境界層積分法による排除厚さの**初期推定器** (初回 NS 用の壁だけに使う)。

計画: plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.1。

方法 = Sivells CONTUR (AEDC-TR-78-63 §5, Eq. 61–90) の軸対称 von Kármán 運動量積分:

$$\frac{d\theta}{dx} + \theta\left[\frac{2-M^2+H}{M(1+\frac{\gamma-1}{2}M^2)}\frac{dM}{dx}
  + \frac{1}{r_w}\frac{dr_w}{dx}\right] = \frac{C_f}{2}\sec\phi_w \qquad (61)$$

- プロファイル: $q/q_e=(z/\delta)^{1/N}$ (67)、$\rho/\rho_e=T_e/T$ (68)、
  $T = T_w + a(T_{aw}-T_w)\,u + [T_e - a(T_{aw}-T_w) - T_w]\,u^2$, $u=q/q_e$ (69; $a$=1 Crocco / 0 放物)。
- $\theta,\delta^*$ は横曲率重み $(1 - z\cos\phi_w/r_w)$ 込み (62, 63)。$H=\delta^*/\theta$。
- 摩擦: $C_f = C_{f_i}/F_c$, $R_{\theta_i}=F_{R_\delta}R_{\theta_c}$ (70–73)、
  $F_c=[\int_0^1(\rho/\rho_e)^{1/2}du]^{-2}$, $F_{R_\delta}=\mu_e/\mu_w$,
  $C_{f_i}=0.0773/[(\log R_{\theta_i}+4.561)(\log R_{\theta_i}-0.546)]$ (75)、$R_{\theta_c}$ は平板形 $\theta_c$ (74)。
- $N(R_\delta)$ は CONTUR Fig. 10 のテーブル (Eq. 77–82 の連鎖の結果)。
- 各ステーションで状態量 $\theta$ から $\delta$ を根探索し $N,\delta^*,H,\theta_c$ を得る。

Sasman–Cresci (ユーザ計画の指定) は閉包の原論文が手元に無いため未実装。`closure` で差し替え可能な構造にし、
既定は `"contur"`。**CFD との一致率は合否条件にしない** (初期値生成専用)。

版 (`closure_version`、plan tooling-nozzle-isothermal-wall-chain §4.7): `contur_v1` (既定、従来どおり) / `contur_v2` (NS と同じ気体で書いた整合版:
温度分布をエンタルピーで $h(v) = h_w + a(h_{aw}-h_w)v + [h_e - a(h_{aw}-h_w) - h_w]v^2$、$T = h^{-1}(h)$; $h_{aw} = h_e + r(h_0 - h_e)$;
粘性は呼び出し側が渡す `mu_fn` (NS と同じ輸送物性); 運動量式の加速の項は $(2 + H - M^2)\,d\ln u_e/dx$)。

熱境界条件: `thermal_bc = {"mode": "adiabatic"}` または
`{"mode": "prescribed_temperature", "Tw": 300.0}` / `{"mode": "prescribed_temperature", "Tw_table": [[x_rt, Tw], ...]}`。
断熱壁は $T_w = T_{aw} = T_e(1 + r_f\frac{\gamma-1}{2}M^2)$, $r_f = Pr^{1/3}$。
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .deltastar import _sutherland

# CONTUR Fig. 10: N vs log10(R_δ)
_N_TAB_LOGRE = np.array([3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
_N_TAB_N = np.array([4.6, 4.9, 5.4, 6.3, 7.3, 8.7, 10.7, 12.5])
_GL_U, _GL_W = np.polynomial.legendre.leggauss(16)
_GL_U = 0.5 * (_GL_U + 1.0); _GL_W = 0.5 * _GL_W          # [0,1] 上の 16 点 Gauss


def N_of_Redelta(Re_delta: float) -> float:
    return float(np.interp(np.log10(max(Re_delta, 10.0)), _N_TAB_LOGRE, _N_TAB_N))


def _profile_integrals(delta: float, N: float, Tw: float, Taw: float, Te: float,
                       rw: float, cos_phi: float, a: float, thermo=None):
    """(θ, δ*, θ_c, F_c) を Gauss 16 点で。z = δ u^N, dz = N δ u^{N-1} du。
    thermo (contur_v2 の EdgeConditions): 温度分布をエンタルピーで書く (h(v) の 2 次式 → T = h⁻¹)。None は従来の温度形。"""
    u = _GL_U; w = _GL_W
    if thermo is None:
        T = Tw + a * (Taw - Tw) * u + (Te - a * (Taw - Tw) - Tw) * u ** 2
    else:
        hw, haw, he = thermo.H(Tw), thermo.H(Taw), thermo.H(Te)
        T = thermo.Hinv(hw + a * (haw - hw) * u + (he - a * (haw - hw) - hw) * u ** 2)
    rho_rel = Te / np.maximum(T, 1e-30)
    z = delta * u ** N
    dz = N * delta * u ** (N - 1.0)
    curv = 1.0 - z * cos_phi / rw
    theta = float(np.sum(w * curv * rho_rel * u * (1.0 - u) * dz))
    dstar = float(np.sum(w * curv * (1.0 - rho_rel * u) * dz))
    theta_c = float(np.sum(w * rho_rel * u * (1.0 - u) * dz))
    Fc = float(np.sum(w * np.sqrt(rho_rel))) ** -2
    return theta, dstar, theta_c, Fc


class EdgeConditions:
    """設計壁に沿う縁条件 M_e(x), T_e, p_e, ρ_e, u_e, μ_e, γ_e と幾何 r_w, dr_w/dx (x, r は r_t 単位)。"""

    def __init__(self, design_wall, wall_tbl, gas, cp: float, Pt: float, Tt: float, rt_m: float,
                 Pr: float = 0.72, version: str = "contur_v1", mu_fn=None, T_range=None, n_grid: int = 4000):
        from ..evaluate.ic import invert_area_ratio
        self.dw = design_wall; self.rt_m = float(rt_m); self.Pt = float(Pt); self.Tt = float(Tt)
        self.Pr = Pr
        if version not in ("contur_v1", "contur_v2"):
            raise ValueError(f"closure_version {version!r} は contur_v1 | contur_v2")
        if version == "contur_v2" and mu_fn is None:
            raise ValueError("contur_v2 は mu_fn (NS と同じ輸送物性の μ(T)) が要る")
        self.version = version
        self._mu_fn = mu_fn
        self.gas_obj = gas if hasattr(gas, "area_ratio") else None
        self.gamma_cpg = None if self.gas_obj is not None else float(gas)
        self.cp = float(cp)
        wall_tbl = np.asarray(wall_tbl, dtype=float)
        self.x_in = float(design_wall.x_in); self.x_F = float(wall_tbl[-1, 0])
        xg = np.linspace(self.x_in + 1e-6, self.x_F, int(n_grid))     # 縁の状態の格子 (既定 4000; 試験で細分の収束を見る)
        rg = design_wall.r(xg)
        M = np.empty_like(xg)
        m_up = xg <= 0.0
        AR_up = np.maximum(rg[m_up], 1.0 + 1e-12) ** 2
        if self.gas_obj is not None:
            from ..evaluate.ic import _invert_area_ratio_gas
            M[m_up] = _invert_area_ratio_gas(AR_up, np.zeros(int(m_up.sum()), bool), self.gas_obj)
        else:
            M[m_up] = invert_area_ratio(AR_up, np.zeros(int(m_up.sum()), bool), self.gamma_cpg)
        M[~m_up] = np.interp(xg[~m_up], wall_tbl[:, 0], wall_tbl[:, 3])
        M = np.maximum(M, 1e-3)
        self._xg = xg; self._M = M
        # dM/dx: 上流 (亜音速枝) は 1D 等エントロピーの解析式
        #   d ln A / d ln M = -(1 - M²)/(1 + (γ-1)/2 M²)  →  dM/dx = -M (1+(γ-1)/2 M²)/(1-M²) · 2 r_w'/r_w
        # (数値微分は M~0.03 の直管で 1/M 倍に増幅されて θ を波打たせる)。下流は MOC 壁 M の勾配。
        dM_num = np.gradient(M, xg)
        self._dMdx = dM_num
        self._m_up = m_up; self._rg = rg; self._rpg = design_wall.r(xg, 1)
        if self.gas_obj is not None:
            g = self.gas_obj
            self.R = float(g.R)
            # T_e(M): 超音速はテーブル、亜音速は上流枝
            Te = np.where(M >= 1.0, g.T_of_M(np.maximum(M, 1.0)), np.interp(M, g._Mu, g._Tu))
            # p_e/Pt = exp(∫_{Tt}^{Te} cp/(R T) dT) を T の細分で数値積分
            Tt_ = self.Tt
            Tgrid = np.linspace(Tt_, float(np.min(Te)) * 0.999, 3000)
            cpg = g.cp_mass(Tgrid)
            s_R = np.concatenate([[0.0], np.cumsum(0.5 * (cpg[1:] / Tgrid[1:] + cpg[:-1] / Tgrid[:-1]) * np.diff(Tgrid) / self.R)])
            lnp = np.interp(Te, Tgrid[::-1], s_R[::-1])
            pe = self.Pt * np.exp(lnp)
            gam = np.asarray(g.gamma(Te), dtype=float)
        else:
            gam0 = self.gamma_cpg
            self.R = self.cp * (gam0 - 1.0) / gam0
            Te = self.Tt / (1.0 + 0.5 * (gam0 - 1.0) * M ** 2)
            pe = self.Pt * (Te / self.Tt) ** (gam0 / (gam0 - 1.0))
            gam = np.full_like(M, gam0)
        self._Te = Te; self._pe = pe; self._gam = gam
        self._rho = pe / (self.R * Te)
        self._ue = M * np.sqrt(gam * self.R * Te)
        up = self._m_up & (M < 0.98)
        dM_an = -M * (1.0 + 0.5 * (gam - 1.0) * M ** 2) / np.maximum(1.0 - M ** 2, 1e-6) \
            * 2.0 * self._rpg / np.maximum(self._rg, 1e-12)
        self._dMdx = np.where(up, dM_an, self._dMdx)
        if self.version == "contur_v2":
            self._init_v2(up, T_range)

    # --- contur_v2 (plan tooling-nozzle-isothermal-wall-chain §4.7) ---
    def _init_v2(self, up, T_range):
        """エンタルピーの表 h(T) (気体の c_p(T) の積分、CPG は c_p T)、粘性の表 μ(T)、d ln u_e/dx。
        d ln u_e/dx: 上流の亜音速枝 (x ≤ 0, M < 0.98) は面積から −(2 r_w'/r_w)/(1 − M²) (等エントロピー・組成一定で γ を使わない)、
        それ以外は u_e(x) の数値微分。M = 1 の近くは数値微分の側 (0/0 を避ける; 上流の解析式と 0.98 で切り替える)。"""
        Te = self._Te
        lo = float(T_range[0]) if T_range else min(150.0, 0.5 * float(np.min(Te)))
        hi = float(T_range[1]) if T_range else 1.25 * self.Tt
        self._Tg = np.linspace(lo, hi, 20001)
        if self.gas_obj is not None:
            cp = np.asarray(self.gas_obj.cp_mass(self._Tg), dtype=float)
            self._hg = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(self._Tg))])
        else:
            self._hg = self.cp * (self._Tg - self._Tg[0])
        self._mug = np.asarray(self._mu_fn(self._Tg), dtype=float)
        if not (np.all(np.isfinite(self._mug)) and np.all(self._mug > 0)):
            raise ValueError("contur_v2: μ(T) の表が不正")
        dlnue_num = np.gradient(np.log(self._ue), self._xg)
        dlnue_an = -2.0 * self._rpg / np.maximum(self._rg, 1e-12) / np.maximum(1.0 - self._M ** 2, 1e-6)
        self._dlnue = np.where(up, dlnue_an, dlnue_num)
        self._mue = self.mu(self._Te)

    def _check_T(self, T):
        T = np.asarray(T, dtype=float)
        if np.any(T < self._Tg[0]) or np.any(T > self._Tg[-1]):
            raise ValueError(f"contur_v2: 温度 {float(np.min(T)):.1f}〜{float(np.max(T)):.1f} K が表の範囲 "
                             f"{self._Tg[0]:.1f}〜{self._Tg[-1]:.1f} K の外 (端値で外挿しない)")
        return T

    def H(self, T):
        return np.interp(self._check_T(T), self._Tg, self._hg)

    def Hinv(self, h):
        h = np.asarray(h, dtype=float)
        if np.any(h < self._hg[0]) or np.any(h > self._hg[-1]):
            raise ValueError("contur_v2: エンタルピーが表の範囲の外 (端値で外挿しない)")
        return np.interp(h, self._hg, self._Tg)

    def mu(self, T):
        if self.version != "contur_v2":
            return _sutherland(T)
        return np.interp(self._check_T(T), self._Tg, self._mug)

    def at(self, x: float) -> dict:
        xg = self._xg
        f = lambda arr: float(np.interp(x, xg, arr))
        M = f(self._M); Te = f(self._Te); gam = f(self._gam)
        rw = float(self.dw.r(np.array([x]))[0]); drw = float(self.dw.r(np.array([x]), 1)[0])
        rf = self.Pr ** (1.0 / 3.0)
        if self.version == "contur_v2":
            he = float(self.H(Te))
            Taw = float(self.Hinv(he + rf * (float(self.H(self.Tt)) - he)))
            return dict(M=M, dMdx=f(self._dMdx), Te=Te, pe=f(self._pe), rho_e=f(self._rho), ue=f(self._ue),
                        gam=gam, mu_e=f(self._mue), rw=rw, drwdx=drw,
                        cos_phi=float(1.0 / np.sqrt(1.0 + drw ** 2)), Taw=Taw, dlnue_dx=f(self._dlnue), thermo=self)
        Taw = Te * (1.0 + rf * 0.5 * (gam - 1.0) * M ** 2)
        return dict(M=M, dMdx=f(self._dMdx), Te=Te, pe=f(self._pe), rho_e=f(self._rho), ue=f(self._ue),
                    gam=gam, mu_e=float(_sutherland(Te)), rw=rw, drwdx=drw,
                    cos_phi=float(1.0 / np.sqrt(1.0 + drw ** 2)), Taw=Taw)


def _wall_temperature(thermal_bc: dict | None, x: float, Taw: float) -> float:
    if not thermal_bc or str(thermal_bc.get("mode", "adiabatic")) == "adiabatic":
        return Taw
    if str(thermal_bc["mode"]) == "prescribed_temperature":
        if "Tw_table" in thermal_bc and thermal_bc["Tw_table"] is not None:
            tb = np.asarray(thermal_bc["Tw_table"], dtype=float)
            return float(np.interp(x, tb[:, 0], tb[:, 1]))
        return float(thermal_bc["Tw"])
    raise ValueError(f"thermal_bc.mode = {thermal_bc.get('mode')} は未対応")


def closure_contur(theta_m: float, e: dict, Tw: float, a: float = 1.0, cf_scale: float = 1.0, n_scale: float = 1.0) -> dict:
    """θ (物理 [m]) と縁条件から δ, N, δ*, H, θ_c, C_f を決める (CONTUR 閉包)。長さは [m]。"""
    rw_m = e["rw_m"]; Te = e["Te"]; Taw = e["Taw"]; cos_phi = e["cos_phi"]
    rho_e, ue, mu_e = e["rho_e"], e["ue"], e["mu_e"]
    th_obj = e.get("thermo")                 # contur_v2 のときだけ (温度形の差し替えを壊さないよう None なら渡さない)
    prof = (lambda *args: _profile_integrals(*args)) if th_obj is None else (lambda *args: _profile_integrals(*args, thermo=th_obj))

    def theta_of_delta(delta):
        N = n_scale * N_of_Redelta(rho_e * ue * delta / mu_e)
        return prof(delta, N, Tw, Taw, Te, rw_m, cos_phi, a)[0] - theta_m

    # 探索区間は壁半径と 1 m の小さい方で切る (平面入口 flat_plate_integral は rw_m=1e30 を渡す)
    L = min(rw_m, 1.0)
    lo, hi = 1e-12 * L, 0.9 * L
    if theta_of_delta(hi) < 0.0:            # θ が壁半径に対して大きすぎる (非物理) → 上限で打ち切り
        delta = hi
    elif theta_of_delta(lo) >= 0.0:         # θ が実質ゼロ → 下限
        delta = lo
    else:
        delta = brentq(theta_of_delta, lo, hi, xtol=1e-12 * L, maxiter=200)
    N = n_scale * N_of_Redelta(rho_e * ue * delta / mu_e)
    th, ds, th_c, Fc = prof(delta, N, Tw, Taw, Te, rw_m, cos_phi, a)
    Re_theta_c = rho_e * ue * th_c / mu_e
    F_Rd = mu_e / float(_sutherland(Tw) if th_obj is None else th_obj.mu(Tw))
    Re_theta_i = max(F_Rd * Re_theta_c, 300.0)                 # Eq. 75 の有効域下限で床
    lg = np.log10(Re_theta_i)
    Cfi = 0.0773 / ((lg + 4.561) * (lg - 0.546))
    Cf = cf_scale * Cfi / Fc
    return dict(delta=delta, N=N, dstar=ds, H=ds / max(th, 1e-30), theta_c=th_c, Cf=Cf, Cfi=Cfi,
                Fc=Fc, F_Rdelta=F_Rd, Re_theta_c=Re_theta_c, Re_delta=rho_e * ue * delta / mu_e)


def integral_bl(design_wall, wall_tbl, gas, cp: float, Pt: float, Tt: float, rt_m: float,
                thermal_bc: dict | None = None, theta0_m: float | None = None,
                x_virtual_m: float | None = None, a_crocco: float = 1.0, closure: str = "contur",
                x_out=None, rtol: float = 1e-6, cf_scale: float = 1.0, n_scale: float = 1.0,
                closure_version: str = "contur_v1", mu_fn=None, edge_n_grid: int = 4000) -> dict:
    r"""入口 $x_{in}$ から $x_F$ まで Eq. (61) を前進積分し、δ*_n(x) と半径方向補正 δ_r = δ*_n/cos φ_w を返す。

    cf_scale / n_scale: CFD 較正用の C_f 倍率と N 倍率 (既定 1 = CONTUR そのまま; plan verification-m6-axis-wave-mesh-su2 §5.1 #8b)。
    theta0_m: 入口の運動量厚さ [m]。None なら仮想発達長 x_virtual_m (None = 入口直管長) の乱流平板
    θ0 = 0.036 x_v Re_{x_v}^{-0.2}。戻り値の x, δ 系は r_t 単位、θ 等の物理量は [m] も併記。
    rtol: RK45 の相対許容差。戻り値の `solve_ivp` に実際に渡した値と評価回数を残す (`settings` には入れない —
    prepare_ns の `delta_r_initial.json` を変えないため; plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c)。"""
    if closure != "contur":
        raise NotImplementedError(f"closure={closure!r} は未実装 (contur のみ)")
    ec = EdgeConditions(design_wall, wall_tbl, gas, cp, Pt, Tt, rt_m, version=closure_version, mu_fn=mu_fn, n_grid=edge_n_grid)
    x0, x1 = ec.x_in + 1e-6, ec.x_F
    e0 = ec.at(x0)
    if theta0_m is None:
        if x_virtual_m is None:
            x_virtual_m = float(getattr(design_wall, "L_pipe", 0.5)) * rt_m
        Re_x = e0["rho_e"] * e0["ue"] * x_virtual_m / e0["mu_e"]
        theta0_m = 0.036 * x_virtual_m * Re_x ** -0.2
        theta0_source = f"flat-plate virtual length {x_virtual_m:.4g} m (Re_x {Re_x:.3g})"
    else:
        theta0_source = "user"

    def station(x, theta_rt):
        e = ec.at(x); e["rw_m"] = e["rw"] * rt_m
        Tw = _wall_temperature(thermal_bc, x, e["Taw"])
        c = closure_contur(theta_rt * rt_m, e, Tw, a=a_crocco, cf_scale=cf_scale, n_scale=n_scale)
        return e, Tw, c

    def rhs(x, y):
        theta_rt = max(float(y[0]), 1e-12)
        e, Tw, c = station(x, theta_rt)
        M = e["M"]; gam = e["gam"]
        if closure_version == "contur_v2":
            # 等エントロピー・組成一定: d ln ρ_e = −M² d ln u_e (γ を使わない、plan §4.7-4)
            term = (2.0 + c["H"] - M ** 2) * e["dlnue_dx"] + e["drwdx"] / e["rw"]
        else:
            term = (2.0 - M ** 2 + c["H"]) / (M * (1.0 + 0.5 * (gam - 1.0) * M ** 2)) * e["dMdx"] \
                + e["drwdx"] / e["rw"]
        return [0.5 * c["Cf"] / e["cos_phi"] - theta_rt * term]

    ivp = dict(method="RK45", rtol=rtol, atol=1e-14, max_step=(x1 - x0) / 400.0)
    sol = solve_ivp(rhs, (x0, x1), [theta0_m / rt_m], dense_output=True, **ivp)
    if not sol.success:
        raise RuntimeError(f"integral_bl: 積分失敗 ({sol.message})")
    xs = np.asarray(x_out, dtype=float) if x_out is not None else np.linspace(x0, x1, 1500)
    th = np.maximum(sol.sol(xs)[0], 1e-12)
    rows = {k: [] for k in ("theta", "dstar_n", "H", "N", "delta", "Cf", "M", "Te", "Tw", "Taw", "cos_phi", "Re_theta_c")}
    for x, t in zip(xs, th):
        e, Tw, c = station(float(x), float(t))
        rows["theta"].append(t); rows["dstar_n"].append(c["dstar"] / rt_m); rows["H"].append(c["H"])
        rows["N"].append(c["N"]); rows["delta"].append(c["delta"] / rt_m); rows["Cf"].append(c["Cf"])
        rows["M"].append(e["M"]); rows["Te"].append(e["Te"]); rows["Tw"].append(Tw); rows["Taw"].append(e["Taw"])
        rows["cos_phi"].append(e["cos_phi"]); rows["Re_theta_c"].append(c["Re_theta_c"])
    out = {k: np.asarray(v) for k, v in rows.items()}
    out["x"] = xs
    out["delta_r"] = out["dstar_n"] / out["cos_phi"]
    settings_v = {} if closure_version == "contur_v1" else {"closure_version": closure_version}
    out["settings"] = dict(model="contur_momentum_integral", closure=closure, a_crocco=a_crocco,
                           cf_scale=cf_scale, n_scale=n_scale, **settings_v,
                           thermal_bc=(thermal_bc or {"mode": "adiabatic"}), theta0_m=float(theta0_m),
                           theta0_source=theta0_source, x_virtual_m=x_virtual_m, rt_m=rt_m, Pt=Pt, Tt=Tt,
                           gas=("semiperfect" if ec.gas_obj is not None else f"cpg gamma={ec.gamma_cpg}"))
    # 積分の実効の設定 (solve_ivp に渡した値) と手間。settings とは分ける (上の docstring の rtol の項)
    out["solve_ivp"] = dict(ivp, nfev=int(sol.nfev), n_steps=int(len(sol.t) - 1))
    return out


def delta_r_function(res: dict):
    """integral_bl の結果から PhysicalNozzleWall(delta_r_x=...) に渡す callable (r_t 単位)。"""
    x, d = res["x"], res["delta_r"]
    return lambda xq, _x=x, _d=d: np.interp(xq, _x, _d)


def flat_plate_integral(x_m, Me: float, Te: float, Pe: float, gam: float, cp: float, theta0_m: float,
                        Tw=None, a_crocco: float = 1.0, r_rec: float | None = None) -> dict:
    r"""平板 (一定外縁条件・平面) 用の入口: $d\theta/dx = C_f/2$ を CONTUR 閉包で積分する
    (plan tooling-nozzle-isothermal-wall-chain §4.3-3: 巨大 $r_w$ の流用でなく明示の平面入口)。

    x_m: 積分する x [m] (昇順、x_m[0] で θ = theta0_m)。Tw: None = 断熱 / float / callable(x)。
    r_rec: 回復係数 (None = Pr^{1/3} = 0.896)。戻り: dict(x, theta, dstar, H, delta, N, Cf, Tw, Taw) 全て [m]/無次元。
    """
    R = cp * (gam - 1.0) / gam
    r_rec = float(0.72 ** (1.0 / 3.0)) if r_rec is None else float(r_rec)
    Taw = Te * (1.0 + r_rec * 0.5 * (gam - 1.0) * Me ** 2)
    rho_e = Pe / (R * Te); ue = Me * np.sqrt(gam * R * Te); mu_e = float(_sutherland(Te))
    e = dict(rw_m=1e30, Te=Te, Taw=Taw, cos_phi=1.0, rho_e=rho_e, ue=ue, mu_e=mu_e)

    def Tw_at(x):
        if Tw is None: return Taw
        return float(Tw(x)) if callable(Tw) else float(Tw)

    def rhs(x, y):
        c = closure_contur(max(float(y[0]), 1e-12), e, Tw_at(x), a=a_crocco)
        return [0.5 * c["Cf"]]

    xs = np.asarray(x_m, dtype=float)
    sol = solve_ivp(rhs, (xs[0], xs[-1]), [theta0_m], method="RK45", rtol=1e-7, atol=1e-14,
                    dense_output=True, max_step=(xs[-1] - xs[0]) / 200.0)
    if not sol.success:
        raise RuntimeError(f"flat_plate_integral: 積分失敗 ({sol.message})")
    th = np.maximum(sol.sol(xs)[0], 1e-12)
    rows = {k: [] for k in ("dstar", "H", "delta", "N", "Cf", "Tw")}
    for x, t in zip(xs, th):
        c = closure_contur(float(t), e, Tw_at(float(x)), a=a_crocco)
        rows["dstar"].append(c["dstar"]); rows["H"].append(c["H"]); rows["delta"].append(c["delta"])
        rows["N"].append(c["N"]); rows["Cf"].append(c["Cf"]); rows["Tw"].append(Tw_at(float(x)))
    out = {k: np.asarray(v) for k, v in rows.items()}
    out.update(x=xs, theta=th, Taw=Taw, rho_e=rho_e, ue=ue, mu_e=mu_e)
    return out
