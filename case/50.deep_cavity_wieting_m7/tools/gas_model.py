#!/usr/bin/env python3
"""case/50 の気体モデル — メタン–空気燃焼生成物 (NASA TN D-5908 の試験気体)。

W70 (TN D-5908) の物性は ref.18 = Leyhe & Howell, NASA TN D-914 (1962)
"Calculation Procedure for Thermodynamic, Transport, and Flow Properties of the
Combustion Products of a Hydrocarbon Fuel Mixture Burned in Air" の手順による。
ここでは同じ**考え方** (燃焼生成物の凍結組成 + 混合則) を、リポジトリの NASA-9 熱力学
(design/forge_design/gas) と Chapman–Enskog + Wilke/Mason–Saxena の輸送で再構成する。

- 熱力学: NASA-9 (ResolvedSpeciesDB) — cp(T), h(T), R, γ(T)
- 輸送:   Chapman–Enskog (Neufeld Ω) + Wilke (μ) / Mason–Saxena (λ)、λ は modified Eucken
- 組成:   CH4 + 空気の断熱定圧燃焼 (完全燃焼・凍結)。Tt から燃空比を逆算する。
          1650–1900 K では解離は無視できる (平衡計算との差は CEA で後日確認する)。

**これは W70 の原手順そのものではない**ので、Pr・cp を報告値 (Pr≈0.75) と突き合わせて使う。
"""
import sys, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.composition import ResolvedSpeciesDB          # noqa: E402

R_UNIV = 8.31446261815324        # J/(mol·K)。DB の MW は kg/mol (N2 = 0.0280134)
LHV_CH4 = 50.03e6                # CH4 の低位発熱量 [J/kg] (298 K, 生成水は気体)
                                 # NASA-9 の h_mass は**顕熱** (基準 ~298 K) なので、反応熱は別に与える

# Lennard-Jones パラメータ (σ [Å], ε/k [K]) — Poling et al., The Properties of Gases and Liquids
LJ = {
    "N2":  (3.798, 71.4),
    "O2":  (3.467, 106.7),
    "CO2": (3.941, 195.2),
    "H2O": (2.641, 809.1),       # 極性。Stockmayer 補正はしない (低温域で数 % の過小評価)
    "AR":  (3.542, 93.3),
}
# 空気 (モル分率)
AIR_X = {"N2": 0.78084, "O2": 0.20946, "AR": 0.00934, "CO2": 0.00036}


def _omega_mu(Tstar):
    """Neufeld ら (1972) の衝突積分 Ω_μ。"""
    A, B, C, D, E, F = 1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787
    return A * Tstar ** (-B) + C * math.exp(-D * Tstar) + E * math.exp(-F * Tstar)


class CombustionProducts:
    """CH4 + 空気の燃焼生成物 (凍結組成)。"""

    def __init__(self, Tt, T_react=300.0, db=None):
        self.db = db if db is not None else ResolvedSpeciesDB.builtin()
        self.Tt = float(Tt)
        self.T_react = float(T_react)
        self.phi = self._solve_phi()
        self.Y = self._products(self.phi)
        self.MW = {s: self.db.MW(s) for s in self.Y}          # kg/mol
        self.R = R_UNIV * sum(y / self.MW[s] for s, y in self.Y.items())   # J/(kg·K)
        self.X = self._mole(self.Y)

    # --- 組成 ---
    def _products(self, phi):
        """当量比 phi の完全燃焼生成物 (質量分率)。CH4 + 2 O2 -> CO2 + 2 H2O。"""
        n_air = 1.0                                            # 空気 1 mol 基準
        x = dict(AIR_X)
        n_o2 = x["O2"] * n_air
        n_ch4 = phi * n_o2 / 2.0                               # 化学量論は CH4:O2 = 1:2
        n = {
            "N2":  x["N2"] * n_air,
            "AR":  x["AR"] * n_air,
            "O2":  n_o2 - 2.0 * n_ch4,
            "CO2": x["CO2"] * n_air + n_ch4,
            "H2O": 2.0 * n_ch4,
        }
        tot_m = sum(v * self.db.MW(s) for s, v in n.items())
        return {s: v * self.db.MW(s) / tot_m for s, v in n.items() if v > 0}

    def _mole(self, Y):
        n = {s: y / self.MW[s] for s, y in Y.items()}
        tot = sum(n.values())
        return {s: v / tot for s, v in n.items()}

    def _fuel_frac_of(self, phi):
        """当量比 phi のときの燃料質量分率 (反応物全体に対する)。"""
        x = dict(AIR_X)
        n_ch4 = phi * x["O2"] / 2.0
        m_air = sum(x[s] * self.db.MW(s) for s in x)
        m_f = n_ch4 * 0.0160425                                  # CH4 の MW [kg/mol]
        return m_f / (m_f + m_air)

    def _solve_phi(self):
        """断熱火炎温度が Tt になる当量比を二分法で解く。

        NASA-9 の `h_mass` は顕熱 (基準 ~298 K) なので、反応熱は LHV で与える:
            f · LHV = h_prod(Tt) − h_prod(T_react)     (定圧・断熱・完全燃焼)
        """
        def f(phi):
            Y = self._products(phi)
            h = self.db.h_mass(Y, np.array([self.Tt, self.T_react]))
            return self._fuel_frac_of(phi) * LHV_CH4 - float(h[0] - h[1])
        lo, hi = 1e-4, 0.999
        flo, fhi = f(lo), f(hi)
        if flo * fhi > 0:
            raise ValueError(f"Tt={self.Tt} K が当量比 {lo}–{hi} の範囲外 (f={flo:.3g},{fhi:.3g})")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    # --- 熱力学 ---
    def cp(self, T):
        return float(self.db.cp_mass(self.Y, np.atleast_1d(np.asarray(T, float)))[0])

    def h(self, T):
        return float(self.db.h_mass(self.Y, np.atleast_1d(np.asarray(T, float)))[0])

    def gamma(self, T):
        cp = self.cp(T)
        return cp / (cp - self.R)

    def a(self, T):
        return math.sqrt(self.gamma(T) * self.R * T)

    # --- 輸送 ---
    def _mu_species(self, s, T):
        sig, eps = LJ[s]
        Ts = T / eps
        Mg = self.MW[s] * 1.0e3                       # g/mol (Chapman–Enskog の式は g/mol)
        return 2.6693e-6 * math.sqrt(Mg * T) / (sig ** 2 * _omega_mu(Ts))   # Pa·s

    def _lam_species(self, s, T):
        """modified Eucken: λ = μ (1.32 cv + 1.77 R_s)。"""
        mu = self._mu_species(s, T)
        Rs = R_UNIV / self.MW[s]
        cps = float(self.db.cp_mass({s: 1.0}, np.array([T]))[0])
        cvs = cps - Rs
        return mu * (1.32 * cvs + 1.77 * Rs)

    def _phi_wilke(self, i, j, mu_i, mu_j):
        Mi, Mj = self.MW[i], self.MW[j]
        num = (1.0 + math.sqrt(mu_i / mu_j) * (Mj / Mi) ** 0.25) ** 2
        den = math.sqrt(8.0 * (1.0 + Mi / Mj))
        return num / den

    def mu(self, T):
        mus = {s: self._mu_species(s, T) for s in self.X}
        out = 0.0
        for i, xi in self.X.items():
            den = sum(self.X[j] * self._phi_wilke(i, j, mus[i], mus[j]) for j in self.X)
            out += xi * mus[i] / den
        return out

    def lam(self, T):
        mus = {s: self._mu_species(s, T) for s in self.X}
        lams = {s: self._lam_species(s, T) for s in self.X}
        out = 0.0
        for i, xi in self.X.items():
            den = sum(self.X[j] * self._phi_wilke(i, j, mus[i], mus[j]) for j in self.X)
            out += xi * lams[i] / den
        return out

    def Pr(self, T):
        return self.mu(T) * self.cp(T) / self.lam(T)

    def summary(self):
        return dict(Tt=self.Tt, phi=self.phi, fuel_mass_frac=self._fuel_mass_frac(),
                    Y={k: round(v, 6) for k, v in self.Y.items()},
                    X={k: round(v, 6) for k, v in self.X.items()},
                    R=self.R)

    def _fuel_mass_frac(self):
        return self._fuel_frac_of(self.phi)


def _self_check():
    """空気で輸送モデルを検算 (CE + Wilke が既知値を再現するか)。"""
    db = ResolvedSpeciesDB.builtin()
    air = CombustionProducts.__new__(CombustionProducts)
    air.db = db
    air.Y = {"N2": 0.75518, "O2": 0.23139, "AR": 0.012885, "CO2": 0.000545}
    air.MW = {s: db.MW(s) for s in air.Y}
    air.R = R_UNIV * sum(y / air.MW[s] for s, y in air.Y.items())
    air.X = air._mole(air.Y)
    rows = []
    for T, mu_ref, pr_ref in ((300.0, 1.846e-5, 0.707), (600.0, 3.017e-5, 0.702), (1000.0, 4.152e-5, 0.726)):
        rows.append((T, air.mu(T), mu_ref, air.Pr(T), pr_ref))
    print("空気での検算 (参考値: Poling/NIST)")
    print(f"{'T [K]':>7} {'mu CE':>11} {'mu ref':>11} {'差 [%]':>8} {'Pr CE':>8} {'Pr ref':>8}")
    ok = True
    for T, mu, mur, pr, prr in rows:
        d = 100 * (mu / mur - 1)
        ok &= abs(d) < 5.0 and abs(pr - prr) < 0.05
        print(f"{T:7.0f} {mu:11.4e} {mur:11.4e} {d:8.2f} {pr:8.3f} {prr:8.3f}")
    print("VERDICT:", "PASS" if ok else "FAIL", "(空気で μ ±5 % / Pr ±0.05)")
    return ok


if __name__ == "__main__":
    _self_check()
    print()
    for Tt in (1650.0, 1880.0):
        g = CombustionProducts(Tt)
        s = g.summary()
        print(f"--- Tt = {Tt:.0f} K ---")
        print(f"  当量比 phi = {s['phi']:.4f}   燃料質量分率 = {s['fuel_mass_frac']*100:.2f} %")
        print(f"  Y = {s['Y']}")
        print(f"  R = {s['R']:.2f} J/(kg·K)")
        for T in (180.0, 294.0, 570.0, 1000.0, Tt):
            print(f"  T={T:7.1f} K : cp={g.cp(T):7.1f}  gamma={g.gamma(T):.4f}  "
                  f"mu={g.mu(T):.4e}  lam={g.lam(T):.4e}  Pr={g.Pr(T):.4f}")
