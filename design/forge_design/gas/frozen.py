r"""凍結組成 (frozen) の thermally-perfect 混合気体 — SERN ⑤ R3 (plan §5.1 R3, codex C2 採用 2026-09-13)。

排気 = CEA 平衡組成 (燃焼器出口 station 3) を**凍結**した多成分 TP、外気 = 空気。forge では両者を
**擬似種 2 種** (`EXH`, `AIR`; `mixture_pseudo_species` で NASA-9 を厳密混合) の `thermalMethod: 2` として解く
(多成分 implicit 不安定 [[wys-tp-divergence-is-cold-not-multispecies]] を避けつつ、せん断層の混合は 2 種の輸送で表す)。
入口状態 (ρ, u)・IC の内部エネルギー・理想推力の正規化を**同じ NASA-9 物性**で計算するための Python 側の実体。

熱力学 (質量基準, NASA-9): $c_p(T)$, $h(T)$ (絶対; 生成エンタルピー込み), $s^0(T)$ (1 bar 基準)。
forge の `thermoHrefTemp` (sensible datum: $h_s(T_{\rm ref})=0$) と同じ基準で $h_{\rm sens}$ を出し、IC の
$\rho e = \rho\,(h_{\rm sens} - R T) + \tfrac12 \rho |u|^2$ を組む。
"""
from __future__ import annotations

import numpy as np

from .semiperfect import LJ_PARAMS, RU, SPECIES_NASA9, T_FLOOR, T_MID, _cp_R_raw, _h_RT_raw, mixture_pseudo_species

# 乾燥空気 (モル分率; CO2 は微量なので N2 に含める)
AIR_MOLE = {"N2": 0.78084, "O2": 0.20946, "AR": 0.00934}
# 1 bar 基準 (NASA-9 の a8 は 1 bar)
P_STD = 1.0e5


def mole_to_mass(x: dict, db=None) -> dict:
    """モル分率 → 質量分率 (`gas.composition.mole_to_mass` へ委譲; 単一ソース)。"""
    from .composition import mole_to_mass as _m2m
    return _m2m(x, db)


def _s0_R_raw(a, T):
    """NASA-9 の標準エントロピー s0/R (1 bar)。"""
    return (-a[0] / (2.0 * T * T) - a[1] / T + a[2] * np.log(T) + a[3] * T + a[4] * T * T / 2.0
            + a[5] * T ** 3 / 3.0 + a[6] * T ** 4 / 4.0 + a[8])


class FrozenGas:
    """凍結組成の TP 混合気体 (質量分率 Y)。`name` は forge 擬似種名。"""

    def __init__(self, Y_mass: dict, name: str = "MIX", href_T: float = 298.15, db=None):
        # db: ResolvedSpeciesDB (省略時は内蔵)。輸送種 (lump) の熱力学は lump 内質量分率でこのクラスを作れば厳密に同じ
        from .composition import ResolvedSpeciesDB
        self._db = db if db is not None else ResolvedSpeciesDB.builtin()
        Y = {k.upper(): float(v) for k, v in Y_mass.items()}
        tot = sum(Y.values())
        self.Y = {k: v / tot for k, v in Y.items() if v > 0.0}
        self._db.require(self.Y)
        self.name = name; self.href_T = float(href_T)
        self.MW = 1.0 / sum(y / self._db.MW(k) for k, y in self.Y.items())
        self.R = RU / self.MW
        self._h_ref = float(self._h_abs(np.array([self.href_T]))[0]) if self.href_T > 0 else 0.0

    @classmethod
    def from_mole(cls, x_mole: dict, name: str = "MIX", href_T: float = 298.15, db=None) -> "FrozenGas":
        return cls(mole_to_mass(x_mole, db), name, href_T, db)

    @classmethod
    def air(cls, href_T: float = 298.15, db=None) -> "FrozenGas":
        return cls.from_mole(AIR_MOLE, "AIR", href_T, db)

    # --- 種ごとの多項式評価 (T_FLOOR 未満は cp 凍結 = semiperfect.py と同じ約束) ---
    def _per_species(self, T, fn_raw, fn_floor):
        T = np.atleast_1d(np.asarray(T, dtype=float)); out = np.zeros_like(T)
        for k, y in self.Y.items():
            e = self._db[k]; w = y * RU / e.MW
            lo, hi = np.asarray(e.low, float), np.asarray(e.high, float)
            v = np.where(T < e.Tmid, fn_raw(lo, np.maximum(T, T_FLOOR)), fn_raw(hi, T))
            if fn_floor is not None:
                v = np.where(T < T_FLOOR, fn_floor(lo, T), v)
            out += w * v
        return out

    def cp_mass(self, T):
        return self._per_species(T, _cp_R_raw, lambda a, T: _cp_R_raw(a, T_FLOOR) * np.ones_like(T))

    def _h_abs(self, T):
        # h = R T (h/RT); T<T_FLOOR: h(TF) − cp_F (TF − T)
        def raw(a, T):
            return _h_RT_raw(a, T) * T
        def floor(a, T):
            return _h_RT_raw(a, T_FLOOR) * T_FLOOR - _cp_R_raw(a, T_FLOOR) * (T_FLOOR - T)
        return self._per_species(T, raw, floor)

    def h_mass(self, T):
        """絶対エンタルピー (生成込み) [J/kg]。"""
        return self._h_abs(T)

    def h_sens(self, T):
        """forge `thermoHrefTemp` と同じ基準 (h_sens(href_T) = 0) の顕エンタルピー [J/kg]。"""
        return self._h_abs(T) - self._h_ref

    def e_sens(self, T):
        return self.h_sens(T) - self.R * np.atleast_1d(np.asarray(T, float))

    def s0_mass(self, T):
        """標準エントロピー (1 bar) [J/kg/K]。T<T_FLOOR は s(TF) − cp_F ln(TF/T)。"""
        return self._per_species(T, _s0_R_raw, lambda a, T: _s0_R_raw(a, T_FLOOR) - _cp_R_raw(a, T_FLOOR) * np.log(T_FLOOR / np.maximum(T, 1.0)))

    def entropy(self, T, P):
        return self.s0_mass(T) - self.R * np.log(np.asarray(P, float) / P_STD)

    def gamma(self, T):
        cp = self.cp_mass(T); return cp / (cp - self.R)

    def a(self, T):
        return np.sqrt(self.gamma(T) * self.R * np.atleast_1d(np.asarray(T, float)))

    def state(self, M: float, P: float, T: float) -> dict:
        ro = P / (self.R * T); a = float(self.a(T)[0]); u = M * a
        return {"M": float(M), "P": float(P), "T": float(T), "ro": float(ro), "u": float(u), "a": a, "R": self.R,
                "gamma_T": float(self.gamma(T)[0]), "cp_T": float(self.cp_mass(T)[0]), "h_sens": float(self.h_sens(T)[0]),
                "gas": self.name}

    def pseudo_species_db(self) -> dict:
        return mixture_pseudo_species(self.Y, self.name, db=self._db)

    def summary(self) -> dict:
        return {"name": self.name, "Y": self.Y, "MW": self.MW, "R": self.R, "href_T": self.href_T}


def isentropic_T(gas: FrozenGas, T1: float, p1: float, p2: float) -> float:
    """(T1, p1) から p2 まで等エントロピー変化した温度 (s0(T) − R ln p = 一定)。"""
    from scipy.optimize import brentq
    target = float(gas.entropy(T1, p1)[0])
    f = lambda T: float(gas.entropy(T, p2)[0]) - target
    lo, hi = (50.0, T1) if p2 < p1 else (T1, 20000.0)
    return float(brentq(f, lo, hi, xtol=1e-9))


def ideal_gross_thrust_frozen(gas: FrozenGas, M_in: float, T_in: float, p_in: float, p_a: float) -> tuple:
    """入口一様流 (M_in, T_in, p_in; H = 1) を p_a まで等エントロピー膨張した理想総推力 F/(p_in H) = ṁ u_e /(p_in H) と出口 M。
    CPG (`rao_planar.ideal_gross_thrust`) と同じ正規化。定数 cp の種 (Ar) では CPG と一致する。"""
    st = gas.state(M_in, p_in, T_in)
    h0 = float(gas.h_mass(T_in)[0]) + 0.5 * st["u"] ** 2
    T_e = isentropic_T(gas, T_in, p_in, p_a)
    u_e = float(np.sqrt(max(2.0 * (h0 - float(gas.h_mass(T_e)[0])), 0.0)))
    M_e = u_e / float(gas.a(T_e)[0])
    return st["ro"] * st["u"] * u_e / p_in, M_e, T_e
