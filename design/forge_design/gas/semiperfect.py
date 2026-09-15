r"""Semi-perfect (thermally perfect, frozen 組成) 気体モデル — NASA-9 多項式 (CEA)。

forge 本体の内蔵 DB (`solver_density_cuda/cuda_forge/thermo_d.cu::builtinDB`, CEA
McBride–Gordon 2002 の 2 区間 200–1000–6000 K) と**同一係数**を Python 側に持ち、
設計 (MOC・遷音速・面積比) と CFD (forge TP, `thermalMethod: 1`) の熱力学を一致させる。

**MOC が γ に依存する箇所** (これだけ差し替えれば特性線法は thermally perfect でも成立):

- Prandtl–Meyer 関数 $\nu(M)$: 等エントロピー膨張で
  $d\nu = \sqrt{M^2-1}\,\frac{dV}{V}$。$V(T)$ は $h_0 = h(T) + V^2/2$ から、
  $M(T) = V/a(T)$、$a^2 = \gamma(T) R T$ ($\gamma = c_p/c_v$ 局所値)。
  $\nu$ を $T$ でパラメトライズして数値積分し、$M \leftrightarrow \nu$ を単調テーブルで引く。
- 質量流束密度 $\rho V/(\rho_0 a_0)$、面積比 $A/A^*$: 同じ $T$ パラメトライズで閉形式なし
  → テーブル。
- 音速 $a(T)$・$\gamma(T)$: 遷音速解 (Hall) の $\gamma$ には**スロート温度での局所 γ** を渡す
  (Hall 級数は定数 γ 前提。スロート近傍の $\gamma$ 変化は $T$ 変化が小さいので 2 次)。

これは「有効 γ の完全気体」より正確で、「反応平衡 (CEA equilibrium)」より粗い —
組成凍結 (frozen) は膨張ノズルで標準的な近似 (Anderson, Zucrow–Hoffman)。
"""
from __future__ import annotations

import numpy as np

RU = 8.314462618  # J/(mol K)

# NASA-9: cp/R = a0 T^-2 + a1 T^-1 + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
#         h/(RT) = -a0 T^-2 + a1 ln T / T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
# 係数は forge thermo_d.cu builtinDB と同一 (CEA)。
SPECIES_NASA9 = {
    "N2": dict(MW=0.0280134,
               low=[2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03,
                    1.384646189e-05, -9.625793620e-09, 2.519705809e-12, 7.108460860e+02,
                    -1.076003744e+01],
               high=[5.877124060e+05, -2.239249073e+03, 6.066949220e+00, -6.139685500e-04,
                     1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 1.283210415e+04,
                     -1.586640027e+01]),
    "O2": dict(MW=0.0319988,
               low=[-3.425563420e+04, 4.847000970e+02, 1.119010961e+00, 4.293889240e-03,
                    -6.836300520e-07, -2.023372700e-09, 1.039040018e-12, -3.391454870e+03,
                    1.849699470e+01],
               high=[-1.037939022e+06, 2.344830282e+03, 1.819732036e+00, 1.267847582e-03,
                     -2.188067988e-07, 2.053719572e-11, -8.193467050e-16, -1.689010929e+04,
                     1.738716506e+01]),
    "CO2": dict(MW=0.0440095,
                low=[4.943650540e+04, -6.264116010e+02, 5.301725240e+00, 2.503813816e-03,
                     -2.127308728e-07, -7.689988780e-10, 2.849677801e-13, -4.528198460e+04,
                     -7.048279440e+00],
                high=[1.176962419e+05, -1.788791477e+03, 8.291523190e+00, -9.223156780e-05,
                      4.863676880e-09, -1.891053312e-12, 6.330036590e-16, -3.908350590e+04,
                      -2.652669281e+01]),
    "H2O": dict(MW=0.0180153,
                low=[-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,
                     -7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04,
                     1.724205775e+01],
                high=[1.034972096e+06, -2.412698562e+03, 4.646110780e+00, 2.291998307e-03,
                      -6.836830480e-07, 9.426468930e-11, -4.822380530e-15, -1.384286509e+04,
                      -7.978148510e+00]),
    "AR": dict(MW=0.039948,
               low=[0, 0, 2.5, 0, 0, 0, 0, -7.453750000e+02, 4.379674910e+00],
               high=[0, 0, 2.5, 0, 0, 0, 0, -7.453750000e+02, 4.379674910e+00]),
    # --- H2-air 燃焼生成物 (SERN ⑤ R3, 2026-09-13): CEA2 thermo.inp (McBride–Gordon 2002) から転記。
    #     凍結組成の擬似種 (mixture_pseudo_species) に畳むための係数で、forge 内蔵 DB には無い種。
    "H2": dict(MW=0.00201588,
               low=[4.078323210e+04, -8.009186040e+02, 8.214702010e+00, -1.269714457e-02, 1.753605076e-05,
                    -1.202860270e-08, 3.368093490e-12, 2.682484665e+03, -3.043788844e+01],
               high=[5.608128010e+05, -8.371504740e+02, 2.975364532e+00, 1.252249124e-03, -3.740716190e-07,
                     5.936625200e-11, -3.606994100e-15, 5.339824410e+03, -2.202774769e+00]),
    "OH": dict(MW=0.01700734,
               low=[-1.998858990e+03, 9.300136160e+01, 3.050854229e+00, 1.529529288e-03, -3.157890998e-06,
                    3.315446180e-09, -1.138762683e-12, 2.991214235e+03, 4.674110790e+00],
               high=[1.017393379e+06, -2.509957276e+03, 5.116547860e+00, 1.305299930e-04, -8.284322260e-08,
                     2.006475941e-11, -1.556993656e-15, 2.019640206e+04, -1.101282337e+01]),
    "H": dict(MW=0.00100794,
              low=[0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 2.547370801e+04, -4.466828530e-01],
              high=[6.078774250e+01, -1.819354417e-01, 2.500211817e+00, -1.226512864e-07, 3.732876330e-11,
                    -5.687744560e-15, 3.410210197e-19, 2.547486398e+04, -4.481917770e-01]),
    "NO": dict(MW=0.0300061,
               low=[-1.143916503e+04, 1.536467592e+02, 3.431468730e+00, -2.668592368e-03, 8.481399120e-06,
                    -7.685111050e-09, 2.386797655e-12, 9.098214410e+03, 6.728725490e+00],
               high=[2.239018716e+05, -1.289651623e+03, 5.433936030e+00, -3.656034900e-04, 9.880966450e-08,
                     -1.416076856e-11, 9.380184620e-16, 1.750317656e+04, -8.501669090e+00]),
    "O": dict(MW=0.0159994,
              low=[-7.953611300e+03, 1.607177787e+02, 1.966226438e+00, 1.013670310e-03, -1.110415423e-06,
                   6.517507500e-10, -1.584779251e-13, 2.840362437e+04, 8.404241820e+00],
              high=[2.619020262e+05, -7.298722030e+02, 3.317177270e+00, -4.281334360e-04, 1.036104594e-07,
                    -9.438304330e-12, 2.725038297e-16, 3.392428060e+04, -6.679585350e-01]),
    "CO": dict(MW=0.0280101,
               low=[1.489045326e+04, -2.922285939e+02, 5.724527170e+00, -8.176235030e-03, 1.456903469e-05,
                    -1.087746302e-08, 3.027941827e-12, -1.303131878e+04, -7.859241350e+00],
               high=[4.619197250e+05, -1.944704863e+03, 5.916714180e+00, -5.664282830e-04, 1.398814540e-07,
                     -1.787680361e-11, 9.620935570e-16, -2.466261084e+03, -1.387413108e+01]),
}
# Lennard-Jones (σ [Å], ε/k_B [K]; Svehla 1962 / Chemkin transport)。擬似種の輸送係数は質量分率加重 (粗い近似で十分)
LJ_PARAMS = {"N2": (3.621, 97.53), "O2": (3.458, 107.4), "CO2": (3.763, 244.0), "H2O": (2.605, 572.4), "AR": (3.330, 136.5),
             "H2": (2.827, 59.7), "OH": (3.147, 79.8), "H": (2.708, 37.0), "NO": (3.492, 116.7), "O": (3.050, 106.7), "CO": (3.690, 91.7)}
T_MID = 1000.0


def _coef(name, T):
    sp = SPECIES_NASA9[name]
    return np.asarray(sp["low"] if T < T_MID else sp["high"], dtype=float)


# NASA-9 の下限 (200 K) 未満は係数外挿になり cp が 100 K 以下で非物理に増大する。
# T_FLOOR 未満は cp を T_FLOOR の値で凍結し、h は T_FLOOR で連続に接続する (2026-08-17、
# M6 設計で出口 137 K が下限を割るため導入。M4.2 [出口 247 K] では不変)。
# forge 側の擬似種 (mixture_pseudo_species) も同じ凍結を 2 区間 DB で表現し、設計と CFD を揃える。
T_FLOOR = 200.0


def _cp_R_raw(a, T):
    return a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4


def _h_RT_raw(a, T):
    return (-a[0] / T**2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3
            + a[5] * T**3 / 4 + a[6] * T**4 / 5 + a[7] / T)


def _cp_R(a, T):
    T = np.asarray(T, dtype=float)
    return np.where(T < T_FLOOR, _cp_R_raw(a, T_FLOOR), _cp_R_raw(a, T))


def _h_RT(a, T):
    T = np.asarray(T, dtype=float)
    # T<T_FLOOR: h(T) = h(T_FLOOR) - cp_F (T_FLOOR - T)  →  h/(RT) = [h_F/R - cpF (T_FLOOR - T)]/T
    hF_R = _h_RT_raw(a, T_FLOOR) * T_FLOOR
    cpF = _cp_R_raw(a, T_FLOOR)
    frozen = (hF_R - cpF * (T_FLOOR - T)) / np.maximum(T, 1e-30)
    return np.where(T < T_FLOOR, frozen, _h_RT_raw(a, T))


class GasCPG:
    """完全気体 (回帰対照)。`GasSemiPerfect` と同じ API。"""

    def __init__(self, gamma: float = 1.4, cp: float = 1004.5):
        self.gamma_ref = float(gamma)
        self.cp_ref = float(cp)
        self.R = cp * (gamma - 1.0) / gamma
        self.kind = "cpg"

    def gamma(self, T):
        return np.full_like(np.asarray(T, dtype=float), self.gamma_ref)

    def nu(self, M):
        from ..geometry.moc_kernel import pm_nu
        return pm_nu(M, self.gamma_ref)

    def mach_of_nu(self, nu):
        from ..geometry.moc_kernel import pm_mach_vec
        return pm_mach_vec(np.atleast_1d(nu), self.gamma_ref)

    def area_ratio(self, M):
        g = self.gamma_ref
        M = np.asarray(M, dtype=float)
        return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) \
            ** (0.5 * (g + 1.0) / (g - 1.0))

    def T_of_M(self, M, Tt):
        return Tt / (1.0 + 0.5 * (self.gamma_ref - 1.0) * np.asarray(M, float) ** 2)

    def gamma_throat(self, Tt):
        return self.gamma_ref


class GasSemiPerfect:
    r"""thermally perfect・frozen 組成 (質量分率 Y) 気体。全温 $T_t$ 固定でテーブル化。

    使い方: `gas = GasSemiPerfect({"N2":0.72,"CO2":0.15,"H2O":0.11,"O2":0.02}, Tt=1000)`
    → `gas.nu(M)`, `gas.mach_of_nu(nu)`, `gas.area_ratio(M)`, `gas.gamma(T)`,
    `gas.gamma_throat()`, `gas.T_of_M(M)`。
    """

    def __init__(self, Y: dict, Tt: float, T_min: float | None = None, n_tab: int = 6000, db=None):
        # db: ResolvedSpeciesDB (省略時は呼び出し時点の内蔵 SPECIES_NASA9)。設計 (MOC)・擬似種・IC・CFD が同じ DB を使う (plan cea-mole-fraction M1)
        from .composition import ResolvedSpeciesDB
        self._db = db if db is not None else ResolvedSpeciesDB.builtin()
        Y = {k.upper(): float(v) for k, v in Y.items()}
        tot = sum(Y.values())
        self.Y = {k: v / tot for k, v in Y.items()}
        self._db.require(self.Y)
        self.Tt = float(Tt)
        self.kind = "semiperfect"
        # 混合の R (質量基準)
        self.R = RU * sum(y / self._db.MW(k) for k, y in self.Y.items())
        # T テーブル (Tt から T_min まで単調減 = 等エントロピー膨張)。
        # T_min は M~8 相当まで届く比 (CPG γ=1.4 で T/Tt=1/13.8) を既定に。NASA-9 の
        # 下限 200 K より下は forge 本体と同じ「係数クランプ」の外挿になる点に注意。
        if T_min is None:
            T_min = max(self.Tt / 14.0, 60.0)
        T = np.linspace(self.Tt, T_min, n_tab)
        cp = self.cp_mass(T)
        h = self.h_mass(T)
        h0 = float(h[0])
        V2 = 2.0 * (h0 - h)
        V = np.sqrt(np.maximum(V2, 0.0))
        gam = cp / (cp - self.R)
        a = np.sqrt(gam * self.R * T)
        M = V / a
        # ν(M): dν = sqrt(M²-1) dV/V for M>1 (積分は T に沿って)
        sup = M >= 1.0
        i_star = int(np.argmax(sup))                # 最初の M≥1 (ソニック点)
        # ソニック点を精密化 (線形補間) — テーブル起点
        if i_star == 0:
            raise ValueError("Tt が低すぎるか組成が不正: 既に M>1")
        f = (1.0 - M[i_star - 1]) / (M[i_star] - M[i_star - 1])
        T_star = T[i_star - 1] + f * (T[i_star] - T[i_star - 1])
        self.T_star = float(T_star)
        self.gamma_star = float(np.interp(T_star, T[::-1], gam[::-1]))
        # 超音速枝のみでテーブル
        Ts = np.concatenate([[T_star], T[i_star:]])
        cps = self.cp_mass(Ts); hs = self.h_mass(Ts)
        Vs = np.sqrt(np.maximum(2.0 * (h0 - hs), 0.0))
        gs = cps / (cps - self.R)
        a_s = np.sqrt(gs * self.R * Ts)
        Ms = np.maximum(Vs / a_s, 1.0)
        Ms[0] = 1.0
        integrand = np.sqrt(np.maximum(Ms**2 - 1.0, 0.0)) / np.maximum(Vs, 1e-30)
        nu_s = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(Vs))])
        # ρ V / (ρ* a*) と A/A*: ρ ∝ P/(RT), P/Pt = exp(-∫ cp/T dT / R) (等エントロピー、cp(T))
        # s(T)/R = ∫ cp/(R T) dT — 数値積分
        s_R = np.concatenate([[0.0], np.cumsum(0.5 * (cps[1:] / Ts[1:] + cps[:-1] / Ts[:-1])
                                               * np.diff(Ts) / self.R)])   # Ts 減少なので負
        # P/P* = exp(s_R) * (T/T*)^0 … 厳密には ln(P/P*) = ∫ (cp/R) dT/T (等エントロピー)
        P_rel = np.exp(s_R)
        rho_rel = P_rel / (Ts / T_star)
        flux_rel = rho_rel * Vs / (rho_rel[0] * Vs[0])
        self._M = Ms; self._nu = nu_s; self._T = Ts; self._gam = gs
        self._AR = 1.0 / np.maximum(flux_rel, 1e-30)
        self._flux = flux_rel
        # 上流 (亜音速) 枝の A/A* — IC 用
        Tu = T[:i_star][::-1]                       # T_star 側 → Tt 側 (M 減少)
        cpu = self.cp_mass(Tu); hu = self.h_mass(Tu)
        Vu = np.sqrt(np.maximum(2.0 * (h0 - hu), 0.0))
        gu = cpu / (cpu - self.R); au = np.sqrt(gu * self.R * Tu)
        Mu = np.minimum(Vu / au, 1.0)
        s_Ru = np.concatenate([[0.0], np.cumsum(0.5 * (cpu[1:] / Tu[1:] + cpu[:-1] / Tu[:-1])
                                                * np.diff(Tu) / self.R)])
        P_u = np.exp(s_Ru); rho_u = P_u / (Tu / T_star)
        flux_u = rho_u * Vu / (rho_rel[0] * Vs[0])
        self._Mu = Mu[::-1]; self._ARu = (1.0 / np.maximum(flux_u, 1e-30))[::-1]
        self._Tu = Tu[::-1]                            # 亜音速枝の T (M 増加順)

    # --- 混合熱力学 (質量基準) ---
    def cp_mass(self, T):
        return self._db.cp_mass(self.Y, T)

    def h_mass(self, T):
        return self._db.h_mass(self.Y, T)

    def _h_mass_legacy(self, T):   # (旧実装; ResolvedSpeciesDB.h_mass と同式。参照用に残す)
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.zeros_like(T)
        for k, y in self.Y.items():
            sp = SPECIES_NASA9[k]
            hRT = np.where(T < T_MID, _h_RT(np.asarray(sp["low"]), T),
                           _h_RT(np.asarray(sp["high"]), T))
            out += y * hRT * RU * T / sp["MW"]
        return out

    def gamma(self, T):
        cp = self.cp_mass(T)
        return cp / (cp - self.R)

    def gamma_throat(self, Tt=None):
        return self.gamma_star

    # --- MOC が使う関数 ---
    def nu(self, M):
        return np.interp(np.asarray(M, dtype=float), self._M, self._nu)

    def mach_of_nu(self, nu):
        return np.interp(np.asarray(nu, dtype=float), self._nu, self._M)

    def area_ratio(self, M):
        M = np.asarray(M, dtype=float)
        return np.where(M >= 1.0, np.interp(M, self._M, self._AR),
                        np.interp(M, self._Mu, self._ARu))

    def T_of_M(self, M, Tt=None):
        return np.interp(np.asarray(M, dtype=float), self._M, self._T)

    def mass_flux_density(self, M):
        """ρV/(ρ* a*) (M≥1)。CPG の `_mass_flux_density` はよどみ量規格化なので
        比だけ使う場面 (壁閉包の相対流束) では同等。"""
        return np.interp(np.asarray(M, dtype=float), self._M, self._flux)

    def summary(self) -> dict:
        return {"kind": self.kind, "Y": self.Y, "Tt": self.Tt, "R": self.R,
                "T_star": self.T_star, "gamma_star": self.gamma_star,
                "gamma_Tt": float(self.gamma(self.Tt)[0]),
                "gamma_M4": float(self.gamma(self.T_of_M(4.0))[0]) if self._M[-1] > 4 else None,
                "T_M4": float(self.T_of_M(4.0)) if self._M[-1] > 4 else None,
                "AR_M4": float(self.area_ratio(4.0)) if self._M[-1] > 4 else None}


def mixture_pseudo_species(Y: dict, name: str = "MIX", freeze_low_T: bool = False, db=None) -> dict:
    r"""frozen 組成の混合物を **単一の擬似種** (NASA-9) にまとめる (forge の
    `physProp.speciesDBFile` 用)。

    質量基準の $c_p^{\rm mix}(T)=\sum_k Y_k c_{p,k}(T)$ は各種の多項式の線形和なので、
    モル基準 NASA-9 に戻すと $a_i^{\rm mix} = \sum_k Y_k \frac{MW_{\rm mix}}{MW_k}\,a_{i,k}$
    ($i=0..7$、$a_8$ [エントロピー定数] も同式で可) と**厳密に**混合できる。
    $MW_{\rm mix} = 1/\sum_k (Y_k/MW_k)$。

    多成分 TP (`thermalMethod: 2` + 複数 species) の既知の implicit 不安定性
    ([[wys-tp-divergence-is-cold-not-multispecies]]) を避け、**1 種の TP** として
    forge を回すための道具。組成が凍結している設計 (膨張ノズル) ではこれで正確。
    戻り: forge の `speciesDBFile` yaml にそのまま書ける dict。
    """
    from .composition import ResolvedSpeciesDB
    db = db if db is not None else ResolvedSpeciesDB.builtin()
    Y = {k.upper(): float(v) for k, v in Y.items()}
    tot = sum(Y.values()); Y = {k: v / tot for k, v in Y.items()}
    db.require(Y)
    MW_mix = 1.0 / sum(y / db.MW(k) for k, y in Y.items())
    low = np.zeros(9); high = np.zeros(9)
    for k, y in Y.items():
        e = db[k]; w = y * MW_mix / e.MW
        low += w * np.asarray(e.low, dtype=float)
        high += w * np.asarray(e.high, dtype=float)
    # LJ は質量分率加重 (輸送は粗い近似で十分 — 粘性は Sutherland 側で扱う)
    sig = sum(y * db[k].LJ_sigma for k, y in Y.items())
    eps = sum(y * db[k].LJ_eps_kB for k, y in Y.items())
    # forge DB は 2 区間しか持てない。freeze_low_T=True なら設計側の T_FLOOR 凍結と揃えて
    #   Tmid = T_FLOOR: 下 = 定数 cp(T_FLOOR) (h を連続接続), 上 = 元 low 係数 (200–1000 K)
    # とするが、これは元 high (>1000 K) を捨てるので **Tt ≤ 1000 K 専用**。既定 (False) は
    # 元の 200–1000–6000 形式 (Tt=1550 K の案件で必要)。場が 200 K を割る設計 (Tt 1000 K で M6 等)
    # のときだけ True にすること。
    if freeze_low_T:
        cF = float(_cp_R_raw(low, T_FLOOR))
        a7F = float(T_FLOOR * (_h_RT_raw(low, T_FLOOR) - cF))
        const = [0.0, 0.0, cF, 0.0, 0.0, 0.0, 0.0, a7F, float(low[8])]
        return {name: {"MW": float(MW_mix), "LJ_sigma": float(sig), "LJ_eps_kB": float(eps),
                       "Tlo": 50.0, "Tmid": float(T_FLOOR), "Thi": 6000.0,
                       "nasa9_low": const, "nasa9_high": [float(v) for v in low]}}
    return {name: {"MW": float(MW_mix), "LJ_sigma": float(sig), "LJ_eps_kB": float(eps),
                   "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
                   "nasa9_low": [float(v) for v in low],
                   "nasa9_high": [float(v) for v in high]}}


def mixture_pseudo_species_split(Y: dict, keep=("H2O",), name_dry: str = "MIXDRY", db=None) -> tuple:
    r"""**凝縮向け分割**: `keep` の種 (既定 H₂O) を独立種のまま残し、**それ以外を 1 つの擬似種**
    `name_dry` に畳む (計画 plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)。

    戻り: (db, Y_split, order)。`db` は forge `speciesDBFile` 用 dict (擬似種 + 残した種を
    **内蔵と同じ係数で明示的に書き出す** = 自己完結)、`Y_split` は `{name_dry: 1-ΣY_keep, keep...}`、
    `order` は `species:` に書く順序 (`[name_dry, *keep]`; 凝縮種 index は `order.index("H2O")`)。
    混合則は `mixture_pseudo_species` と同じ質量分率線形混合 (dry 部は全種独立と厳密に同じ熱力学)。
    """
    Y = {k.upper(): float(v) for k, v in Y.items()}
    tot = sum(Y.values()); Y = {k: v / tot for k, v in Y.items()}
    keep = tuple(k.upper() for k in keep)
    for k in keep:
        if k not in Y:
            raise KeyError(f"mixture_pseudo_species_split: {k} が組成に無い ({list(Y)})")
    Y_dry = {k: v for k, v in Y.items() if k not in keep}
    y_dry_tot = sum(Y_dry.values())
    if y_dry_tot <= 0.0:
        raise ValueError("mixture_pseudo_species_split: 乾き成分が無い")
    from .composition import ResolvedSpeciesDB
    rdb = db if db is not None else ResolvedSpeciesDB.builtin()
    db = mixture_pseudo_species({k: v / y_dry_tot for k, v in Y_dry.items()}, name_dry, db=rdb)
    for k in keep:
        db[k] = rdb[k].to_db_dict()
    Y_split = {name_dry: y_dry_tot, **{k: Y[k] for k in keep}}
    return db, Y_split, [name_dry, *keep]
