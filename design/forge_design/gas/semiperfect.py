r"""Semi-perfect (thermally perfect, frozen 組成) 気体モデル — NASA-9 多項式 (CEA)。

forge 本体の内蔵 DB (`solver_density_cuda/input/speciesDB.cpp::speciesDB_builtin`, CEA
McBride–Gordon 2002 thermo.inp そのもの; 設計側はその先頭 2 区間 200–1000–6000 K だけを持つ) と**同じ共通データ** (`solver_density_cuda/data/species/forge_species_v1.yaml`) を読み、
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

from pathlib import Path

import numpy as np

RU = 8.314462618  # J/(mol K)

# NASA-9: cp/R = a0 T^-2 + a1 T^-1 + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
#         h/(RT) = -a0 T^-2 + a1 ln T / T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
# 係数・MW・LJ・原子組成は共通 species データ (solver_density_cuda/data/species/forge_species_v1.yaml; forge 本体の内蔵 DB と
# 同じファイル) の legacy_builtin: design の種から作る (plans/active/thermophysics-solver-owned-species-db.md §5.1 #4)。
# 名前は従来どおり大文字キー (canonical ID `Ar` → `AR`)。canonical ID + 別名表への移行は plan #8。
# 注: cea_thermo_to_species_db.py はこのファイルをパスで単独 import するので、ここでは相対 import をしない。
SPECIES_DATA_FILE = Path(__file__).resolve().parents[3] / "solver_density_cuda" / "data" / "species" / "forge_species_v1.yaml"
SPECIES_DATA_SCHEMA = "forge_species_data_v1"


# 設計側の温度域 (plan thermophysics-solver-owned-species-db §5.1 #13-3, 2026-10-01 決定 案 B)。共通データの種は
# CEA thermo.inp そのもの (200–1000–6000 K の 2 区間、または 6000–20000 K を足した 3 区間) で、設計側は**先頭 2 区間だけ**を持つ。
# 6000 K 超を第 2 区間の外挿で黙って評価しない: 内蔵種の cp/h/s° を T > DESIGN_T_MAX で評価すると例外
# (composition.ResolvedSpeciesDB.species_* と evaluate/ic.py)。ソルバの温度反転も 6000 K でクランプする
# (`cuda_forge/dependentVariables_d.cu` DEPVAR_TMAX)。T < 200 K の扱い (端で cp 固定・h 線形) は従来どおり。
DESIGN_T_BOUNDS = (200.0, 1000.0, 6000.0)
DESIGN_T_MAX = DESIGN_T_BOUNDS[-1]
_DESIGN_T_BOUNDS_3 = DESIGN_T_BOUNDS + (20000.0,)


# LJ パラメータの集合 (共通データの LJ_sets; plan thermophysics-solver-owned-species-db §4.10, #14)。ソルバの physProp.ljSource と
# 同じ規則 (順序付きの集合名リストの先頭から探す) の Python 鏡像。既定はソルバと同じ [gri30, svehla1962] (2026-10-01 ユーザ決定)。
LJ_SET_NAMES = ("gri30", "svehla1962", "legacy_v1")
LJ_SOURCE_DEFAULT = ("gri30", "svehla1962")


def check_lj_source(lj_source=None) -> tuple:
    """ljSource の検査 (None = 既定)。空・未知の集合名・重複は ValueError (ソルバ speciesDB_checkLjSource と同じ)。"""
    if lj_source is None:
        return LJ_SOURCE_DEFAULT
    src = tuple(str(s) for s in lj_source)
    if not src:
        raise ValueError("ljSource が空 (集合名を 1 つ以上: gri30, svehla1962, legacy_v1)")
    for k, s in enumerate(src):
        if s not in LJ_SET_NAMES:
            raise ValueError(f"ljSource: 未知の LJ 集合 {s!r} (gri30, svehla1962, legacy_v1)")
        if s in src[:k]:
            raise ValueError(f"ljSource: LJ 集合 {s!r} が 2 回ある")
    return src


def _load_design_species(path=SPECIES_DATA_FILE):
    """共通データから (SPECIES_NASA9, LJ_SETS, 原子組成) を従来の形・順序・大文字キーで返す。
    LJ_SETS は {従来キー: {集合名: (σ, ε/k_B)}} (どの集合にも無い種は空 dict)。
    区間は [200,1000],[1000,6000] (+ 任意の [6000,20000]) だけを受け、先頭 2 区間を low/high に取る (第 3 区間は捨てる)。
    それ以外の区間構成 (1 区間・非標準の区切り) は ValueError。"""
    import yaml
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("schema") != SPECIES_DATA_SCHEMA:
        raise ValueError(f"{path}: schema が {SPECIES_DATA_SCHEMA} でない")
    nasa9, lj, atoms = {}, {}, {}
    for e in raw["species"]:
        if "design" not in (e.get("legacy_builtin") or []):
            continue
        key = str(e["id"]).upper()   # 従来キー (Ar → AR)
        if key in nasa9:
            raise ValueError(f"{path}: 従来キー {key} が重複 ({e['id']})")
        if e.get("phase") != "gas":
            raise ValueError(f"{path}: {e['id']} の phase {e.get('phase')} は内蔵種に使えない (gas のみ)")
        iv = e["intervals"]
        # 設計側は 2 区間 200–1000–6000 K (T_MID, ResolvedSpeciesDB の既定区切り)。CEA の第 3 区間 6000–20000 K は持たない。
        bounds = [float(iv[0]["Tlo"])] + [float(v["Thi"]) for v in iv]
        contiguous = all(float(iv[k]["Tlo"]) == float(iv[k - 1]["Thi"]) for k in range(1, len(iv)))
        if not contiguous or tuple(bounds) not in (DESIGN_T_BOUNDS, _DESIGN_T_BOUNDS_3):
            raise ValueError(f"{path}: {e['id']} の温度区間 {bounds} が 200–1000–6000 K (+ 任意の 6000–20000 K) でない")
        low, high = [float(v) for v in iv[0]["coeffs"]], [float(v) for v in iv[1]["coeffs"]]
        if len(low) != 9 or len(high) != 9:
            raise ValueError(f"{path}: {e['id']} の係数は 9 個ずつ必要")
        nasa9[key] = dict(MW=float(e["MW"]), low=low, high=high)
        # LJ は出典別の集合 LJ_sets (plan §4.10, #14)。解決は lj_params (ljSource の先頭から探す)
        lj[key] = {s: (float(v["sigma"]), float(v["eps_kB"])) for s, v in (e.get("LJ_sets") or {}).items()}
        if e.get("atoms") is not None:
            atoms[key] = dict(e["atoms"])
    return nasa9, lj, atoms


def check_design_T(T, what="") -> None:
    """内蔵種 (共通データの先頭 2 区間) を T > DESIGN_T_MAX で評価しようとしたら例外 (第 2 区間を外挿しない; §5.1 #13-3)。"""
    Tmax = float(np.nanmax(np.asarray(T, dtype=float))) if np.size(T) else -np.inf
    if Tmax > DESIGN_T_MAX:
        raise ValueError(f"{what}: T = {Tmax!r} K は設計側の温度域 (≤ {DESIGN_T_MAX} K) を超える。"
                         "内蔵種は CEA の先頭 2 区間だけを持ち、6000 K 超を外挿しない (plan thermophysics-solver-owned-species-db #13-3)")


SPECIES_NASA9, LJ_SETS, SPECIES_ATOMS = _load_design_species()


def lj_params(lj_source=None) -> dict:
    """{従来キー: (σ [Å], ε/k_B [K])}: 各種の LJ を ljSource (None = 既定 [gri30, svehla1962]) の先頭から探して最初にある集合の値。
    どの集合にも無い種は含めない (ソルバは LJ なしとして LJ を読む使い方で拒否する)。"""
    src = check_lj_source(lj_source)
    out = {}
    for k, sets in LJ_SETS.items():
        for s in src:
            if s in sets:
                out[k] = sets[s]
                break
    return out


# Lennard-Jones (σ [Å], ε/k_B [K]) は LJ_PARAMS (既定の ljSource [gri30, svehla1962] で解決; ソルバの既定と同じ)。
# 擬似種の輸送係数は質量分率加重 (粗い近似で十分)
LJ_PARAMS = lj_params()
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
