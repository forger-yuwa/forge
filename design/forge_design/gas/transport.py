"""NS (`physProp.viscMethod: 2`) と同じ混合気の分子粘性 μ(T) — 積分法 (CONTUR) の `contur_v2` 用
(plan tooling-nozzle-isothermal-wall-chain §4.7-3、§5.1 #10c)。

データはソルバと同じ共通ファイル (`solver_density_cuda/data/species/forge_transport_v1.yaml` の CEA の種ごとのフィットと相互作用フィット、
`forge_species_v1.yaml` の MW) を読む。モデル:
  - `cea`: ln η = A ln T + B/T + C/T² + D [µP]、区間は「T ≤ Thi の最初の区間、無ければ最後」(cea2.f TRANIN)。
  - `custom:h2o_iapws_cea_v1`: T ≤ 500 K は IAPWS 2008 の希薄気体 μ₀、700 K 以上は CEA、その間は log の smoothstep でつなぐ。
    253.15 K 未満は 253.15 K の両対数勾配で外挿 (ソルバ・独立参照 `tests/unit/transport_reference.py` と同じ)。
  - 混合則: CEA (cea2.f TRANP) μ = Σ_i η_i X_i / Σ_j φ_ij X_j、φ_ij = 2 W_j η_i / (η_ij (W_i + W_j)) (i ≠ j)、φ_ii = 1。
    η_ij は相互作用フィット、無い組は剛体球 (cea2.f 5565–5570)。
組成は凍結 (境界層の中で一定) とし、質量分率 Y からモル分率 X を作る。
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import yaml

DATA = Path(__file__).resolve().parents[3] / "solver_density_cuda" / "data" / "species"
_H_IAPWS = [1.67752, 2.20462, 0.6366564, -0.241605]
_TC = 647.096


def _fit(rows, T: float) -> float:
    use = next((r for r in rows if T <= r[1]), rows[-1])
    A, B, C, D = use[2:6]
    return math.exp(A * math.log(T) + B / T + C / T ** 2 + D)


def _iapws_mu(T: float) -> float:
    Tb = T / _TC
    return 100.0 * math.sqrt(Tb) / sum(_H_IAPWS[k] / Tb ** k for k in range(4)) * 1e-6


def _iapws_slope(T: float) -> float:
    Tb = T / _TC
    S = sum(c * Tb ** (-k) for k, c in enumerate(_H_IAPWS))
    dS = sum(-k * c * Tb ** (-k - 1) for k, c in enumerate(_H_IAPWS))
    return 0.5 - Tb * dS / S


class MixtureViscosity:
    """μ(T) [Pa·s]。Y: {種: 質量分率} (和で正規化)、models: {種: "cea" | "custom:h2o_iapws_cea_v1"} (problem の gas.transport)。"""

    def __init__(self, Y: dict, models: dict, transport_yaml: Path | None = None, species_yaml: Path | None = None):
        tr = yaml.safe_load(open(transport_yaml or DATA / "forge_transport_v1.yaml"))
        sp = yaml.safe_load(open(species_yaml or DATA / "forge_species_v1.yaml"))
        self._V = {e["id"]: e["V"] for e in tr["species"] if e.get("V")}
        self._P = {}
        for e in tr.get("interactions", []):
            if e.get("V"):
                a, b = e["pair"]
                self._P[(a, b)] = self._P[(b, a)] = e["V"]
        mw = {e["id"]: float(e["MW"]) for e in sp["species"] if "MW" in e}
        names = [k for k, v in Y.items() if float(v) > 0.0]
        missing = [k for k in names if k not in models]
        if missing:
            raise ValueError(f"MixtureViscosity: 輸送モデルの無い種 {missing} (gas.transport に書く)")
        for k in names:
            if models[k] not in ("cea", "custom:h2o_iapws_cea_v1"):
                raise ValueError(f"MixtureViscosity: 種 {k} のモデル {models[k]!r} は未対応 (cea | custom:h2o_iapws_cea_v1)")
            if models[k] == "cea" and k not in self._V:
                raise ValueError(f"MixtureViscosity: 種 {k} の CEA 粘性フィットが無い")
            if models[k] == "custom:h2o_iapws_cea_v1" and k != "H2O":
                raise ValueError("custom:h2o_iapws_cea_v1 は H2O 専用")
            if k not in mw:
                raise ValueError(f"MixtureViscosity: 種 {k} の MW が forge_species_v1.yaml に無い")
        self.names = names
        self.models = {k: models[k] for k in names}
        self.W = np.array([mw[k] * 1e3 for k in names])                    # [g/mol] (φ では比だけが効く)
        y = np.array([float(Y[k]) for k in names]); y = y / y.sum()
        n = y / self.W
        self.X = n / n.sum()

    def species_mu(self, k: str, T: float) -> float:
        if self.models[k] == "cea":
            return _fit(self._V[k], T) * 1e-7
        # custom:h2o_iapws_cea_v1
        cm = lambda t: _fit(self._V["H2O"], t) * 1e-7  # noqa: E731
        if T >= 700.0:
            return cm(T)
        if T < 253.15:
            return _iapws_mu(253.15) * (T / 253.15) ** _iapws_slope(253.15)
        if T <= 500.0:
            return _iapws_mu(T)
        s = (T - 500.0) / 200.0
        w = 3.0 * s * s - 2.0 * s ** 3
        return math.exp((1.0 - w) * math.log(_iapws_mu(T)) + w * math.log(cm(T)))

    def _pair(self, i: int, j: int, T: float, eta) -> float:
        a, b = self.names[i], self.names[j]
        if (a, b) in self._P:
            return _fit(self._P[(a, b)], T) * 1e-7
        Wi, Wj = self.W[i], self.W[j]
        ratio = math.sqrt(Wj / Wi)
        e = 5.656854 * eta[i] * math.sqrt(Wj / (Wi + Wj))
        return e / (1.0 + math.sqrt(ratio * eta[i] / eta[j])) ** 2

    def mu_scalar(self, T: float) -> float:
        T = float(T)
        if not (T > 0.0 and math.isfinite(T)):
            raise ValueError(f"MixtureViscosity: 温度が不正 ({T})")
        n = len(self.names)
        eta = [self.species_mu(k, T) for k in self.names]
        if n == 1:
            return eta[0]
        mu = 0.0
        for i in range(n):
            sv = 0.0
            for j in range(n):
                if i == j:
                    phi = 1.0
                else:
                    phi = 2.0 / (self._pair(i, j, T, eta) * (self.W[i] + self.W[j])) * self.W[j] * eta[i]
                sv += phi * self.X[j]
            mu += eta[i] * self.X[i] / sv
        return mu

    def __call__(self, T):
        if np.ndim(T) == 0:
            return self.mu_scalar(float(T))
        return np.array([self.mu_scalar(float(t)) for t in np.ravel(T)]).reshape(np.shape(T))
