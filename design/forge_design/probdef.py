"""問題定義 YAML の読込と検証。

1 案件 = 1 YAML。パラメータは spec (仕様固定) / derived (派生・閉ループ) /
dv (設計変数) の 3 区分で宣言する (親計画 §4.1, §4.6(a))。
Phase 0 は type: thruster_bell のみ実装。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import yaml

KNOWN_TYPES = ("thruster_bell", "wind_tunnel_axisym", "wind_tunnel_axisym_walldriven",
               "wind_tunnel_axisym_axismach", "sern_2d")


@dataclass
class Problem:
    name: str
    type: str
    gamma: float
    cp: float
    spec: dict
    dv: dict
    geometry: dict
    mesh: dict
    evaluate: dict
    raw: dict = field(repr=False, default_factory=dict)
    path: str | None = field(repr=False, default=None)   # 問題 YAML の所在 (gas.species_db の相対パス基準)

    @property
    def R_gas(self) -> float:
        return self.cp * (self.gamma - 1.0) / self.gamma

    # --- ガスモデル (2026-08-17): gas.model = cpg (既定) | semiperfect ---
    # semiperfect: gas.species {名: 質量分率}、Tt は spec.Tt。NASA-9 (CEA) の
    # thermally-perfect・frozen 組成。設計 (MOC) と CFD (forge TP 擬似種) で同一係数。
    @property
    def gas_model(self):
        gs = self.raw.get("gas", {})
        kind = str(gs.get("model", "cpg"))
        if kind == "cpg":
            from .gas import GasCPG
            return GasCPG(self.gamma, self.cp)
        if kind == "semiperfect":
            from .gas import GasSemiPerfect
            Y, _, _ = self.gas_composition
            return GasSemiPerfect(Y, Tt=float(self.spec["Tt"]), db=self.species_db)
        if kind == "frozen_tp":
            # SERN ⑤ R3: 逆設計 (平面 MOC) は設計点 γ の CPG のまま (形状パラメータ化)。CFD・入口状態・正規化は
            # runner 側で FrozenGas (排気 = CEA 凍結組成, 外気 = 空気) を使う
            from .gas import GasCPG
            return GasCPG(self.gamma, self.cp)
        raise ValueError(f"gas.model '{kind}' は未知 (cpg | semiperfect | frozen_tp)")

    # --- 壁の熱境界条件 (2026-09-12, plan tooling-nozzle-isothermal-wall-chain §4.1) ---
    # spec.wall_thermal: {mode: adiabatic} (既定) | {mode: isothermal, Tw: <K>}。
    # bcond (wall / wall_isothermal+Ts)・積分法初期壁 (thermal_bc)・帳簿の 3 箇所が全てここを読む (単一ソース)。
    @property
    def wall_thermal(self) -> dict:
        wt = self.spec.get("wall_thermal") or {"mode": "adiabatic"}
        mode = str(wt.get("mode", "adiabatic"))
        if mode == "adiabatic":
            return {"mode": "adiabatic"}
        if mode == "isothermal":
            Tw = float(wt["Tw"])
            if not Tw > 0.0:
                raise ValueError("spec.wall_thermal.Tw は正の温度 [K]")
            return {"mode": "isothermal", "Tw": Tw}
        raise ValueError(f"spec.wall_thermal.mode '{mode}' は未知 (adiabatic | isothermal)")

    @property
    def wall_thermal_bc_integral(self) -> dict:
        """積分法初期壁 (`feedback/deltastar_integral.integral_bl`) に渡す thermal_bc。"""
        wt = self.wall_thermal
        if wt["mode"] == "isothermal":
            return {"mode": "prescribed_temperature", "Tw": wt["Tw"]}
        return {"mode": "adiabatic"}

    def wall_bcond_line(self, euler: bool, phys_id: int = 3, output: int = 1) -> str:
        """forge bcondConfig の壁 1 行。Euler は slip、NS は wall_thermal に従い wall / wall_isothermal (Ts)。"""
        if euler:
            return f"{{physID: {phys_id}, kind: slip,             outputHDFflg: {output}, ints: , floats: }}"
        wt = self.wall_thermal
        if wt["mode"] == "isothermal":
            return (f"{{physID: {phys_id}, kind: wall_isothermal,  outputHDFflg: {output}, ints: , "
                    f"floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {wt['Tw']}}}}}")
        return f"{{physID: {phys_id}, kind: wall,             outputHDFflg: {output}, ints: , floats: }}"

    # --- 組成の単一ソース (2026-09-16, plan thermophysics-cea-mole-fraction-species §4.1) ---
    @property
    def species_db(self):
        """解決済み種 DB: 内蔵 SPECIES_NASA9 に `gas.species_db` (問題 YAML の所在基準の相対パス可) を上書きしたもの。
        設計 (MOC)・擬似種・IC・`species_db.yaml` 出力の全経路がこれを使う。"""
        from .gas.composition import ResolvedSpeciesDB
        path = self.raw.get("gas", {}).get("species_db")
        if not path:
            return ResolvedSpeciesDB.builtin()
        from pathlib import Path
        pth = Path(str(path))
        if not pth.is_absolute() and self.path:
            pth = Path(self.path).resolve().parent / pth
        if not pth.exists():
            raise FileNotFoundError(f"gas.species_db '{path}' が無い (解決: {pth})")
        return ResolvedSpeciesDB.from_file(pth)

    @property
    def composition_basis(self) -> str:
        return str(self.raw.get("gas", {}).get("composition_basis", "mass")).lower()

    @property
    def gas_composition(self) -> tuple:
        """(Y 正規化, X 正規化, 入力の総和) — `gas.species` を `gas.composition_basis` (mass | mole) で解釈。"""
        from .gas.composition import composition_to_mass
        gs = self.raw.get("gas", {})
        if "species" not in gs:
            raise ValueError("gas.species が無い")
        return composition_to_mass(dict(gs["species"]), self.composition_basis, self.species_db)

    @property
    def condensing_species(self) -> str | None:
        """凝縮種の名前 (正本)。`gas.condensing_species`、無ければ旧 `evaluate.tp_keep_species`、それも無ければ組成に H2O があれば H2O。"""
        gs = self.raw.get("gas", {})
        name = gs.get("condensing_species") or self.evaluate.get("tp_keep_species")
        if name:
            return str(name).upper()
        if "species" in gs and any(str(k).upper() == "H2O" for k in gs["species"]):
            return "H2O"
        return None

    def species_layout(self, streams: dict | None = None):
        """統一 `tp_species` スキーマを解決した輸送種配置 (ノズル: 単一流れ `inflow`)。SERN は runner が流れ辞書を渡す。"""
        from .gas.composition import parse_tp_species, resolve_species_layout
        if streams is None:
            Y, _, _ = self.gas_composition
            streams = {"inflow": Y}
        cond_on = bool(self.evaluate.get("condensation"))
        return resolve_species_layout(parse_tp_species(self.evaluate), streams, self.species_db,
                                      self.condensing_species, condensation=cond_on)

    @property
    def is_semiperfect(self) -> bool:
        return str(self.raw.get("gas", {}).get("model", "cpg")) == "semiperfect"

    @property
    def is_frozen_tp(self) -> bool:
        """SERN ⑤ R3: 排気 = 凍結組成 TP 擬似種 (`gas.exhaust_composition` [モル分率] / 作動点 `gas.composition`)、外気 = 空気。"""
        return str(self.raw.get("gas", {}).get("model", "cpg")) == "frozen_tp"


def load_problem(path) -> Problem:
    with open(path) as f:
        raw = yaml.safe_load(f)
    errs = []
    for key in ("name", "type", "gas", "spec", "dv", "geometry", "mesh", "evaluate"):
        if key not in raw:
            errs.append(f"必須セクション '{key}' が無い")
    if errs:
        raise ValueError("問題定義 YAML: " + "; ".join(errs))
    if raw["type"] not in KNOWN_TYPES:
        raise ValueError(f"type '{raw['type']}' は未実装 (Phase 0 は {KNOWN_TYPES})")

    prob = Problem(
        name=raw["name"],
        type=raw["type"],
        gamma=float(raw["gas"].get("gamma", 1.4)),
        cp=float(raw["gas"].get("cp", 1004.5)),
        spec=raw["spec"],
        dv=raw["dv"],
        geometry=raw["geometry"],
        mesh=raw["mesh"],
        evaluate=raw["evaluate"],
        raw=raw,
        path=str(path),
    )
    _validate(prob)
    return prob


def _validate(p: Problem) -> None:
    errs = []
    if p.type == "thruster_bell":
        for k in ("Pt", "Tt", "p_ambient", "r_throat"):
            if k not in p.spec:
                errs.append(f"spec.{k} が必要")
        # 過拘束チェック: 面積比 (or 出口径) と L は dv、出口マッハは指定不可
        if "M_exit" in p.spec:
            errs.append("thruster_bell で spec.M_exit は指定不可 (面積比から従属)")
        for k in ("eps",):
            if k not in p.spec:
                errs.append("spec.eps (面積比 Ae/At) が必要")
    elif p.type in ("wind_tunnel_axisym", "wind_tunnel_axisym_walldriven",
                    "wind_tunnel_axisym_axismach"):
        for k in ("Pt", "Tt"):
            if k not in p.spec:
                errs.append(f"spec.{k} が必要")
        # (D_e, r_throat, M_design) は独立指定 2 つまで (過拘束チェック — 親計画 §4.6(a))
        given = [k for k in ("D_e", "r_throat", "M_design") if k in p.spec]
        if len(given) > 2:
            errs.append(f"(D_e, r_throat, M_design) は 2 つまで指定可 (過拘束: {given})")
        if "M_design" not in p.spec:
            errs.append("spec.M_design が必要 (現実装は M_design + r_throat の組を要求)")
        if "r_throat" not in p.spec:
            errs.append("spec.r_throat が必要 (D_e 従属モードは未実装 — 閉ループ派生で追加予定)")
    # gas (2026-09-16): 組成の基準・種名・凝縮種・tp_species スキーマを入力段階で検査 (plan cea-mole-fraction §4.1–4.2)
    gs = p.raw.get("gas", {})
    if str(gs.get("model", "cpg")) == "semiperfect":
        try:
            if p.composition_basis not in ("mass", "mole"):
                errs.append(f"gas.composition_basis '{p.composition_basis}' は未知 (mass | mole)")
            p.gas_composition
            if p.evaluate.get("cfd_gas", "same") != "cpg":
                p.species_layout()
        except (ValueError, KeyError, FileNotFoundError) as ex:
            errs.append(f"gas: {ex}")
    elif "composition_basis" in gs or "species_db" in gs or "condensing_species" in gs:
        if str(gs.get("model", "cpg")) not in ("semiperfect", "frozen_tp"):
            errs.append("gas.composition_basis / species_db / condensing_species は gas.model semiperfect | frozen_tp でのみ有効")
    # dv の bound 検査
    for name, d in p.dv.items():
        if isinstance(d, dict) and not d.get("fixed", False):
            if "min" in d and "max" in d:
                v = float(d.get("value", d["min"]))
                if not (d["min"] <= v <= d["max"]):
                    errs.append(f"dv.{name}: value {v} が [min,max] 外")
    if errs:
        raise ValueError("問題定義検証: " + "; ".join(errs))


def dv_value(p: Problem, name: str, default=None):
    d = p.dv.get(name, default)
    if isinstance(d, dict):
        return float(d["value"])
    if d is None:
        raise KeyError(f"dv.{name} が無い")
    return float(d)
