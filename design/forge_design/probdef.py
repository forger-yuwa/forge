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
            return GasSemiPerfect(dict(gs["species"]), Tt=float(self.spec["Tt"]))
        raise ValueError(f"gas.model '{kind}' は未知 (cpg | semiperfect)")

    @property
    def is_semiperfect(self) -> bool:
        return str(self.raw.get("gas", {}).get("model", "cpg")) == "semiperfect"


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
