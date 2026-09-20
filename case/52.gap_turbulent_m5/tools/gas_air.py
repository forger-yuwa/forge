#!/usr/bin/env python3
"""case/52 の気体 — 乾燥空気の semi-perfect (TP)。case/50 の CombustionProducts と同じ API を出す。

熱力学: NASA-9 (ResolvedSpeciesDB)、輸送: Chapman–Enskog + Wilke / Mason–Saxena (forge の viscMethod 2 と同族)。
組成は case/49 と同じ (N2 0.75518 / O2 0.23139 / AR 0.012885 / CO2 0.000545)。
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "case" / "50.deep_cavity_wieting_m7" / "tools"))
sys.path.insert(0, str(HERE.parents[2] / "design"))
from gas_model import CombustionProducts, R_UNIV                 # noqa: E402
from forge_design.gas.composition import ResolvedSpeciesDB       # noqa: E402

DRY_AIR = {"N2": 0.75518, "O2": 0.23139, "AR": 0.012885, "CO2": 0.000545}


def dry_air(Tt: float, Y: dict | None = None) -> CombustionProducts:
    """乾燥空気 TP (燃焼しないので当量比の逆算は行わない)。"""
    db = ResolvedSpeciesDB.builtin()
    g = CombustionProducts.__new__(CombustionProducts)
    g.db = db
    g.Tt = float(Tt)
    g.phi = 0.0
    g.Y = dict(Y or DRY_AIR)
    g.MW = {s: db.MW(s) for s in g.Y}
    g.R = R_UNIV * sum(y / g.MW[s] for s, y in g.Y.items())
    g.X = g._mole(g.Y)
    return g


if __name__ == "__main__":
    g = dry_air(1222.73)
    print(f"R = {g.R:.3f} J/(kg·K)   (case/49 derived: 287.048)")
    for T in (216.65, 500.0, 700.0, 1126.3):
        print(f"  T={T:7.2f} : cp={g.cp(T):7.2f}  gamma={g.gamma(T):.4f}  "
              f"mu={g.mu(T):.4e}  lam={g.lam(T):.4e}  Pr={g.Pr(T):.4f}")
