"""NASA TM X-71972 TABLE 1 + CEA2 → SERN 作動点 (plan §4.10) の導出。

TABLE 1 は station 3 (燃焼器出口 = ノズル入口) の p, T, V と当量比 φ・作動モードを与えるが、
M_3 と γ を与えない。そこを CEA2 (`tmx1.inp` / `tmxair.inp`, tp 問題) で埋める。

  p_inf = q_inf / (0.7 M_inf^2)      (q_inf = 71850 N/m^2 = 1500 psf, γ_air 1.4)
  T_inf = US Standard Atmosphere 1976 (p_inf から逆算)
  R     = 8314.462 / Mbar_CEA,  cp = γ R/(γ-1)   ← γ と R を合わせる CPG (音速と EOS が CEA と一致)
  M_3   = V_3 / sqrt(γ R T_3)

実行: python3 case/46.sern_design/cea/tmx_operating_points.py
"""
import re
from pathlib import Path

CEA = Path(__file__).resolve().parent
Q_INF, G_AIR, RU = 71850.0, 1.4, 8314.462

# TABLE 1 (p.31): M_inf, φ, mode, p3 [N/m^2], T3 [K], V3 [m/s], CEA ケース名
TABLE1 = [
    (4,  1.0, "ramjet",   157717, 2343, 1093, "m4_on"),
    (4,  0.0, "off",       17888,  338, 1072, "m4_off"),
    (6,  1.0, "scramjet", 101027, 2328, 1621, "m6_on"),
    (6,  0.0, "off",       16327,  448, 1655, "m6_off"),
    (10, 1.5, "scramjet",  57935, 2222, 2837, "m10on"),
    (10, 0.0, "off",       14172,  772, 2831, "m10off"),
]


def cea_case(case: str) -> tuple[float, float, dict]:
    """CEA の tp 出力から (GAMMAs, Mbar, モル分率 dict) を拾う (R3: 凍結組成の入力)。"""
    for f in ("tmx1.out", "tmxair.out"):
        txt = (CEA / f).read_text().splitlines()
        for i, L in enumerate(txt):
            if L.strip().startswith("CASE =") and L.split("=")[1].strip() == case:
                blk = txt[i:i + 90]
                mw = gam = None; x = {}; in_x = False
                for s_ in blk:
                    t = s_.strip()
                    if t.startswith("M, (1/n)"):
                        mw = float(t.split()[-1])
                    elif t.startswith("GAMMAs"):
                        gam = float(t.split()[-1])
                    elif t.startswith("MOLE FRACTIONS"):
                        in_x = True; continue
                    elif in_x:
                        if not t:
                            continue                      # MOLE FRACTIONS 直後の空行
                        if "THERMODYNAMIC" in t:
                            in_x = False; continue
                        parts = t.split()
                        if len(parts) >= 2 and len(parts) % 2 == 0:
                            for a, b in zip(parts[::2], parts[1::2]):
                                x[a.lstrip("*")] = float(b)
                if mw and gam:
                    return gam, mw, x
    raise KeyError(f"CEA ケース '{case}' が見つからない ({CEA} の .out を作り直すこと)")


def cea_props(case: str) -> tuple[float, float]:
    """CEA の tp 出力から (GAMMAs, Mbar) を拾う。"""
    gam, mw, _ = cea_case(case)
    return gam, mw


def std_atm_T(p: float) -> float:
    """US Standard Atmosphere 1976 の 11-32 km から p [Pa] に対応する T [K]。"""
    if p >= 5474.89:                       # 11-20 km: 等温層
        return 216.65
    return 216.65 * (p / 5474.89) ** (-1.0 / 34.1632)   # 20-32 km: L = +1 K/km


def frozen_props(x_mole: dict, T3: float):
    """凍結組成の (R, cp(T3), γ(T3), a(T3)) — forge_design.gas.frozen.FrozenGas (NASA-9, CEA と同一係数)。"""
    import sys
    sys.path.insert(0, str(CEA.parents[2] / "design"))
    from forge_design.gas.frozen import FrozenGas
    g = FrozenGas.from_mole(x_mole, "EXH")
    return g.R, float(g.cp_mass(T3)[0]), float(g.gamma(T3)[0]), float(g.a(T3)[0]), g


def main():
    print(f"{'M_inf':>5} {'phi':>4} {'mode':>9} {'p_inf':>7} {'T_inf':>6} {'NPR':>6} "
          f"{'gamma':>6} {'R':>6} {'cp':>7} {'a3':>7} {'M_in':>6} | {'gam_fr':>6} {'cp_fr':>7} {'a3_fr':>7} {'M_in_fr':>7}")
    rows = []
    for M_inf, phi, mode, p3, T3, V3, case in TABLE1:
        gam, mw, x = cea_case(case)
        R = RU / mw
        cp = gam * R / (gam - 1.0)
        a3 = (gam * R * T3) ** 0.5
        p_inf = Q_INF / (0.5 * G_AIR * M_inf**2)
        T_inf = std_atm_T(p_inf)
        R_fr, cp_fr, gam_fr, a_fr, gas = frozen_props(x, T3)
        rows.append(dict(M_inf=M_inf, phi=phi, mode=mode, p_inf=p_inf, T_inf=T_inf, x=x,
                         p3=p3, T3=T3, gam=gam, cp=cp, M_in=V3 / a3, npr=p3 / p_inf,
                         gam_fr=gam_fr, cp_fr=cp_fr, M_in_fr=V3 / a_fr, R_fr=R_fr, MW_fr=gas.MW))
        print(f"{M_inf:5} {phi:4.1f} {mode:>9} {p_inf:7.0f} {T_inf:6.1f} {p3/p_inf:6.1f} "
              f"{gam:6.4f} {R:6.1f} {cp:7.1f} {a3:7.1f} {V3/a3:6.4f} | {gam_fr:6.4f} {cp_fr:7.1f} {a_fr:7.1f} {V3/a_fr:7.4f}")

    print("\n--- problem YAML に貼る形 (plan §4.10 の生産構成: m6_on / m10on / m4_off) ---")
    for r, (name, w) in zip([rows[2], rows[4], rows[1]], [("m6_on", 0.5), ("m10_on", 0.3), ("m4_off", 0.2)]):
        print(f"  - name: {name}\n    weight: {w}\n"
              f"    external: {{M_inf: {r['M_inf']}.0, p_inf: {r['p_inf']:.0f}.0, T_inf: {r['T_inf']:.1f}}}\n"
              f"    inflow: {{M_in: {r['M_in']:.4f}, p_in: {r['p3']}.0, T_in: {r['T3']}.0}}\n"
              f"    gas: {{gamma: {r['gam']:.4f}, cp: {r['cp']:.1f}}}")
    print("\n--- R3 (gas.model: frozen_tp): 凍結組成の擬似種。M_in は凍結音速で、gamma/cp は凍結 γ(T3)/cp(T3) (設計 kernel の CPG 値) ---")
    for r, (name, w) in zip([rows[2], rows[4], rows[1]], [("m6_on", 0.5), ("m10_on", 0.3), ("m4_off", 0.2)]):
        # 種名はクォートする: YAML は NO を真偽値 False に読む (2026-09-13 に踏んだ)
        comp = ", ".join(f"'{k}': {v:.5f}" for k, v in sorted(r["x"].items(), key=lambda kv: -kv[1]) if v >= 1e-5)
        print(f"  - name: {name}\n    weight: {w}\n"
              f"    external: {{M_inf: {r['M_inf']}.0, p_inf: {r['p_inf']:.0f}.0, T_inf: {r['T_inf']:.1f}}}\n"
              f"    inflow: {{M_in: {r['M_in_fr']:.4f}, p_in: {r['p3']}.0, T_in: {r['T3']}.0}}\n"
              f"    gas: {{gamma: {r['gam_fr']:.4f}, cp: {r['cp_fr']:.1f}, composition: {{{comp}}}}}   # MW {r['MW_fr']*1e3:.3f} (CEA {RU/r['R_fr']:.3f})")


if __name__ == "__main__":
    main()
