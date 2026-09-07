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


def cea_props(case: str) -> tuple[float, float]:
    """CEA の tp 出力から (GAMMAs, Mbar) を拾う。"""
    for f in ("tmx1.out", "tmxair.out"):
        txt = (CEA / f).read_text().splitlines()
        for i, L in enumerate(txt):
            if L.strip().startswith("CASE =") and L.split("=")[1].strip() == case:
                blk = txt[i:i + 60]
                mw = gam = None
                for s in blk:
                    if s.lstrip().startswith("M, (1/n)"):
                        mw = float(s.split()[-1])
                    elif s.lstrip().startswith("GAMMAs"):
                        gam = float(s.split()[-1])
                if mw and gam:
                    return gam, mw
    raise KeyError(f"CEA ケース '{case}' が見つからない ({CEA} の .out を作り直すこと)")


def std_atm_T(p: float) -> float:
    """US Standard Atmosphere 1976 の 11-32 km から p [Pa] に対応する T [K]。"""
    if p >= 5474.89:                       # 11-20 km: 等温層
        return 216.65
    return 216.65 * (p / 5474.89) ** (-1.0 / 34.1632)   # 20-32 km: L = +1 K/km


def main():
    print(f"{'M_inf':>5} {'phi':>4} {'mode':>9} {'p_inf':>7} {'T_inf':>6} {'NPR':>6} "
          f"{'gamma':>6} {'R':>6} {'cp':>7} {'a3':>7} {'M_in':>6}")
    rows = []
    for M_inf, phi, mode, p3, T3, V3, case in TABLE1:
        gam, mw = cea_props(case)
        R = RU / mw
        cp = gam * R / (gam - 1.0)
        a3 = (gam * R * T3) ** 0.5
        p_inf = Q_INF / (0.5 * G_AIR * M_inf**2)
        T_inf = std_atm_T(p_inf)
        rows.append(dict(M_inf=M_inf, phi=phi, mode=mode, p_inf=p_inf, T_inf=T_inf,
                         p3=p3, T3=T3, gam=gam, cp=cp, M_in=V3 / a3, npr=p3 / p_inf))
        print(f"{M_inf:5} {phi:4.1f} {mode:>9} {p_inf:7.0f} {T_inf:6.1f} {p3/p_inf:6.1f} "
              f"{gam:6.4f} {R:6.1f} {cp:7.1f} {a3:7.1f} {V3/a3:6.4f}")

    print("\n--- problem YAML に貼る形 (plan §4.10 の生産構成: m6_on / m10on / m4_off) ---")
    for r, (name, w) in zip([rows[2], rows[4], rows[1]], [("m6_on", 0.5), ("m10_on", 0.3), ("m4_off", 0.2)]):
        print(f"  - name: {name}\n    weight: {w}\n"
              f"    external: {{M_inf: {r['M_inf']}.0, p_inf: {r['p_inf']:.0f}.0, T_inf: {r['T_inf']:.1f}}}\n"
              f"    inflow: {{M_in: {r['M_in']:.4f}, p_in: {r['p3']}.0, T_in: {r['T3']}.0}}\n"
              f"    gas: {{gamma: {r['gam']:.4f}, cp: {r['cp']:.1f}}}")


if __name__ == "__main__":
    main()
