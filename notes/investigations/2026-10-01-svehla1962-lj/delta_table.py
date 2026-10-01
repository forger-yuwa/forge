#!/usr/bin/env python3
"""#14-L1 (vi): LJ の既定変更 (legacy_v1 → [gri30, svehla1962]) による D_ij と kinetic μ・λ の Δ 表 (合否なし・記録)。

  python3 notes/investigations/2026-10-01-svehla1962-lj/delta_table.py > notes/investigations/2026-10-01-svehla1962-lj/delta_table.md

- LJ は共通データ (solver_density_cuda/data/species/forge_species_v1.yaml) の LJ_sets を、ソルバと同じ規則 (ljSource の先頭から探す)
  で独立参照 tests/unit/transport_reference.py が解決したもの。
- kinetic μ・λ: physProp.transport で全実種を kinetic にしたときの混合物 (独立参照 Reference.state; 実装とは試験で ≤1e-12)。
- D_ij: 輸送種 (lump は 1 種) どうしの Chapman-Enskog 二元拡散係数 (cuda_forge/thermo_d.cuh thermo_Dbinary と同じ式を double で;
  σ_ij = (σ_i+σ_j)/2, ε_ij = √(ε_i ε_j), Ω(1,1)* は Neufeld)。lump の LJ はソルバの合成と同じ質量分率平均 (暫定規約, plan #7)。
  D_ij ∝ 1/P なので相対差は圧力に依らない (101325 Pa で評価)。
- 組成: case/44 va3 (run_0509/run_0532 型: MIXDRY lump + H2O、入口 Y_H2O 0.03769536)、case/16 (N2 + H2O、Y_H2O 0.01095)、
  参考に SERN 代表 (m6_on のノズル入口 11 種; plan の指定外)。
"""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "solver_density_cuda", "tests", "unit"))
import transport_reference as TR   # noqa: E402

OLD, NEW = ("legacy_v1",), ("gri30", "svehla1962")
TS = (300.0, 600.0, 1000.0, 2000.0)
MIXDRY = {"name": "MIXDRY", "lump": {"N2": 0.7088730183520131, "O2": 0.23037573447245457, "AR": 0.008503874694843205,
                                     "CO2": 0.05224737248068922}, "basis": "mole"}
SERN = ["N2", "H2O", "H2", "AR", "OH", "O2", "NO", "H", "O", "CO2", "CO"]
SERN_Y = [0.732638, 0.241109, 0.00103232, 0.0125424, 0.00420498, 0.00521324, 0.00253027, 5.15746e-05, 0.000242324,
          0.000378318, 5.73292e-05]
CASES = [("case/44 va3 (MIXDRY lump + H2O)", [MIXDRY, "H2O"], [1.0 - 0.03769536, 0.03769536]),
         ("case/16 (N2 + H2O)", ["N2", "H2O"], [1.0 - 0.01095, 0.01095]),
         ("参考: SERN 代表 m6_on ノズル入口 (11 種)", SERN, SERN_Y)]


def omega11(t):
    t = min(max(t, 0.3), 100.0)
    return 1.06036 * t**-0.15610 + 0.19300 * math.exp(-0.47635 * t) + 1.03587 * math.exp(-1.52996 * t) + 1.76474 * math.exp(-3.89411 * t)


def transport_lj(sp, lj_source):
    """輸送種 (lump は質量分率平均) の (MW, σ, ε)。"""
    def real(n):
        cid = TR.ALIAS[n]
        s, e, _ = TR.builtin_lj(cid, lj_source)
        return TR.BUILTIN[cid]["MW"], s, e
    if isinstance(sp, str):
        return real(sp)
    mem = list(sp["lump"].items())
    tot = sum(v for _, v in mem)
    x = [v / tot for _, v in mem]
    r = [real(k) for k, _ in mem]
    MW = sum(xi * m for xi, (m, _, _) in zip(x, r))
    Y = [xi * m / MW for xi, (m, _, _) in zip(x, r)]
    return MW, sum(y * s for y, (_, s, _) in zip(Y, r)), sum(y * e for y, (_, _, e) in zip(Y, r))


def dbin(a, b, T, P=101325.0):
    Mi, si, ei = a
    Mj, sj, ej = b
    sig, eps = 0.5 * (si + sj), math.sqrt(ei * ej)
    return 1.8583e-3 * math.sqrt(T**3 * (1 / (Mi * 1e3) + 1 / (Mj * 1e3))) / ((P / 101325.0) * sig * sig * omega11(T / eps)) * 1e-4


def name(sp):
    return sp if isinstance(sp, str) else sp["name"]


def main():
    print("# LJ 既定変更の Δ 表 (#14-L1 (vi); 合否なし・記録)\n")
    print("生成: `python3 notes/investigations/2026-10-01-svehla1962-lj/delta_table.py` (式と入力は同スクリプトの docstring)。"
          "Δ = 新既定 [gri30, svehla1962] / legacy_v1 − 1。\n")
    for title, species, Y in CASES:
        print(f"## {title}\n")
        # 解決した LJ (実種)
        reals = []
        for s in species:
            reals += [s] if isinstance(s, str) else list(s["lump"])
        print("| 実種 | legacy_v1 σ/ε | 新既定 σ/ε (集合) |")
        print("| --- | --- | --- |")
        for r in reals:
            cid = TR.ALIAS[r]
            so, eo, _ = TR.builtin_lj(cid, OLD)
            sn, en, setn = TR.builtin_lj(cid, NEW)
            mark = "" if (so, eo) == (sn, en) else " **変化**"
            print(f"| {r} | {so}/{eo} | {sn}/{en} ({setn}){mark} |")
        print()
        tr = {r: "kinetic" for r in reals}
        ro, rn = TR.Reference(species, tr, lj_source=OLD), TR.Reference(species, tr, lj_source=NEW)
        Wm = [y / m for y, m in zip(Y, ro.mw)]
        X = [w / sum(Wm) for w in Wm]
        lo = [transport_lj(s, OLD) for s in species]
        ln = [transport_lj(s, NEW) for s in species]
        pairs = [(i, j) for i in range(len(species)) for j in range(i + 1, len(species))]
        print("| T [K] | Δμ (kinetic 混合) | Δλ (kinetic 混合) | 絶対値最大の ΔD_ij (組) | ΔD_ij (" + ", ".join(
            f"{name(species[i])}-{name(species[j])}" for i, j in pairs[:3]) + (" …" if len(pairs) > 3 else "") + ") |")
        print("| --- | --- | --- | --- | --- |")
        for T in TS:
            so, sn = ro.state(T, X), rn.state(T, X)
            dmu, dla = sn["mu"] / so["mu"] - 1, sn["lam"] / so["lam"] - 1
            dd = [(dbin(ln[i], ln[j], T) / dbin(lo[i], lo[j], T) - 1, i, j) for i, j in pairs]
            w = max(dd, key=lambda v: abs(v[0]))
            print(f"| {T:.0f} | {dmu:+.3e} | {dla:+.3e} | {w[0]:+.3e} ({name(species[w[1]])}-{name(species[w[2]])}) | "
                  + ", ".join(f"{d:+.3e}" for d, _, _ in dd[:3]) + " |")
        print(f"\n組成 (輸送種のモル分率): " + ", ".join(f"{name(s)} {x:.6g}" for s, x in zip(species, X)) + "\n")


if __name__ == "__main__":
    main()
