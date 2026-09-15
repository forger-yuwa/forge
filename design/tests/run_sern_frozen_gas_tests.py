#!/usr/bin/env python3
"""⑤ SERN R3: 凍結組成 TP 擬似種 (gas/frozen.py) と runner の物性配管の単体テスト。"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.gas.frozen import AIR_MOLE, FrozenGas, ideal_gross_thrust_frozen, isentropic_T, mole_to_mass  # noqa: E402
from forge_design.gas.semiperfect import SPECIES_NASA9, GasCPG  # noqa: E402
from forge_design.geometry.rao_planar import ideal_gross_thrust  # noqa: E402

FAIL = 0
def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1

# --- 純 N2: NIST 値 (R=296.8, cp(300)=1040, cp(2000)=1284 J/kgK, h(298.15)≈0) ---
n2 = FrozenGas({"N2": 1.0}, "N2")
check("N2: R = 296.8", abs(n2.R - 296.80) < 0.1, f"{n2.R:.2f}")
check("N2: cp(300) = 1040", abs(n2.cp_mass(300.0)[0] - 1040.0) < 2.0, f"{n2.cp_mass(300.0)[0]:.1f}")
check("N2: cp(2000) = 1284", abs(n2.cp_mass(2000.0)[0] - 1284.0) < 2.0, f"{n2.cp_mass(2000.0)[0]:.1f}")
check("N2: 絶対 h(298.15) ≈ 0 (基準元素)", abs(n2.h_mass(298.15)[0]) < 50.0, f"{n2.h_mass(298.15)[0]:.1f}")
check("N2: h_sens(298.15) = 0 (thermoHrefTemp 基準)", abs(n2.h_sens(298.15)[0]) < 1e-6)
check("N2: dh/dT = cp (数値微分)", abs((n2.h_mass(1000.5)[0] - n2.h_mass(999.5)[0]) - n2.cp_mass(1000.0)[0]) < 0.05)
check("N2: T ds = dh − v dp (等圧で ds = cp dT/T)", abs((n2.s0_mass(1000.5)[0] - n2.s0_mass(999.5)[0]) - n2.cp_mass(1000.0)[0] / 1000.0) < 1e-4)
check("N2: Tmid 1000 K で h 連続", abs(n2.h_mass(1000.0 + 1e-6)[0] - n2.h_mass(1000.0 - 1e-6)[0]) < 5.0)

# --- 空気 ---
air = FrozenGas.air()
check("空気: R = 287.05 ± 0.1", abs(air.R - 287.05) < 0.1, f"{air.R:.2f}")
check("空気: γ(300) = 1.400 ± 0.002", abs(air.gamma(300.0)[0] - 1.400) < 0.002, f"{air.gamma(300.0)[0]:.4f}")
check("空気: a(288.15) = 340.3 ± 0.3", abs(air.a(288.15)[0] - 340.3) < 0.3, f"{air.a(288.15)[0]:.2f}")
check("mole_to_mass: 空気 Y_N2 = 0.7552", abs(mole_to_mass(AIR_MOLE)["N2"] - 0.7552) < 5e-4)

# --- Ar (定数 cp): CPG と一致するはず ---
ar = FrozenGas({"AR": 1.0}, "AR")
g_ar = float(ar.gamma(500.0)[0])
check("Ar: γ = 5/3", abs(g_ar - 5.0 / 3.0) < 1e-6)
F_fr, Me_fr, Te = ideal_gross_thrust_frozen(ar, 2.0, 1500.0, 1e5, 5e3)
F_cpg, Me_cpg = ideal_gross_thrust(2.0, 0.05, g_ar)
check("Ar: 理想推力 (凍結 TP) = CPG 解析式", abs(F_fr / F_cpg - 1.0) < 1e-6 and abs(Me_fr - Me_cpg) < 1e-6, f"{F_fr:.6f} vs {F_cpg:.6f}, M_e {Me_fr:.4f}/{Me_cpg:.4f}")
T2 = isentropic_T(ar, 1500.0, 1e5, 5e3)
check("Ar: 等エントロピー T2 = T1 (p2/p1)^((γ−1)/γ)", abs(T2 / (1500.0 * 0.05 ** (1 - 1 / g_ar)) - 1.0) < 1e-6)

# --- 排気 (m6_on, CEA tp 平衡組成を凍結): MW = CEA の M 24.430 ---
x_m6 = {"N2": 0.63890, "H2O": 0.32695, "H2": 0.01251, "AR": 0.00767, "OH": 0.00604, "O2": 0.00398, "NO": 0.00206, "H": 0.00125, "O": 0.00037, "CO2": 0.00021, "CO": 0.00005}
exh = FrozenGas.from_mole(x_m6, "EXH")
check("排気 m6_on: MW = 24.430 (CEA)", abs(exh.MW * 1e3 - 24.430) < 0.01, f"{exh.MW*1e3:.3f}")
cp_fr = float(exh.cp_mass(2328.0)[0]); g_fr = float(exh.gamma(2328.0)[0])
check("排気 m6_on: 凍結 cp(2328) < 平衡 Cp 2551 (反応寄与なし)、> 1500", 1500.0 < cp_fr < 2551.0, f"cp_fr {cp_fr:.1f}, γ_fr {g_fr:.4f}")
check("排気 m6_on: 凍結 γ(2328) > 平衡 GAMMAs 1.1828", g_fr > 1.1828, f"{g_fr:.4f}")
F, Me, Te = ideal_gross_thrust_frozen(exh, 1.67, 2328.0, 101027.0, 2851.0)
check("排気 m6_on: 理想膨張 NPR 35.4 → T_e < T_in, M_e > M_in, F_nd > 0", 0 < Te < 2328.0 and Me > 1.67 and F > 0, f"T_e {Te:.0f} K, M_e {Me:.3f}, F/(p_in H) {F:.3f}")
db = exh.pseudo_species_db()
check("擬似種 DB: EXH 1 種、nasa9 9 係数 × 2 区間、MW 一致", set(db) == {"EXH"} and len(db["EXH"]["nasa9_low"]) == 9 and abs(db["EXH"]["MW"] - exh.MW) < 1e-12)
# 擬似種の係数で cp を再評価すると混合 cp と一致 (mixture_pseudo_species の厳密性)
a_lo = np.array(db["EXH"]["nasa9_low"]); T = 800.0
cp_pseudo = (a_lo[0] / T**2 + a_lo[1] / T + a_lo[2] + a_lo[3] * T + a_lo[4] * T**2 + a_lo[5] * T**3 + a_lo[6] * T**4) * 8.314462618 / db["EXH"]["MW"]
check("擬似種 DB: 係数から cp(800) が混合 cp と一致", abs(cp_pseudo / exh.cp_mass(T)[0] - 1.0) < 1e-10)

# --- runner: 生産 YAML の外部動圧 = 71850 Pa (codex C2 の表を潰す) ---
from forge_design.evaluate import runner_sern as R  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
yml = Path(__file__).resolve().parents[2] / "case" / "46.sern_design" / "problem_moo_frozen_tp_cycle3op.yaml"
if yml.exists():
    for op in ("m6_on", "m10_on", "m4_off"):
        p = load_problem(yml); R.select_operating_point(p, op); st = R.gas_states(p)
        q = st["q_inf"]
        check(f"frozen_tp {op}: 外部動圧 ½ρu² = 71850 Pa ± 0.2 %", abs(q / 71850.0 - 1.0) < 2e-3, f"{q:.0f} Pa, ρ∞ {st['ext']['ro']:.4f}, u∞ {st['ext']['u']:.1f}")
        check(f"frozen_tp {op}: 入口 Y=[1,0], 外気 Y=[0,1], 擬似種 [EXH, AIR]", st["exhaust"]["Y"] == [1.0, 0.0] and st["ext"]["Y"] == [0.0, 1.0] and st["species"] == ["EXH", "AIR"])
        F_nd, M_e = R.ideal_thrust(p, st)
        check(f"frozen_tp {op}: 理想推力 F/(p_in H) 有限・出口 M > 入口 M", np.isfinite(F_nd) and M_e > st["exhaust"]["M"], f"F {F_nd:.3f}, M_e {M_e:.3f}")
    p = load_problem(yml); R.select_operating_point(p, "m6_on"); st = R.gas_states(p)
    up = np.array([True, False]); ic = R.region_ic_arrays(up, st, p.gamma)
    gx = R.frozen_gases(p)["exhaust"]
    check("region_ic_arrays (frozen): roe = ρ(h_sens − RT) + ½ρu², roY0 = ρ (排気側)", abs(ic["roe"][0] - st["exhaust"]["ro"] * (gx.e_sens(st["exhaust"]["T"])[0] + 0.5 * st["exhaust"]["u"] ** 2)) < 1e-6 * abs(ic["roe"][0])
          and ic["roY0"][0] == st["exhaust"]["ro"] and ic["roY1"][0] == 0.0 and ic["roY0"][1] == 0.0 and ic["roY1"][1] == st["ext"]["ro"])
    cfg = R._solver_config(p, 100, 10, 0.5, 1000.0)
    check("solverConfig (frozen): thermalMethod 2 + species [\"EXH\", \"AIR\"] (引用符付き) + thermoHrefTemp", "thermalMethod: 2" in cfg and 'species: ["EXH", "AIR"]' in cfg and "thermoHrefTemp: 298.15" in cfg)
    bc = R._bcond_config(p, st)
    check("bcondConfig (frozen): 入口 Y0/Y1 が排気 (1,0)・外気 (0,1)", "Y0: 1, Y1: 0" in bc.split("inlet_nozzle")[1].split("\n")[0] and "Y0: 0, Y1: 1" in bc.split("inlet_ext")[1].split("\n")[0])
    # cpg 側は無変更 (回帰)
    p0 = load_problem(Path(__file__).resolve().parents[2] / "case" / "46.sern_design" / "problem_moo_sst_node_cycle3op.yaml"); R.select_operating_point(p0, "m6_on"); st0 = R.gas_states(p0)
    check("cpg 回帰: gas_states に Y 無し・thermalMethod 0", "Y" not in st0["exhaust"] and "thermalMethod: 0" in R._solver_config(p0, 100, 10, 0.5, 1000.0) and "Y0" not in R._bcond_config(p0, st0))
    check("cpg 回帰: 外部動圧が −15 % ずれる (codex C2 の再現; m6_on)", abs(st0["q_inf"] / 71850.0 - 0.845) < 0.01, f"{st0['q_inf']:.0f} Pa")
else:
    print("skip: problem_moo_frozen_tp_cycle3op.yaml が無い")


# --- 統一 tp_species スキーマ (plan thermophysics-cea-mole-fraction-species §4.5 / §6 SERN, 2026-09-16) ---
if yml.exists():
    import h5py, tempfile
    # full: 輸送種 = 排気 ∪ 外気、トレーサ roXi、IC に roY0..N-1 + roXi、bcond に Xi
    p = load_problem(yml); p.evaluate["tp_species"] = {"mode": "full"}; R.select_operating_point(p, "m6_on"); st = R.gas_states(p)
    check("full m6_on: 輸送種 11 種 (排気 ∪ 外気), tracer, 入口 Y 和 1, Xi 1/0", len(st["species"]) == 11 and st["tracer"] and abs(sum(st["exhaust"]["Y"]) - 1) < 1e-12
          and abs(sum(st["ext"]["Y"]) - 1) < 1e-12 and st["exhaust"]["Xi"] == 1.0 and st["ext"]["Xi"] == 0.0)
    cfg = R._solver_config(p, 100, 10, 0.5, 1000.0)
    check("full m6_on: solverConfig に 11 種 (引用符付き) と tracer: exhaust", f"species: [{', '.join(chr(34)+k+chr(34) for k in st['species'])}]" in cfg and "tracer: exhaust" in cfg)
    import yaml as _yaml
    _cfg = _yaml.safe_load(cfg.replace('"{', '{'))
    check("full m6_on: 生成 config を yaml.safe_load しても NO/N/... が文字列のまま (codex result M1)", [str(x) for x in _cfg["physProp"]["species"]] == st["species"] and all(isinstance(x, str) for x in _cfg["physProp"]["species"]))
    ic = R.region_ic_arrays(np.array([True, False]), st, p.gamma)
    check("full m6_on: IC に roY0..roY10 と roXi (排気側 ρ, 外気側 0)", all(f"roY{i}" in ic for i in range(11)) and ic["roXi"][0] == st["exhaust"]["ro"] and ic["roXi"][1] == 0.0
          and abs(sum(ic[f"roY{i}"][0] for i in range(11)) - st["exhaust"]["ro"]) < 1e-9 * st["exhaust"]["ro"])
    bc = R._bcond_config(p, st)
    check("full m6_on: bcond の排気入口に Y0..Y10 と Xi: 1", "Y10:" in bc.split("inlet_nozzle")[1].split("\n")[0] and "Xi: 1" in bc.split("inlet_nozzle")[1].split("\n")[0] and "Xi: 0" in bc.split("inlet_ext")[1].split("\n")[0])
    # 輸送種ごとの FrozenGas の和 = 排気ガスの熱力学 (lump 線形混合の厳密性)
    g = R.frozen_gases(p); Tq = 2000.0
    e_sum = sum(y * gg.e_sens(Tq)[0] for y, gg in zip(st["exhaust"]["Y"], g["transported"])); e_ref = g["exhaust"].e_sens(Tq)[0]
    check("full m6_on: Σ Y_s e_sens,s(T) = 排気 e_sens(T)", abs(e_sum / e_ref - 1) < 1e-12, f"{e_sum:.6e} vs {e_ref:.6e}")
    # lumped + keep [H2O]: EXH = 1 − Y_H2O、m4_off (H2O 無し) でも配置が同じ
    tp = {"mode": "lumped", "lumps": {"EXH": {"from": "stream", "stream": "inflow"}, "AIR": {"from": "stream", "stream": "external"}}, "keep": ["H2O"]}
    p = load_problem(yml); p.evaluate["tp_species"] = tp; R.select_operating_point(p, "m6_on"); st = R.gas_states(p)
    check("lumped+keep m6_on: [EXH, AIR, H2O], 排気 [0.7589, 0, 0.2411], **tracer 有り** (Y_EXH<1 で流入元ラベルにならない; codex result M8)", st["species"] == ["EXH", "AIR", "H2O"] and abs(st["exhaust"]["Y"][2] - 0.2411091186) < 1e-9 and st["tracer"] and st["exhaust"]["Xi"] == 1.0)
    from forge_design.gas.composition import species_meta as _smeta, _exhaust_fraction_spec
    g2 = R.frozen_gases(p)
    check("lumped+keep: species_meta.exhaust_fraction = tracer Xi", _smeta(g2["layout"])["exhaust_fraction"] == {"kind": "tracer", "array": "Xi", "conserved": "roXi"})
    p_al = load_problem(yml); R.select_operating_point(p_al, "m6_on"); g_al = R.frozen_gases(p_al)
    check("別名 [EXH, AIR]: exhaust_fraction = species Y0 (EXH), tracer 無し", _exhaust_fraction_spec(g_al["layout"]) == {"kind": "species", "array": "Y0", "conserved": "roY0", "species": "EXH"} and not g_al["layout"].tracer)
    p = load_problem(yml); p.evaluate["tp_species"] = tp; R.select_operating_point(p, "m4_off"); st4 = R.gas_states(p)
    check("lumped+keep m4_off: 同じ配置 [EXH, AIR, H2O] で Y_H2O = 0", st4["species"] == ["EXH", "AIR", "H2O"] and st4["exhaust"]["Y"] == [1.0, 0.0, 0.0])
    # restart_by_index / warm_from_same_mesh: 全 roY + roXi を引き継ぐ (codex M4 の既存バグ修正)
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "res.h5"; dst = Path(td) / "sern.h5"
        with h5py.File(src, "w") as f:
            v = f.create_group("VALUE"); v["ro"] = np.array([2.0, 2.0]); v["roUx"] = np.array([1.0, 1.0]); v["roUy"] = np.zeros(2); v["roUz"] = np.zeros(2); v["roe"] = np.array([5.0, 5.0])
            v["roY0"] = np.array([0.8, 0.8]); v["roY1"] = np.array([1.2, 1.2]); v["roXi"] = np.array([0.4, 0.4]); v["wall_dist"] = np.array([9.0, 9.0])
        with h5py.File(dst, "w") as f:
            v = f.create_group("VALUE"); v["ro"] = np.array([1.0, 1.0]); v["roUx"] = np.zeros(2); v["roUy"] = np.zeros(2); v["roUz"] = np.zeros(2); v["roe"] = np.ones(2)
            v["roY0"] = np.array([1.0, 1.0]); v["roY1"] = np.zeros(2); v["wall_dist"] = np.array([1.0, 1.0])
        R.restart_by_index(src, dst)
        with h5py.File(dst) as f:
            ok = (f["VALUE/roY0"][0] == 0.8 and f["VALUE/roY1"][0] == 1.2 and "roXi" in f["VALUE"] and abs(f["VALUE/roXi"][0] - 0.4) < 1e-7
                  and f["VALUE/wall_dist"][0] == 1.0 and abs((f["VALUE/roY0"][0] + f["VALUE/roY1"][0]) / f["VALUE/ro"][0] - 1.0) < 1e-12)
        check("restart_by_index: roY0/roY1/roXi を引き継ぎ ΣρY = ρ、wall_dist は触らない", bool(ok))
        from forge_design.evaluate import runner_sern3d as R3
        run3 = Path(td) / "run3"; run3.mkdir(); (run3 / R3.MESH).write_bytes(Path(dst).read_bytes())
        with h5py.File(run3 / R3.MESH, "r+") as f:
            f["VALUE/roY0"][:] = 1.0; f["VALUE/roY1"][:] = 0.0
        R3.warm_from_same_mesh(run3, src)
        with h5py.File(run3 / R3.MESH) as f:
            check("warm_from_same_mesh (3D): roY/roXi を引き継ぐ", f["VALUE/roY1"][0] == 1.2 and "roXi" in f["VALUE"])

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
