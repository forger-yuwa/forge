#!/usr/bin/env python3
"""semi-perfect (NASA-9, frozen) ガスモデルの単体検証。

1. 一定 cp 擬似種で CPG に機械精度退化 (ν, A/A*, T(M))
2. NASA-9 の cp が文献値 (N2/O2/CO2/H2O) に一致
3. 混合擬似種 (mixture_pseudo_species) が多成分混合と一致
4. MOC カーネルの pm_nu/pm_mach がガスモデルで往復
5. 有効 γ 一定の CPG が燃焼ガスの A/A* を系統的に誤る (semi-perfect の存在意義)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forge_design.gas.semiperfect as spm  # noqa: E402
from forge_design.gas import GasCPG, GasSemiPerfect  # noqa: E402
from forge_design.geometry.moc_kernel import pm_mach, pm_mach_vec, pm_nu  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# 1. CPG 退化
spm.SPECIES_NASA9["AIRX"] = dict(MW=0.0289647, low=[0, 0, 3.5, 0, 0, 0, 0, -1043.1, 3.0],
                                 high=[0, 0, 3.5, 0, 0, 0, 0, -1043.1, 3.0])
g = GasSemiPerfect({"AIRX": 1.0}, Tt=1000.0, n_tab=20000)
cpg = GasCPG(gamma=1.4, cp=float(g.cp_mass(500.0)[0]))
Ms = np.array([1.5, 2.0, 3.0, 4.0, 5.0])
check(f"CPG 退化: ν 一致 (max|Δ| {np.max(np.abs(g.nu(Ms)-cpg.nu(Ms))):.1e})",
      float(np.max(np.abs(g.nu(Ms) - cpg.nu(Ms)))) < 1e-6)
check(f"CPG 退化: A/A* 一致 (max rel {np.max(np.abs(g.area_ratio(Ms)/cpg.area_ratio(Ms)-1)):.1e})",
      float(np.max(np.abs(g.area_ratio(Ms) / cpg.area_ratio(Ms) - 1))) < 1e-6)
check(f"CPG 退化: γ* = {g.gamma_star:.6f}", abs(g.gamma_star - 1.4) < 1e-6)
check(f"CPG 退化: T* = {g.T_star:.3f} (理論 833.333)", abs(g.T_star - 1000 / 1.2) < 0.01)

# 2. NASA-9 文献値
for sp, T, ref, tol in (("N2", 300, 1040, 2), ("CO2", 1000, 1234, 2), ("H2O", 1000, 2290, 5), ("O2", 300, 918, 2)):
    v = float(GasSemiPerfect({sp: 1.0}, Tt=1500.0).cp_mass(T)[0])
    check(f"NASA-9 cp {sp}@{T}K = {v:.1f} (文献 ≈{ref})", abs(v - ref) < tol)

# 3. 混合擬似種
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
mix = spm.mixture_pseudo_species(Y, "MIXT")
spm.SPECIES_NASA9["MIXT"] = dict(MW=mix["MIXT"]["MW"], low=mix["MIXT"]["nasa9_low"],
                                 high=mix["MIXT"]["nasa9_high"])
gm = GasSemiPerfect(Y, Tt=1000.0); g1 = GasSemiPerfect({"MIXT": 1.0}, Tt=1000.0)
Ts = np.array([300., 900., 1000., 1500.])
check(f"混合擬似種: cp 一致 (max|Δ| {np.max(np.abs(gm.cp_mass(Ts)-g1.cp_mass(Ts))):.1e})",
      float(np.max(np.abs(gm.cp_mass(Ts) - g1.cp_mass(Ts)))) < 1e-8)
check(f"混合擬似種: h 一致 (max|Δ| {np.max(np.abs(gm.h_mass(Ts)-g1.h_mass(Ts))):.1e})",
      float(np.max(np.abs(gm.h_mass(Ts) - g1.h_mass(Ts)))) < 1e-4)
check(f"混合擬似種: R 一致 ({gm.R:.4f} vs {g1.R:.4f})", abs(gm.R - g1.R) < 1e-9)

# 4. MOC カーネル往復
for M in (1.2, 2.5, 4.0):
    nu = float(pm_nu(M, gm)); Mb = pm_mach(nu, gm)
    check(f"pm_nu/pm_mach ガス往復 M={M}: {Mb:.6f}", abs(Mb - M) < 1e-4)
Mv = pm_mach_vec(pm_nu(Ms, gm), gm)
check(f"pm_mach_vec ガス往復 (max|Δ| {np.max(np.abs(Mv-Ms)):.1e})", float(np.max(np.abs(Mv - Ms))) < 1e-4)

# 5. 有効 γ 一定 CPG の系統誤差
ar_sp = float(gm.area_ratio(4.0))
ar_star = float(GasCPG(gamma=gm.gamma_star).area_ratio(4.0))
ar_m4 = float(GasCPG(gamma=float(gm.gamma(gm.T_of_M(4.0))[0])).area_ratio(4.0))
print(f"info 燃焼ガス A/A*(M4): semi-perfect {ar_sp:.3f} / CPG(γ*={gm.gamma_star:.3f}) {ar_star:.3f} / CPG(γ_M4) {ar_m4:.3f}")
check("semi-perfect の A/A* は γ* と γ_M4 の CPG の間にある (単一 γ では表現不能)",
      min(ar_star, ar_m4) < ar_sp < max(ar_star, ar_m4))
check(f"CPG(γ*) は出口半径を >5% 誤る ({100*(np.sqrt(ar_star/ar_sp)-1):+.1f}%)",
      abs(np.sqrt(ar_star / ar_sp) - 1) > 0.05)


# --- split 擬似種 (MIXDRY + H2O) の等価性 (2026-08-17, tp-split-h2o plan) ---
from forge_design.gas.semiperfect import RU, SPECIES_NASA9, mixture_pseudo_species_split, _cp_R, _h_RT
Ysp = {"N2": 0.75277, "CO2": 0.17426, "O2": 0.02286, "H2O": 0.05011}
dbs, Ysplit, order = mixture_pseudo_species_split(Ysp)
gfull = GasSemiPerfect(Ysp, Tt=1000.0)
def _cp_db(e, T):
    a = np.asarray(e["nasa9_low"] if T < e["Tmid"] else e["nasa9_high"]); return float(_cp_R(a, T)) * RU / e["MW"]
def _h_db(e, T):
    a = np.asarray(e["nasa9_low"] if T < e["Tmid"] else e["nasa9_high"]); return float(_h_RT(a, T)) * RU * T / e["MW"]
errcp = max(abs(sum(Ysplit[k] * _cp_db(dbs[k], T) for k in order) - float(gfull.cp_mass(T)[0])) for T in (250., 400., 800., 1200.))
errh = max(abs(sum(Ysplit[k] * _h_db(dbs[k], T) for k in order) - float(gfull.h_mass(T)[0])) for T in (250., 400., 800., 1200.))
check(f"split: order={order}, Σ Y={sum(Ysplit.values()):.6f}", order == ["MIXDRY", "H2O"] and abs(sum(Ysplit.values()) - 1) < 1e-12)
check(f"split: 混合 cp が全種混合と一致 (max|Δ| {errcp:.1e})", errcp < 1e-9)
check(f"split: 混合 h が全種混合と一致 (max|Δ| {errh:.1e})", errh < 1e-6)
check("split: H2O エントリは内蔵 NASA-9 と同係数", dbs["H2O"]["nasa9_low"] == [float(v) for v in SPECIES_NASA9["H2O"]["low"]])


# --- 組成の単一ソース (plan thermophysics-cea-mole-fraction-species §6 単体 (a)–(e), 2026-09-16) ---
from forge_design.gas import composition as C  # noqa: E402


def _chk(name, cond, info=""):
    check(name + (f"  [{info}]" if info else ""), cond)

import tempfile, pathlib  # noqa: E402
db = C.ResolvedSpeciesDB.builtin()
# (a) mole↔mass 往復
x_va3 = {"H2O": 6.09135e-2, "N2": 6.64860e-1, "O2": 2.16072e-1, "AR": 7.97588e-3, "CO2": 4.90034e-2}
Y = C.mole_to_mass(x_va3, db); X2 = C.mass_to_mole(Y, db); x_norm, tot = C.normalize_fractions(x_va3)
_chk("composition (a): mole→mass→mole 往復 rtol 1e-12", all(abs(X2[k] / x_norm[k] - 1.0) < 1e-12 for k in x_norm), f"Σ入力 {tot:.6f}")
# (b) va3 のモル分率 (Σ 0.998825) → Y_H2O = 0.03769539643
_chk("composition (b): va3 X (Σ 0.998825) → Y_H2O = 0.03769539643 (rtol 1e-9)", abs(Y["H2O"] / 0.03769539643 - 1.0) < 1e-9, f"{Y['H2O']:.11f}")
Yb, Xb, totb = C.composition_to_mass(x_va3, "mole", db)
_chk("composition (b'): composition_to_mass(mole) = mole_to_mass (rtol 1e-12), 入力総和を返す", all(abs(Yb[k] / Y[k] - 1) < 1e-12 for k in Y) and abs(totb - 0.998825) < 1e-6)
# (c) full と lumped の混合 cp(T)/h(T) が 200–3000 K で一致 (rtol 1e-12 + atol cp 1e-9 / h 1e-6)
tp_l = C.parse_tp_species({"tp_species": "split_h2o"}); tp_f = C.parse_tp_species({"tp_species": "full"})
Ll = C.resolve_species_layout(tp_l, {"inflow": Y}, db, "H2O", condensation=True)
Lf = C.resolve_species_layout(tp_f, {"inflow": Y}, db, "H2O", condensation=True)
T = np.linspace(200.0, 3000.0, 57)
def mix(L, fn):
    Yt = dict(zip(L.species, L.Y_transport("inflow")))
    out = np.zeros_like(T)
    for sname, y in Yt.items():
        e = L.entries[sname]; tmp = C.ResolvedSpeciesDB({sname: e})
        out += y * getattr(tmp, fn)({sname: 1.0}, T)
    return out
cp_l, cp_f = mix(Ll, "cp_mass"), mix(Lf, "cp_mass"); h_l, h_f = mix(Ll, "h_mass"), mix(Lf, "h_mass")
ok_cp = np.all(np.abs(cp_l - cp_f) <= 1e-12 * np.abs(cp_f) + 1e-9); ok_h = np.all(np.abs(h_l - h_f) <= 1e-12 * np.abs(h_f) + 1e-6)
_chk("composition (c): lumped [MIXDRY,H2O] と full の cp(T) 一致 (200–3000 K)", bool(ok_cp), f"max|Δcp| {np.max(np.abs(cp_l-cp_f)):.2e}")
_chk("composition (c): lumped と full の h(T) 一致 (h(T_ref)=0 近傍含む)", bool(ok_h), f"max|Δh| {np.max(np.abs(h_l-h_f)):.2e}")
_chk("composition (c'): 旧 split_h2o の順序 [MIXDRY, H2O]・凝縮 index 1・入口ベクトル和 1", Ll.species == ["MIXDRY", "H2O"] and Ll.cond_index == 1 and abs(sum(Ll.Y_transport("inflow")) - 1) < 1e-12)
_chk("composition (c''): full の順序 = 組成順, 凝縮 index = H2O の位置", Lf.species == ["H2O", "N2", "O2", "AR", "CO2"] and Lf.cond_index == 0)
# (d) species_db.yaml: 解析後の値が一致、由来コメントあり
txt = C.species_db_yaml(Ll); import yaml as _y; parsed = _y.safe_load(txt)
_chk("composition (d): species_db.yaml を再解析すると MW/係数が一致し、lumped/source コメントがある",
      abs(parsed["MIXDRY"]["MW"] - Ll.entries["MIXDRY"].MW) < 1e-15 and parsed["H2O"]["nasa9_low"] == list(db["H2O"].low)
      and "# lumped: {" in txt and "# source: CEA thermo.inp" in txt)
meta = C.species_meta(Ll)
_chk("composition (d'): species_meta の展開行列 (MIXDRY→実種) の和が 1、streams.Y_transport が一致", abs(sum(meta["expansion"]["MIXDRY"].values()) - 1) < 1e-12 and meta["streams"]["inflow"]["Y_transport"] == Ll.Y_transport("inflow"))
# (e) 拒否条件
def rejects(name, fn):
    try:
        fn(); _chk(name, False, "例外が出ない")
    except (ValueError, KeyError) as ex:
        _chk(name, True, str(ex)[:60])
rejects("reject: 凝縮 ON なのに凝縮種が keep に無い (pseudo)", lambda: C.resolve_species_layout(C.parse_tp_species({"tp_species": "pseudo"}), {"inflow": Y}, db, "H2O", condensation=True))
rejects("reject: 擬似種名が実種名と衝突", lambda: C.resolve_species_layout(C.parse_tp_species({"tp_species": {"mode": "lumped", "lumps": {"N2": {"from": "composition"}}, "keep": ["H2O"]}}), {"inflow": Y}, db, "H2O"))
rejects("reject: 空 lump", lambda: C.resolve_species_layout(C.parse_tp_species({"tp_species": {"mode": "lumped", "lumps": {"MIXDRY": {"from": "composition", "exclude": ["N2", "O2", "AR", "CO2"]}}, "keep": ["H2O"]}}), {"inflow": Y}, db, "H2O"))
rejects("reject: 未配分 (exclude したのに keep に無い)", lambda: C.resolve_species_layout(C.parse_tp_species({"tp_species": {"mode": "lumped", "lumps": {"MIXDRY": {"from": "composition", "exclude": ["CO2"]}}, "keep": ["H2O"]}}), {"inflow": Y}, db, "H2O"))
rejects("reject: 二重配分 (2 つの composition lump)", lambda: C.resolve_species_layout(C.parse_tp_species({"tp_species": {"mode": "lumped", "lumps": {"A": {"from": "composition"}, "B": {"from": "composition"}}, "keep": ["H2O"]}}), {"inflow": Y}, db, "H2O"))
rejects("reject: mapping と tp_keep_species の併用", lambda: C.parse_tp_species({"tp_species": {"mode": "full"}, "tp_keep_species": "H2O"}))
rejects("reject: 未知の composition_basis", lambda: C.composition_to_mass(x_va3, "volume", db))
rejects("reject: 負の分率", lambda: C.mole_to_mass({"N2": 1.0, "O2": -0.1}, db))
rejects("reject: DB に無い種", lambda: C.mole_to_mass({"N2": 1.0, "XE": 0.1}, db))
rejects("reject: 真偽値キー (NO 無引用)", lambda: C.mole_to_mass({"N2": 1.0, False: 0.1}, db))
# 外部 DB (CEA 直読み形式) の上書きと温度区切りの検査
with tempfile.TemporaryDirectory() as td:
    f = pathlib.Path(td) / "db.yaml"
    f.write_text('"XE":\n  MW: 0.131293\n  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n  nasa9_low: [0,0,2.5,0,0,0,0,-745.375,6.16]\n  nasa9_high: [0,0,2.5,0,0,0,0,-745.375,6.16]\n  atoms: {XE: 1}\n'
                 '"HE":\n  MW: 0.004002602\n  Tlo: 200.0\n  Tmid: 6000.0\n  Thi: 20000.0\n  nasa9_low: [0,0,2.5,0,0,0,0,-745.375,0.93]\n  nasa9_high: [0,0,2.5,0,0,0,0,-745.375,0.93]\n')
    db2 = C.ResolvedSpeciesDB.from_file(f)
    _chk("外部 DB: 内蔵に上書き追加 (XE, HE) され N2 は内蔵のまま", "XE" in db2 and "HE" in db2 and db2["N2"].source == C.BUILTIN_SOURCE and db2["XE"].atoms == {"XE": 1.0})
    rejects("reject: 温度区切りが標準でない種 (HE 200/6000/20000) を lump に畳む", lambda: C.lump_entry("MIXX", {"N2": 0.5, "HE": 0.5}, db2))
    Lx = C.resolve_species_layout(C.parse_tp_species({"tp_species": "full"}), {"inflow": {"N2": 0.5, "HE": 0.5}}, db2)
    _chk("外部 DB: full なら非標準区切りの種も輸送種にできる", Lx.species == ["N2", "HE"] and Lx.entries["HE"].Thi == 20000.0)
# SERN 型: 流れ lump + keep の質量配分 (codex M3): 純排気でも Y_EXH = 1 − Y_H2O
x_m6 = {"N2": 0.63890, "H2O": 0.32695, "H2": 0.01251, "AR": 0.00767, "OH": 0.00604, "O2": 0.00398, "NO": 0.00206, "H": 0.00125, "O": 0.00037, "CO2": 0.00021, "CO": 0.00005}
Y6 = C.mole_to_mass(x_m6, db); Yair = C.mole_to_mass({"N2": 0.78084, "O2": 0.20946, "AR": 0.00934}, db)
Ls = C.resolve_species_layout(C.parse_tp_species({"tp_species": {"mode": "lumped", "lumps": {"EXH": {"from": "stream", "stream": "inflow"}, "AIR": {"from": "stream", "stream": "external"}}, "keep": ["H2O"]}}),
                              {"inflow": Y6, "external": Yair}, db)
_chk("SERN lumped+keep: 輸送種 [EXH, AIR, H2O], 排気入口 [1−Y_H2O, 0, Y_H2O] (0.2411)", Ls.species == ["EXH", "AIR", "H2O"] and abs(Ls.Y_transport("inflow")[2] - 0.2411091186) < 1e-9 and abs(Ls.Y_transport("inflow")[0] - (1 - 0.2411091186)) < 1e-9 and Ls.Y_transport("external") == [0.0, 1.0, 0.0])
_chk("SERN lumped+keep: EXH の lump 内組成に H2O が無い (二重計上なし)", "H2O" not in Ls.lumps["EXH"]["members"])
La = C.resolve_species_layout(C.parse_tp_species({"tp_species": ["EXH", "AIR"]}), {"inflow": Y6, "external": Yair}, db)
_chk("SERN 別名 [EXH, AIR]: 輸送種 [EXH, AIR], 入口 [1,0]/[0,1], トレーサ無し", La.species == ["EXH", "AIR"] and La.Y_transport("inflow") == [1.0, 0.0] and La.Y_transport("external") == [0.0, 1.0] and not La.tracer)
Lf2 = C.resolve_species_layout(C.parse_tp_species({"tp_species": "full"}), {"inflow": Y6, "external": Yair}, db)
_chk("SERN full: 輸送種 = 排気 ∪ 外気 (11 種), トレーサ有り, 各流れの入口ベクトル和 1", Lf2.n == 11 and Lf2.tracer and all(abs(sum(Lf2.Y_transport(s)) - 1) < 1e-12 for s in ("inflow", "external")))
# 元素質量分率の診断
Z = C.element_mass_fractions(Y6, db)
_chk("元素質量分率: Σ_e Z_e = 1", abs(sum(Z.values()) - 1.0) < 1e-6, f"{sum(Z.values()):.8f}")

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
