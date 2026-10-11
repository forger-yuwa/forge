#!/usr/bin/env python3
"""③ベル 2 モードの共通評価器 — TP のガスと 1D 理想値・runner の配線 (plan tooling-nozzle-bell-two-mode §5.1 #2a)。

合格条件 (plan §5.1 #2a): 1D の凍結 TP の出口温度・ṁ_1D・出口速度・出口圧力・F_1D が、同じ条件の CEA 凍結計算
(燃焼室で凍結、面積比 50) と相対 0.1 % 以内。反応物の状態 (LOX/LH2 の温度) と凍結位置 (燃焼室) を固定して記録する。
加えて: MOC 側 (GasSemiPerfect)・CFD 側 (lump EXH)・IC/理想値 (FrozenGas) が同じ組成・DB・基準から作られること、
問題 YAML の単一ソース (gas.exhaust) の検査、TP の solverConfig/bcondConfig が推奨レシピに沿うこと、CPG 経路の config が
変更前と文字単位で同じこと、TP の壁の平滑化 (P 不変・TP の内部エネルギー)。

CEA2 (FCEA2) が見つかれば排気ファイルを作り直して照合し (保存ファイルとの一致も確かめる)、無ければ case/68 の保存ファイルの
CEA 参照表で照合する。実行: design/.venv-opt の Python で `python design/tests/run_bell_tp_tests.py`。
"""
import hashlib
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.gas.composition import ResolvedSpeciesDB  # noqa: E402
from forge_design.gas.rocket_exhaust import (DEFAULT_PRODUCTS, ChamberSpec, RocketExhaust, cea_dir,  # noqa: E402
                                             cea_numbers, cea_reference_at, load_exhaust_file, resolve_exhaust)
from forge_design.gas.semiperfect import RU, GasSemiPerfect  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CASE68 = ROOT / "case" / "68.bell_two_mode"
EXH_FILE = CASE68 / "exhaust_lox_lh2_pc100_of6.yaml"
PROBLEM = CASE68 / "problem_bell_tp_top.yaml"
TOL_1D = 1.0e-3          # 合格条件: 相対 0.1 %

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:
        return str(e) or type(e).__name__
    return None


# ---------------------------------------------------------------- CEA の数値欄の読み方
check("cea_numbers: RHO の '4.6118 0' と '2.1329-2'", cea_numbers(" 4.6118 0 2.8605 0 2.1329-2") == [4.6118, 2.8605, 0.021329])
check("cea_numbers: trace の '2.854 -3' と '1.883-12'", cea_numbers("2.854 -3 1.811 -3 1.883-12") == [0.002854, 0.001811, 1.883e-12])
check("cea_numbers: 負の値・普通の小数", cea_numbers("-986.31 -2168.77 17.6132") == [-986.31, -2168.77, 17.6132])

# ---------------------------------------------------------------- 排気ファイル (CEA2 で作り直せればそれ、無ければ保存物)
stored = load_exhaust_file(EXH_FILE)
try:
    cea_dir()
    live = resolve_exhaust(ChamberSpec(), supar=stored["cea_input"]["supar"])
except FileNotFoundError as e:
    live = None
    print(f"note: CEA2 が見つからないので保存ファイル {EXH_FILE.name} の CEA 参照表で照合する ({e})")
data = live if live is not None else stored
src = "CEA2 を今回実行" if live is not None else "保存ファイル"
print(f"note: 照合の参照 = {src}; {data['cea_provenance']['program']} / thermo.inp {data['cea_provenance']['thermo_inp_update']}")
if live is not None:
    check("保存ファイルが CEA2 の作り直しと一致 (組成・Tc・Pc・参照表)",
          live["chamber"] == stored["chamber"] and live["cea_frozen_reference"] == stored["cea_frozen_reference"]
          and live["cea_input"] == stored["cea_input"])

# 反応物の状態と凍結位置を固定して記録する (plan §4.1 / §5.1 #2a)
ci = data["cea_input"]
check("反応物: 燃料 H2(L) 20.27 K", ci["fuel"] == "H2(L)" and ci["fuel_T"] == 20.27, f"{ci['fuel']} {ci['fuel_T']} K")
check("反応物: 酸化剤 O2(L) 90.17 K", ci["oxid"] == "O2(L)" and ci["oxid_T"] == 90.17, f"{ci['oxid']} {ci['oxid_T']} K")
check("燃焼室圧 10 MPa・O/F 6 (質量比)", ci["pc"] == 1.0e7 and ci["of_ratio"] == 6.0 and ci["of_ratio_basis"] == "mass")
check("凍結位置 = 燃焼室 (rocket frozen nfz=1)", ci["freeze"].startswith("chamber") and "nfz=1" in ci["freeze"], ci["freeze"])
check("生成物の候補 = forge の熱物性・輸送データにある 6 種", tuple(ci["products"]) == DEFAULT_PRODUCTS, str(ci["products"]))
ch = data["chamber"]
check("燃焼室温度 3523.96 K (CEA2)", ch["Tc_K"] == 3523.96, f"{ch['Tc_K']}")
om = data["omitted_species"]
check("除いた生成物 (HO2・H2O2・O3) のモル分率の和 < 1e-4、全温の差 < 0.5 K",
      om["sum"] < 1e-4 and abs(om["dTc_K"]) < 0.5, f"sum {om['sum']:.3e}, dTc {om['dTc_K']:+.2f} K, {sorted(om['mole_fractions'])}")

# ---------------------------------------------------------------- 1D 理想値 vs CEA 凍結 (合格条件: ε 50 で相対 0.1 %)
db = ResolvedSpeciesDB.builtin()
ex = RocketExhaust(data, db, 298.15)
for eps in ci["supar"]:
    i1 = ex.ideal_1d(float(eps))
    c = cea_reference_at(data, float(eps))
    tag = "合格条件" if float(eps) == 50.0 else "追加"
    for k, lab in (("T_e", "出口温度"), ("mdot_per_At", "ṁ_1D/A* (= Pc/c*)"), ("u_e", "出口速度 (= Isp)"),
                   ("p_e", "出口圧力"), ("F_vac_per_At", "F_1D,vac/A* (= ṁ·Ivac)")):
        rel = i1[k] / c[k] - 1.0
        check(f"[{tag}] ε {eps:g}: {lab} が CEA 凍結と相対 0.1 % 以内", abs(rel) <= TOL_1D,
              f"1D {i1[k]:.6g} / CEA {c[k]:.6g} → {rel:+.2e}")
i50 = ex.ideal_1d(50.0)
c50 = cea_reference_at(data, 50.0)
check("スロート: T*・p* が CEA 凍結と相対 1e-4 以内",
      abs(i50["T_star"] / c50["T_star"] - 1) < 1e-4 and abs(i50["p_star"] / c50["p_star"] - 1) < 1e-4,
      f"T* {i50['T_star']:.2f}/{c50['T_star']:.2f}, p* {i50['p_star']:.5g}/{c50['p_star']:.5g}")
check("出口 M が CEA 凍結と相対 1e-3 以内", abs(i50["M_e"] / c50["M_e"] - 1) < 1e-3, f"{i50['M_e']:.4f}/{c50['M_e']:.3f}")
check("1D: u* = a* (チョーク条件の解)", abs(i50["u_star"] / i50["a_star"] - 1) < 1e-10)
check("1D: p_a を与えると F = F_vac − p_a ε A*",
      abs(ex.ideal_1d(50.0, p_a=1.0e4)["F_per_At"] - (i50["F_vac_per_At"] - 1.0e4 * 50.0)) < 1e-6 * i50["F_vac_per_At"])

# ---------------------------------------------------------------- 契約: MOC・CFD・IC/理想値が同じ組成・DB・基準
fg = ex.frozen_gas()
mg = ex.moc_gas()
L = ex.layout()
check("契約: MOC (GasSemiPerfect)・FrozenGas・ex の質量分率が同一", fg.Y == ex.Y and mg.Y == ex.Y)
check("契約: CFD の輸送種 = 単一の lump EXH、lump の質量分率が同一",
      L.species == ["EXH"] and max(abs(L.lumps["EXH"]["mass"][k] - ex.Y[k]) for k in ex.Y) < 1e-15)
check("契約: 気体定数 R が一致 (MOC / FrozenGas / lump の MW)",
      abs(mg.R / fg.R - 1) < 1e-13 and abs(RU / L.entries["EXH"].MW / fg.R - 1) < 1e-13, f"R {fg.R:.6f}")
Tg = np.array([300.0, 1000.0, 2000.0, 3523.96])
check("契約: h(T) が MOC と FrozenGas でビット一致 (同じ DB の関数)", np.array_equal(mg.h_mass(Tg), fg.h_mass(Tg)))
cp_lump = ResolvedSpeciesDB.species_cp_R(L.entries["EXH"], Tg) * RU / L.entries["EXH"].MW
check("契約: lump EXH の NASA-9 cp が構成種の質量分率混合と一致 (相対 1e-12)",
      float(np.max(np.abs(cp_lump / fg.cp_mass(Tg) - 1))) < 1e-12)
check("契約: MOC の全温 = 燃焼室温度", mg.Tt == ex.Tc)
check("契約: MOC のスロート温度が 1D と相対 2e-4 以内 (表の分解能)",
      abs(mg.T_star / i50["T_star"] - 1) < 2e-4, f"{mg.T_star:.3f}/{i50['T_star']:.3f}")
check("契約: MOC の面積比 A/A*(M_e,1D) = 50 (相対 1e-3)", abs(float(mg.area_ratio(i50["M_e"])) / 50.0 - 1) < 1e-3,
      f"{float(mg.area_ratio(i50['M_e'])):.4f}")
check("契約: FrozenGas の顕エンタルピー基準 h_sens(298.15) = 0", abs(float(fg.h_sens(298.15)[0])) < 1e-6)

# ---------------------------------------------------------------- 問題 YAML (gas.exhaust の単一ソース)
p = load_problem(PROBLEM)
check("問題: spec.Pt = Pc・spec.Tt = Tc (排気ファイルから)", p.spec["Pt"] == ch["Pc_Pa"] and p.spec["Tt"] == ch["Tc_K"])
check("問題: gas_model は契約の MOC 気体 (Tt = Tc、同じ Y)",
      isinstance(p.gas_model, GasSemiPerfect) and p.gas_model.Tt == ch["Tc_K"] and p.gas_model.Y == p.exhaust.Y)
check("問題: 輸送種 = [EXH] (tp_species の既定)", list(p.species_layout().species) == ["EXH"])
check("問題: transport は lump の構成実種 6 種", sorted(p.transport_for_ns()) == sorted(DEFAULT_PRODUCTS))


def _bad(mut):
    d = yaml.safe_load(PROBLEM.read_text())
    mut(d)
    with tempfile.TemporaryDirectory() as td:
        pth = Path(td) / "bad.yaml"
        (Path(td) / EXH_FILE.name).write_text(EXH_FILE.read_text())
        pth.write_text(yaml.safe_dump(d, allow_unicode=True))
        return raises(lambda: load_problem(pth), (ValueError, FileNotFoundError))


check("拒否: gas.exhaust と spec.Tt の併用", _bad(lambda d: d["spec"].update(Tt=3500.0)) is not None)
check("拒否: gas.exhaust と spec.Pt の併用", _bad(lambda d: d["spec"].update(Pt=1.0e7)) is not None)
check("拒否: gas.exhaust と gas.species の併用", _bad(lambda d: d["gas"].update(species={"H2O": 1.0})) is not None)
check("拒否: gas.model cpg で gas.exhaust", _bad(lambda d: d["gas"].update(model="cpg")) is not None)
check("拒否: 排気ファイルが無い", _bad(lambda d: d["gas"].update(exhaust="missing.yaml")) is not None)
check("拒否: transport の漏れ", _bad(lambda d: d["gas"]["transport"].pop("OH")) is not None)

# ---------------------------------------------------------------- runner: TP の config / BC (forge は起動しない)
from forge_design.evaluate import runner as R  # noqa: E402

BANNED = ("bndFirstOrder", "isCompressible", "meshFormat", "nodeAxisDirichlet", "kInf", "omegaInf", "lineDtWallRelief",
          "implicitRelaxSST", "detectNaNInterval")
for stage, want in (("soft", (0, 0, 10, 0.5)), ("mid", (0, 0, 10, 1.0)), ("main", (1, 2, 4, 4.0))):
    txt = R._solver_config_tp(p, 1000, 250, stage)
    c = yaml.safe_load(txt)
    pp, dt, tb = c["physProp"], c["time"]["deltaT"], c["turbulence"]
    check(f"TP config [{stage}]: 段の値 (convMethod, limiter, nStepInner, cfl_pseudo)",
          (c["space"]["convMethod"], c["space"]["limiter"], c["time"]["nStepInner"], dt["cfl_pseudo"]) == want)
    check(f"TP config [{stage}]: thermalMethod 2・viscMethod 2・thermoHrefTemp 298.15・lump EXH",
          pp["thermalMethod"] == 2 and pp["viscMethod"] == 2 and pp["thermoHrefTemp"] == 298.15
          and pp["species"][0]["name"] == "EXH" and pp["species"][0]["basis"] == "mole")
    check(f"TP config [{stage}]: physProp.transport が構成実種 6 種", sorted(pp["transport"]) == sorted(DEFAULT_PRODUCTS))
    check(f"TP config [{stage}]: SST 低 Re (wallTreatmentSST 0)・implicitRelax 0.7・output level 1・blockDPLUR",
          tb["model"] == "sst" and tb["wallTreatmentSST"] == 0 and dt["implicitRelax"] == 0.7 and dt["blockDPLUR"] == 1
          and c["output"]["level"] == 1 and c["time"]["timeIntegration"] == 11)
    check(f"TP config [{stage}]: node + nodeWallDirichlet + nodeInletCornerWall + 軸対称",
          c["mesh"]["discretization"] == "node" and c["mesh"]["nodeWallDirichlet"] == 1
          and c["mesh"]["nodeInletCornerWall"] == 1 and c["mesh"]["isAxisymmetric"] == 1)
    check(f"TP config [{stage}]: 廃止・禁止キーを含まない", not any(b in txt for b in BANNED) and "\n  ro:" not in txt)
    check(f"TP config [{stage}]: 参照 cp・gamma が混合気の R を与える",
          abs(pp["cp"] * (pp["gamma"] - 1) / pp["gamma"] / fg.R - 1) < 1e-12)
check("TP config: 品質検査用の cell 変換は node 専用キーを落とす",
      "nodeInletCornerWall" not in R._solver_config_tp(p, 10, 10, "soft", discretization="cell"))
bc = yaml.safe_load(R._bcond_config(p))
check("TP bcond: 入口 Pt = Pc・Tt = Tc、出口 outflow (真空)、壁 wall (断熱)",
      bc["inlet"]["floats"]["Pt"] == ch["Pc_Pa"] and bc["inlet"]["floats"]["Tt"] == ch["Tc_K"]
      and bc["outlet"]["kind"] == "outflow" and bc["wall"]["kind"] == "wall")
p.evaluate["outlet_kind"] = "outlet_statPress"
check("TP bcond: evaluate.outlet_kind outlet_statPress を明示すれば静圧出口",
      yaml.safe_load(R._bcond_config(p))["outlet"]["kind"] == "outlet_statPress")
p.evaluate["outlet_kind"] = "bogus"
check("TP bcond: 未知の outlet_kind は拒否", raises(lambda: R._bcond_config(p), ValueError) is not None)
p.evaluate.pop("outlet_kind")
p.evaluate["wall_treatment_sst"] = 1
check("TP config: 壁関数 (wall_treatment_sst 1) は拒否", raises(lambda: R._solver_config_tp(p, 10, 10), ValueError) is not None)
p.evaluate.pop("wall_treatment_sst")
check("TP 段: 冷間 IC から main は拒否 (soft から)", raises(lambda: R._tp_stage("main", None), ValueError) is not None)
check("TP 段: 既定は冷間 soft / warm main", R._tp_stage(None, None) == "soft" and R._tp_stage(None, "x.h5") == "main")
id1 = R.ideal_1d_tp(p)
check("ideal_1d_tp: A_t = π r_t²、F = F_vac (p_a 0)、ṁ ≈ 271 kg/s",
      abs(id1["A_t"] - np.pi * 0.14 ** 2) < 1e-15 and id1["F"] == id1["F_vac"] and abs(id1["mdot"] - 270.87) < 0.1,
      f"ṁ {id1['mdot']:.3f} kg/s, F_vac {id1['F_vac'] / 1e3:.2f} kN")

# ---------------------------------------------------------------- CPG 経路の config が変更前 (HEAD 32779703) と文字単位で同じ
CPG_YAML = """name: bell_top_smoke
type: thruster_bell
gas: {gamma: 1.4, cp: 1004.5}
spec: {Pt: 4.0e6, Tt: 1500.0, p_ambient: 20000.0, r_throat: 0.010, eps: 9.0}
dv:
  theta_a_deg: {value: 20.0, min: 12.0, max: 30.0}
  theta_e_deg: {value: 8.0,  min: 4.0,  max: 14.0}
  L_over_rt:   {value: 7.0,  min: 5.0,  max: 12.0}
geometry: {bell_type: top, r_inlet: 2.5, L_pipe: 0.5, L_contract: 3.0, Ru: 1.5, Rd: 0.4, phi_u_deg: 25.0}
mesh: {ni: 221, nj: 65, wall_first_frac: 5.0e-3, throat_refine: 3.0}
evaluate: {turbulence: sst, nStepOuter: 12000, outStepInterval: 4000, objectives: [eta_cf, length]}
"""
with tempfile.TemporaryDirectory() as td:
    (Path(td) / "cpg.yaml").write_text(CPG_YAML)
    pc = load_problem(Path(td) / "cpg.yaml")
    h_cfg = hashlib.sha256(R._solver_config(pc, 12000, 4000).encode()).hexdigest()
    h_bc = hashlib.sha256(R._bcond_config(pc).encode()).hexdigest()
    msg = raises(lambda: R.prepare(Path(td) / "cpg.yaml", Path(td) / "x", stage="soft"), ValueError)
    check("CPG: TP 用の stage を渡すと拒否 (run dir を作る前)", msg is not None and "stage" in msg
          and not (Path(td) / "x").exists(), msg or "")
check("CPG: solverConfig が変更前と同一 (sha256; case/40 problem_bell_top.yaml と同じ入力)",
      h_cfg == "d86128d33d76eafebe8c6e1f2b0367c9ec82ccf3ebccc9acbb47c2f29d2e885c", h_cfg[:16])
check("CPG: bcondConfig が変更前と同一 (sha256; outlet_statPress のまま)",
      h_bc == "b3ed2129024be5c52295b0d21558d9bc6f61ca3e6d5dd3875f060a4dae2c8d8b", h_bc[:16])

# ---------------------------------------------------------------- TP の壁の平滑化 (P 不変、roe = ρ e_TP(T))
import h5py  # noqa: E402

Te = np.linspace(800.0, 3500.0, 40)
check("TP: T → e → T の往復 (表の逆引き) が 0.01 K 以内",
      float(np.max(np.abs(R._tp_T_of_e(fg, fg.e_sens(Te)) - Te))) < 0.01)
with tempfile.TemporaryDirectory() as td:
    h5 = Path(td) / "w.h5"
    n, nw = 30, 20
    Tw = 2500.0 + 300.0 * (-1.0) ** np.arange(nw)          # 壁列の T 市松 (等圧)
    P0 = 2.0e6
    T = np.r_[Tw, np.full(n - nw, 2000.0)]
    ro = P0 / (fg.R * T)
    ux = np.r_[np.full(nw, 5.0), np.full(n - nw, 800.0)]
    roe = ro * (fg.e_sens(T) + 0.5 * ux * ux)
    with h5py.File(h5, "w") as f:
        g = f.create_group("BCONDS/b3")
        g.attrs["bcondKind"] = np.bytes_("wall")
        g.create_dataset("iCells", data=np.arange(nw))
        f.create_dataset("CELLS/centCoords", data=np.c_[np.arange(n, dtype=float), np.zeros(n), np.zeros(n)].ravel())
        for k, v in (("ro", ro), ("roe", roe), ("roUx", ro * ux), ("roUy", 0 * ro), ("roUz", 0 * ro)):
            f.create_dataset(f"VALUE/{k}", data=v.astype(np.float32))
    R._smooth_wall_entropy(h5, 1.2, 3000.0, gas=fg)
    with h5py.File(h5) as f:
        ro2, roe2, rux2 = f["VALUE/ro"][:].astype(float), f["VALUE/roe"][:].astype(float), f["VALUE/roUx"][:].astype(float)
    T2 = R._tp_T_of_e(fg, roe2[:nw] / ro2[:nw])
    P2 = ro2[:nw] * fg.R * T2
    check("TP 平滑化: 壁列の P が不変 (float32 の丸め以内)", float(np.max(np.abs(P2 / P0 - 1))) < 1e-5,
          f"{float(np.max(np.abs(P2 / P0 - 1))):.1e}")
    check("TP 平滑化: 壁列の T の市松が 1/10 以下 (内側の点)",
          float(np.max(np.abs(np.diff(T2[3:-3], 2)))) < 0.1 * 4 * 300.0, f"{float(np.max(np.abs(np.diff(T2[3:-3], 2)))):.1f} K")
    check("TP 平滑化: 壁列の運動量 0・壁外は不変",
          float(np.max(np.abs(rux2[:nw]))) == 0.0 and np.array_equal(ro2[nw:], ro[nw:].astype(np.float32).astype(float)))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
