#!/usr/bin/env python3
"""湿り場の変換と H2O 潜熱の気液ペア (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10; §6 V7 (e)(e′)) の試験。

  python3 solver_density_cuda/tests/unit/test_convert_latent_pair.py [--forge BIN]     (GPU 不要; BIN は --resolve-species 対応)

(1) Python の気液ペア (convert_species_field.H2OLatentPair) が、試験内で独立に書いた L = h_v − h_l (種 DB の外挿規約 + 液相の延長規約)
    と 120–400 K で相対 1e-12、ソルバ (tests/unit/test_cond_latent_pair.cu の出力値) と 1e-6 J/kg で一致する
(2) V7(e′) 0 step A/B (codex diagnose 2026-09-27): 150 K, N2/H2O 0.95/0.05 (総質量分率), g 0.01, 気相同一・L だけ異なる
    (SRC = #10 以前の h2o_latent, DST = 気液ペア)。A = 気相差だけの旧式 → 補正 0、宛先 EOS で反転した T − 150 = −0.005802 K を再現。
    B = 気液を含む全差 (roe_delta_per_mass) → 補正 +4.6592 J/kg。**B の合格: 独立な二相エネルギー評価との補正誤差 ≤1e-6 J/kg、T 相対誤差 ≤1e-8**
    段 3 (#13-3) で H2O の MW が CEA の 0.01801528 に変わり (旧 0.0180153; #10 以前の h2o_latent も旧 MW)、新 L が 150 K で +3.101 J/kg
    (相対 +1.11e-6) 動いたので、アンカーは codex の値をその分だけ動かした値 (補正 +4.6282 J/kg = 0.01 × 462.819 J/kg,
    T − 150 = −0.005802 × 4.6282/4.6592 K) にした。MW を戻した差 (datum・延長規約だけの差) は test_cond_latent_pair.cu が −465.920 J/kg で見る
(3) V7(e): 同じ潜熱モデルで datum (thermoHrefTemp 0 → 298.15) と種の順序を変える湿り場の変換で T が保たれる (相対 ≤1e-8)
(4) 潜熱モデル不明 (記録に液相なし) の湿りセルは拒否、乾いたセルは通す。液相モデルだけの違いで再構成が発火する (_db_differs)
(5) 端から端 (forge --resolve-species を使う; 無ければ SKIP): #10 以前の記録 (液相なし) で印を付けた H2O 湿り場を、同じ config の宛先へ変換する
    - 既定 → 拒否 (移行手順 --src-latent legacy-v0 を案内)
    - --src-latent legacy-v0 → 通り、書かれた roe/ρ が旧 → 新の潜熱差ぶん (g ΔL) だけ動き (float32 量子化の範囲で)、宛先の属性が付く
    - 記録に液相がある場で --src-latent legacy-v0 → 拒否
規約: [PASS]/[FAIL]/[SKIP]、失敗があれば非ゼロ終了。
"""
import argparse, hashlib, os, re, shutil, subprocess, sys, tempfile

import numpy as np
import h5py
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, "..", ".."))
TOOLS = os.path.join(SOLVER, "tools")
sys.path.insert(0, TOOLS)
import forge_species as fsp  # noqa: E402
import convert_species_field as csf  # noqa: E402
from total_quantities import _TPGas, RU  # noqa: E402

DATA = os.path.join(SOLVER, "data", "species", "forge_species_v1.yaml")
# 150 K・g 0.01 の補正 g (L_old − L_new) のアンカー [J/kg]: codex diagnose 2026-09-27 の 4.6592 (旧 MW) を、段 3 (#13-3) の H2O MW 変更で
# 新 L が動いた分 (test_cond_latent_pair.cu の表: 150 K の L_new − L_old = −462.819076 J/kg) に合わせた値
CORR_150 = 0.01 * 462.819076
g_fail = 0


def check(ok, what, val=None, tol=None):
    global g_fail
    extra = "" if val is None else f"  (value {val:.6e}" + ("" if tol is None else f", tol {tol:.1e}") + ")"
    print(("[PASS] " if ok else "[FAIL] ") + what + extra)
    if not ok:
        g_fail += 1


# ----------------------------------------------------------------------------- 共通データから記録と同じ形のエントリを作る
def common_entries():
    raw = yaml.safe_load(open(DATA))
    sp = {str(e["id"]): e for e in raw["species"]}

    def gas(i):
        e = sp[i]; iv = e["intervals"]
        return {"MW": float(e["MW"]), "Tlo": float(iv[0]["Tlo"]), "Tmid": float(iv[0]["Thi"]), "Thi": float(iv[1]["Thi"]),
                "nasa9_low": [float(x) for x in iv[0]["coeffs"]], "nasa9_high": [float(x) for x in iv[1]["coeffs"]],
                "LJ_sigma": float(e["LJ"]["sigma"]), "LJ_eps_kB": float(e["LJ"]["eps_kB"])}
    L = sp["H2O(L)"]
    cond = {"name": "H2O(L)", "phase": "condensed", "pair_of": str(L["pair_of"]), "gas_index": 1, "gas_name": "H2O",
            "MW": float(L["MW"]), "Tlo": float(L["intervals"][0]["Tlo"]), "Thi": float(L["intervals"][0]["Thi"]),
            "nasa9": [float(x) for x in L["intervals"][0]["coeffs"]],
            "extension": {"below": L["extension"]["below"], "above": L["extension"]["above"]},
            "rule": fsp.SPECIES_CONDENSED_EXTENSION, "latent": fsp.SPECIES_CONDENSED_LATENT, "datum": fsp.SPECIES_CONDENSED_DATUM}
    return {"N2": gas("N2"), "H2O": gas("H2O")}, cond


# 独立参照: 種 DB と同じ気相の評価 (T<Tlo は c_p(Tlo) 一定の線形、T<Tmid で low) と液相の延長規約を素直に書いたもの
def ref_latent(gas, cond, T):
    def h_mol(a, T):
        return RU * T * (-a[0] / T**2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4
                         + a[6] * T**4 / 5 + a[7] / T)

    def cp_mol(a, T):
        return RU * (a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4)

    def hv(T):
        a = gas["nasa9_low"] if T < gas["Tmid"] else gas["nasa9_high"]
        if T < gas["Tlo"]:
            return (h_mol(gas["nasa9_low"], gas["Tlo"]) + cp_mol(gas["nasa9_low"], gas["Tlo"]) * (T - gas["Tlo"])) / gas["MW"]
        return h_mol(a, T) / gas["MW"]
    al, Tlo, Thi, MW = cond["nasa9"], cond["Tlo"], cond["Thi"], gas["MW"]

    def hl_mass(a, T):
        # 液相は a0=1.3e9, a1=−2.4e7 の項が強く打ち消すので、ソルバ (cond_liquid_h_abs_poly) と同じ演算順 hRT·(R_u/MW)·T で評価する
        # (順序が違うと丸めだけで ~1e-11 相対ずれる; 規約の独立性は区間・延長の書き方で保つ)
        return (-a[0] / (T * T) + a[1] * np.log(T) / T + a[2] + a[3] * T / 2.0 + a[4] * T * T / 3.0
                + a[5] * T * T * T / 4.0 + a[6] * T * T * T * T / 5.0 + a[7] / T) * (RU / MW) * T

    def hl(T):
        if T < Tlo:
            cpl = (hl_mass(al, Tlo + 0.5) - hl_mass(al, Tlo - 0.5 + 1.0e-9)) / 1.0
            return hl_mass(al, Tlo) - cpl * (Tlo - T)
        return hl_mass(al, min(T, Thi))
    Tg = max(T, 45.0)
    return min(max(hv(Tg) - hl(Tg), 1.5e6), 3.5e6)


def unit_tests():
    DB, cond = common_entries()
    pair = csf.H2OLatentPair(DB["H2O"], cond)
    Ts = np.arange(120.0, 400.0 + 1e-9, 0.25)
    worst = max(abs(pair(np.array([T]))[0] - ref_latent(DB["H2O"], cond, T)) / ref_latent(DB["H2O"], cond, T) for T in Ts)
    check(worst <= 1e-12, "(1) Python pair vs independent L = h_v - h_l, 120-400 K every 0.25 K (rel)", worst, 1e-12)
    # ソルバ (C++) の値 (tests/unit/test_cond_latent_pair.cu の表 L_new, 小数 6 桁)。段 3 (#13-3) の H2O MW 変更後の値
    # (変更前: 120 K 2864570.940375, 150 K 2793256.307606, 250 K 2555665.244644, 298.15 K 2442581.404111; 相対 +1.11e-6)
    cxx = {120.0: 2864574.120533, 150.0: 2793259.408592, 250.0: 2555668.081863, 298.15: 2442584.115788}
    wc = max(abs(pair(np.array([T]))[0] - v) for T, v in cxx.items())
    check(wc <= 1.0e-6, "(1) Python pair vs solver h2o_latent (C++ table, 6 decimals) [J/kg]", wc, 1e-6)

    # ---- (2) V7(e′) A/B ----
    names = ["N2", "H2O"]; Tref = 298.15
    gas = _TPGas(DB, names, Tref)
    Y = [np.array([0.95]), np.array([0.05])]
    T0 = np.array([150.0]); g = np.array([0.01]); ro = 1.0
    eos_s = csf.CondEOS(1, carrier=True, h2o_latent=csf.h2o_latent_legacy_v0, h2o_latent_key=csf.LEGACY_LATENT_KEY)
    eos_d = csf.CondEOS(1, carrier=True, h2o_latent=pair, h2o_latent_key=pair.key)
    Rm = gas.Rmix(Y)
    e_src = gas.h(Y, T0) - Rm * T0 + eos_s.e_liquid_term(T0, g, Rm)      # 旧モデルで作った湿り場の内部エネルギー
    # A: 気相差だけ (旧式) → 補正 0
    dA = (gas.h(Y, T0) - Rm * T0) - (gas.h(Y, T0) - Rm * T0)
    TA = csf.T_from_e(gas, Y, e_src + dA, T0, g=g, eos=eos_d)
    # 期待値: codex の −0.005802 K (旧 MW) を、補正量の比 4.6282/4.6592 (段 3 の MW 由来の ΔL) で動かしたもの
    TA_exp = -0.005802 * CORR_150 / 4.6592
    check(abs(dA[0]) == 0.0 and abs((TA[0] - 150.0) - TA_exp) <= 5e-7,
          f"(2) A (gas-only formula): correction {dA[0]:.1f} J/kg, T - 150 = {TA[0] - 150.0:.7f} K (expect {TA_exp:.7f})",
          abs((TA[0] - 150.0) - TA_exp), 5e-7)
    # B: 気液を含む全差
    dB = csf.roe_delta_per_mass(gas, gas, Y, Y, T0, g, g, eos_s, eos_d)
    Lold = csf.h2o_latent_legacy_v0(150.0); Lnew = ref_latent(DB["H2O"], cond, 150.0)
    dRef = (0.01 * (461.5 * 150.0 - Lnew)) - (0.01 * (461.5 * 150.0 - Lold))   # 独立な二相エネルギー評価の差 (気相は同一)
    check(abs(dB[0] - dRef) <= 1e-6 and abs(dRef - CORR_150) <= 5e-5,
          f"(2) B (gas + liquid): correction {dB[0]:.7f} J/kg vs independent {dRef:.7f} (anchor {CORR_150:.4f} = codex 4.6592 moved by the #13-3 MW)",
          abs(dB[0] - dRef), 1e-6)
    TB = csf.T_from_e(gas, Y, e_src + dB, T0, g=g, eos=eos_d)
    check(abs(TB[0] - 150.0) / 150.0 <= 1e-8, f"(2) B: T after conversion = {TB[0]:.10f} K (relative error)", abs(TB[0] - 150.0) / 150.0, 1e-8)
    check(abs(dB[0] - dA[0]) > 1.0, "(2) A and B differ (the liquid term is what the old formula dropped)", abs(dB[0] - dA[0]))

    # ---- (3) V7(e): 同じ潜熱モデルで datum と種順を変える ----
    gs0 = _TPGas(DB, ["N2", "H2O"], 0.0); gd1 = _TPGas(DB, ["H2O", "N2"], 298.15)
    Tw = np.array([150.0, 199.9, 230.0, 270.0, 300.0]); gw = np.array([1e-4, 0.01, 0.02, 0.04, 0.001])
    Ys = [np.full(5, 0.95), np.full(5, 0.05)]; Yd = [Ys[1], Ys[0]]
    e_s = gs0.h(Ys, Tw) - gs0.Rmix(Ys) * Tw + eos_d.e_liquid_term(Tw, gw, gs0.Rmix(Ys))
    e_d = e_s + csf.roe_delta_per_mass(gs0, gd1, Ys, Yd, Tw, gw, gw, eos_d, eos_d)
    Tc = csf.T_from_e(gd1, Yd, e_d, Tw * 1.05, g=gw, eos=eos_d)
    wT = float(np.max(np.abs(Tc - Tw) / Tw))
    check(wT <= 1e-8, "(3) wet field, datum 0 -> 298.15 and species order swapped: T preserved (rel)", wT, 1e-8)

    # ---- (4) 不明な潜熱モデル ----
    eos_u = csf.CondEOS(1, carrier=True, h2o_latent=None, why_unknown="test: record without liquid phase")
    z = eos_u.e_liquid_term(np.array([200.0, 250.0]), np.zeros(2), 280.0)
    check(np.all(z == 0.0), "(4) unknown latent model: dry cells (g=0) pass without evaluating L")
    try:
        eos_u.e_liquid_term(np.array([200.0]), np.array([0.01]), 280.0); ok = False; msg = ""
    except SystemExit as e:
        ok = True; msg = str(e)
    check(ok and "--src-latent legacy-v0" in msg, "(4) unknown latent model: wet cell refused with the migration hint")
    lay = {"Tref": 298.15, "db": {"N2": DB["N2"], "H2O": DB["H2O"]}, "names": ["N2", "H2O"]}
    d1, why1 = csf._db_differs(lay, lay, eos_s, eos_d)
    d2, _ = csf._db_differs(lay, lay, eos_d, eos_d)
    check(d1 and not d2 and "latent" in why1, f"(4) liquid-model-only difference triggers roe reconstruction ({why1})")


# ----------------------------------------------------------------------------- (5) 端から端
# case/44 run_0509 由来の有効な config (tests/unit/test_forge_species_record.py と同じ) に凝縮 (H2O carrier) を足したもの
CFG = ('mesh: {discretization: "node", meshFileName: "in.h5", valueFileName: "in.h5"}\n'
       'gpu: 1\nsolver: "SLAU"\n'
       'physProp: {thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1220.7, gamma: 1.31526, '
       'species: ["N2", "H2O"], thermoHrefTemp: 298.15}\n'
       'time:\n  unsteady: 0\n  dualTime: 0\n  last: {nStepOuter: 1}\n'
       '  deltaT: {control: 1, dt: 1e-8, cfl: 6.0, cfl_pseudo: 6.0, dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, '
       'implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}\n'
       '  outStepStart: 0\n  outStepInterval: 1\n  timeIntegration: 11\n  nStepInner: 4\n'
       'space: {convMethod: 1, limiter: 2}\nturbulence: {model: "none"}\ninitial: "uniform_p101325_u10"\n'
       'condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condGasSpecies: 1}\n')
N = 6


def make_run(d):
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "solverConfig.yaml"), "w").write(CFG)
    meta = {"mode": "full", "species": ["N2", "H2O"], "keep": [], "condensing_species": "H2O", "condensing_index": 1,
            "tracer": {"enabled": False, "name": None, "definition": None}, "lumps": {},
            "expansion": {"N2": {"N2": 1.0}, "H2O": {"H2O": 1.0}}, "streams": {}, "MW": {"N2": 0.0280134, "H2O": 0.01801528}}
    yaml.safe_dump(meta, open(os.path.join(d, "species_meta.yaml"), "w"), sort_keys=False)


def write_h5(path, DB, T, g, latent, ro=0.8, u=50.0):
    gas = _TPGas(DB, ["N2", "H2O"], 298.15)
    Y = [np.full(N, 0.95), np.full(N, 0.05)]
    e = gas.h(Y, T) - gas.Rmix(Y) * T + g * (461.5 * T - latent(T))
    ro = np.full(N, ro)
    V = {"ro": ro, "roUx": ro * u, "roUy": np.zeros(N), "roUz": np.zeros(N), "roe": ro * (e + 0.5 * u * u),
         "roY0": ro * Y[0], "roY1": ro * Y[1], "rog_0": ro * g, "roQ0_0": ro * 1e12, "roQ1_0": ro * 1e4, "roQ2_0": ro * 1e-4}
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.zeros(3 * N))
        for k, v in V.items():
            f.create_dataset("VALUE/" + k, data=np.asarray(v, np.float64 if k == "roe" else np.float32))
    return ro, e


def stamp_pre10(run_dir, h5path):
    """run の今の記録から液相ブロックを除いた「#10 以前の記録」を作り、h5 に印を付ける (新しい記録は消す)。"""
    r = fsp.resolve_species(run_dir, inplace=True)
    rec = r["record"]
    text = open(r["record_path"]).read()
    i = text.index("condensed:\n"); j = text.index("provenance:\n")
    text = text[:i] + text[j:]
    h_old = hashlib.sha256(fsp.compat_text(rec["schema"], rec["datum"], rec["thermoHrefTemp"], rec["extrapolation"], rec["species"],
                                           rec["transport_compat"], None).encode()).hexdigest()
    text = re.sub(r'compat_hash: "[0-9a-f]{64}"', f'compat_hash: "{h_old}"', text)
    os.remove(r["record_path"])
    name = f"resolved_species_{h_old[:16]}.yaml"
    open(os.path.join(run_dir, name), "w").write(text)
    old = fsp.load_record(os.path.join(run_dir, name))
    assert old["consistent"] and old["condensed"] is None, old["problems"]
    fsp.write_species_attrs(h5path, {"species_hash": h_old, "species_record_sha256": old["integrity"],
                                     "species_record_file": name, "species_input_unverified": 0})
    return r["hash"], h_old


def run_convert(src, dst, extra=()):
    cmd = [sys.executable, os.path.join(TOOLS, "convert_species_field.py"), os.path.join(src, "in.h5"), os.path.join(dst, "in.h5"),
           "--meta", os.path.join(dst, "species_meta.yaml"), *extra]
    env = dict(os.environ); env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def e2e_tests(forge):
    try:
        exe = fsp.find_forge(forge)
    except fsp.SpeciesResolveUnavailable as e:
        exe = None; print(f"[SKIP] (5) end-to-end: {e}")
    if exe is None:
        if forge is None:
            print("[SKIP] (5) end-to-end: no forge with --resolve-species (--forge / FORGE_BIN)")
        return
    os.environ["FORGE_BIN"] = exe
    DB, cond = common_entries()
    pair = csf.H2OLatentPair(DB["H2O"], cond)
    root = tempfile.mkdtemp(prefix="test_convert_latent_pair_")
    try:
        src, dst = os.path.join(root, "src"), os.path.join(root, "dst")
        make_run(src); make_run(dst)
        T = np.array([150.0, 180.0, 199.0, 230.0, 260.0, 290.0]); g = np.array([0.01, 0.02, 0.005, 0.03, 0.01, 0.001])
        ro, e_old = write_h5(os.path.join(src, "in.h5"), DB, T, g, csf.h2o_latent_legacy_v0)
        write_h5(os.path.join(dst, "in.h5"), DB, T, np.zeros(N), pair)
        h_new, h_old = stamp_pre10(src, os.path.join(src, "in.h5"))
        check(h_new != h_old, f"(5) pre-#10 record stamped on SRC ({h_old[:16]}; the config now resolves to {h_new[:16]} with the liquid phase)")
        rc, out = run_convert(src, dst)
        check(rc != 0 and "--src-latent legacy-v0" in out, f"(5) default: pre-#10 wet field refused with the migration hint (rc={rc})")
        if not (rc != 0 and "--src-latent legacy-v0" in out):
            print("      " + "\n      ".join(out.strip().splitlines()[-8:]))
        rc, out = run_convert(src, dst, ("--src-latent", "legacy-v0"))
        ok = rc == 0 and "all checks passed" in out
        check(ok, f"(5) --src-latent legacy-v0: converted (rc={rc})")
        if not ok:
            print("      " + "\n      ".join(out.strip().splitlines()[-12:]))
            return
        with h5py.File(os.path.join(dst, "in.h5"), "r") as f:
            roe_w = f["VALUE/roe"][:].astype(np.float64); ro_w = f["VALUE/ro"][:].astype(np.float64)
            at = dict(f.attrs)
        de = roe_w / ro_w - (e_old + 0.5 * 50.0 ** 2)
        expect = g * (csf.h2o_latent_legacy_v0(T) - pair(T))
        q = 2.0 * np.spacing(np.abs(roe_w).astype(np.float32)).astype(np.float64) / ro_w + 1e-9   # 書き込み dtype の量子化
        werr = float(np.max(np.abs(de - expect) / q))
        print("      written Δ(roe)/ρ [J/kg]: " + ", ".join(f"{x:.4f}" for x in de) + "  expected g(L_old−L_new): " + ", ".join(f"{x:.4f}" for x in expect))
        check(werr <= 1.0, "(5) written roe moved by g(L_old - L_new) within the write-dtype quantisation (max |err|/quantum)", werr, 1.0)
        check(abs(expect[0] - CORR_150) < 1e-4, f"(5) 150 K, g 0.01 cell: expected correction {expect[0]:.4f} J/kg (anchor {CORR_150:.4f}; codex 4.6592 before #13-3)")
        hs = at.get("species_hash"); hs = hs.decode() if isinstance(hs, bytes) else hs
        check(hs == h_new, f"(5) destination stamped with the current (liquid-phase) species hash {str(hs)[:16]}")
        # 記録に液相がある場で legacy-v0 → 拒否
        src2 = os.path.join(root, "src2"); make_run(src2)
        write_h5(os.path.join(src2, "in.h5"), DB, T, g, pair)
        r2 = fsp.resolve_species(src2, inplace=True)
        fsp.write_species_attrs(os.path.join(src2, "in.h5"), {"species_hash": r2["hash"], "species_record_sha256": r2["record"]["integrity"],
                                                              "species_record_file": r2["record_file"], "species_input_unverified": 0})
        rc, out = run_convert(src2, dst, ("--src-latent", "legacy-v0", "--dry-run"))
        check(rc != 0 and "矛盾" in out, f"(5) record with the liquid phase + --src-latent legacy-v0 -> refused (rc={rc})")
        rc, out = run_convert(src2, dst, ("--dry-run",))
        check(rc == 0 and "all checks passed" in out, f"(5) record with the liquid phase, same model on both sides -> passes (rc={rc})")
        if rc != 0:
            print("      " + "\n      ".join(out.strip().splitlines()[-8:]))
    finally:
        shutil.rmtree(root, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=None, help="--resolve-species 対応の forge (既定: FORGE_BIN)")
    a = ap.parse_args()
    unit_tests()
    e2e_tests(a.forge)
    print(("ALL PASS" if g_fail == 0 else "FAILED") + f" ({g_fail} failures)")
    return 0 if g_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
