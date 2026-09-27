"""CEA (thermo.inp / trans.inp) と forge (共通データ + LJ/Chapman–Enskog/修正 Eucken) の単成分物性の比較データを作る。
usage (リポジトリルートから): python3 notes/investigations/2026-09-27-cea-vs-forge-properties/compute.py > .../data.json
forge 側: solver_density_cuda/data/species/forge_species_v1.yaml の MW・NASA-9 (2 区間)・LJ。
  μ = 2.6693e-6 √(M[g/mol] T) / (σ² Ω22*(T/ε)) [Pa s] (Neufeld Ω22; cuda_forge/thermo_d.cuh thermo_mu_species と同式)
  λ = μ (cp + 1.25 R)  (修正 Eucken; thermo_lambda_species)。cp は forge の NASA-9 (範囲外は端でクランプ + 定 cp; 本比較は 200–6000 K 内)。
CEA 側: thermo.inp の cp・MW、trans.inp の ln η[μP] = A lnT + B/T + C/T² + D, ln λ[μW/(cm K)] 同形 (区間外は最寄り区間で外挿 = CEA 本体と同じ)。
低温 H2O の代替: IAPWS 2008 希薄粘性 μ0、IAPWS 2011 希薄熱伝導 λ0 (公式適用域 253.15 K 以上)。
"""
import importlib.util, json, math, re
import yaml
R_U = 8.314462618
spec = importlib.util.spec_from_file_location("ceat", "solver_density_cuda/tools/cea_thermo_to_species_db.py")
ceat = importlib.util.module_from_spec(spec); spec.loader.exec_module(ceat)
TH = ceat.parse_thermo_inp(".venv-cea/nasa_cea/thermo.inp")
FD = {e["id"]: e for e in yaml.safe_load(open("solver_density_cuda/data/species/forge_species_v1.yaml"))["species"]}

# ---- trans.inp ----
lines = open(".venv-cea/nasa_cea/trans.inp").read().splitlines()
num = re.compile(r"-?\d\.\d+E[ +-]\d+")
TR = {}; i = 1
while i < len(lines):
    h = lines[i]; m = re.search(r"V(\d)C(\d)", h)
    if not m: i += 1; continue
    a, b = h[0:15].strip(), h[15:30].strip(); nv, nc = int(m.group(1)), int(m.group(2)); rec = {"V": [], "C": [], "src": h[40:].strip()}
    for k in range(nv + nc):
        L = lines[i + 1 + k]
        rec[L[1]].append((float(L[2:11]), float(L[11:20]), [float(x.replace("E ", "E+")) for x in num.findall(L[20:])]))
    if not b: TR[a] = rec
    i += 1 + nv + nc

def fit(ivs, T):
    lo = min(ivs, key=lambda v: v[0]); hi = max(ivs, key=lambda v: v[1])
    use = next((v for v in ivs if v[0] <= T <= v[1]), lo if T < lo[0] else hi)
    A, B, C, D = use[2]
    return math.exp(A * math.log(T) + B / T + C / T**2 + D)
def cea_mu(sp, T): return fit(TR[sp]["V"], T) * 1e-7          # μP -> Pa s
def cea_lam(sp, T): return fit(TR[sp]["C"], T) * 1e-4         # μW/(cm K) -> W/(m K)
def cea_range(sp, k): v = TR[sp][k]; return (min(x[0] for x in v), max(x[1] for x in v))

def nasa_cp(a): return lambda T: (a[0]/T**2 + a[1]/T + a[2] + a[3]*T + a[4]*T**2 + a[5]*T**3 + a[6]*T**4)
def forge_cp(e, T):
    iv = e["intervals"]; Tc = min(max(T, iv[0]["Tlo"]), iv[-1]["Thi"])
    a = next(v["a"] if "a" in v else v["coeffs"] for v in iv if v["Tlo"] <= Tc <= v["Thi"])
    return nasa_cp(a)(Tc) * R_U / e["MW"]
def cea_cp(name, T):
    rec = TH[name]; iv = rec["intervals"]
    a = next(v[2] for v in iv if v[0] <= T <= v[1])
    return nasa_cp(a)(T) * R_U / (rec["MW"] * 1e-3)

def om22(Ts): return 1.16145/Ts**0.14874 + 0.52487*math.exp(-0.77320*Ts) + 2.16178*math.exp(-2.43787*Ts)
def forge_mu(e, T):
    M = e["MW"] * 1e3; s = e["LJ"]["sigma"]; eps = e["LJ"]["eps_kB"]
    return 2.6693e-6 * math.sqrt(M * T) / (s * s * om22(T / eps))
def forge_lam(e, T): return forge_mu(e, T) * (forge_cp(e, T) + 1.25 * R_U / e["MW"])

def iapws_mu0(T):
    Tb = T / 647.096; H = [1.67752, 2.20462, 0.6366564, -0.241605]
    return 100 * math.sqrt(Tb) / sum(H[k] / Tb**k for k in range(4)) * 1e-6
def iapws_lam0(T):
    Tb = T / 647.096; L = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]
    return math.sqrt(Tb) / sum(L[k] / Tb**k for k in range(5)) * 1e-3

SPECIES = ["N2", "O2", "Ar", "CO2", "H2O", "He", "H2", "CO", "NO", "OH", "H", "O"]
TCOLS = [200, 300, 600, 1000, 2000, 3000]
out = {"species": [], "h2o": {}, "meta": {}}
for sp in SPECIES:
    e = FD.get(sp)
    if e is None: continue
    row = {"id": sp, "MW_forge": e["MW"], "MW_cea": TH[sp]["MW"] * 1e-3, "LJ": e.get("LJ"), "has_trans": sp in TR}
    row["MW_rel"] = e["MW"] / (TH[sp]["MW"] * 1e-3) - 1
    row["cp"], row["mu"], row["lam"] = {}, {}, {}
    for T in TCOLS:
        row["cp"][T] = forge_cp(e, T) / cea_cp(sp, T) - 1
        if sp in TR and e.get("LJ"):
            vlo, vhi = cea_range(sp, "V"); clo, chi = cea_range(sp, "C") if TR[sp]["C"] else (None, None)
            row["mu"][T] = {"rel": forge_mu(e, T) / cea_mu(sp, T) - 1, "in": vlo <= T <= vhi}
            if TR[sp]["C"]:
                row["lam"][T] = {"rel": forge_lam(e, T) / cea_lam(sp, T) - 1, "in": clo <= T <= chi}
    if sp in TR:
        row["trans_src"] = TR[sp]["src"]; row["V_range"] = cea_range(sp, "V")
        row["C_range"] = cea_range(sp, "C") if TR[sp]["C"] else None
    out["species"].append(row)
# 曲線 (μ・λ の forge/CEA 比、8 種)
Tg = [200 + 10 * k for k in range(0, 281)]  # 200..3000
out["curves"] = {}
for sp in ["N2", "O2", "Ar", "CO2", "H2O", "He", "H2", "CO"]:
    e = FD[sp]; vlo, vhi = cea_range(sp, "V"); clo, chi = cea_range(sp, "C")
    out["curves"][sp] = {
        "mu": [[T, 100 * (forge_mu(e, T) / cea_mu(sp, T) - 1), vlo <= T <= vhi] for T in Tg],
        "lam": [[T, 100 * (forge_lam(e, T) / cea_lam(sp, T) - 1), clo <= T <= chi] for T in Tg]}
# H2O 低温の作戦比較
e = FD["H2O"]; Th = [150 + 5 * k for k in range(0, 191)]  # 150..1100
def blend_mu(T):
    if T >= 373.2: return cea_mu("H2O", T)
    if T >= 253.15: return iapws_mu0(T)
    n = math.log(iapws_mu0(263.15) / iapws_mu0(253.15)) / math.log(263.15 / 253.15)
    return iapws_mu0(253.15) * (T / 253.15)**n
def blend_lam(T):
    if T >= 373.2: return cea_lam("H2O", T)
    if T >= 253.15: return iapws_lam0(T)
    n = math.log(iapws_lam0(263.15) / iapws_lam0(253.15)) / math.log(263.15 / 253.15)
    return iapws_lam0(253.15) * (T / 253.15)**n
out["h2o"] = {"T": Th,
    "mu": {"forge_LJ": [forge_mu(e, T) * 1e6 for T in Th], "CEA": [cea_mu("H2O", T) * 1e6 for T in Th],
           "IAPWS": [iapws_mu0(T) * 1e6 for T in Th], "blend": [blend_mu(T) * 1e6 for T in Th]},
    "lam": {"forge_LJ": [forge_lam(e, T) * 1e3 for T in Th], "CEA": [cea_lam("H2O", T) * 1e3 for T in Th],
            "IAPWS": [iapws_lam0(T) * 1e3 for T in Th], "blend": [blend_lam(T) * 1e3 for T in Th]},
    "cea_lo": 373.2, "iapws_lo": 253.15}
out["meta"] = {"TCOLS": TCOLS, "forge_data": "solver_density_cuda/data/species/forge_species_v1.yaml",
               "cea": ".venv-cea/nasa_cea/{thermo,trans}.inp"}
print(json.dumps(out))
