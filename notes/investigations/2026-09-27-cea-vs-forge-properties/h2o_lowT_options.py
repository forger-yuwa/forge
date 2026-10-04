"""373 K 未満の H2O の μ・λ の案を、単成分と N2 との混合物で比べる (plan #5t, §10)。
案: CEA 外挿 / 案 B (373.2 K 以上 CEA, 253.15–373.2 K IAPWS 希薄気体, 253 K 未満は IAPWS の傾きで冪外挿) /
    LJ (forge 現行, 双極子補正なし) / LJ + Stockmayer (Brokaw 1969 の近似 Ω22 = Ω22_LJ + 0.2 δ*²/T*, δ* = μd²/(2εσ³), μd 1.844 D; Chemkin/Cantera 系の扱いの近似)。
λ は LJ 系では修正 Eucken λ = μ(cp + 1.25R)、cp は forge の NASA-9。
混合は CEA frozen 混合則 (mixing_ab.py の mixB と同じ; FCEA2 と ≤0.013 % で一致を確認済み)。N2 は CEA、H2O–N2 相互作用は trans.inp (有効 300–1000 K; 300 K 未満は外挿)。
"""
import math, re, yaml, importlib.util, sys
sys.path.insert(0, "notes/investigations/2026-09-27-cea-vs-forge-properties")
spec = importlib.util.spec_from_file_location("mab", "notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py")
src = open(spec.origin).read().split("def read_out")[0]      # 関数定義だけ借りる
ns = {}; exec(src, ns)
fit, TR, MW = ns["fit"], ns["TR"], ns["MW"]
R_U = 8.314462618
FD = {e["id"]: e for e in yaml.safe_load(open("solver_density_cuda/data/species/forge_species_v1.yaml"))["species"]}
def cp_h2o(T):
    iv = FD["H2O"]["intervals"]; Tc = min(max(T, iv[0]["Tlo"]), iv[-1]["Thi"]); a = next(v["coeffs"] for v in iv if v["Tlo"] <= Tc <= v["Thi"])
    return (a[0]/Tc**2 + a[1]/Tc + a[2] + a[3]*Tc + a[4]*Tc**2 + a[5]*Tc**3 + a[6]*Tc**4) * R_U / FD["H2O"]["MW"]
def om22(Ts): return 1.16145/Ts**0.14874 + 0.52487*math.exp(-0.77320*Ts) + 2.16178*math.exp(-2.43787*Ts)
SIG, EPS, MUD = 2.605, 572.4, 1.844
dstar = (MUD*1e-18)**2 / (2 * EPS*1.380649e-16 * (SIG*1e-8)**3)
def mu_lj(T, polar=False):
    Ts = T / EPS; om = om22(Ts) + (0.2 * dstar**2 / Ts if polar else 0.0)
    return 2.6693e-6 * math.sqrt(18.0153 * T) / (SIG**2 * om)      # Pa s
def lam_lj(T, polar=False): return mu_lj(T, polar) * (cp_h2o(T) + 1.25 * R_U / 0.0180153)
def iap_mu(T): Tb=T/647.096; H=[1.67752,2.20462,0.6366564,-0.241605]; return 100*math.sqrt(Tb)/sum(H[k]/Tb**k for k in range(4))*1e-6
def iap_lam(T): Tb=T/647.096; L=[2.443221e-3,1.323095e-2,6.770357e-3,-3.454586e-3,4.096266e-4]; return math.sqrt(Tb)/sum(L[k]/Tb**k for k in range(5))*1e-3
def pw(f, T):
    if T >= 253.15: return f(T)
    n = math.log(f(263.15)/f(253.15))/math.log(263.15/253.15); return f(253.15)*(T/253.15)**n
OPT = {
 "CEA外挿":      (lambda T: fit(TR[("H2O","")]["V"],T)*1e-7, lambda T: fit(TR[("H2O","")]["C"],T)*1e-4),
 "案B":          (lambda T: fit(TR[("H2O","")]["V"],T)*1e-7 if T>=373.2 else pw(iap_mu,T), lambda T: fit(TR[("H2O","")]["C"],T)*1e-4 if T>=373.2 else pw(iap_lam,T)),
 "LJ(現行)":     (lambda T: mu_lj(T), lambda T: lam_lj(T)),
 "LJ+双極子":    (lambda T: mu_lj(T,True), lambda T: lam_lj(T,True)),
}
def mix(X, T, opt):
    mu_w, la_w = OPT[opt]
    e = {"N2": fit(TR[("N2","")]["V"],T)*1e-7, "H2O": mu_w(T)}
    c = {"N2": fit(TR[("N2","")]["C"],T)*1e-4, "H2O": la_w(T)}
    eij = fit(TR[("H2O","N2")]["V"], T)*1e-7
    sp = [s for s in X if X[s] > 0]; mu = lam = 0.0
    for i in sp:
        sv = sc = 0.0
        for j in sp:
            if i == j: phi = psi = 1.0
            else:
                phi = 2.0/(eij*(MW[i]+MW[j]))*MW[j]*e[i]
                psi = phi*(1+2.41*(MW[i]-MW[j])*(MW[i]-0.142*MW[j])/(MW[i]+MW[j])**2)
            sv += phi*X[j]; sc += psi*X[j]
        mu += e[i]*X[i]/sv; lam += c[i]*X[i]/sc
    return mu, lam
print(f"δ* (H2O, μd 1.844 D) = {dstar:.3f}")
print("\n[単成分 H2O] 各案 / 案B − 1 [%]   (案B の値: μ μPa·s, λ mW/mK)")
print("  T[K]  案B μ    案B λ  | CEA外挿 μ/λ     | LJ(現行) μ/λ    | LJ+双極子 μ/λ   | 対CEA(373K以上)")
for T in (200,220,250,273,300,350,373.2,600,1000):
    bm, bl = OPT["案B"][0](T), OPT["案B"][1](T)
    s = f"{T:6.1f} {bm*1e6:7.3f} {bl*1e3:7.2f} |"
    for o in ("CEA外挿","LJ(現行)","LJ+双極子"):
        s += f" {100*(OPT[o][0](T)/bm-1):+6.1f}/{100*(OPT[o][1](T)/bl-1):+6.1f} |"
    print(s)
print("\n[混合物 N2 + H2O 蒸気, CEA 混合則] 各案 / 案B − 1 [%] (μ/λ)")
for xw in (0.017, 0.06, 0.2):
    print(f" X_H2O = {xw}")
    for T in (200,250,300,350):
        X={"N2":1-xw,"H2O":xw}; bm,bl = mix(X,T,"案B"); s=f"   {T:4d} K |"
        for o in ("CEA外挿","LJ(現行)","LJ+双極子"):
            m,l = mix(X,T,o); s += f" {o} {100*(m/bm-1):+5.2f}/{100*(l/bl-1):+5.2f} |"
        print(s)
