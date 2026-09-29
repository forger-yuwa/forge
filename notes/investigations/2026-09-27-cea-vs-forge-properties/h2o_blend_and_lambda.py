"""(1) CEA と IAPWS の H2O μ・λ の差の温度分布 (つなぎ位置を決める)
   (2) kinetic 側の λ: 修正 Eucken と Chemkin/Cantera 系 (Warnatz の内部自由度分離 + Parker の回転緩和数) を CEA と比較。
Warnatz: λ = μ/W (f_tr Cv_tr + f_rot Cv_rot + f_vib Cv_vib), f_tr = 5/2 (1 − 2/π · Cv_rot/Cv_tr · A/B), f_rot = (ρD/μ)(1 + 2/π · A/B), f_vib = ρD/μ,
  A = 5/2 − ρD/μ, B = Z_rot + 2/π (5/3 · Cv_rot/R + ρD/μ), ρD/μ = 6/5 · Ω11/Ω22 の逆数… (ρD_kk/μ_k = 6 A*/5, A* = Ω22/Ω11)。
  Z_rot(T) = Z_rot(298) F(298)/F(T), F(T) = 1 + π^{3/2}/2 (ε/kT)^{1/2} + (π²/4 + 2)(ε/kT) + π^{3/2} (ε/kT)^{3/2}。
  極性 (H2O): Ω22 += 0.2 δ*²/T*, Ω11 += 0.19 δ*²/T* (Brokaw)。Z_rot(298): N2 4.0, O2 3.8, CO2 2.1, H2O 4.0 (Chemkin transport data)。
"""
import math, yaml, re
R_U = 8.314462618
FD = {e["id"]: e for e in yaml.safe_load(open("solver_density_cuda/data/species/forge_species_v1.yaml"))["species"]}
lines = open(".venv-cea/nasa_cea/trans.inp").read().splitlines(); num = re.compile(r"-?\d\.\d+E[ +-]\d+"); TR={}; i=1
while i < len(lines):
    h=lines[i]; m=re.search(r"V(\d)C(\d)",h)
    if not m: i+=1; continue
    a,b=h[0:15].strip(),h[15:30].strip(); nv,nc=int(m.group(1)),int(m.group(2)); rec={"V":[],"C":[]}
    for k in range(nv+nc):
        L=lines[i+1+k]; rec[L[1]].append((float(L[2:11]),float(L[11:20]),[float(x.replace("E ","E+")) for x in num.findall(L[20:])]))
    if not b: TR[a]=rec
    i+=1+nv+nc
def fit(ivs,T):
    lo=min(ivs,key=lambda v:v[0]); hi=max(ivs,key=lambda v:v[1]); use=next((v for v in ivs if v[0]<=T<=v[1]), lo if T<lo[0] else hi)
    A,B,C,D=use[2]; return math.exp(A*math.log(T)+B/T+C/T**2+D)
cmu=lambda s,T: fit(TR[s]["V"],T)*1e-7; clam=lambda s,T: fit(TR[s]["C"],T)*1e-4
def iap_mu(T): Tb=T/647.096; H=[1.67752,2.20462,0.6366564,-0.241605]; return 100*math.sqrt(Tb)/sum(H[k]/Tb**k for k in range(4))*1e-6
def iap_lam(T): Tb=T/647.096; L=[2.443221e-3,1.323095e-2,6.770357e-3,-3.454586e-3,4.096266e-4]; return math.sqrt(Tb)/sum(L[k]/Tb**k for k in range(5))*1e-3
print("(1) H2O: CEA / IAPWS − 1 [%]")
for T in (373.2,400,450,500,600,700,800,900,1000,1073.2):
    print(f"  {T:7.1f} K  μ {100*(cmu('H2O',T)/iap_mu(T)-1):+6.2f}   λ {100*(clam('H2O',T)/iap_lam(T)-1):+6.2f}")
# IAPWS の単調性 (下端)
for f,nm in ((iap_mu,'μ'),(iap_lam,'λ')):
    Tm=min(range(100,400), key=lambda T: f(T)); print(f"  IAPWS 希薄 {nm} の最小点 ≈ {Tm} K (これより下で T を下げると増える)")
def cpR(s,T):
    iv=FD[s]["intervals"]; Tc=min(max(T,iv[0]["Tlo"]),iv[-1]["Thi"]); a=next(v["coeffs"] for v in iv if v["Tlo"]<=Tc<=v["Thi"])
    return a[0]/Tc**2+a[1]/Tc+a[2]+a[3]*Tc+a[4]*Tc**2+a[5]*Tc**3+a[6]*Tc**4
def om22(x): return 1.16145/x**0.14874+0.52487*math.exp(-0.77320*x)+2.16178*math.exp(-2.43787*x)
def om11(x): return 1.06036/x**0.15610+0.19300*math.exp(-0.47635*x)+1.03587*math.exp(-1.52996*x)+1.76474*math.exp(-3.89411*x)
GEOM={"N2":"lin","O2":"lin","CO2":"lin","H2O":"nonlin"}; ZROT={"N2":4.0,"O2":3.8,"CO2":2.1,"H2O":4.0}; DIP={"H2O":1.844}
def props(s,T,model):
    e=FD[s]; M=e["MW"]; sig=e["LJ"]["sigma"]; eps=e["LJ"]["eps_kB"]; Ts=T/eps
    d=0.0
    if s in DIP: d=(DIP[s]*1e-18)**2/(2*eps*1.380649e-16*(sig*1e-8)**3)
    o22=om22(Ts)+0.2*d*d/Ts; o11=om11(Ts)+0.19*d*d/Ts
    mu=2.6693e-6*math.sqrt(M*1e3*T)/(sig**2*o22)
    R=R_U/M; cp=cpR(s,T)*R
    if model=="eucken": return mu, mu*(cp+1.25*R)
    cv=cp-R; cvt=1.5*R; cvr=R if GEOM[s]=="lin" else 1.5*R; cvv=max(cv-cvt-cvr,0.0)
    rDmu=6.0/5.0*o22/o11
    F=lambda TT: 1+math.pi**1.5/2*(eps/TT)**0.5+(math.pi**2/4+2)*(eps/TT)+math.pi**1.5*(eps/TT)**1.5
    Z=ZROT[s]*F(298.0)/F(T)
    A=2.5-rDmu; B=Z+2/math.pi*(5/3*cvr/R+rDmu)
    ftr=2.5*(1-2/math.pi*cvr/cvt*A/B); frot=rDmu*(1+2/math.pi*A/B); fvib=rDmu
    return mu, mu*(ftr*cvt+frot*cvr+fvib*cvv)
print("\n(2) kinetic 側 λ の CEA 比 [%] (μ は双極子補正込みの Chapman–Enskog)")
print("  種     T[K]  |  μ     | λ 修正Eucken | λ Warnatz(Chemkin系)")
for s in ("N2","O2","CO2","H2O"):
    for T in (300,600,1000,2000):
        if s=="H2O" and T<373: Tref="(IAPWS比)"; mc,lc=iap_mu(T),iap_lam(T)
        else: Tref=""; mc,lc=cmu(s,T),clam(s,T)
        mu,le=props(s,T,"eucken"); _,lw=props(s,T,"warnatz")
        print(f"  {s:4s} {T:6d}  | {100*(mu/mc-1):+6.1f} | {100*(le/lc-1):+8.1f}    | {100*(lw/lc-1):+8.1f} {Tref}")
