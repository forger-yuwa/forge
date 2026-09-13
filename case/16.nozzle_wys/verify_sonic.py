#!/usr/bin/env python3
"""res_*.h5 の sonic / gamma を host 側で (P, ρ, T, Y, g) から再計算して照合する (実カーネル照合; plan §6 判定 4)。
  legacy (condSonicModel 0): c=√(γ_mix R_mix T), γ=γ_mix (全蒸気 NASA-9)
  two-phase (1):           c²=γ_2φ R_eff T, R_eff=R_mix−g R_w, c_p,2φ=c_p^全蒸気−g L'(T), γ_2φ=c_p,2φ/(c_p,2φ−R_eff)
usage: verify_sonic.py RUN_DIR RES.h5 --model 0|1"""
import sys, os, argparse, h5py, yaml, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("res"); ap.add_argument("--model", type=int, required=True)
a = ap.parse_args()
Ru = 8.314462618
db = yaml.safe_load(open(os.path.join(a.run, "species_db.yaml")))
sp = yaml.safe_load(open(os.path.join(a.run, "solverConfig.yaml")))["physProp"]["species"]
def cp_nasa9(T, c):  # cp/R = a0/T² + a1/T + a2 + a3 T + a4 T² + a5 T³ + a6 T⁴
    return c[0]/T**2 + c[1]/T + c[2] + c[3]*T + c[4]*T**2 + c[5]*T**3 + c[6]*T**4
def cp_species(name, T):
    s = db[name]; T = np.clip(T, s["Tlo"], s["Thi"])   # thermo_d.cuh と同じく範囲外は端でクランプ
    c = np.where(T < s["Tmid"], cp_nasa9(T, s["nasa9_low"]), cp_nasa9(T, s["nasa9_high"]))
    return c*Ru/s["MW"]
# h2o_latent (condensationProperties_d.cuh) の移植
ag = [-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03, -7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04]
al = [1.326371304e+09, -2.448295388e+07, 1.879428776e+05, -7.678995050e+02, 1.761556813e+00, -2.151167128e-03, 1.092570813e-06, 1.101760476e+08]
Rw = 461.5
def h_mass(c, T):
    hRT = -c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T
    return hRT*Rw*T
def latent(T):
    Tg = np.clip(T, 45.0, 1000.0); Tf = 273.15
    hv = h_mass(ag, Tg)
    h0 = h_mass(al, Tf); cpl = (h_mass(al, Tf + 0.5) - h_mass(al, Tf - 0.5 + 1e-9))/1.0
    hl = np.where(Tg >= Tf, h_mass(al, np.minimum(Tg, 373.15)), h0 - cpl*(Tf - Tg))
    return np.clip(hv - hl, 1.5e6, 3.5e6)
f = h5py.File(a.res, "r"); V = f["VALUE"]
T = np.array(V["T"], dtype=np.float64); P = np.array(V["P"], dtype=np.float64); ro = np.array(V["ro"], dtype=np.float64)
Y = [np.array(V[f"Y{i}"], dtype=np.float64) for i in range(len(sp))]; ysum = sum(Y); Y = [y/ysum for y in Y]
g = np.array(V["g_0"], dtype=np.float64) if "g_0" in V else np.zeros_like(T)
cp = sum(Y[i]*cp_species(sp[i], T) for i in range(len(sp))); Rmix = sum(Y[i]*Ru/db[sp[i]]["MW"] for i in range(len(sp)))
gmix = cp/(cp - Rmix); c2 = gmix*Rmix*T; gam = gmix.copy()
if a.model == 1:
    wet = g > 1e-12
    dL = (latent(T + 0.1) - latent(T - 0.1))/0.2
    Reff = Rmix - g*Rw; cp2 = cp - g*dL; cv2 = cp2 - Reff
    ok = wet & (cv2 > 0.05*cp2)
    g2 = np.where(ok, cp2/np.maximum(cv2, 1e-30), gmix); c2 = np.where(ok, g2*Reff*T, c2); gam = g2
c_host = np.sqrt(c2); c_res = np.array(V["sonic"], dtype=np.float64)
rel = np.abs(c_host - c_res)/c_res
wetmask = g > 1e-12
print(f"sonic: nodes={len(T)} wet={int(wetmask.sum())}  max rel |c_host−c_res| all={rel.max():.2e}  wet={rel[wetmask].max() if wetmask.any() else 0:.2e}  dry={rel[~wetmask].max():.2e}")
okg = True
if "gamma" in V:
    gr = np.array(V["gamma"], dtype=np.float64); relg = np.abs(gam - gr)/gr
    print(f"gamma: max rel |γ_host−γ_res| all={relg.max():.2e}  wet={relg[wetmask].max() if wetmask.any() else 0:.2e}  dry={relg[~wetmask].max():.2e}")
    okg = np.all(np.isfinite(relg)) and relg.max() < 1e-5
elif a.model == 1:
    print("gamma: NOT IN OUTPUT (output.level 2 で gamma を書く必要あり) -> FAIL"); okg = False
if wetmask.any():
    j = np.argmax(g); print(f"at g_max={g[j]:.5f}: T={T[j]:.2f} c_res={c_res[j]:.3f} c_host={c_host[j]:.3f} c_allvap={np.sqrt(gmix[j]*Rmix[j]*T[j]):.3f} ({(c_res[j]/np.sqrt(gmix[j]*Rmix[j]*T[j])-1)*100:+.2f} %)")
okc = np.all(np.isfinite(rel)) and rel.max() < 1e-5
print("VERDICT(verify_sonic):", "PASS" if (okc and okg) else "FAIL")
sys.exit(0 if (okc and okg) else 1)
