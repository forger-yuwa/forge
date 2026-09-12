#!/usr/bin/env python3
"""res_*.h5 の condTheta_0 (非等温補正 θ) を host 側で (T, P, ρ, Y, g) と種 DB から再計算して照合 (plan condensation-kantrowitz-carrier §6 判定 5)。
mode は solverConfig の condKantrowitz。S>1 のセルのみ比較。condLim_0 の統計も出す。
usage: verify_theta.py RUN_DIR RES.h5"""
import sys, os, h5py, yaml, numpy as np
run, res = sys.argv[1], sys.argv[2]
cfg = yaml.safe_load(open(os.path.join(run, "solverConfig.yaml"))); cond = cfg["condensation"]; mode = int(cond.get("condKantrowitz", 0))
sp = cfg["physProp"]["species"]; db = yaml.safe_load(open(os.path.join(run, "species_db.yaml"))); iv = int(cond["condGasSpecies"])
Ru = 8.314462618; kB = 1.380649e-23; NA = 6.02214076e23; Rw = 461.5; Mw = 0.0180153; cvv_const = 1393.5; cpv_const = 1855.0
def cp_nasa9(T, c): return c[0]/T**2 + c[1]/T + c[2] + c[3]*T + c[4]*T**2 + c[5]*T**3 + c[6]*T**4
def cp_mass(name, T):
    s = db[name]; T = np.clip(T, s["Tlo"], s["Thi"]); c = np.where(T < s["Tmid"], cp_nasa9(T, s["nasa9_low"]), cp_nasa9(T, s["nasa9_high"])); return c*Ru/s["MW"]
# h2o_latent / h2o_psat (condensationProperties_d.cuh の移植)
ag = [-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03, -7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04]
al = [1.326371304e+09, -2.448295388e+07, 1.879428776e+05, -7.678995050e+02, 1.761556813e+00, -2.151167128e-03, 1.092570813e-06, 1.101760476e+08]
def h_mass(c, T):
    return (-c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T)*Rw*T
def latent(T):
    Tg = np.clip(T, 45.0, 1000.0); Tf = 273.15; hv = h_mass(ag, Tg); h0 = h_mass(al, Tf); cpl = (h_mass(al, Tf + 0.5) - h_mass(al, Tf - 0.5 + 1e-9))
    hl = np.where(Tg >= Tf, h_mass(al, np.minimum(Tg, 373.15)), h0 - cpl*(Tf - Tg)); return np.clip(hv - hl, 1.5e6, 3.5e6)
def psat_mk(T):  # Murphy & Koop 2005 liquid
    return np.exp(54.842763 - 6763.22/T - 4.210*np.log(T) + 0.000367*T + np.tanh(0.0415*(T - 218.8))*(53.878 - 1331.22/T - 9.44523*np.log(T) + 0.014025*T))
f = h5py.File(res, "r"); V = f["VALUE"]
T = np.array(V["T"], dtype=float); P = np.array(V["P"], dtype=float); ro = np.array(V["ro"], dtype=float)
Y = [np.array(V[f"Y{i}"], dtype=float) for i in range(len(sp))]; g = np.array(V["g_0"], dtype=float)
pv = ro*np.maximum(Y[iv] - g, 0)*Rw*T; ps = psat_mk(T); S = pv/ps
b = latent(T)/(Rw*T)
Mv = db[sp[iv]]["MW"]
if mode == 1:
    gam = cpv_const/cvv_const; th = np.maximum(2*(gam - 1)/(gam + 1)*b*(b - 0.5), 0.0)
elif mode >= 2:
    cvv = (cp_mass(sp[iv], T) - Ru/Mv)*Mv/Ru; av = np.maximum(Y[iv] - g, 0)/Mv; cs = np.zeros_like(T)
    for i in range(len(sp)):
        if i == iv: continue
        Mi = db[sp[i]]["MW"]; cvi = (cp_mass(sp[i], T) - Ru/Mi)*Mi/Ru; cs += np.where(Y[i] > 0, Y[i]/Mi*np.sqrt(Mv/Mi)*(cvi + 0.5), 0.0)
    qhat = b - 0.5 - (np.log(np.maximum(S, 1e-300)) if mode == 3 else 0.0); den = av*(cvv + 0.5) + cs
    th = np.where((den > 0) & (av > 0), av*qhat**2/np.maximum(den, 1e-300), 0.0)
else: th = np.zeros_like(T)
thr = np.array(V["condTheta_0"], dtype=float); lim = np.array(V["condLim_0"], dtype=float) if "condLim_0" in V else None
ms = (S > 1.0)                 # 過飽和セル (condLim 統計用)
m = ms & (thr > 0)             # θ 照合は θ>0 のセル (mode 0 では空)
# 相対差は θ_res>0.1 のセルで評価 (蒸気枯渇 a_v→0 で θ→0 のセルは float 保存の Y,g の丸めが相対差を増幅するので絶対差で見る)
mc = m & (thr > 0.1)
rel = np.abs(th[mc] - thr[mc])/np.maximum(thr[mc], 1e-300) if mc.any() else np.array([0.0])
absd = np.abs(th[m] - thr[m]).max() if m.any() else 0.0
# psat の実装差 (Murphy-Koop 移植) は ln S にしか入らないので mode 3 だけ影響: 参考に S の再現も出す
Sres = np.array(V["condS_0"], dtype=float); relS = np.abs(S[m] - Sres[m])/np.maximum(Sres[m], 1e-300) if m.any() else np.array([0.0])
print(f"mode={mode}  S>1 cells={int(m.sum())}  theta_res range [{thr[m].min() if m.any() else 0:.3g}, {thr[m].max() if m.any() else 0:.3g}]  max rel |theta_host-theta_res| (theta>0.1) = {rel.max():.2e}, max abs = {absd:.2e}  (S recon rel {relS.max():.1e})")
if lim is not None and ms.any():
    nuc = ms & (g < 1e-4); print(f"condLim_0: min over S>1 = {lim[ms].min():.3f}; over nucleation zone (S>1, g<1e-4) min={lim[nuc].min() if nuc.any() else float('nan'):.3f} mean={lim[nuc].mean() if nuc.any() else float('nan'):.3f}; cells with lim<0.9: {int((lim[ms] < 0.9).sum())}")
ok = rel.max() < 2e-4 and absd < 1e-3
print("VERDICT(verify_theta):", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
