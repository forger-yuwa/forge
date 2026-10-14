"""plan tooling-nozzle-core-grid §4.90 の材料 (探索的・登録外): step 20000 の残差を圧力の成分 R_p/c² とエントロピーの成分に分ける。usage: python3 cg_res_split.py <run> ... (各 run の res_20000)"""
import json, h5py, numpy as np, sys
# 残差を圧力の成分 R_p = r·R (r = ∂p/∂Q、TP の一般 EOS: (χ + κ e_k, −κu, −κv, −κw, κ)、χ = c² − κ(h − e_k)) と
# エントロピーの成分 R_s = R_ρ − R_p/c² に分け、領域ごとに体積重み付き RMS (単位体積あたり) を出す。κ = γ − 1、γ = c²ρ/P。
runs = sys.argv[1:]
for run in runs:
    info = json.load(open(run + "/prepare_info.json")); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    h = h5py.File(run + "/res_20000.h5", "r")
    v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("res_ro", "res_roUx", "res_roUy", "res_roe", "volume", "ro", "P", "Ux", "Uy", "sonic", "h0")}
    x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S; nj = x.size // ni; jw = nj - 1 - np.arange(x.size) % nj
    c2 = v["sonic"] ** 2; gam = c2 * v["ro"] / v["P"]; kap = gam - 1; ek = 0.5 * (v["Ux"] ** 2 + v["Uy"] ** 2); hs = v["h0"] - ek
    chi = c2 - kap * hs
    Rp = (chi + kap * ek) * v["res_ro"] - kap * v["Ux"] * v["res_roUx"] - kap * v["Uy"] * v["res_roUy"] + kap * v["res_roe"]
    Rs = v["res_ro"] - Rp / c2
    V = v["volume"]
    out = []
    for lab, m in (("縮流部", (x >= -8) & (x <= -2)), ("壁際", (x >= -2) & (x <= 0) & (jw <= 10)), ("島", (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100))):
        rms = lambda z: float(np.sqrt(np.sum(V[m] * (z[m] / V[m]) ** 2) / np.sum(V[m])))
        a, b = rms(Rp / c2), rms(Rs)
        out.append("%s 圧力の成分 %.3e エントロピーの成分 %.3e 比 %.3g" % (lab, a, b, a / b))
    print(run, " | ".join(out))
