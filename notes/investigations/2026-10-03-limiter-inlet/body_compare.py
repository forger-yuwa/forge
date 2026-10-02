# limiter-inlet-column-oscillation §4.3: control (cfl 2) vs cfl 1 — R_inlet = sqrt(sum_{i in S} res_ro^2), global rms_* floors (tail 10 % / 20 %)
import h5py, numpy as np, pandas as pd, glob
C = "/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
R = {"ctrl_cfl2": "run_0527_liminlet_ctrl_cfl2", "cfl1": "run_0528_liminlet_cfl1"}
cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
out = {}
for k, r in R.items():
    d = pd.read_csv(C + r + "/residual_history.csv").groupby("step").last().reset_index(); n = len(d)
    for tl in (0.1, 0.2):
        t = d.iloc[int((1 - tl) * n):]
        for c in cols: out.setdefault(k, {})[f"{c}@{int(tl*100)}%"] = float(np.median(t[c]))
    fs = sorted(glob.glob(C + r + "/res_*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))
    fs = [f for f in fs if int(f.split("_")[-1][:-3]) >= 0.8 * 48000]
    f0 = h5py.File(fs[0]); x = np.array(f0["MESH/COORD"][...]).reshape(-1, 3)[:, 0] * 1e3
    S = (np.abs(x + 60) < 0.05) | (np.abs(x + 59.37) < 0.05)
    Ri, Rg, sh = [], [], []
    for f in fs:
        res = np.array(h5py.File(f)["VALUE"]["res_ro"][...]).ravel().astype(np.float64)
        Ri.append(np.sqrt(np.sum(res[S] ** 2))); Rg.append(np.sqrt(np.sum(res ** 2))); sh.append(np.sum(res[S] ** 2) / np.sum(res ** 2))
    out[k]["R_inlet med(tail20%)"] = float(np.median(Ri)); out[k]["R_global med"] = float(np.median(Rg)); out[k]["inlet share med"] = float(np.median(sh))
    out[k]["n_snap"] = len(fs)
F = pd.DataFrame(out); F["cfl1/ctrl"] = F["cfl1"] / F["ctrl_cfl2"]
print(F.to_string(float_format=lambda v: f"{v:.3e}"))
