# limiter-inlet-column-oscillation §4.3: control (cfl 2) vs cfl 1 — R_inlet = sqrt(sum_{i in S} res_ro^2), global rms_* floors (tail 10 % / 20 %)
import h5py, numpy as np, pandas as pd, glob
C = "/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
import sys
R = {"ctrl": sys.argv[1], "test": sys.argv[2]} if len(sys.argv) > 2 else {"ctrl": "run_0527_liminlet_ctrl_cfl2", "test": "run_0528_liminlet_cfl1"}
cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
out = {}
for k, r in R.items():
    d = pd.read_csv(C + r + "/residual_history.csv").groupby("step").last().reset_index(); n = len(d)
    for tl in (0.1, 0.2):
        t = d.iloc[int((1 - tl) * n):]
        for c in cols: out.setdefault(k, {})[f"{c}@{int(tl*100)}%"] = float(np.median(t[c]))
    fs = sorted(glob.glob(C + r + "/res_*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))
    nmax = max(int(f.split("_")[-1][:-3]) for f in fs)
    for tl in (0.1, 0.2):   # codex 2026-10-03 (policy) m: R_入口 も末尾 10 %・20 % の両方で集計する
        fz = [f for f in fs if int(f.split("_")[-1][:-3]) >= (1 - tl) * nmax]
        f0 = h5py.File(fz[0]); x = np.array(f0["MESH/COORD"][...]).reshape(-1, 3)[:, 0] * 1e3
        S = (np.abs(x + 60) < 0.05) | (np.abs(x + 59.37) < 0.05)
        Ri, Rg, sh = [], [], []
        for f in fz:
            res = np.array(h5py.File(f)["VALUE"]["res_ro"][...]).ravel().astype(np.float64)
            Ri.append(np.sqrt(np.sum(res[S] ** 2))); Rg.append(np.sqrt(np.sum(res ** 2))); sh.append(np.sum(res[S] ** 2) / np.sum(res ** 2))
        out[k][f"R_inlet@{int(tl*100)}%"] = float(np.median(Ri)); out[k][f"R_global@{int(tl*100)}%"] = float(np.median(Rg))
        out[k][f"inlet share@{int(tl*100)}%"] = float(np.median(sh)); out[k][f"n_snap@{int(tl*100)}%"] = len(fz)
F = pd.DataFrame(out); F["test/ctrl"] = F["test"] / F["ctrl"]
print(F.to_string(float_format=lambda v: f"{v:.3e}"))
