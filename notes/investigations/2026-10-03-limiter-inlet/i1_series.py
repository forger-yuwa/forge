# limiter-inlet-column-oscillation §5.1 #1 (I1): psi / res_ro time series at the residual-dominant nodes of run_0526 (5-step output, 200 steps)
import h5py, numpy as np, glob, pandas as pd, sys
C = "/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
D = C + "run_0526_liminlet_i1_ts/"
fs = sorted(glob.glob(D + "res_*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))
# res_0 は IC の書き出しで limiter/res が未計算 (0) なので除く
fs = [f for f in fs if not f.endswith("/res_0.h5")]
steps = [int(s.split("_")[-1][:-3]) for s in fs]
f0 = h5py.File(fs[0]); xyz = np.array(f0["MESH/COORD"][...]).reshape(-1, 3) * 1e3
V = ["limiter_ro", "limiter_P", "limiter_Ux", "limiter_Uy", "res_ro"]
A = {k: np.array([np.array(h5py.File(f)["VALUE"][k][...]).ravel() for f in fs]) for k in V}
res = np.abs(A["res_ro"]); tot = (res ** 2).sum(1)
x = xyz[:, 0]
def grp(m): return (res[:, m] ** 2).sum(1) / tot
print("snapshots", len(fs), "steps", steps[0], "..", steps[-1])
for nm, m in (("x=-60", np.abs(x + 60) < 0.05), ("x=-59.37", np.abs(x + 59.37) < 0.05), ("x=-32 wall", (np.abs(x + 32) < 0.05)), ("x=95", np.abs(x - 95) < 0.05)):
    s = grp(m); print(f"share {nm:11s} mean {s.mean():.3f} min {s.min():.3f} max {s.max():.3f}")
# nodes: top 12 by time-mean res^2
o = np.argsort((res ** 2).mean(0))[::-1][:12]
rows = []
for i in o:
    r = {"x": round(x[i], 2), "y": round(xyz[i, 1], 2), "res_mean": res[:, i].mean(), "res_cv": res[:, i].std() / res[:, i].mean()}
    for k in V[:4]:
        p = A[k][:, i]; d = np.diff(p)
        r[k[8:] + "_mean"] = p.mean(); r[k[8:] + "_std"] = p.std()
        nz = d[np.abs(d) > 1e-4]
        r[k[8:] + "_flipfrac"] = np.mean(np.sign(nz[1:]) != np.sign(nz[:-1])) if nz.size > 2 else np.nan
    # correlation of |res| with the psi of the variable that moves most
    kmax = max(V[:4], key=lambda k: A[k][:, i].std())
    r["most_moving"] = kmax[8:]; r["corr_res_psi"] = np.corrcoef(res[:, i], A[kmax][:, i])[0, 1] if A[kmax][:, i].std() > 0 else np.nan
    rows.append(r)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(pd.DataFrame(rows).to_string(float_format=lambda v: f"{v:.3g}"))
# residual history periodicity (outer steps)
h = pd.read_csv(D + "residual_history.csv").groupby("step").last()
r = np.log10(h["rms_ro"].values); r = r - r.mean()
sp = np.abs(np.fft.rfft(r)) ** 2; fr = np.fft.rfftfreq(r.size, 1)
k = np.argsort(sp[1:])[::-1][:3] + 1
print("rms_ro range", h["rms_ro"].min(), h["rms_ro"].max(), "dominant periods [steps]", [round(1 / fr[j], 1) for j in k], "power share", [round(sp[j] / sp[1:].sum(), 3) for j in k])
print("rms_ro step-to-step sign-flip fraction", np.mean(np.sign(np.diff(r)[1:]) != np.sign(np.diff(r)[:-1])))
