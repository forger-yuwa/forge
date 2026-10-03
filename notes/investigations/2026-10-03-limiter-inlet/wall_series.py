# limiter-inlet-column-oscillation §4.3 report quantities (dry): upper-wall p/p0 (mean over the Fig.3 x>=10 mm points, x=42, 52 mm) and
# upper-wall T mean, per res_*.h5 -> CSV for check_quasisteady --series-csv. Same wall-node set (wall_dist == 0, y > 0) for every run.
import h5py, numpy as np, glob, sys, csv
C = "/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
run = sys.argv[1]; P0 = 59070.0
xe = []
for r in csv.DictReader(open(C + "wyslouzil_fig3_pp0.csv")):
    k = [c for c in r if c.lower().startswith("x")][0]; xv = float(r[k]); xe.append(xv)
xe = np.array(xe) * 10.0; xe = xe[xe >= 10.0]   # CSV の x は cm
fs = sorted(glob.glob(C + run + "/res_*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))
fs = [f for f in fs if not f.endswith("/res_0.h5")]
w = csv.writer(open(C + run + "/wall_series.csv", "w")); w.writerow(["step", "pw_mean_x10", "pw42", "pw52", "Tw_mean_x10_K"])
for f in fs:
    h = h5py.File(f); V = h["VALUE"]; xyz = np.array(h["MESH/COORD"][...]).reshape(-1, 3)
    wd = np.array(V["wall_dist"][...]).ravel(); m = (wd == 0) & (xyz[:, 1] > 0)
    o = np.argsort(xyz[m, 0]); xw = xyz[m, 0][o] * 1e3
    pw = np.array(V["P"][...]).ravel()[m][o] / P0; tw = np.array(V["T"][...]).ravel()[m][o]
    w.writerow([int(f.split("_")[-1][:-3]), np.interp(xe, xw, pw).mean(), np.interp(42, xw, pw), np.interp(52, xw, pw), np.interp(xe, xw, tw).mean()])
print(run, len(fs), "rows, x points", xe.size)
