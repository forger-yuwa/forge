import json, sys, os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from forge_design.probdef import load_problem
from forge_design.evaluate import runner_axismach as RA
from forge_design.geometry import moc_inverse as MI
rec = {}
orig = MI._axis_grid
def spy(x0, x_end, n, dx0=None, *a, **k):
    rec.update(x0=float(x0), x_end=float(x_end), n=int(n), dx0=dx0)
    return orig(x0, x_end, n, dx0, *a, **k)
MI._axis_grid = spy
p = load_problem(sys.argv[1]); p.geometry["Lc_mode"] = "explicit"; p.geometry["n_axis_inv"] = 2400
L0 = 37.2020873623039; p.dv["L_c"] = {"value": L0}
d = RA.design_chain(p)
xE0 = float(d["law"].x_E); off = rec["x_end"] - xE0; xA = xE0 - L0
print("x0", rec["x0"], "x_end", rec["x_end"], "x_E", xE0, "x_end - x_E", off, "dx0", rec["dx0"])
fine = json.load(open(sys.argv[2]))["rows"]
L = np.array([r["L_c"] for r in fine]); f = np.array([r["f"] for r in fine])
# 各 L_c で x_E より下流で最初の軸の節点の番号と、x_E − (直前の節点)
idx, frac = [], []
for Lc in L:
    xE = xA + Lc
    g = orig(rec["x0"], xE + off, rec["n"], rec["dx0"])
    j = int(np.searchsorted(g, xE)); idx.append(j); frac.append((xE - g[j - 1]) / (g[j] - g[j - 1]))
idx = np.array(idx); frac = np.array(frac)
jf = np.flatnonzero(np.abs(np.diff(f) - 0.292 * np.diff(L)) > 1e-3)
ji = np.flatnonzero(np.diff(idx) != 0)
print("f の段の位置 (L_c の区間):", [(round(L[i], 5), round(L[i + 1], 5), round(f[i + 1] - f[i], 5)) for i in jf])
print("x_E が軸の節点をまたぐ区間:", [(round(L[i], 5), round(L[i + 1], 5)) for i in ji])
g = orig(rec["x0"], rec["x_end"], rec["n"], rec["dx0"]); j = int(np.searchsorted(g, xE0))
print("x_E 付近の軸の間隔", g[j] - g[j - 1], "節点の動き d(x_node)/dL_c (数値)",
      (orig(rec["x0"], rec["x_end"] + 1e-3, rec["n"], rec["dx0"])[j] - g[j]) / 1e-3)
