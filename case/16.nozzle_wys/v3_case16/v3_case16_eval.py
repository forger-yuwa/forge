"""#10 case/16 V3 評価 (plan §6 V3 case/16 条件, 2026-09-25 固定): run_0506 (flag0) / run_0507 (flag1) 分岐。
全 res_*.h5 に extract_wall_pp0.py を当て、3 点 p/p0・出口中心線 M・flag 差 L∞ (x∈[10,94] mm) の系列 CSV を書く。"""
import subprocess, sys, glob, re, numpy as np, h5py
from pathlib import Path
C = Path(__file__).resolve().parents[1]; OUT = Path(__file__).resolve().parent
R0, R1 = C/"run_0506_node_lam_flag0_br", C/"run_0507_node_lam_flag1_br"
XS = (16.4, 45.6, 85.0)
def steps(r): return sorted(int(re.search(r"res_(\d+)\.h5", f).group(1)) for f in glob.glob(str(r/"res_*.h5")) if re.search(r"res_\d+\.h5$", f))
def wall(r, n):
    o = OUT/f"{r.name}_w{n}.csv"
    if not o.exists(): subprocess.run([sys.executable, str(C/"extract_wall_pp0.py"), str(r/f"res_{n}.h5"), str(o)], check=True, capture_output=True)
    return np.genfromtxt(o, delimiter=",", names=True)
def mexit(r, n):
    with h5py.File(r/f"res_{n}.h5") as f:
        V = f["VALUE"]; P, ro, ux, uy = (V[k][:] for k in ("P","ro","Ux","Uy"))
    with h5py.File(r/"nozzle_user_2d.h5") as g:
        X = g["MESH/COORD"][:].reshape(-1, 3) if g["MESH/COORD"].ndim == 1 else g["MESH/COORD"][:]
    x, y = X[:len(P),0], X[:len(P),1]; xm = x.max(); col = x > xm - 1e-6
    i = np.where(col)[0][np.argmin(np.abs(y[col] - 0.5*(y[col].min()+y[col].max())))]
    return float(np.hypot(ux[i], uy[i]) / np.sqrt(1.4*P[i]/ro[i]))   # γ=1.4 近似 (判定閾 1.7 に対し十分)
S = sorted(set(steps(R0)) & set(steps(R1)))
rows = []
for n in S:
    a, b = wall(R0, n), wall(R1, n)
    m = (a["x_mm"] >= 10) & (a["x_mm"] <= 94); assert np.allclose(a["x_mm"], b["x_mm"])
    d = np.abs(b["contour_wall"][m] - a["contour_wall"][m]) / a["contour_wall"][m]
    rows.append([n] + [np.interp(x, a["x_mm"], a["contour_wall"]) for x in XS] + [np.interp(x, b["x_mm"], b["contour_wall"]) for x in XS]
                + [d.max(), mexit(R0, n), mexit(R1, n)])
hdr = "step,f0_x16,f0_x46,f0_x85,f1_x16,f1_x46,f1_x85,dLinf,Mexit0,Mexit1"
np.savetxt(OUT/"v3_series.csv", np.array(rows), delimiter=",", header=hdr, comments="", fmt="%.8g")
A = np.array(rows); tail = A[A[:,0] >= A[-1,0]*0.6]
print("dumps", len(S), "tail (last 40%) steps", int(tail[0,0]), "-", int(tail[-1,0]))
for j, nm in enumerate(hdr.split(",")[1:], 1): print(f"{nm:8s} tail mean {tail[:,j].mean():.6g}  last {A[-1,j]:.6g}")
# (b) 末尾平均の場で L∞ も出す
Wa = np.mean([wall(R0, int(n))["contour_wall"] for n in tail[:,0]], axis=0); Wb = np.mean([wall(R1, int(n))["contour_wall"] for n in tail[:,0]], axis=0)
x = wall(R0, int(tail[0,0]))["x_mm"]; m = (x >= 10) & (x <= 94)
dd = np.abs(Wb[m]-Wa[m])/Wa[m]; print(f"(b) tail-mean field Linf = {100*dd.max():.4f} % at x={x[m][dd.argmax()]:.2f} mm")
