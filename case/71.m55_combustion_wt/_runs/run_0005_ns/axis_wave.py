# 軸 (j = 0) と j = 9 の M(x) を試験部で比べる (構造格子 index = i·nj + j)
import sys, json, h5py, numpy as np
C = "/home/ubuntu/forge-cardmesh/case/71.m55_combustion_wt/"
Md = 5.5
for run, res in (("run_0003_ns_cont", "res_200000.h5"), ("run_0005_ns", "res_200000.h5"), ("run_0005_ns", "res_100000.h5"),
                 ("run_0001_euler", "res_12000.h5"), ("run_0004_euler", "res_12000.h5")):
    info = json.load(open(C + run + "/prepare_info.json")); ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]); S = float(info["scale_m"])
    with h5py.File(C + run + "/nozzle.h5") as g: X = g["MESH/COORD"][:].reshape(-1, 3)
    try:
        with h5py.File(C + run + "/" + res) as f:
            V = f["VALUE"]; M = np.hypot(V["Ux"][:], V["Uy"][:]) / V["sonic"][:]
    except OSError as e:
        print(run, res, "なし", e); continue
    out = []
    for j in (0, 9):
        idx = np.arange(ni) * nj + j; x = X[idx, 0] / S; m = M[idx]
        sel = (x >= 40) & (x <= 79.5)
        dev = (m[sel] / Md - 1) * 100
        k = np.argmax(np.abs(dev))
        out.append(f"j{j} r/r_w≈{X[idx[-1],1]/X[(ni-1)*nj+nj-1,1]:.2f}: dev [{dev.min():+.3f}, {dev.max():+.3f}] % (max@x {x[sel][k]:.1f})")
    # 波の位置: x 40–80 での極値
    idx = np.arange(ni) * nj; x = X[idx, 0] / S; m = M[idx]; sel = (x >= 38) & (x <= 79.8)
    xs, ms = x[sel], m[sel]; d = np.diff(ms); ext = np.where(np.sign(d[1:]) != np.sign(d[:-1]))[0] + 1
    print(f"{run}/{res}: " + " | ".join(out) + " | 軸の極値 x " + ", ".join(f"{xs[e]:.1f}({(ms[e]/Md-1)*100:+.2f}%)" for e in ext[:8]))
