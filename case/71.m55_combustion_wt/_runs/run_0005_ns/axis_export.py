import json, h5py, numpy as np
C = "/home/ubuntu/forge-cardmesh/case/71.m55_combustion_wt/"
for run, res in (("run_0003_ns_cont", "res_200000.h5"), ("run_0005_ns", "res_200000.h5"), ("run_0004_euler", "res_12000.h5"), ("run_0001_euler", "res_12000.h5")):
    info = json.load(open(C + run + "/prepare_info.json")); ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    with h5py.File(C + run + "/nozzle.h5") as g: X = g["MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(C + run + "/" + res) as f:
        V = f["VALUE"]; M = np.hypot(V["Ux"][:], V["Uy"][:]) / V["sonic"][:]
    cols = []
    for j in (0, 9, 20):
        idx = np.arange(ni) * nj + j; cols += [X[idx, 0], X[idx, 1], M[idx]]
    np.savetxt(f"/home/ubuntu/m55_diag_0005/axisM_{run}.csv", np.c_[tuple(cols)], delimiter=",", header="x0,r0,M0,x9,r9,M9,x20,r20,M20", comments="")
print("ok")
