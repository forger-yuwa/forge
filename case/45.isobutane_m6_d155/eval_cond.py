"""凝縮 ON run の評価: 軸上 onset (g_0 > 1e-4 になる x)、出口の g_0 (軸・最大)、S 最大、出口面の軸 M / コア M と dry 比較。
usage: design/.venv-opt/bin/python eval_cond.py <cond_run> <dry_run> [--euler run]"""
import sys, json, argparse
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.probdef import load_problem
CASE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(); ap.add_argument("cond_run"); ap.add_argument("dry_run"); ap.add_argument("--problem", default="problem_d155_ns.yaml")
a = ap.parse_args(); gas = load_problem(CASE / a.problem).gas_model; Md = 6.0
def load(run):
    rd = CASE / run; info = json.loads((rd / "prepare_info.json").read_text()); ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]; S = info["scale_m"]
    with h5py.File(rd / "nozzle.h5") as f: nc = f["/MESH/COORD"][:].reshape(-1, 3)
    ress = sorted(rd.glob("res_[0-9]*.h5"), key=lambda q: int(q.stem.split("_")[1]))
    out = []
    for res in ress:
        with h5py.File(res) as f:
            ro = f["/VALUE/ro"][:].reshape(ni, nj); Ux = f["/VALUE/Ux"][:].reshape(ni, nj); P = f["/VALUE/P"][:].reshape(ni, nj); T = f["/VALUE/T"][:].reshape(ni, nj)
            g = f["/VALUE/g_0"][:].reshape(ni, nj) if "/VALUE/g_0" in f else np.zeros((ni, nj))
            Smax = f["/VALUE/condS_0"][:].reshape(ni, nj) if "/VALUE/condS_0" in f else np.zeros((ni, nj))
        gam = np.asarray(gas.gamma(T.ravel())).reshape(ni, nj); M = Ux / np.sqrt(gam * P / ro)
        out.append(dict(res=res.name, M=M, g=g, S=Smax, ro=ro, Ux=Ux, T=T))
    return dict(x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S, S=S, snaps=out)
Cn = load(a.cond_run); Dr = load(a.dry_run)
x = Cn["x"][:, 0]; xF = x[-1]
def exit_vals(D, snap, xq):
    i = int(np.argmin(abs(D["x"][:, 0] - xq))); rw = D["r"][i, -1]; mc = D["r"][i] <= 0.6 * rw
    w = snap["ro"][i] * snap["Ux"][i] * D["r"][i]
    return float(snap["M"][i, 0]), float(np.sum(snap["M"][i][mc] * w[mc]) / np.sum(w[mc])), float(snap["g"][i, 0]), float(snap["g"][i].max()), float(snap["T"][i, 0])
print(f"### {a.cond_run} (dry ref {a.dry_run})  r_t {Cn['S']*1e3:.2f} mm")
for s in Cn["snaps"]:
    gax = s["g"][:, 0]; on = np.where(gax > 1e-4)[0]; x_on = x[on[0]] if len(on) else float("nan")
    ax, core, gax_e, gmax_e, Tax = exit_vals(Cn, s, xF - 1.0)
    print(f"  {s['res']:14s} onset(axis g>1e-4) x={x_on:6.1f} r_t | 出口 x={xF-1:.1f}: axis M {ax:.4f} ({(ax-Md)/Md*100:+.2f} %) core M {core:.4f} ({(core-Md)/Md*100:+.2f} %) g_axis {gax_e*100:.3f} % g_max {gmax_e*100:.3f} % T_axis {Tax:.1f} K | S_max(domain) {s['S'].max():.2f}")
d = Dr["snaps"][-1]; ax, core, *_ = exit_vals(Dr, d, xF - 1.0)
print(f"  dry ({d['res']}): 出口 axis M {ax:.4f} core M {core:.4f}")
# 試験部の軸 M 分布 (dry vs cond)
print("  軸 M (cond − dry)/Md [%] along test section: " + " ".join(f"{xq}:{(np.interp(xq, x, Cn['snaps'][-1]['M'][:,0]) - np.interp(xq, Dr['x'][:,0], d['M'][:,0]))/Md*100:+.2f}" for xq in [40, 50, 60, 70, 80, 90, 94]))
