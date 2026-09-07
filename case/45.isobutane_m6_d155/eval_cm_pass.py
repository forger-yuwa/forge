"""固定 Euler 基準 δ_r 反復 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) の pass 評価。

指定 NS run について: NS/Euler 質量流量比、出口面 (x_F 近傍) の軸 M とコア質量流束重み平均 M、
軸 M − Euler 軸 M の x 分布、固定点差 (抽出 δ_r vs 与えた δ_r) を表にする。
使い方: design/.venv-opt/bin/python case/45.isobutane_m6_d155/eval_cm_pass.py run_0004_ns_v3 run_0019_ns_cm_pass1 ...
他 case: --case <case_dir> --euler <run名> --problem <yaml> を指定 (Md は prepare_info から)。
"""
import argparse, json, sys
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.probdef import load_problem
from forge_design.metrics.deltastar import deltastar_from_core_matched_euler

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--case", default=str(Path(__file__).resolve().parent))
ap.add_argument("--euler", default="run_0001_euler_shortest_dry")
ap.add_argument("--problem", default="problem_d155_ns.yaml")
ap.add_argument("--core-frac", type=float, default=0.6, help="出口コア平均の半径比")
args = ap.parse_args()
CASE = Path(args.case)
EULER = CASE / args.euler
prob = load_problem(CASE / args.problem); gas = prob.gas_model
Md = float(json.loads((EULER / "prepare_info.json").read_text())["Md"])


def load(run):
    rd = CASE / run; info = json.loads((rd / "prepare_info.json").read_text())
    ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]; S = info["scale_m"]
    with h5py.File(rd / "nozzle.h5") as f: nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda q: int(q.stem.split("_")[1]))[-1]
    with h5py.File(res) as f:
        ro = f["/VALUE/ro"][:].reshape(ni, nj); Ux = f["/VALUE/Ux"][:].reshape(ni, nj); Uy = f["/VALUE/Uy"][:].reshape(ni, nj)
        P = f["/VALUE/P"][:].reshape(ni, nj); T = f["/VALUE/T"][:].reshape(ni, nj)
    x = nc[:, 0].reshape(ni, nj) / S; r = nc[:, 1].reshape(ni, nj) / S
    gam = np.asarray(gas.gamma(T.ravel())).reshape(ni, nj)
    M = np.sqrt(Ux ** 2 + Uy ** 2) / np.sqrt(gam * P / ro)
    return dict(run=run, res=res.name, info=info, x=x, r=r, M=M, q=ro * Ux)


E = load(args.euler)
x_E = float(E["info"]["x_E"]); x_F = float(E["x"][-1, 0])
xs_axis = [round(v, 1) for v in np.concatenate([np.linspace(0.5 * x_E, x_E, 5), np.linspace(x_E, x_F - 1.0, 10)[1:]])]
for run in args.runs:
    N = load(run)
    ext = deltastar_from_core_matched_euler(CASE / run, EULER)
    mf = ext["massflow"]
    print(f"\n### {run} ({N['res']})  mdot_NS/mdot_E = {mf['mdot_ratio']:.4f}  δr_t,eff {mf['delta_r_throat_eff']:.4f} / applied {mf['delta_r_throat_given']:.4f}")
    # 出口面: 軸・コア (r ≤ 0.6 r_w, 質量流束重み)
    for xq in [round(v, 1) for v in np.linspace(x_E + 0.5 * (x_F - x_E), x_F - 1.0, 4)]:
        i = int(np.argmin(abs(N["x"][:, 0] - xq))); rw = N["r"][i, -1]; mc = N["r"][i] <= args.core_frac * rw
        w = N["q"][i] * N["r"][i]
        Mc = float(np.sum(N["M"][i][mc] * w[mc]) / np.sum(w[mc]))
        print(f"  x={N['x'][i,0]:5.1f}: axis M {N['M'][i,0]:.4f} ({(N['M'][i,0]-Md)/Md*100:+.2f} %)  core M {Mc:.4f} ({(Mc-Md)/Md*100:+.2f} %)")
    dev = [(np.interp(xq, N["x"][:, 0], N["M"][:, 0]) - np.interp(xq, E["x"][:, 0], E["M"][:, 0])) / Md * 100 for xq in xs_axis]
    print("  (axis − Euler)/Md [%]: " + " ".join(f"{xq}:{d:+.2f}" for xq, d in zip(xs_axis, dev)))
    print(f"  max|axis−Euler| in test section (x>x_E): {max(abs(d) for xq, d in zip(xs_axis, dev) if xq >= x_E):.2f} %")
    fin = np.isfinite(ext["delta_r_use"]) & (ext["delta_in"] > 1e-3)
    ratio = ext["delta_r_use"][fin] / ext["delta_in"][fin]
    print(f"  fixed point: use/in median {np.median(ratio):.3f}, p10/p90 {np.percentile(ratio,10):.3f}/{np.percentile(ratio,90):.3f}, "
          f"max|use−in| {np.nanmax(np.abs(ext['delta_r_use']-ext['delta_in'])):.4f} r_t; ok {ext['ok'].sum()}/{len(ext['ok'])} hard_ok {ext['hard_ok'].sum()}")
    xq_list = [0, 1, 2, 4, 8] + [round(v, 1) for v in np.linspace(0.5 * x_E, x_F - 1.0, 5)]
    print("  x: in/use: " + " ".join(f"{xq}:{np.interp(xq, ext['x'], ext['delta_in']):.4f}/{np.interp(xq, ext['x'], np.where(np.isfinite(ext['delta_r_use']), ext['delta_r_use'], np.nan)):.4f}" for xq in xq_list))
