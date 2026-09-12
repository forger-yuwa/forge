#!/usr/bin/env python3
"""forge (node, CPG, 素 SST) vs SU2 (axisym SST-V2003m) の同一メッシュ比較 — 等温壁チェーン S3 (plan §4.4)。

SU2 restart_flow.csv (RESTART_ASCII) を座標対応で forge のノード順に戻し、構造格子 (ni,nj) で同じ関数から
  ṁ (2π∫ρu r dr の中央値), 出口コア M, 壁 q_w(x) (3 点片側差分, λ_w=μ_w cp/Pr), Q_w=∫q_w 2πr ds, 壁 y1+,
  ステーション x/rt の δ99 (0.995 ρu), δ* (軸対称質量収支), θ
を出し forge/SU2 の比を表にする。
usage: design/.venv-opt/bin/python compare_nozzle_su2.py FORGE_RUN SU2_DIR [--stations 4,8,12,16,20] [--json out]
"""
import argparse, csv, json, sys
from pathlib import Path
import h5py, numpy as np
from scipy.spatial import cKDTree
CASE = Path(__file__).resolve().parent
GAM, CP, PR = 1.32752, 1190.2, 0.72; R = CP * (GAM - 1) / GAM


def mu_suth(T): return 1.716e-5 * (np.asarray(T, float) / 273.0) ** 1.5 * 384.0 / (np.asarray(T, float) + 111.0)


def d1(y, f):
    h1, h2 = y[1] - y[0], y[2] - y[0]
    return ((f[1] - f[0]) * h2 ** 2 - (f[2] - f[0]) * h1 ** 2) / (h1 * h2 * (h2 - h1))


def load_forge(rd):
    info = json.loads((rd / "prepare_info.json").read_text()); ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]; S = info["scale_m"]
    with h5py.File(rd / "nozzle.h5") as f: nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda q: int(q.stem.split("_")[1]))[-1]
    with h5py.File(res) as f:
        V = f["VALUE"]; g = lambda k: V[k][:].astype(float)
        ro, Ux, Uy, P, T = g("ro"), g("Ux"), g("Uy"), g("P"), g("T")
    return dict(nc=nc, ni=ni, nj=nj, S=S, xE=float(info["x_E"]), ro=ro, u=Ux, v=Uy, P=P, T=T, res=res.name)


def load_su2(su2_dir, nc):
    with open(Path(su2_dir) / "restart_flow.csv") as f:
        rdr = csv.reader(f); hdr = [h.strip().strip('"') for h in next(rdr)]
        data = np.array([[float(v) for v in row] for row in rdr if row])
    col = {h: i for i, h in enumerate(hdr)}
    x, y = data[:, col["x"]], data[:, col["y"]]
    tree = cKDTree(np.c_[x, y]); dist, idx = tree.query(nc[:, :2])   # forge ノード → SU2 行
    assert dist.max() < 1e-6 and len(set(idx.tolist())) == len(idx), dist.max()
    d = data[idx]
    ro = d[:, col["Density"]]; u = d[:, col["Momentum_x"]] / ro; v = d[:, col["Momentum_y"]] / ro
    if "Pressure" in col: P = d[:, col["Pressure"]]
    else: P = (GAM - 1) * (d[:, col["Energy"]] - 0.5 * ro * (u ** 2 + v ** 2))
    T = d[:, col["Temperature"]] if "Temperature" in col else P / (ro * R)
    return dict(ro=ro, u=u, v=v, P=P, T=T)


def metrics(F, D, stations):
    ni, nj, S = F["ni"], F["nj"], F["S"]
    x = F["nc"][:, 0].reshape(ni, nj); r = F["nc"][:, 1].reshape(ni, nj)
    ro, u, v, P, T = (D[k].reshape(ni, nj) for k in ("ro", "u", "v", "P", "T"))
    mu = mu_suth(T); a = np.sqrt(GAM * P / ro); M = np.hypot(u, v) / a
    xs = x[:, 0] / S
    md = 2 * np.pi * np.trapezoid(ro * u * r, r, axis=1); sel = (xs > -2.5) & (xs < F["xE"])
    xF = float(xs[-1]); iE = int(np.argmin(np.abs(xs - (xF - 1.0))))
    rr = r[iE]; q = (ro * u)[iE]; core = rr <= 0.6 * rr[-1]
    Mc = float(np.trapezoid(M[iE, core] * q[core] * rr[core], rr[core]) / np.trapezoid(q[core] * rr[core], rr[core]))
    xw, rw = x[:, -1], r[:, -1]; phi = np.arctan(np.gradient(rw, xw)); c = np.cos(phi)
    y1 = (rw - r[:, -2]) * c; y2 = (rw - r[:, -3]) * c
    ut = lambda j: u[:, j] * c + v[:, j] * np.sin(phi)
    dudn = np.array([d1([0, y1[i], y2[i]], [ut(-1)[i], ut(-2)[i], ut(-3)[i]]) for i in range(ni)])
    dTdn = np.array([d1([0, y1[i], y2[i]], [T[i, -1], T[i, -2], T[i, -3]]) for i in range(ni)])
    tau = mu[:, -1] * np.abs(dudn); yp = y1 * np.sqrt(ro[:, -1] * tau) / mu[:, -1]
    qw = mu[:, -1] * CP / PR * dTdn; ds = np.hypot(np.gradient(xw), np.gradient(rw))
    out = dict(mdot=float(np.median(md[sel])), M_core_exit=Mc, M_axis_exit=float(M[iE, 0]), Tw_exit=float(T[iE, -1]),
               Qw=float(np.sum(qw * 2 * np.pi * rw * ds)), qw_peak=float(qw.max()), qw_peak_x=float(xs[np.argmax(qw)]),
               y1p_throat=float(yp[np.argmin(np.abs(xs))]), y1p_max=float(yp.max()), tau_exit=float(tau[iE]), qw_exit=float(qw[iE]))
    st = {}
    for xq in stations:
        i = int(np.argmin(np.abs(xs - xq))); rq = r[i]; qq = (ro * u)[i]; yy = (rq[-1] - rq)[::-1]; qy = qq[::-1]   # 壁 → 軸
        win = yy <= 0.4 * rq[-1]; jm = int(np.argmax(qy[win])); je = int(np.argmax(qy[win] >= 0.995 * qy[win][jm]))
        qe = qy[je]; sl = slice(0, je + 1); rwq = rq[-1]
        I = np.trapezoid((1 - qy[sl] / qe) * (1 - yy[sl] / rwq), yy[sl])
        dstar = rwq * (1 - np.sqrt(max(1 - 2 * I / rwq, 0)))
        uu = np.hypot(u, v)[i][::-1]; ue = uu[je]
        theta = np.trapezoid(qy[sl] / qe * (1 - uu[sl] / ue) * (1 - yy[sl] / rwq), yy[sl])
        st[xq] = dict(d99=float(yy[je] / S), dstar=float(dstar / S), theta=float(theta / S), qw=float(qw[i]), Cf_tau=float(tau[i]), yp=float(yp[i]))
    out["stations"] = st
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("forge_run"); ap.add_argument("su2_dir")
    ap.add_argument("--stations", default="4,8,12,16,20"); ap.add_argument("--json", default=None)
    a = ap.parse_args(); stations = [float(s) for s in a.stations.split(",")]
    F = load_forge(CASE / a.forge_run); mf = metrics(F, F, stations)
    S = load_su2(CASE / a.su2_dir, F["nc"]); ms = metrics(F, S, stations)
    print(f"=== {a.forge_run} ({F['res']}) vs {a.su2_dir} ===")
    for k in ("mdot", "M_core_exit", "M_axis_exit", "Tw_exit", "Qw", "qw_peak", "qw_peak_x", "y1p_throat", "y1p_max"):
        f_, s_ = mf[k], ms[k]; print(f"  {k:12s} forge {f_:12.6g}  SU2 {s_:12.6g}  ratio {f_/s_ if s_ else float('nan'):.4f}")
    print("  station x/rt: d99 f/s | dstar f/s | theta f/s | tau_w f/s | qw f/s | y1+ f / s")
    for xq in stations:
        f_, s_ = mf["stations"][xq], ms["stations"][xq]
        rat = lambda k: f_[k] / s_[k] if s_[k] else float("nan")
        print(f"   {xq:5.1f}: {f_['d99']:.4f} {rat('d99'):.3f} | {f_['dstar']:.5f} {rat('dstar'):.3f} | {f_['theta']:.5f} {rat('theta'):.3f} | {rat('Cf_tau'):.3f} | {f_['qw']/1e3:8.1f} kW/m2 {rat('qw'):.3f} | {f_['yp']:.2f} / {s_['yp']:.2f}")
    if a.json: Path(a.json).write_text(json.dumps(dict(forge=mf, su2=ms), indent=1, default=float))


if __name__ == "__main__":
    main()
