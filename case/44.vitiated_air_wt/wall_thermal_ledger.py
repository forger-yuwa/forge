#!/usr/bin/env python3
"""ノズル NS run の壁温感度台帳 (plan tooling-nozzle-isothermal-wall-chain §4.5 / §4.6-5)。

各 run (node 構造格子 ni×nj, 壁 = j=nj-1) について:
  - 実測 y1+(x) (壁物性: ρ_w, μ_w, τ_w = μ_w ∂u_t/∂n, 2 次片側差分)、max とスロート/出口の値
  - q_w(x) = λ_w ∂T/∂n cos φ_w (壁へ入る向き正), Q_w = ∫ q_w 2π r_w ds [W], q_w ピーク位置
  - NS/Euler 質量流量比 (metrics.deltastar.massflow_ratio), 出口面 (x_E 近傍) のコア質量流束重み M・軸 M
  - 出口面の全温・全圧分布 (tools/total_quantities.py の API; TP 対応) の質量流束重み平均と壁近傍値
usage: design/.venv-opt/bin/python wall_thermal_ledger.py RUN [RUN ...] [--euler run_0005_va_R2_LU6_Lc8] [--core-frac 0.6] [--json out.json]
"""
import argparse, json, sys
from pathlib import Path
import h5py, numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "design")); sys.path.insert(0, str(ROOT / "solver_density_cuda/tools"))
from forge_design.metrics.deltastar import massflow_ratio
from total_quantities import total_state
CASE = Path(__file__).resolve().parent


def mu_suth(T):
    return 1.716e-5 * (np.asarray(T, float) / 273.0) ** 1.5 * 384.0 / (np.asarray(T, float) + 111.0)


def d1(y, f):
    """y[0] = 壁, 非一様 3 点片側差分 (壁向き +)。"""
    h1, h2 = y[1] - y[0], y[2] - y[0]
    return ((f[1] - f[0]) * h2 ** 2 - (f[2] - f[0]) * h1 ** 2) / (h1 * h2 * (h2 - h1))


def ledger(run, euler, core_frac):
    rd = CASE / run; info = json.loads((rd / "prepare_info.json").read_text())
    ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]; S = info["scale_m"]
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda q: int(q.stem.split("_")[1]))[-1]
    with h5py.File(rd / "nozzle.h5") as f: nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res) as f:
        V = f["VALUE"]; g = lambda k: V[k][:].astype(float).reshape(ni, nj)
        ro, Ux, Uy, P, T = g("ro"), g("Ux"), g("Uy"), g("P"), g("T")
        mu = g("vis_lam") if "vis_lam" in V else mu_suth(T); son = g("sonic") if "sonic" in V else None
        cp_arr = g("cp") if "cp" in V else None
    x = nc[:, 0].reshape(ni, nj); r = nc[:, 1].reshape(ni, nj)
    # 壁 (j = nj-1) と法線: 壁角 φ_w = atan(dr_w/dx)
    xw, rw = x[:, -1], r[:, -1]
    phi = np.arctan(np.gradient(rw, xw)); cphi = np.cos(phi)
    # 壁法線に沿う 3 点: j = nj-1, nj-2, nj-3 (r 方向)。距離は r 差 × cos φ (法線成分)
    y0 = np.zeros(ni); y1 = (rw - r[:, -2]) * cphi; y2 = (rw - r[:, -3]) * cphi
    ut = lambda j: (Ux[:, j] * cphi + Uy[:, j] * np.sin(phi))       # 壁接線速度
    dudn = np.array([d1([y0[i], y1[i], y2[i]], [ut(-1)[i], ut(-2)[i], ut(-3)[i]]) for i in range(ni)])
    dTdn = np.array([d1([y0[i], y1[i], y2[i]], [T[i, -1], T[i, -2], T[i, -3]]) for i in range(ni)])
    muw = mu[:, -1]; row = ro[:, -1]; tau = muw * np.abs(dudn); utau = np.sqrt(np.maximum(tau, 0) / row)
    yplus1 = y1 * row * utau / muw
    # λ_w = μ_w cp / Pr (定 Pr 0.72; cp は TP なら場の cp、CPG は問題定義から)
    cpw = cp_arr[:, -1] if cp_arr is not None else float(json.loads((rd / "prepare_info.json").read_text()).get("cp", 0) or 0) or None
    if cpw is None:
        import yaml
        cfg = yaml.safe_load((rd / "solverConfig.yaml").read_text()); cpw = float(cfg["physProp"]["cp"])
    lam_w = muw * cpw / 0.72
    qw = lam_w * dTdn                                   # 壁へ入る向き正
    ds = np.hypot(np.gradient(xw), np.gradient(rw))
    Qw = float(np.sum(qw * 2 * np.pi * rw * ds))
    xs = xw / S
    # 出口面 (x_E 近傍): コア質量流束重み M・軸 M・全温/全圧
    # 出口面 = メッシュ末端 x_F の 1 r_t 手前 (試験部の一様域; eval_cm_pass.py と同じ)。x_E は膨張終端
    xF = float(xs[-1]); iE = int(np.argmin(np.abs(xs - (xF - 1.0))))
    if son is None:
        gam = 1.4; a = np.sqrt(gam * P / ro)
    else:
        a = son
    M = np.hypot(Ux, Uy) / a
    rr = r[iE, :]; q = (ro * Ux)[iE, :]; core = rr <= core_frac * rr[-1]
    Mcore = float(np.trapezoid(M[iE, core] * q[core] * rr[core], rr[core]) / np.trapezoid(q[core] * rr[core], rr[core]))
    try:
        ts = total_state(str(rd), str(res)); T0 = np.asarray(ts["T0"]).reshape(ni, nj); P0 = np.asarray(ts["P0"]).reshape(ni, nj); t0src = ts.get("method", "h0")
    except BaseException as e:  # 旧バイナリの res (h0 無し): CPG 近似 T + u²/2cp で代用 (台帳に明記)
        cpf = cp_arr if cp_arr is not None else cpw
        T0 = T + 0.5 * (Ux ** 2 + Uy ** 2) / cpf; gam_l = (a ** 2 * ro / P); P0 = P * (T0 / T) ** (gam_l / (gam_l - 1)); t0src = f"fallback T+u2/2cp ({str(e)[:40]})"
    T0core = float(np.trapezoid(T0[iE, core] * q[core] * rr[core], rr[core]) / np.trapezoid(q[core] * rr[core], rr[core]))
    P0core = float(np.trapezoid(P0[iE, core] * q[core] * rr[core], rr[core]) / np.trapezoid(q[core] * rr[core], rr[core]))
    T0all = float(np.trapezoid(T0[iE, :] * q * rr, rr) / np.trapezoid(q * rr, rr))
    mr = massflow_ratio(rd, CASE / euler)
    it = int(np.argmin(np.abs(xs)))
    out = dict(run=run, res=res.name, wall_thermal=info.get("wall_thermal"), mesh=info["mesh"],
               y1p_max=float(yplus1.max()), y1p_at=float(xs[np.argmax(yplus1)]), y1p_throat=float(yplus1[it]), y1p_exit=float(yplus1[iE]),
               Tw_throat=float(T[it, -1]), Tw_exit=float(T[iE, -1]),
               Qw_W=Qw, qw_peak_W_m2=float(qw.max()), qw_peak_x=float(xs[np.argmax(qw)]), qw_throat=float(qw[it]), qw_exit=float(qw[iE]),
               mdot_ratio=mr["mdot_ratio"], M_axis_exit=float(M[iE, 0]), M_core_exit=Mcore, x_exit_eval=float(xs[iE]),
               T0_core_exit=T0core, P0_core_exit=P0core, T0_massavg_exit=T0all, T0_wall_exit=float(T0[iE, -1]), T0_source=t0src,
               delta_r_exit_given=float(info.get("initializer", {}).get("delta_r_exit", np.nan)) if info.get("initializer") else None)
    np.savetxt(rd / "wall_thermal_profile.csv", np.c_[xs, rw / S, T[:, -1], yplus1, tau, qw], delimiter=",", comments="",
               header="x_rt,r_w_rt,Tw,y1plus,tau_w,q_w_W_m2")
    return out


def series(run, core_frac):
    """全 res_*.h5 で出口コア M・軸 M・ṁ (NS 自身の中央値)・Q_w を出し、末尾 50 % の drift を判定。"""
    rd = CASE / run; info = json.loads((rd / "prepare_info.json").read_text())
    ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]; S = info["scale_m"]
    with h5py.File(rd / "nozzle.h5") as f: nc = f["/MESH/COORD"][:].reshape(-1, 3)
    x = nc[:, 0].reshape(ni, nj) / S; r = nc[:, 1].reshape(ni, nj)
    fs = sorted(rd.glob("res_[0-9]*.h5"), key=lambda q: int(q.stem.split("_")[1]))
    fs = [f for f in fs if int(f.stem.split("_")[1]) > 0]
    xF = float(x[-1, 0]); iE = int(np.argmin(np.abs(x[:, 0] - (xF - 1.0))))
    rows = []
    for fp in fs:
        with h5py.File(fp) as f:
            V = f["VALUE"]; g = lambda k: V[k][:].astype(float).reshape(ni, nj)
            ro, Ux, Uy, T = g("ro"), g("Ux"), g("Uy"), g("T"); a = g("sonic") if "sonic" in V else np.sqrt(1.4 * g("P") / ro)
            mu = g("vis_lam") if "vis_lam" in V else mu_suth(T); cpa = g("cp") if "cp" in V else None
        M = np.hypot(Ux, Uy) / a; rr = r[iE, :]; q = (ro * Ux)[iE, :]; core = rr <= core_frac * rr[-1]
        Mc = float(np.trapezoid(M[iE, core] * q[core] * rr[core], rr[core]) / np.trapezoid(q[core] * rr[core], rr[core]))
        md = 2 * np.pi * np.trapezoid(ro * Ux * r, r, axis=1); sel = (x[:, 0] > -2.5) & (x[:, 0] < float(info["x_E"]))
        xw, rw = x[:, -1] * S, r[:, -1]; phi = np.arctan(np.gradient(rw, xw)); cphi = np.cos(phi)
        y1 = (rw - r[:, -2]) * cphi; y2 = (rw - r[:, -3]) * cphi
        dTdn = np.array([d1([0.0, y1[i], y2[i]], [T[i, -1], T[i, -2], T[i, -3]]) for i in range(ni)])
        cpw = cpa[:, -1] if cpa is not None else 1190.2
        qw = mu[:, -1] * cpw / 0.72 * dTdn; ds = np.hypot(np.gradient(xw), np.gradient(rw))
        rows.append((int(fp.stem.split("_")[1]), Mc, float(M[iE, 0]), float(np.median(md[sel])), float(np.sum(qw * 2 * np.pi * rw * ds))))
    rows = np.array(rows); print(f"\n=== series {run}: step / M_core_exit / M_axis_exit / mdot / Q_w ===")
    for rw_ in rows: print(f"  {int(rw_[0]):7d} {rw_[1]:.5f} {rw_[2]:.5f} {rw_[3]:.6g} {rw_[4]:.6g}")
    tol = dict(M_core_exit=0.0005, M_axis_exit=0.001, mdot=0.001, Q_w=0.01); ok = True
    for j, k in enumerate(("M_core_exit", "M_axis_exit", "mdot", "Q_w"), start=1):
        t = rows[len(rows) // 2:, j]; drift = (t.max() - t.min()) / max(abs(t.mean()), 1e-30); ok &= drift <= tol[k]
        print(f"  {k}: tail drift {drift*100:.3f} % ({'STEADY' if drift <= tol[k] else 'DRIFTING'}, tol {tol[k]*100:.2f} %)")
    print("  SERIES VERDICT:", "STEADY" if ok else "NOT STEADY")
    return bool(ok)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--euler", default="run_0005_va_R2_LU6_Lc8")
    ap.add_argument("--series", action="store_true")
    ap.add_argument("--core-frac", type=float, default=0.6); ap.add_argument("--json", default=None)
    a = ap.parse_args()
    rows = [ledger(r, a.euler, a.core_frac) for r in a.runs]
    keys = ("y1p_max", "y1p_throat", "y1p_exit", "Tw_throat", "Tw_exit", "Qw_W", "qw_peak_W_m2", "qw_peak_x", "mdot_ratio", "M_axis_exit", "M_core_exit", "T0_core_exit", "P0_core_exit", "T0_massavg_exit")
    print("quantity".ljust(18) + "".join(r["run"][:34].rjust(36) for r in rows))
    for k in keys:
        print(k.ljust(18) + "".join(f"{r[k]:36.6g}" for r in rows))
    if a.series:
        for r in a.runs: series(r, a.core_frac)
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=1, default=float))


if __name__ == "__main__":
    main()
