#!/usr/bin/env python3
"""case/48 超音速冷却壁平板の評価 (plan tooling-nozzle-isothermal-wall-chain §4.3)。

forge node run (res_*.h5) と SU2 (restart_flow.csv) から、各 x ステーションで
  C_f (壁勾配 2 次片側差分), q_w = λ_w dT/dn (壁へ入る向き正), δ*, θ, H (ρu 0.995 縁), T_w, M_e, Re_θ, y1+
を取り、理論/経験式と比較する:
  - van Driest II (非圧縮基準 Kármán–Schoenherr) の C_f(Re_θ; M_e, T_w/T_aw)
  - Reynolds アナロジー係数 2St/C_f (St = q_w/[ρ_e u_e c_p (T_aw−T_w)])
  - CONTUR 平面積分 (deltastar_integral.flat_plate_integral, θ0 = 最初のステーションの forge θ) の δ*, θ, H
  - 温度–速度関係 (Walz 形 / Duan–Martín 形 C_T=0.8259) との RMS 差 (診断)
  - van Driest 変換速度分布と対数則 (診断)
  - 平板のエネルギー閉合: ∫q_w dx vs 入口・出口の全エンタルピー流束差 (単位幅)
usage:
  cooled_plate_eval.py RUN_DIR [--su2 SU2_DIR] [--stations 0.3,0.6,0.9] [--out PREFIX] [--series] [--json]
  --series: 全 res_*.h5 の時系列で C_f, q_w, δ*, θ, T_w の末尾 50 % drift を判定 (STEADY 公差: δ*/θ 1 %, q_w/C_f 1.5 %)
"""
import argparse, csv, glob, json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

GAM, CP, PR, PRT = 1.4, 1004.5, 0.72, 0.9
R = CP * (GAM - 1) / GAM
M_INF, P_INF, T_INF = 4.19, 5037.4, 283.0
RO_INF = P_INF / (R * T_INF); U_INF = M_INF * math.sqrt(GAM * R * T_INF)
R_REC = PR ** (1.0 / 3.0)


def mu_suth(T):
    T = np.asarray(T, float)
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


# ---------------------------------------------------------------- 読み込み
def load_forge(res_h5):
    import h5py
    with h5py.File(res_h5, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        V = f["VALUE"]
        d = dict(x=c[:, 0], y=c[:, 1], ro=V["ro"][:].astype(float), u=V["Ux"][:].astype(float), v=V["Uy"][:].astype(float),
                 T=V["T"][:].astype(float), P=V["P"][:].astype(float), wd=V["wall_dist"][:].astype(float))
        d["mu"] = V["vis_lam"][:].astype(float) if "vis_lam" in V else mu_suth(d["T"])
        for k in ("k", "omega", "vis_turb"):
            if k in V: d[k] = V[k][:].astype(float)
    return d


def load_su2(su2_dir):
    path = Path(su2_dir) / "restart_flow.csv"
    with open(path) as f:
        rdr = csv.reader(f); hdr = [h.strip().strip('"') for h in next(rdr)]
        data = np.array([[float(v) for v in row] for row in rdr if row])
    col = {h: i for i, h in enumerate(hdr)}
    x, y = data[:, col["x"]], data[:, col["y"]]
    ro = data[:, col["Density"]]
    if "Velocity_x" in col:
        u, v = data[:, col["Velocity_x"]], data[:, col["Velocity_y"]]
    else:
        u, v = data[:, col["Momentum_x"]] / ro, data[:, col["Momentum_y"]] / ro
    T = data[:, col["Temperature"]] if "Temperature" in col else None
    P = data[:, col["Pressure"]] if "Pressure" in col else None
    if P is None:
        E = data[:, col["Energy"]] / ro; P = (GAM - 1) * ro * (E - 0.5 * (u ** 2 + v ** 2))
    if T is None: T = P / (ro * R)
    mu = data[:, col["Laminar_Viscosity"]] if "Laminar_Viscosity" in col else mu_suth(T)
    d = dict(x=x, y=y, ro=ro, u=u, v=v, T=T, P=P, mu=mu, wd=np.full(len(x), np.nan))
    # 壁距離: 平板 (y=0, x>=0) からの距離で代用 (ステーション抽出だけに使う)
    d["wd"] = np.where(x >= -1e-9, y, np.hypot(x, y))
    return d


# ---------------------------------------------------------------- ステーション量
def deriv_wall(y, f):
    """y[0]=壁, 非一様 3 点片側差分で df/dy|_wall。"""
    h1, h2 = y[1] - y[0], y[2] - y[0]
    return ((f[1] - f[0]) * h2 ** 2 - (f[2] - f[0]) * h1 ** 2) / (h1 * h2 * (h2 - h1))


def station(D, xs, ywin=None, edge_frac=0.995, tol=2e-4):
    """x=xs の壁法線カラム (構造格子: x が一致するノード列) から壁量と積分厚さ。"""
    x, y = D["x"], D["y"]
    # 壁ノード列の x の中で最も近いものに丸める
    wallx = np.unique(np.round(x[(np.abs(y) < 1e-9) & (x >= -1e-9)], 7))
    xv = wallx[np.argmin(np.abs(wallx - xs))]
    col = np.where(np.abs(x - xv) < 1e-7)[0]
    col = col[np.argsort(y[col])]
    yy = y[col]; u = D["u"][col]; ro = D["ro"][col]; T = D["T"][col]; P = D["P"][col]; mu = D["mu"][col]
    if ywin is None:
        ywin = 0.6 * (xv + 0.1) * math.tan(math.asin(1.0 / M_INF))  # 前縁マッハ波の下側
    m = yy <= ywin
    rou = ro * u
    imax = int(np.argmax(rou[m]))
    je = int(np.argmax(rou[m] >= edge_frac * rou[m][imax]))   # 壁側から最初に 99.5 % に達する点
    ue, roe, mue, Te, Pe = u[je], ro[je], mu[je], T[je], P[je]
    Me = ue / math.sqrt(GAM * R * Te)
    sl = slice(0, je + 1)
    dstar = float(np.trapezoid(1.0 - rou[sl] / (roe * ue), yy[sl]))
    theta = float(np.trapezoid(rou[sl] / (roe * ue) * (1.0 - u[sl] / ue), yy[sl]))
    Tw = T[0]; muw = mu[0]; row = ro[0]
    dudy = deriv_wall(yy[:3], u[:3]); dTdy = deriv_wall(yy[:3], T[:3])
    tau = muw * dudy
    lam_w = muw * CP / PR
    qw = lam_w * dTdy                      # 壁へ入る向き正 (T が壁から離れて増えるとき正)
    ut = math.sqrt(max(tau, 0.0) / row)
    Taw = Te * (1.0 + R_REC * 0.5 * (GAM - 1) * Me ** 2)
    St = qw / (roe * ue * CP * (Taw - Tw)) if abs(Taw - Tw) > 1.0 else float("nan")
    return dict(x=float(xv), Tw=float(Tw), Te=float(Te), Me=float(Me), Pe=float(Pe), ue=float(ue), roe=float(roe), mue=float(mue),
                Taw=float(Taw), Cf=float(tau / (0.5 * RO_INF * U_INF ** 2)), Cf_e=float(tau / (0.5 * roe * ue ** 2)),
                tau=float(tau), qw=float(qw), St=float(St), dstar=dstar, theta=theta, H=dstar / theta,
                y_edge=float(yy[je]), Re_theta=float(roe * ue * theta / mue), y1p=float(yy[1] * row * ut / muw),
                utau=ut, col=col, y=yy, u=u, ro=ro, T=T, rou=rou, je=je)


# ---------------------------------------------------------------- 理論
def cf_karman_schoenherr(re_theta):
    lg = np.log10(np.asarray(re_theta, float))
    return 1.0 / (17.08 * lg ** 2 + 25.11 * lg + 6.012)


def vd2_factors(Me, Te, Tw, Taw):
    """van Driest II: F_c, F_θ (=μ_e/μ_w)。C_f = C_f,i(F_θ Re_θ)/F_c。"""
    r = (Taw / Te - 1.0) / (0.5 * (GAM - 1) * Me ** 2)   # 実効回復係数
    A2 = r * 0.5 * (GAM - 1) * Me ** 2 * Te / Tw
    B = (1.0 + r * 0.5 * (GAM - 1) * Me ** 2) * Te / Tw - 1.0
    den = math.sqrt(B ** 2 + 4.0 * A2)
    alpha = (2.0 * A2 - B) / den; beta = B / den
    Fc = (Taw / Te - 1.0) / (math.asin(alpha) + math.asin(beta)) ** 2
    Ftheta = float(mu_suth(Te) / mu_suth(Tw))
    return Fc, Ftheta


def cf_vd2(re_theta, Me, Te, Tw, Taw):
    Fc, Ft = vd2_factors(Me, Te, Tw, Taw)
    return cf_karman_schoenherr(re_theta * Ft) / Fc


def temp_velocity(st, alpha=1.0):
    """T(u) の Walz 形 (alpha=1) / Duan–Martín 形 (alpha=0.8259) を forge の T と比較 (0.05<u/ue<0.95 の RMS 相対差)。"""
    ue, Te, Tw, Taw, Me = st["ue"], st["Te"], st["Tw"], st["Taw"], st["Me"]
    sl = slice(0, st["je"] + 1)
    uu = st["u"][sl] / ue; T = st["T"][sl]
    Tr = Taw
    Tm = Tw + (Tr - Tw) * ((1 - alpha) * uu ** 2 + alpha * uu) - R_REC * 0.5 * (GAM - 1) * Me ** 2 * Te * uu ** 2
    m = (uu > 0.05) & (uu < 0.95)
    return float(np.sqrt(np.mean(((Tm[m] - T[m]) / T[m]) ** 2)))


def vd_profile(st):
    """van Driest 変換 u_VD+ (y+ 基準は壁物性)。対数則との差 (30<y+<0.2 δ+) を返す。"""
    sl = slice(0, st["je"] + 1)
    y, u, ro = st["y"][sl], st["u"][sl], st["ro"][sl]
    ut = st["utau"]; row = ro[0]; muw = mu_suth(st["Tw"])
    yp = y * row * ut / muw
    up = u / ut
    uvd = np.concatenate([[0.0], np.cumsum(0.5 * (np.sqrt(ro[1:] / row) + np.sqrt(ro[:-1] / row)) * np.diff(up))])
    dp = st["y_edge"] * row * ut / muw
    m = (yp > 30) & (yp < 0.2 * dp)
    dev = float(np.mean(uvd[m] - (np.log(yp[m]) / 0.41 + 5.0))) if m.sum() > 3 else float("nan")
    return yp, up, uvd, dev


def plate_integrals(D, x_from=0.0, x_to=1.0, n=200):
    """SU2 の CD (REF_AREA=1, q_inf 基準) と HF (∫q_w dx, 単位幅) に対応する平板積分量。"""
    xs = np.linspace(x_from + 2e-3, x_to - 2e-3, n)
    cf = []; qw = []
    for xv in xs:
        st = station(D, xv); cf.append(st["Cf"]); qw.append(st["qw"])
    return dict(CD=float(np.trapezoid(cf, xs)), HF=float(np.trapezoid(qw, xs)))


def load_su2_surface(su2_dir):
    """surface_flow.csv (MARKER_PLOTTING wall): x, Skin_Friction_Coefficient_x, Heat_Flux, Y_Plus。"""
    path = Path(su2_dir) / "surface_flow.csv"
    with open(path) as f:
        rdr = csv.reader(f); hdr = [h.strip().strip('"') for h in next(rdr)]
        data = np.array([[float(v) for v in row] for row in rdr if row])
    col = {h: i for i, h in enumerate(hdr)}
    x = data[:, col["x"]]; o = np.argsort(x)
    out = dict(x=x[o])
    for k, names in (("Cf", ("Skin_Friction_Coefficient_x", "Skin_Friction_Coefficient")), ("qw", ("Heat_Flux",)), ("yplus", ("Y_Plus",))):
        for nm in names:
            if nm in col: out[k] = data[o, col[nm]]; break
    return out


def energy_closure(D):
    """単位幅: 入口 (x=-0.1) と出口 (x=1.0) の全エンタルピー流束差 vs ∫ q_w dx。"""
    def flux(xv):
        col = np.where(np.abs(D["x"] - xv) < 1e-6)[0]; col = col[np.argsort(D["y"][col])]
        y = D["y"][col]; ro = D["ro"][col]; u = D["u"][col]; v = D["v"][col]; T = D["T"][col]
        h0 = CP * T + 0.5 * (u ** 2 + v ** 2)
        return float(np.trapezoid(ro * u * h0, y)), float(np.trapezoid(ro * u, y))
    Hin, min_ = flux(D["x"].min()); Hout, mout = flux(D["x"].max())
    wallx = np.unique(np.round(D["x"][(np.abs(D["y"]) < 1e-9) & (D["x"] >= -1e-9)], 7))
    qs = np.array([station(D, xw)["qw"] for xw in wallx[1:-1:4]])
    Q = float(np.trapezoid(qs, wallx[1:-1:4]))
    return dict(H_in=Hin, H_out=Hout, dH=Hin - Hout, Q_wall=Q, mass_in=min_, mass_out=mout,
                closure=(Hin - Hout) / Q if abs(Q) > 0 else float("nan"))


# ---------------------------------------------------------------- 本体
def evaluate(D, stations, label, contur=True, theta0=None):
    rows = []
    for xs in stations:
        st = station(D, xs)
        Fc, Ft = vd2_factors(st["Me"], st["Te"], st["Tw"], st["Taw"])
        st["Cf_vd2"] = float(cf_vd2(st["Re_theta"], st["Me"], st["Te"], st["Tw"], st["Taw"]))
        st["Fc"] = Fc; st["Ftheta"] = Ft
        st["RAF"] = 2.0 * st["St"] / st["Cf"] if st["Cf"] > 0 else float("nan")
        st["r_rec_meas"] = (st["Tw"] - st["Te"]) / (st["Te"] * 0.5 * (GAM - 1) * st["Me"] ** 2)
        st["dT_walz"] = temp_velocity(st, 1.0); st["dT_dm"] = temp_velocity(st, 0.8259)
        yp, up, uvd, st["vd_logdev"] = vd_profile(st)
        st["_vd"] = (yp, up, uvd)
        rows.append(st)
    if contur:
        from forge_design.feedback.deltastar_integral import flat_plate_integral
        xs = np.linspace(rows[0]["x"], rows[-1]["x"], 40)
        Tw = None if abs(rows[0]["Tw"] - rows[0]["Taw"]) < 0.06 * rows[0]["Taw"] else float(np.mean([r["Tw"] for r in rows]))
        Me = float(np.mean([r["Me"] for r in rows])); Te = float(np.mean([r["Te"] for r in rows])); Pe = float(np.mean([r["Pe"] for r in rows]))
        ci = flat_plate_integral(xs, Me, Te, Pe, GAM, CP, theta0_m=(theta0 if theta0 else rows[0]["theta"]), Tw=Tw)
        for r in rows:
            r["theta_contur"] = float(np.interp(r["x"], xs, ci["theta"])); r["dstar_contur"] = float(np.interp(r["x"], xs, ci["dstar"]))
            r["H_contur"] = float(np.interp(r["x"], xs, ci["H"])); r["Cf_contur"] = float(np.interp(r["x"], xs, ci["Cf"]))
    return rows


def table(rows, label):
    print(f"\n=== {label} ===")
    print(" x     Tw[K]  Me    Re_th  y1+   Cf*1e3 VD2*1e3 ratio | qw[kW/m2] St*1e3 2St/Cf | d*[mm] th[mm]  H   | CONTUR d* th  H  Cf | dT walz/DM  VDlog")
    for r in rows:
        c = f"{r.get('dstar_contur',float('nan'))*1e3:5.2f} {r.get('theta_contur',float('nan'))*1e3:5.3f} {r.get('H_contur',float('nan')):4.2f} {r.get('Cf_contur',float('nan'))*1e3:5.3f}"
        print(f"{r['x']:.2f} {r['Tw']:6.0f} {r['Me']:5.2f} {r['Re_theta']:6.0f} {r['y1p']:5.2f} {r['Cf']*1e3:6.3f} {r['Cf_vd2']*1e3:6.3f} {r['Cf']/r['Cf_vd2']:5.3f} | "
              f"{r['qw']/1e3:8.2f} {r['St']*1e3:6.3f} {r['RAF']:5.2f} | {r['dstar']*1e3:6.3f} {r['theta']*1e3:5.3f} {r['H']:5.2f} | {c} | "
              f"{r['dT_walz']*100:4.1f}% {r['dT_dm']*100:4.1f}% {r['vd_logdev']:5.2f}")


def series(run_dir, stations):
    fs = sorted(glob.glob(str(Path(run_dir) / "res_[0-9]*.h5")), key=lambda s: int(Path(s).stem.split("_")[1]))
    fs = [f for f in fs if int(Path(f).stem.split("_")[1]) > 0]
    keys = ("Cf", "qw", "dstar", "theta", "Tw")
    tol = dict(Cf=0.015, qw=0.015, dstar=0.01, theta=0.01, Tw=0.005)
    hist = {xs: {k: [] for k in keys} for xs in stations}; steps = []
    for f in fs:
        D = load_forge(f); steps.append(int(Path(f).stem.split("_")[1]))
        for xs in stations:
            st = station(D, xs)
            for k in keys: hist[xs][k].append(st[k])
    print(f"\n=== series {run_dir} ({len(fs)} snapshots {steps[0]}..{steps[-1]}) — tail 50 % drift (max-min)/|mean| ===")
    allok = True
    for xs in stations:
        line = []
        for k in keys:
            v = np.array(hist[xs][k]); t = v[len(v) // 2:]
            if k == "qw" and abs(t.mean()) < 1e-3 * RO_INF * U_INF * CP * 100.0:   # 断熱 (q_w≈0): 相対 drift は無意味 → 絶対床で判定
                drift = (t.max() - t.min()) / (RO_INF * U_INF * CP * 100.0)
            else:
                drift = (t.max() - t.min()) / max(abs(t.mean()), 1e-30)
            ok = drift <= tol[k]; allok &= ok
            line.append(f"{k} {t[-1]:.4g} drift {drift*100:.2f}% {'STEADY' if ok else 'DRIFTING'}")
        print(f"x={xs:.2f}: " + " | ".join(line))
    print("SERIES VERDICT:", "STEADY" if allok else "NOT STEADY")
    return allok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--su2", default=None); ap.add_argument("--stations", default="0.3,0.45,0.6,0.75,0.9")
    ap.add_argument("--out", default=None); ap.add_argument("--series", action="store_true"); ap.add_argument("--json", action="store_true")
    ap.add_argument("--ref", default=None, help="比の基準 run (断熱 A): Cf/δ*/θ の比を VD-II/CONTUR の比と比べる")
    ap.add_argument("--closure", action="store_true"); ap.add_argument("--integrals", action="store_true")
    a = ap.parse_args()
    stations = [float(v) for v in a.stations.split(",")]
    fs = sorted(glob.glob(str(Path(a.run) / "res_[0-9]*.h5")), key=lambda s: int(Path(s).stem.split("_")[1]))
    D = load_forge(fs[-1]); rows = evaluate(D, stations, a.run)
    table(rows, f"forge {a.run} ({Path(fs[-1]).name})")
    out = {"forge": [{k: v for k, v in r.items() if not k.startswith("_") and not isinstance(v, np.ndarray)} for r in rows]}
    if a.closure:
        ec = energy_closure(D); out["closure"] = ec
        print(f"energy closure: dH(in-out) {ec['dH']:.1f} W/m  Q_wall {ec['Q_wall']:.1f} W/m  ratio {ec['closure']:.3f}  (mass in/out {ec['mass_in']:.4f}/{ec['mass_out']:.4f})")
    if a.integrals:
        pi = plate_integrals(D); out["integrals"] = pi
        print(f"plate integrals (x 0-1): CD=int Cf dx = {pi['CD']:.5g}   HF=int q_w dx = {pi['HF']:.6g} W/m")
    if a.su2 and (Path(a.su2) / "surface_flow.csv").exists():
        sf = load_su2_surface(a.su2)
        print("\nSU2 surface_flow.csv vs forge (wall):  x   Cf_su2  Cf_forge  ratio |  qw_su2  qw_forge  ratio | y+_su2")
        for r in rows:
            k = int(np.argmin(np.abs(sf["x"] - r["x"])))
            cs = sf.get("Cf", np.full_like(sf["x"], np.nan))[k]; qs = sf.get("qw", np.full_like(sf["x"], np.nan))[k]; yp = sf.get("yplus", np.full_like(sf["x"], np.nan))[k]
            print(f"   {r['x']:.2f}  {cs*1e3:6.3f}  {r['Cf']*1e3:6.3f}  {r['Cf']/cs if cs else float('nan'):.4f} | {qs/1e3:8.2f} {r['qw']/1e3:8.2f} {r['qw']/qs if qs else float('nan'):.4f} | {yp:.2f}")
        out["su2_surface"] = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in sf.items()}
    if a.su2 and (Path(a.su2) / "restart_flow.csv").exists():
        S = load_su2(a.su2); srows = evaluate(S, stations, a.su2)
        table(srows, f"SU2 {a.su2}")
        print("\nforge/SU2:  x    Cf     qw     dstar  theta")
        for r, s in zip(rows, srows):
            print(f"          {r['x']:.2f} {r['Cf']/s['Cf']:.4f} {r['qw']/s['qw'] if s['qw'] else float('nan'):.4f} {r['dstar']/s['dstar']:.4f} {r['theta']/s['theta']:.4f}")
        out["su2"] = [{k: v for k, v in r.items() if not k.startswith("_") and not isinstance(v, np.ndarray)} for r in srows]
    if a.ref:
        rf = sorted(glob.glob(str(Path(a.ref) / "res_[0-9]*.h5")), key=lambda s: int(Path(s).stem.split("_")[1]))
        Rr = evaluate(load_forge(rf[-1]), stations, a.ref)
        print(f"\nratio to {a.ref}:  x   Cf_ratio  VD2_ratio  | dstar_ratio CONTUR_ratio | theta_ratio CONTUR_ratio")
        for r, q in zip(rows, Rr):
            vd = cf_vd2(r["Re_theta"], r["Me"], r["Te"], r["Tw"], r["Taw"]) / cf_vd2(q["Re_theta"], q["Me"], q["Te"], q["Tw"], q["Taw"])
            print(f"                 {r['x']:.2f}  {r['Cf']/q['Cf']:.4f}   {vd:.4f}   | {r['dstar']/q['dstar']:.4f}  {r['dstar_contur']/q['dstar_contur']:.4f} | {r['theta']/q['theta']:.4f}  {r['theta_contur']/q['theta_contur']:.4f}")
    if a.series:
        out["series_steady"] = series(a.run, [0.3, 0.6, 0.9])
    if a.out:
        plot(rows, out.get("su2"), a.out)
    if a.json or a.out:
        Path((a.out or str(Path(a.run) / "cooled_plate_eval")) + ".json").write_text(json.dumps(out, indent=1, default=float))


def plot(rows, srows, prefix):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    x = [r["x"] for r in rows]
    ax[0, 0].plot([r["Re_theta"] for r in rows], [r["Cf"] * 1e3 for r in rows], "o-", label="forge")
    ax[0, 0].plot([r["Re_theta"] for r in rows], [r["Cf_vd2"] * 1e3 for r in rows], "k--", label="van Driest II (K–S)")
    if srows: ax[0, 0].plot([r["Re_theta"] for r in srows], [r["Cf"] * 1e3 for r in srows], "s-", label="SU2")
    ax[0, 0].set_xlabel("Re_theta"); ax[0, 0].set_ylabel("Cf x1e3"); ax[0, 0].legend()
    ax[0, 1].plot(x, [r["qw"] / 1e3 for r in rows], "o-", label="forge q_w"); 
    if srows: ax[0, 1].plot(x, [r["qw"] / 1e3 for r in srows], "s-", label="SU2")
    ax[0, 1].set_xlabel("x [m]"); ax[0, 1].set_ylabel("q_w [kW/m2]"); ax[0, 1].legend()
    ax[0, 2].plot(x, [r["RAF"] for r in rows], "o-"); ax[0, 2].axhspan(1.0, 1.2, alpha=0.2); ax[0, 2].set_ylabel("2St/Cf"); ax[0, 2].set_xlabel("x [m]")
    ax[1, 0].plot(x, [r["dstar"] * 1e3 for r in rows], "o-", label="δ* forge"); ax[1, 0].plot(x, [r["dstar_contur"] * 1e3 for r in rows], "k--", label="δ* CONTUR")
    ax[1, 0].plot(x, [r["theta"] * 1e3 for r in rows], "o-", label="θ forge"); ax[1, 0].plot(x, [r["theta_contur"] * 1e3 for r in rows], "k:", label="θ CONTUR")
    if srows: ax[1, 0].plot(x, [r["dstar"] * 1e3 for r in srows], "s", label="δ* SU2"); ax[1, 0].plot(x, [r["theta"] * 1e3 for r in srows], "s", label="θ SU2")
    ax[1, 0].set_xlabel("x [m]"); ax[1, 0].set_ylabel("[mm]"); ax[1, 0].legend(fontsize=8)
    r = rows[len(rows) // 2]; yp, up, uvd = r["_vd"]
    m = yp > 0.5
    ax[1, 1].semilogx(yp[m], uvd[m], "-", label="u_VD+ forge"); ax[1, 1].semilogx(yp[m], up[m], ":", label="u+ raw")
    yl = np.logspace(1, 3.5, 50); ax[1, 1].semilogx(yl, np.log(yl) / 0.41 + 5.0, "k--", label="log law"); ax[1, 1].set_xlabel("y+"); ax[1, 1].legend()
    sl = slice(0, r["je"] + 1); uu = r["u"][sl] / r["ue"]
    ax[1, 2].plot(uu, r["T"][sl] / r["Te"], "-", label="forge")
    for al, lb in ((1.0, "Walz"), (0.8259, "Duan–Martín")):
        Tm = r["Tw"] + (r["Taw"] - r["Tw"]) * ((1 - al) * uu ** 2 + al * uu) - R_REC * 0.5 * (GAM - 1) * r["Me"] ** 2 * r["Te"] * uu ** 2
        ax[1, 2].plot(uu, Tm / r["Te"], "--", label=lb)
    ax[1, 2].set_xlabel("u/ue"); ax[1, 2].set_ylabel("T/Te"); ax[1, 2].legend()
    fig.tight_layout(); fig.savefig(prefix + ".png", dpi=110)


if __name__ == "__main__":
    main()
