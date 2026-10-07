"""MOC 壁の同時当てはめ: 形状だけの A/B (CFD 0 step)。plan verification-m6-axis-wave-mesh-su2 §5.1 #13 (事前登録済み, commit 72cabeaf)。
同一 MOC 点群 (n_axis_inv 2400)、ノット h0 0.0125 → h1 0.5、変えるのは λ だけ (A=0, B=1e-9)。比較基準に現行の補間壁 (同じ点群) も並べる。
usage: design/.venv-opt/bin/python moc_wall_fit_ab.py → _band_ab/moc_wall_fit_ab.json
滑らかさ (r″ の局所 3 次多項式残差) は #14 で指標を作り直すまで参考値 (閾値なし)。
"""
import copy, json, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from scipy.signal import savgol_filter
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas, delta_r_from_table
from forge_design.geometry.wall_axismach import PhysicalNozzleWall

C = Path(__file__).resolve().parent
N_AXIS, H0, H1, XG = 2400, 0.0125, 0.5, 6.0
GATE_DR, GATE_DTH = 5e-6, 0.005


def joint_fit(tb, R, lam, k=5, sig_r=1e-6, sig_th=1e-4):
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (XG - x0), 1.0); xs.append(xs[-1] + H0 + (H1 - H0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * H1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1
    E = np.eye(nc)
    D = lambda xq, d: np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1)
    w = np.gradient(x); w = w / w.mean()
    A = (B0.T * (w / sig_r ** 2)) @ B0 + (B1.T * (w / sig_th ** 2)) @ B1
    b = B0.T @ (w * r / sig_r ** 2) + B1.T @ (w * np.tan(th) / sig_th ** 2)
    if lam > 0:
        xq = np.linspace(x0, xe, 8000); B3 = D(xq, 3)
        A = A + lam * (B3.T * np.gradient(xq)) @ B3 / sig_r ** 2
    Cm = np.array([D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0], D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]])
    dv = np.array([r[0], x0 / R, 1.0 / R, r[-1], np.tan(th[-1])])
    K = np.block([[A, Cm.T], [Cm, np.zeros((5, 5))]])
    c = np.linalg.solve(K, np.r_[b, dv])[:nc]
    return BSpline(t, c, k), nc


def n_extrema(v):
    dv = np.diff(v); dv = dv[np.abs(dv) > 1e-14]
    return int(np.count_nonzero(np.diff(np.sign(dv)) != 0))


def sg_resid(xg, f2, win):
    n = int(round(win / (xg[1] - xg[0]))) | 1
    res = f2 - savgol_filter(f2, n, 3)
    m = (xg > xg[0] + win) & (xg < xg[-1] - win)
    out = {}
    for a, b in ((0, 2), (2, 20), (20, 60), (60, 96)):
        mm = m & (xg >= a) & (xg < b)
        out[f"[{a},{b})"] = float("%.3e" % np.abs(res[mm]).max()) if mm.any() else None
    return out


p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = N_AXIS
d = design_chain(p); tb = d["wall_inv"]; R = d["R"]; S = float(p.spec["r_throat"])
dtab = np.loadtxt(C / "_band_ab/delta_r_c2final_run0051.csv", delimiter=",", skiprows=1)  # C2 の δ_r (run_0051 の delta_r_initial.csv の写し; 両壁で共通)
f_dr = delta_r_from_table(dtab[:, 0], dtab[:, 1])
walls = {"interp": (d["wall"], None)}
for name, lam in (("A_lam0", 0.0), ("B_lam1e-9", 1e-9)):
    s, nc = joint_fit(tb, R, lam); w = copy.copy(d["wall"]); w._spl = s; walls[name] = (w, nc)
x0, xe = float(tb[0, 0]), float(tb[-1, 0])
xg = np.arange(x0, xe, 1e-3)
out = {"n_axis_inv": N_AXIS, "n_pts": len(tb), "x0": x0, "x_e": xe, "knots": dict(h0=H0, h1=H1, xg=XG), "gates": dict(dr=GATE_DR, dtheta_deg=GATE_DTH)}
for name, (w, nc) in walls.items():
    dr = w.r(tb[:, 0]) - tb[:, 1]; dth = np.degrees(np.arctan(w.r(tb[:, 0], 1)) - tb[:, 2])
    i_r, i_t = int(np.argmax(np.abs(dr))), int(np.argmax(np.abs(dth)))
    r2 = w.r(xg, 2); rr = w.r(xg)
    eps = 1e-7
    jL = [float(w.r(np.r_[x0 - eps], k)[0]) for k in range(3)]; jR = [float(w.r(np.r_[x0 + eps], k)[0]) for k in range(3)]
    PW = PhysicalNozzleWall(w, tb, S, float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp, offset="radial", delta_r_x=f_dr)
    xp = np.arange(PW.x_throat + 0.01, PW.x_e, 1e-3)
    o = dict(n_cp=nc, max_dr_pts=float(np.abs(dr).max()), x_max_dr=float(tb[i_r, 0]),
             max_dtheta_pts_deg=float(np.abs(dth).max()), x_max_dtheta=float(tb[i_t, 0]),
             gate_fidelity=bool(np.abs(dr).max() <= GATE_DR and np.abs(dth).max() <= GATE_DTH),
             monotone_r=bool(np.all(np.diff(rr) > 0)), r2_sign_changes=n_extrema(np.r_[0, np.diff(np.sign(r2))]) if False else int(np.count_nonzero(np.diff(np.sign(r2)) != 0)),
             r2_local_extrema=n_extrema(r2), r3_sign_changes=int(np.count_nonzero(np.diff(np.sign(w.r(xg, 3))) != 0)),
             start_jump_r_r1_r2=[jR[k] - jL[k] for k in range(3)], start_r1_r2=jR[1:], exit_dr=float(w.r(np.r_[xe])[0] - tb[-1, 1]),
             exit_r1=float(w.r(np.r_[xe], 1)[0]),
             r2_sg_resid_ref={f"win{wn}": sg_resid(xg, r2, wn) for wn in (0.25, 0.5, 1.0, 2.0)},
             physical=dict(validate=PW.validate(), x_throat=PW.x_throat, r_throat=PW.r_throat, kappa_throat=PW.kappa_throat,
                           r_exit=float(PW.r(np.r_[PW.x_e])[0]), r1_exit=float(PW.r(np.r_[PW.x_e], 1)[0]),
                           r2_local_extrema=n_extrema(PW.r(xp, 2)),
                           r2_sg_resid_ref={f"win{wn}": sg_resid(xp, PW.r(xp, 2), wn) for wn in (0.5, 1.0)}))
    out[name] = o
PI = None
for name in ("A_lam0", "B_lam1e-9"):
    w = walls[name][0]; wi = walls["interp"][0]
    out[name]["vs_interp_design_max_dr"] = float(np.abs(w.r(xg) - wi.r(xg)).max())
    out[name]["vs_interp_design_max_dtheta_deg"] = float(np.degrees(np.abs(np.arctan(w.r(xg, 1)) - np.arctan(wi.r(xg, 1))).max()))
(C / "_band_ab/moc_wall_fit_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
