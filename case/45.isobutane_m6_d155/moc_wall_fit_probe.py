"""MOC 壁点の表現 (補間 vs 位置+壁角の同時当てはめ) と MOC 解像度 n_axis_inv の効果を測る (形状だけ, CFD なし)。
plan: plans/active/verification-m6-axis-wave-mesh-su2.md §5.1 #12 (設計壁の r'' リップル)。
usage: design/.venv-opt/bin/python moc_wall_fit_probe.py → _band_ab/moc_wall_fit_probe.json
r'' の高周波 = 0.5 r_t 移動平均からの残差 (報告ツールと同じ定義)。帯 [0.5,20) は滑らかな曲線でも同じ値になる (定義の癖) ことも出す。
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem

C = Path(__file__).resolve().parent
def joint_fit(tb, R, h0=0.05, h1=0.5, xg=6.0, sig_r=1e-6, sig_th=1e-4, lam=1e-3, k=5):
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    x0, xe = x[0], x[-1]
    # ノット間隔: スロート近傍 h0 → xg 以降 h1 (滑らかに)
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (xg - x0), 1.0); xs.append(xs[-1] + h0 + (h1 - h0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * h1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1
    B0 = BSpline.design_matrix(x, t, k).toarray()
    E = np.eye(nc)
    D = lambda xq, d: np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B1 = D(x, 1)
    ds = np.gradient(x); w = ds / ds.mean()
    xq = np.linspace(x0, xe, 6000); B3 = D(xq, 3); wq = np.gradient(xq)
    A = (B0.T * (w / sig_r**2)) @ B0 + (B1.T * (w / sig_th**2)) @ B1 + lam * (B3.T * wq) @ B3 / sig_r**2
    b = B0.T @ (w * r / sig_r**2) + B1.T @ (w * np.tan(th) / sig_th**2)
    C = np.array([D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0], D(np.r_[xe], 1)[0]])
    dv = np.array([r[0], x0 / R, 1.0 / R, np.tan(th[-1])])
    K = np.block([[A, C.T], [C, np.zeros((4, 4))]])
    c = np.linalg.solve(K, np.r_[b, dv])[:nc]
    return BSpline(t, c, k), nc


BANDS = ((0.5, 2), (2, 6), (6, 20), (20, 60), (60, 95))
xx = np.arange(0.53, 94.9, 0.001); K = 500


def band(x, v):
    return {f"[{a},{b})": float("%.3e" % np.abs(v[(x >= a) & (x < b)]).max()) for a, b in BANDS}


def hf(v):
    h = v - np.convolve(v, np.ones(K) / K, "same"); h[:K] = 0; h[-K:] = 0
    return h


W, out = {}, {}
for n in (600, 1200, 2400, 9600):
    p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = n
    d = design_chain(p); tb = d["wall_inv"]; s, nc = joint_fit(tb, d["R"], lam=1e-7); W[n] = (d["wall"], s)
    out[n] = {"n_pts": len(tb), "n_cp_fit": nc,
              "interp_r2_hf": band(xx, hf(d["wall"].r(xx, 2))), "fit_r2_hf": band(xx, hf(s(xx, 2))),
              "interp_dtheta_pts_deg": band(tb[:, 0], np.degrees(np.arctan(d["wall"].r(tb[:, 0], 1)) - tb[:, 2])),
              "fit_dr_pts": band(tb[:, 0], s(tb[:, 0]) - tb[:, 1]),
              "fit_dtheta_pts_deg": band(tb[:, 0], np.degrees(np.arctan(s(tb[:, 0], 1)) - tb[:, 2]))}
rf = W[9600][1]
for n in (600, 1200, 2400):
    s = W[n][1]
    out[n]["fit_dr_vs_fit9600"] = band(xx, s(xx) - rf(xx))
    out[n]["fit_dtheta_vs_fit9600_deg"] = band(xx, np.degrees(np.arctan(s(xx, 1)) - np.arctan(rf(xx, 1))))
(C / "_band_ab/moc_wall_fit_probe.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
