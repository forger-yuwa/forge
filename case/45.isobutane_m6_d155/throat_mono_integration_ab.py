"""plan tooling-nozzle-throat-monotone-r2 の諮問 (notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md) の判別 A/B (CFD 0 step)。
単調拘束付き当てはめ C の平滑化項 ∫(r‴)² の積分方法だけを変える:
  A = 現行 (全長 8000 点の Σ(r‴)²·gradient(x)) / B = 各非零ノット区間で 3 点 Gauss (r‴ は区間 2 次式なので厳密)。
MOC 点群・λ 1e-9・ノット・始終点拘束・単調区間 [0,1.5]・QP 解法は固定。
判定 (事前閾値): D_k = max_[0,0.3] |r_B^(k) − r_A^(k)| / max_[0,0.3] |r_A^(k)| (k = 3, 4)。D3 または D4 > 0.10 → 積分依存あり。両方 ≤ 0.10 → 積分誤差主因説を棄却。
参考 (判定外): 単調拘束なしの現行 V0 にも同じ A/B をかけ、r″ の山の高さが積分方法で変わるかを記録する。
usage: [CASE_RUNS=<run_* のある case dir>] python3 throat_mono_integration_ab.py [OUT_JSON]
"""
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/throat_mono_integration_ab.json"
p = load_problem(C / "problem_d155_ns_finemesh_recal_final.yaml")
p.geometry["initial_line_run"] = str((Path(os.environ.get("CASE_RUNS", C)) / p.geometry["initial_line_run"]).resolve())
d = design_chain(p); tb = d["wall_inv"]; R = float(d["R"])
x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]


def fit(tb, integ="grid", mono=None, lam=1e-9, k=5, sig_r=1e-6, sig_th=1e-4, h0=0.0125, h1=0.5, x_g=6.0):
    """joint_fit_wall と同じ目的関数・ノット・拘束 (始点 r, r′=x0/R, r″=1/R; 出口 r, r′)。integ: 'grid' (現行) / 'gauss' (区間 3 点 Gauss)。"""
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]; x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (x_g - x0), 1.0); xs.append(xs[-1] + h0 + (h1 - h0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * h1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1; E = np.eye(nc)
    D = lambda xq, dd: np.array([BSpline(t, E[i], k)(xq, dd) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1); w = np.gradient(x); w = w / w.mean()
    if integ == "grid":
        xq = np.linspace(x0, xe, 8000); wq = np.gradient(xq)
    else:
        br = np.unique(t); gp, gw = np.polynomial.legendre.leggauss(3)
        a, b = br[:-1], br[1:]; xq = (0.5 * (b - a)[:, None] * gp[None, :] + 0.5 * (a + b)[:, None]).ravel()
        wq = (0.5 * (b - a)[:, None] * gw[None, :]).ravel()
    B3 = D(xq, 3)
    A = (B0.T * (w / sig_r ** 2)) @ B0 + (B1.T * (w / sig_th ** 2)) @ B1 + lam * (B3.T * wq) @ B3 / sig_r ** 2
    b = B0.T @ (w * r / sig_r ** 2) + B1.T @ (w * np.tan(th) / sig_th ** 2)
    sc = float(np.abs(np.diag(A)).max()); A = A / sc; b = b / sc
    Ce = np.array([D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0], D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]])
    de = np.array([r[0], x0 / R, 1.0 / R, r[-1], np.tan(th[-1])]); ne = len(de)
    if mono:
        M3 = np.array([BSpline(t, E[i], k).derivative(3).c[:nc - 3] for i in range(nc)]).T; t3 = t[3:len(t) - 3]
        G = M3[[j for j in range(nc - 3) if t3[j + 3] > mono[0] + 1e-12 and t3[j] < mono[1]]]
        G = G / np.abs(G).max(axis=1, keepdims=True)
    else:
        G = np.zeros((0, nc))
    act = []
    for it in range(200):
        Cm = np.vstack([Ce, G[act]]) if act else Ce; dd = np.r_[de, np.zeros(len(act))]; m = len(dd)
        sol = np.linalg.solve(np.block([[A, Cm.T], [Cm, np.zeros((m, m))]]), np.r_[b, dd]); c, mu = sol[:nc], sol[nc + ne:]
        if len(act) and mu.min() < -1e-12 * max(1.0, np.abs(mu).max()):
            act.pop(int(np.argmin(mu))); continue
        if not len(G) or (G @ c).max() <= 1e-10 * max(1.0, float(np.abs(G @ c).max())):
            break
        act.append(int(np.argmax(G @ c)))
    else:
        raise RuntimeError("有効制約法が 200 回で収束しない")
    s = BSpline(t, c, k)
    kkt = dict(eq_resid=float(np.abs(Ce @ c - de).max()), ineq_max=float((G @ c).max()) if len(G) else None,
               mu_min=float(mu.min()) if len(act) else None, n_active=len(act), iters=it + 1,
               active_t3=[[float(t[3 + j]), float(t[3 + j + 3])] for j in []])
    if mono:
        sel = [j for j in range(nc - 3) if t3[j + 3] > mono[0] + 1e-12 and t3[j] < mono[1]]
        kkt["active_support"] = [[float(t3[sel[a]]), float(t3[sel[a] + 3])] for a in act]
    # 区間ごとの厳密な ∫(r‴)² (3 点 Gauss)
    br = np.unique(t); gp, gw = np.polynomial.legendre.leggauss(3)
    return s, kkt, br, gp, gw


def exact_int_r3sq(s, br, gp, gw, a=None, bnd=None):
    lo, hi = br[:-1], br[1:]
    if a is not None:
        m = (hi > a) & (lo < bnd); lo, hi = np.maximum(lo[m], a), np.minimum(hi[m], bnd)
    xq = 0.5 * (hi - lo)[:, None] * gp[None, :] + 0.5 * (hi + lo)[:, None]; wq = 0.5 * (hi - lo)[:, None] * gw[None, :]
    return float((s(xq.ravel(), 3) ** 2 * wq.ravel()).sum())


xg = np.linspace(0, 0.3, 30001)
out = {"note": "判別 A/B (積分方法のみ)。D_k の閾値 0.10 は診断の事前閾値", "x_points_lt0p3": [float(v) for v in x[x < 0.3]]}
for tag, mono in (("C_mono", (0.0, 1.5)), ("V0_ref", None)):
    S = {}
    for integ in ("grid", "gauss"):
        s, kkt, br, gp, gw = fit(tb, integ, mono)
        nr = x > 0
        S[integ] = (s, dict(kkt=kkt,
                            r2_max_0_0p3=float(s(xg, 2).max()), x_r2_max=float(xg[np.argmax(s(xg, 2))]),
                            r3_max_0_0p3=float(s(xg, 3).max()), r3_min_0_0p3=float(s(xg, 3).min()),
                            r4_absmax_0_0p3=float(np.abs(s(xg, 4)).max()),
                            int_r3sq_0_0p3=exact_int_r3sq(s, br, gp, gw, 0.0, 0.3), int_r3sq_0_1p5=exact_int_r3sq(s, br, gp, gw, 0.0, 1.5),
                            int_r3sq_all=exact_int_r3sq(s, br, gp, gw),
                            dth_first_deg=float(np.degrees(np.arctan(s(x[1], 1)) - th[1])),
                            dth_max_ge2_deg=float(np.abs(np.degrees(np.arctan(s(x[2:], 1)) - th[2:])).max()),
                            dr_max_ge1_rt=float(np.abs(s(x[1:]) - r[1:]).max())))
    sA, sB = S["grid"][0], S["gauss"][0]
    Dk = {f"D{k}": float(np.abs(sB(xg, k) - sA(xg, k)).max() / np.abs(sA(xg, k)).max()) for k in (2, 3, 4)}
    xa = np.linspace(0, x[-1], 400001)
    out[tag] = dict(grid=S["grid"][1], gauss=S["gauss"][1], **Dk,
                    dr_all_max_rt=float(np.abs(sB(xa) - sA(xa)).max()),
                    dth_all_max_deg=float(np.abs(np.degrees(np.arctan(sB(xa, 1)) - np.arctan(sA(xa, 1)))).max()))
c = out["C_mono"]
out["verdict"] = ("積分依存あり (D3 または D4 > 0.10): 厳密積分で形状ゲートを組み直す" if max(c["D3"], c["D4"]) > 0.10
                  else "積分誤差主因説を棄却 (D3, D4 ≤ 0.10): 急変は拘束に伴う形状上の代償として評価")
OUT.parent.mkdir(parents=True, exist_ok=True); OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
