"""初期線整合化の A/B (形状のみ, CFD 0 step)。plan verification-m6-axis-wave-mesh-su2 §5.1 #13 ⑦ (事前登録 commit 3088b853)。
A: 現行 Hall 初期線 (throat_characteristic の M・θ をそのまま)。
B: r 配列・Hall の θ(r)・軸端の x, M, θ を固定し、軸→壁の各区間で
   Δr = tan[½(θᵢ₋₁+θᵢ−μᵢ₋₁−μᵢ)]·Δx と Δ(θ+ν) = S₋ (moc_kernel と同じ源項・軸極限・ν↔M) を xᵢ・νᵢ について連立。
前提ゲート: B の適合残差・幾何角残差がともに最大 1e-8 rad 以下、全点 M>1・有限、壁抽出の欠落なし。
判定: 0<ξ<0.08 (ξ = x − x_w) の |a_θ| ≤ 0.01932° → 第 1 仮説を支持 / 減少 80 % 未満 → 棄却。
usage: design/.venv-opt/bin/python initial_line_ab.py → _band_ab/initial_line_ab.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.evaluate import runner_axismach as RA
from forge_design.geometry import moc_inverse as MI
from forge_design.geometry.moc_kernel import pm_nu, pm_mach, _sin_over_r_vec

C = Path(__file__).resolve().parent
N_START = 161


def line_residuals(x, r, M, th, g):
    nu = np.array([float(pm_nu(float(m), g)) for m in M]); mu = np.arcsin(1.0 / M)
    Gs, Es = [], []
    for i in range(1, len(x)):
        ang = 0.5 * (th[i - 1] + th[i] - mu[i - 1] - mu[i])
        Gs.append(np.arctan2(r[i] - r[i - 1], x[i] - x[i - 1]) - (ang if ang > -np.pi / 2 else ang) if x[i] != x[i - 1] else 0.0)
        Gs[-1] = np.arctan((r[i] - r[i - 1]) / (x[i] - x[i - 1])) - ang
        fB = 0.5 * (np.sin(mu[i]) * float(_sin_over_r_vec(np.r_[r[i]], np.r_[th[i]], np.r_[r[i - 1]], np.r_[th[i - 1]])[0])
                    + np.sin(mu[i - 1]) * float(_sin_over_r_vec(np.r_[r[i - 1]], np.r_[th[i - 1]], np.r_[r[i]], np.r_[th[i]])[0]))
        Es.append((th[i - 1] + nu[i - 1]) - (th[i] + nu[i]) - fB / np.cos(ang) * (x[i - 1] - x[i]))
    return np.array(Gs), np.array(Es)


def consistent_line(xq, rq, Mq, tq, g, iters=200):
    x = np.array(xq, float); r = np.array(rq, float); th = np.array(tq, float)
    nu = np.array([float(pm_nu(float(m), g)) for m in Mq]); M = np.array(Mq, float)
    for i in range(1, len(r)):
        for _ in range(iters):
            mu_p, mu_i = np.arcsin(1 / M[i - 1]), np.arcsin(1 / M[i])
            ang = 0.5 * (th[i - 1] + th[i] - mu_p - mu_i)
            x_new = x[i - 1] + (r[i] - r[i - 1]) / np.tan(ang)
            fB = 0.5 * (np.sin(mu_i) * float(_sin_over_r_vec(np.r_[r[i]], np.r_[th[i]], np.r_[r[i - 1]], np.r_[th[i - 1]])[0])
                        + np.sin(mu_p) * float(_sin_over_r_vec(np.r_[r[i - 1]], np.r_[th[i - 1]], np.r_[r[i]], np.r_[th[i]])[0]))
            nu_new = (th[i - 1] + nu[i - 1]) - th[i] - fB / np.cos(ang) * (x[i - 1] - x_new)
            dx, dn = abs(x_new - x[i]), abs(nu_new - nu[i])
            x[i], nu[i] = x_new, nu_new; M[i] = float(pm_mach(nu_new, g))
            if dx < 1e-15 and dn < 1e-15:
                break
    return x, r, M, th


class _Proxy:
    def __init__(self, ht, g):
        self._ht, self._g = ht, g

    def __getattr__(self, k):
        return getattr(self._ht, k)

    def throat_characteristic(self, n=61, **kw):
        xq, rq, Mq, tq = self._ht.throat_characteristic(n=n, **kw)
        return consistent_line(xq, rq, Mq, tq, self._g)


_inv = RA.inverse_design
_cpl = MI.cplus_flux_wall
res = {}
for arm in ("A_hall", "B_consistent"):
    CAP = {}

    def wrap(throat, target, **kw):
        g = kw.get("gamma")
        th_obj = _Proxy(throat, g) if arm == "B_consistent" else throat
        xq, rq, Mq, tq = th_obj.throat_characteristic(n=N_START)
        CAP["line"] = (np.array(xq), np.array(rq), np.array(Mq), np.array(tq)); CAP["g"] = g
        out = _inv(th_obj, target, **kw); CAP["wall"] = out["wall"]
        return out

    def cap_cpl(levels, init_cum, mdot_star, gamma=1.4):
        CAP["mstar"] = float(mdot_star)
        return _cpl(levels, init_cum, mdot_star, gamma)

    RA.inverse_design = wrap; MI.cplus_flux_wall = cap_cpl
    p = RA.load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = 2400; p.geometry["n_start"] = N_START
    err = None
    try:
        RA.design_chain(p)
    except Exception as e:  # noqa: BLE001  (B は壁足 x<0 で設計壁の構築が拒否される — MOC 壁表は捕捉済み)
        err = f"{type(e).__name__}: {e}"
    RA.inverse_design = _inv; MI.cplus_flux_wall = _cpl
    if "wall" not in CAP:
        res[arm] = {"error": err}; continue
    x, r, M, th = CAP["line"]; G, E = line_residuals(x, r, M, th, CAP["g"])
    w = CAP["wall"]; xw = float(w[0, 0]); xi = w[:, 0] - xw
    m = (xi > 0) & (xi < 0.08)
    ct = np.polyfit(xi[m], np.tan(w[m, 2]), 1); cr = np.polyfit(xi[m], (w[m, 1] - 1) / xi[m], 1)
    res[arm] = dict(design_chain_error=err, line_max_geom_resid_rad=float(np.abs(G).max()), line_max_compat_resid_rad=float(np.abs(E).max()),
                    line_sum_compat_resid_deg=float(np.degrees(E.sum())), line_M_min=float(M.min()), line_finite=bool(np.all(np.isfinite([x, M]))),
                    foot_x=xw, foot_M=float(w[0, 3]), foot_theta_deg=float(np.degrees(w[0, 2])), axis_x=float(x[0]), mdot_star=CAP.get("mstar"),
                    n_wall=len(w), wall_end_xi=float(xi[-1]), wall_end_r=float(w[-1, 1]),
                    start_n_pts=int(m.sum()), a_theta_deg=float(np.degrees(ct[1])), kappa_theta=float(ct[0]),
                    a_theta_fit_resid_max_deg=float(np.degrees(np.abs(np.tan(w[m, 2]) - np.polyval(ct, xi[m])).max())),
                    a_r_deg=float(np.degrees(cr[1])), kappa_r=float(2 * cr[0]), _wall=w)
A, B = res["A_hall"], res["B_consistent"]
gate = (B["line_max_geom_resid_rad"] <= 1e-8 and B["line_max_compat_resid_rad"] <= 1e-8 and B["line_M_min"] > 1 and B["line_finite"])
wa, wb = A.pop("_wall"), B.pop("_wall")
xi_a, xi_b = wa[:, 0] - wa[0, 0], wb[:, 0] - wb[0, 0]
mm = xi_b > 40
out = dict(A_hall=A, B_consistent=B, premise_gate=bool(gate),
           d_foot_M=B["foot_M"] - A["foot_M"], d_mstar_rel=(B["mdot_star"] - A["mdot_star"]) / A["mdot_star"],
           d_wall_end_xi=B["wall_end_xi"] - A["wall_end_xi"],
           dr_xi_gt40_max=float(np.abs(wb[mm, 1] - np.interp(xi_b[mm], xi_a, wa[:, 1])).max()),
           a_theta_reduction=1 - abs(B["a_theta_deg"]) / abs(A["a_theta_deg"]))
if not gate:
    out["verdict"] = "前提ゲート不成立 → 診断不能"
elif abs(B["a_theta_deg"]) <= 0.01932:
    out["verdict"] = "|a_θ| ≤ 0.01932° → 第 1 仮説 (Hall 初期線の C⁻ 不整合が主因) を支持"
elif out["a_theta_reduction"] < 0.8:
    out["verdict"] = "減少 80 % 未満 → この整合化で大半を除去できるという主因仮説を棄却"
else:
    out["verdict"] = "判定境界 (要確認)"
(C / "_band_ab/initial_line_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
