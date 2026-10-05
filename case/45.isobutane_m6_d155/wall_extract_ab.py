"""壁抽出の A/B (MOC 網を固定、CFD 0 step)。plan verification-m6-axis-wave-mesh-su2 §5.1 #13 ⑥ (事前登録 commit 828e12aa)。
A: 現行 cplus_flux_wall (C⁺ 線上の累積流束が m* に達する点)。
B: 始点 (0,1,θ_Hall) から各 C⁺ 線との交点を流線の台形式 rᵢ−rᵢ₋₁ = ½(tanθᵢ₋₁+tanθᵢ)(xᵢ−xᵢ₋₁) で順に決める (状態は網の線分から線形補間)。
主指標: 0<x<0.08 の a_θ (tanθ ≈ a+κx)。判定: B で |a_θ| が 80 % 以上減る → 流束閉包側を優先 / 減らない → 「壁抽出の変更だけで解消」を棄却。
usage: design/.venv-opt/bin/python wall_extract_ab.py → _band_ab/wall_extract_ab.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.geometry import moc_inverse as MI
from forge_design.evaluate.runner_axismach import design_chain, load_problem

C = Path(__file__).resolve().parent
CAP = {}
_orig = MI.cplus_flux_wall


def _capture(levels, init_cum, mdot_star, gamma=1.4):
    CAP.update(levels=levels, init_cum=np.array(init_cum), mdot_star=mdot_star, gamma=gamma)
    return _orig(levels, init_cum, mdot_star, gamma)


MI.cplus_flux_wall = _capture
p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = 2400; p.geometry["n_start"] = 161
d = design_chain(p)
MI.cplus_flux_wall = _orig
lev, cum0, mstar, g = CAP["levels"], CAP["init_cum"], CAP["mdot_star"], CAP["gamma"]
n = lev.shape[1]
wallA = _orig(lev, cum0, mstar, g)


def seg_flux(x, r, th, M):
    Mm = 0.5 * (M[1:] + M[:-1]); thm = 0.5 * (th[1:] + th[:-1]); rm = 0.5 * (r[1:] + r[:-1])
    return MI._mass_flux_density(Mm, g) * (np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x)) * 2 * np.pi * rm


start = MI.cplus_lines(lev, n - 1)[0]                     # 初期値線の壁足 (0,1)
pts = [(float(start[0]), float(start[1]), float(start[2]), float(start[4]), 0.0)]
stop_reason = None
for m in range(n - 2, -1, -1):
    line = MI.cplus_lines(lev, m)
    if len(line) < 2:
        continue
    x, r, th, nu, M = (line[:, i] for i in range(5))
    x0, r0, t0 = pts[-1][0], pts[-1][1], pts[-1][2]
    f = r - r0 - 0.5 * (np.tan(t0) + np.tan(th)) * (x - x0)
    ok = x > x0 - 1e-12
    idx = [j for j in range(len(f) - 1) if ok[j + 1] and f[j] < 0 <= f[j + 1]]
    if not idx:
        stop_reason = f"C+ 線 m={m} と交わらない"
        break
    j = idx[0]
    lo, hi = 0.0, 1.0
    for _ in range(80):                                    # 線分内の二分法 (θ も線形補間するので f は非線形)
        s_ = 0.5 * (lo + hi)
        xs, rs, ts = x[j] + s_ * (x[j + 1] - x[j]), r[j] + s_ * (r[j + 1] - r[j]), th[j] + s_ * (th[j + 1] - th[j])
        if rs - r0 - 0.5 * (np.tan(t0) + np.tan(ts)) * (xs - x0) < 0:
            lo = s_
        else:
            hi = s_
    s_ = 0.5 * (lo + hi)
    P = [v[j] + s_ * (v[j + 1] - v[j]) for v in (x, r, th, M)]
    # 流束残差: 現行の流束式で B の交点までの累積
    dd = seg_flux(x[:j + 1], r[:j + 1], th[:j + 1], M[:j + 1])
    xe, re_, te, Me = np.r_[x[j], P[0]], np.r_[r[j], P[1]], np.r_[th[j], P[2]], np.r_[M[j], P[3]]
    F = cum0[m] + dd.sum() + seg_flux(xe, re_, te, Me).sum()
    pts.append((P[0], P[1], P[2], P[3], float(F - mstar)))
    if P[0] > 1.0:
        break
wallB = np.array(pts)


def start_fit(w):
    m = (w[:, 0] > 0) & (w[:, 0] < 0.08)
    x, r, th = w[m, 0], w[m, 1], w[m, 2]
    cr = np.polyfit(x, (r - 1) / x, 1); ct = np.polyfit(x, np.tan(th), 1)
    sec = np.degrees(np.arctan(np.diff(w[:3, 1]) / np.diff(w[:3, 0])) - np.arctan(np.tan(0.5 * (w[1:3, 2] + w[:2, 2]))))
    return dict(n_pts=int(m.sum()), a_theta_deg=float(np.degrees(ct[1])), kappa_theta=float(ct[0]), a_r_deg=float(np.degrees(cr[1])),
                kappa_r=float(2 * cr[0]), first_secant_minus_mean_deg=[float(v) for v in sec])


fa, fb = start_fit(wallA), start_fit(wallB)
xc = wallB[(wallB[:, 0] > 0) & (wallB[:, 0] <= min(1.0, wallA[-1, 0])), 0]
out = dict(settings=dict(n_axis_inv=2400, n_start=161, mdot_star=float(mstar), n_levels=list(lev.shape)),
           A_flux_closure=fa, B_streamline=fb, stop_reason=stop_reason, B_n_pts=len(wallB), B_x_end=float(wallB[-1, 0]),
           B_flux_resid_max_abs=float(np.abs(wallB[1:, 4]).max()), B_flux_resid_rel_max=float(np.abs(wallB[1:, 4]).max() / mstar),
           B_flux_resid_first5=[float(v) for v in wallB[1:6, 4]],
           pos_diff_B_minus_A_max=float(np.abs(wallB[np.isin(wallB[:, 0], xc), 1] - np.interp(xc, wallA[:, 0], wallA[:, 1])).max()) if len(xc) else None,
           theta_diff_B_minus_A_max_deg=float(np.degrees(np.abs(wallB[np.isin(wallB[:, 0], xc), 2] - np.interp(xc, wallA[:, 0], wallA[:, 2])).max())) if len(xc) else None)
red = 1 - abs(fb["a_theta_deg"]) / abs(fa["a_theta_deg"])
out["a_theta_reduction"] = red
out["verdict"] = ("B で |a_θ| が 80 % 以上減少 → 第 1 仮説 (角度場の不整合が支配) を撤回し流束閉包側を優先" if red >= 0.8
                  else "減少 80 % 未満 → 第 2 仮説 (壁抽出の変更だけで解消) を棄却 (第 1 仮説の確定ではない)")
(C / "_band_ab/wall_extract_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
