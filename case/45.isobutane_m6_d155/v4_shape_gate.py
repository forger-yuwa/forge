"""V4 壁の段 1 形状ゲート (事前登録: plan verification-m6-axis-wave-mesh-su2 §5.1 #16, commit 7556edae)。CFD 0 step。
F1 MOC 出力点 i≥1 で |Δr| ≤ 5e-6 r_t・|Δθ| ≤ 0.005° / F2 r″ が [0,0.3] に内部極値なし・r‴ max ≤ 1.0 / F3 λ ∈ {1e-10,1e-9,1e-8} で r″(0) 幅 ≤ 0.05・r′(0) 幅 ≤ 0.04°
F4 縮流部 μ ≤ 20、x_min ∈ [−0.01, 0]、1 − r_min ≤ 1e-5、x<x_min で r′ ≤ 0、x=0 で C⁰/C¹/C² / F5 (併記) PhysicalNozzleWall (δ_r run_0051) の κ_t・r″ の山・r‴ max、第 1 区間の自己不整合 c₀。
usage: design/.venv-opt/bin/python v4_shape_gate.py → _band_ab/v4_shape_gate.json
"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas, delta_r_from_table  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402
from wall_v4 import build_v4  # noqa: E402

p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); d = design_chain(p); tb = d["wall_inv"]; x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
W = {lam: build_v4(d, lam) for lam in (1e-10, 1e-9, 1e-8)}
w = W[1e-9]; out = {"fit": w.fit_info}
dr = np.abs(w.r(x) - r)[1:]; dth = np.abs(np.degrees(np.arctan(w.r(x, 1)) - th))[1:]
out["F1"] = dict(dr_max=float(dr.max()), dth_max_deg=float(dth.max()), n_fail=int((dr > 5e-6).sum() + (dth > 0.005).sum()), ok=bool(dr.max() <= 5e-6 and dth.max() <= 0.005))
xg = np.linspace(0, 0.3, 30001); r2 = w.r(xg, 2); r3 = w.r(xg, 3)
n_int = int(np.count_nonzero(np.diff(np.sign(np.diff(r2))) != 0))
out["F2"] = dict(r2_interior_extrema=n_int, r3_max=float(np.abs(r3).max()), ok=bool(n_int == 0 and np.abs(r3).max() <= 1.0))
r2s = [W[l].fit_info["r2_0"] for l in W]; r1s = [np.degrees(np.arctan(W[l].fit_info["r1_0"])) for l in W]
out["F3"] = dict(r2_0=r2s, r1_0_deg=r1s, r2_width=float(np.ptp(r2s)), r1_width_deg=float(np.ptp(r1s)),
                 ok=bool(np.ptp(r2s) <= 0.05 and np.ptp(r1s) <= 0.04))
L = w.up.L_U; xu = np.linspace(-L, 0.0, 200001); ru = w.r(xu); i = int(np.argmin(ru)); e = 1e-7
jumps = [float(w.r(np.r_[e], k)[0] - w.r(np.r_[-e], k)[0]) for k in range(3)]
out["F4"] = dict(mu=w.up.mu, x_min=float(xu[i]), one_minus_rmin=float(1 - ru[i]), r1_max_before_xmin=float(w.r(xu[:i], 1).max()),
                 jump_r_r1_r2_at0=jumps,
                 ok=bool(w.up.mu <= 20 and -0.01 <= xu[i] <= 0 and 1 - ru[i] <= 1e-5 and w.r(xu[:i], 1).max() <= 0 and max(abs(j) for j in jumps) < 1e-5))
# 全区間の r″ の極値数 (参考) と補間壁
xa = np.linspace(0, x[-1], 400001)
out["ref"] = dict(r2_extrema_0_xe_v4=int(np.count_nonzero(np.diff(np.sign(np.diff(w.r(xa, 2)))) != 0)),
                  r2_extrema_0_xe_interp=int(np.count_nonzero(np.diff(np.sign(np.diff(d["wall"].r(xa, 2)))) != 0)))
c0 = abs(np.degrees(np.arctan((r[1] - r[0]) / (x[1] - x[0])) - 0.5 * (th[0] + th[1])))
dtab = np.loadtxt(C / "_band_ab/delta_r_c2final_run0051.csv", delimiter=",", skiprows=1); f = delta_r_from_table(dtab[:, 0], dtab[:, 1])
F5 = {"c0_first_segment_deg": float(c0)}
for nm, wall in (("interp", d["wall"]), ("v4", w)):
    PW = PhysicalNozzleWall(wall, tb, float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp, offset="radial", delta_r_x=f)
    xp = np.linspace(PW.x_throat + 1e-4, PW.x_throat + 0.3, 20001)
    F5[nm] = dict(kappa_t=PW.kappa_throat, r2_bump=float(PW.r(xp, 2).max() - PW.r(xp[:1], 2)[0]), r3_max=float(np.abs(PW.r(xp, 3)).max()), validate=PW.validate())
out["F5"] = F5
out["all_ok"] = all(out[k]["ok"] for k in ("F1", "F2", "F3", "F4"))
(C / "_band_ab/v4_shape_gate.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
