"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11d (codex diagnose 2026-10-05, 事前登録): ② の計算の 0 step A/B。
変えるのは solve_rt と dC(k_f) に渡す未緩和 δ(x_F) だけ (A = 粗 run_0105 0.725276、B = 細 run_0103 0.733096)。forge は起動しない。
共通: r_t,prev 76.7531 mm、problem_d155_ns_finemesh_pin.yaml (run_0092 の壁)、CONTUR 積分法・平滑化は c2pin_solve.py と同じ。"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.feedback.deltastar_integral import integral_bl  # noqa: E402
from forge_design.metrics.deltastar import smooth_delta_quintic  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas  # noqa: E402
R_EXIT, S_PREV = 0.775, 0.0767531
p0 = load_problem(C / "problem_d155_ns_finemesh_pin.yaml"); d = design_chain(p0)
rF, xF = float(d["wall_inv"][-1, 1]), float(d["wall_inv"][-1, 0])
out = {"x_F": xF, "r_F": rF}
for arm, dm in (("A_coarse_run0105", 0.725276), ("B_fine_run0103", 0.733096)):
    rt = S_PREV; it = []
    for _ in range(6):
        de = dm * (rt / S_PREV) ** -0.2; rn = R_EXIT / (rF + de); it.append((rt, de, rn)); rt = rn
    target = dm * (rt / S_PREV) ** -0.2

    def dC(kf):
        res = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p0), p0.cp, float(p0.spec["Pt"]), float(p0.spec["Tt"]), rt,
                          thermal_bc=p0.wall_thermal_bc_integral, cf_scale=kf)
        f_s, _ = smooth_delta_quintic(res["x"], res["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)
        return float(f_s(np.r_[xF])[0])
    k0, k1 = 1.02573, 1.05; f0, f1 = dC(k0) - target, dC(k1) - target; hist = [(k0, f0), (k1, f1)]
    for _ in range(20):
        k2 = k1 - f1 * (k1 - k0) / (f1 - f0); f2 = dC(k2) - target; hist.append((k2, f2))
        k0, f0, k1, f1 = k1, f1, k2, f2
        if abs(f2) < 1e-9:
            break
    out[arm] = dict(delta_in=dm, r_t_mm=rt * 1e3, delta_target=target, k_f=k1, k_f_resid_rel=f1 / target, k_f_hist=hist,
                    R_exit_check_m=rt * (rF + target))
    print(arm, f"r_t {rt*1e3:.5f} mm  k_f {k1:.6f}  resid_rel {f1/target:.2e}  R {rt*(rF+target):.7f}")
(C / "_band_ab" / "c2pin_ab_delta.json").write_text(json.dumps(out, indent=1))
