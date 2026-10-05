"""P5 ②③ (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #7、C2 方式 = verification-m6 §5.1 #9 と同じ手順):
pass 1 の NS から E 法 (band_select="edge") で出口 δ_E(x_F) を測り、出口半径 0.775 m から r_t を解き (solve_rt, Re^−0.2 換算)、
新 r_t で CONTUR の k_f を δ_C(x_F; k_f) = δ_E(x_F)·(r_t,new/r_t,prev)^−0.2 に合わせて解く。結果を最終 problem に書く。
usage: python3 c2pin_solve.py PASS1_RUN EULER_REF_RUN → c2pin_solve.json, problem_d155_ns_c2pin_final{,_cond}.yaml"""
import json, re, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.feedback.deltastar_loop import extract_and_merge, solve_rt  # noqa: E402
from forge_design.feedback.deltastar_integral import integral_bl  # noqa: E402
from forge_design.metrics.deltastar import smooth_delta_quintic  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas  # noqa: E402
R_EXIT = 0.775
p1, eu = C / sys.argv[1], C / sys.argv[2]
summ = extract_and_merge(p1, eu, band_select="edge")
rs = solve_rt(C / "problem_d155_ns_c2pin.yaml", R_EXIT, prev_run=p1, euler_run=eu)
rt_new, rt_prev = rs["r_t_m"], rs["r_t_prev_m"]
nx = np.loadtxt(p1 / "delta_r_next.csv", delimiter=",", skiprows=1)
# 新 r_t の problem (形は r_t 無次元で不変)
src = (C / "problem_d155_ns_c2pin.yaml").read_text()
fin = re.sub(r"^(  r_throat:)\s*[0-9.eE+-]+", rf"\1 {rt_new:.7f}", src, count=1, flags=re.M).replace("name: isobutane_m6_d155_ns_c2pin", "name: isobutane_m6_d155_ns_c2pin_final")
(C / "problem_d155_ns_c2pin_final.yaml").write_text(fin)
p = load_problem(C / "problem_d155_ns_c2pin_final.yaml"); d = design_chain(p); xF = float(d["wall_inv"][-1, 0])
dE = float(np.interp(xF, nx[:, 0], nx[:, 1])); target = dE * (rt_new / rt_prev) ** -0.2


def dC(kf):
    res = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt_new,
                      thermal_bc=p.wall_thermal_bc_integral, cf_scale=kf)
    f_s, _ = smooth_delta_quintic(res["x"], res["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)
    return float(f_s(np.r_[xF])[0])


k0, k1 = 1.02573, 1.05; f0, f1 = dC(k0) - target, dC(k1) - target; hist = [(k0, f0), (k1, f1)]
for _ in range(12):
    k2 = k1 - f1 * (k1 - k0) / (f1 - f0); f2 = dC(k2) - target; hist.append((k2, f2))
    k0, f0, k1, f1 = k1, f1, k2, f2
    if abs(f2) < 1e-7:
        break
cond = (C / "problem_d155_ns_c2final_cond.yaml").read_text()
# cond 版: c2pin_final に cond 行を足し、step 数を c2final_cond と同じに
fc = fin.replace("name: isobutane_m6_d155_ns_c2pin_final", "name: isobutane_m6_d155_ns_c2pin_final_cond")
m_cond = re.search(r"^  condensation: .*$", cond, re.M).group(0)
fc = re.sub(r"^  tp_species: .*$", "  tp_species: split_h2o\n" + m_cond, fc, count=1, flags=re.M)
fc = re.sub(r"^(  nStepOuter:) \d+", r"\1 12000", fc, count=1, flags=re.M); fc = re.sub(r"^(  outStepInterval:) \d+", r"\1 4000", fc, count=1, flags=re.M)
(C / "problem_d155_ns_c2pin_final_cond.yaml").write_text(fc)
out = dict(pass1=str(p1.name), euler_ref=str(eu.name), extract={k: v for k, v in summ.items() if k != "massflow"}, solve_rt=rs,
           delta_E_xF=dE, delta_target_new_rt=target, k_f=float(k1), k_f_resid=float(f1), k_f_hist=hist, x_F=xF,
           R_exit_pred_m=float(rt_new * (float(d["wall_inv"][-1, 1]) + target)))
(C / "c2pin_solve.json").write_text(json.dumps(out, indent=1, default=str)); print(json.dumps({k: out[k] for k in ("delta_E_xF", "k_f", "k_f_resid", "R_exit_pred_m")}, default=str), rt_new)
