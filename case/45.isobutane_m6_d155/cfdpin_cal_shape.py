"""CFD ピン V3′ の投入前チェック (CFD 0 step): M_design 5.999584 の壁が pin 壁 (M_design 6) から r·(≤ 1e-4) の一様な縮小 + x_E/x_F の変化 ≤ 0.01 r_t か。
usage: design/.venv-opt/bin/python cfdpin_cal_shape.py [凍結元 run] → _band_ab/cfdpin_cal_shape.json"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from cfd_initial_line import pinned_factory  # noqa: E402
from moc_wall_fit_ab import joint_fit  # noqa: E402
src = sys.argv[1] if len(sys.argv) > 1 else "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k"
H = RA.HallThroat; RA.HallThroat = pinned_factory(src, ["res_6000.h5"])
W = {}
for nm, prob in (("pin", "problem_d155_euler_c2final_n2400.yaml"), ("pincal", "problem_d155_euler_c2final_n2400_pincal.yaml")):
    d = RA.design_chain(RA.load_problem(C / prob)); s, _ = joint_fit(d["wall_inv"], d["R"], 1e-9); W[nm] = (d, s)
RA.HallThroat = H
(dA, sA), (dB, sB) = W["pin"], W["pincal"]
xx = np.linspace(0, min(dA["wall_inv"][-1, 0], dB["wall_inv"][-1, 0]) - 0.01, 20001)
rel = (sB(xx) - sA(xx)) / sA(xx)
out = dict(x_E=[dA["x_E"], dB["x_E"]], x_F=[float(dA["wall_inv"][-1, 0]), float(dB["wall_inv"][-1, 0])], r_F=[float(dA["wall_inv"][-1, 1]), float(dB["wall_inv"][-1, 1])],
           rel_dr_min=float(rel.min()), rel_dr_max=float(rel.max()), rel_dr_exit=float(rel[-1]))
out["ok"] = bool(np.abs(rel).max() <= 1e-4 and abs(out["x_F"][1] - out["x_F"][0]) <= 0.01 and abs(out["x_E"][1] - out["x_E"][0]) <= 0.01)
(C / "_band_ab/cfdpin_cal_shape.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
