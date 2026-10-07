"""スロート始点の MOC 壁の分解: (r−1)/x ≈ a + κx/2、tanθ ≈ a + κx (a = 折れ角, κ = 始点曲率)。
Hall 場の流れ角を R の円弧壁上で評価して壁傾きと比べる。plan verification-m6-axis-wave-mesh-su2 §5.1 #13 (原因の切り分け、形状のみ)。
usage: design/.venv-opt/bin/python throat_kink_probe.py → _band_ab/throat_kink_probe.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.geometry.transonic import HallThroat
from forge_design.evaluate.runner_axismach import design_chain, load_problem

C = Path(__file__).resolve().parent
out = {"start_fit": {}, "hall_on_circle": []}
for n_axis in (1200, 2400):
    for ns in (21, 41, 81, 161):
        p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = n_axis; p.geometry["n_start"] = ns
        d = design_chain(p); tb = d["wall_inv"]; m = (tb[:, 0] > 0) & (tb[:, 0] < 0.08)
        x, r, th = tb[m, 0], tb[m, 1], tb[m, 2]
        cr = np.polyfit(x, (r - 1) / x, 1); ct = np.polyfit(x, np.tan(th), 1)
        out["start_fit"][f"n{n_axis}_ns{ns}"] = dict(n_pts=int(m.sum()), a_from_r_deg=float(np.degrees(cr[1])), kappa_from_r=float(2 * cr[0]),
                                                    a_from_theta_deg=float(np.degrees(ct[1])), kappa_from_theta=float(ct[0]))
R = d["R"]; g = d["gamma_hall"]; ht = HallThroat(R=R, gamma=g)
for x in (0.0129, 0.0258, 0.052, 0.1057, 0.2):
    rc = 1 + R - np.sqrt(R ** 2 - x ** 2); th_h = float(ht.theta(np.array([x]), np.array([rc]))[0]); th_w = float(np.arcsin(x / R))
    out["hall_on_circle"].append(dict(x=x, hall_theta_deg=float(np.degrees(th_h)), circle_slope_deg=float(np.degrees(th_w)), ratio=th_h / th_w))
out["R"], out["gamma_hall"] = R, g
(C / "_band_ab/throat_kink_probe.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
