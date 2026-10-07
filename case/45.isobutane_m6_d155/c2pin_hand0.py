"""P5 手 0 (plan tooling-nozzle-cfd-pinned-initial-line §9 2026-10-05 事前登録、CFD 0 step): run_0090 の M 不足が δ 超過の 1D 換算で説明できるか。
ΔM(x) = M_NS − M_Euler(ピン run_0086) at r/r_w = 0.1、ΔM_pred = −M·(2Δδ/r_w)/[(M²−1)/(1+(γ−1)M²/2)]、Δδ = δ_E − δ_C [r_t]。
usage: design/.venv-opt/bin/python c2pin_hand0.py [case_dir] → _band_ab/c2pin_hand0.json"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402
from forge_design.feedback.deltastar_loop import read_delta_r_next  # noqa: E402
D = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
NS, EU = D / "run_0090_ns_c2pin_final", D / "run_0086_euler_wallfit_pincal_r1_ext6k"
# δ_E = 未緩和の平滑化抽出 (delta_E 列)。2026-10-05 訂正 (codex result M1): 以前は第 2 列 (ω 緩和後の δ_target) を読んでいた
nx = read_delta_r_next(NS / "_extract_edge/delta_r_next.csv")
di = np.loadtxt(NS / "delta_r_initial.csv", delimiter=",", skiprows=1)
FN, FE = load_field(NS, "res_12000.h5"), load_field(EU, "res_6000.h5")
G = 1.274
xF = float(json.load(open(D / "c2pin_solve.json"))["x_F"])
rows = []
for x in (20.0, 40.0, 60.0, 80.0, xF - 0.5):
    dd = float(np.interp(x, nx["x_rt"], nx["delta_E"]) - np.interp(x, di[:, 0], di[:, 1]))
    rw = float(np.interp(x, FE["X"][:, -1], FE["R"][:, -1]))
    Mn = float(eta_line(FN, "M", 0.1, np.r_[x])[0]); Me = float(eta_line(FE, "M", 0.1, np.r_[x])[0])
    am = (Me ** 2 - 1) / (1 + 0.5 * (G - 1) * Me ** 2)
    pred = -Me * (2 * dd / rw) / am; dm = Mn - Me
    rows.append(dict(x=x, d_delta_rt=dd, d_delta_rel=dd / float(np.interp(x, di[:, 0], di[:, 1])), r_w=rw, M_euler=Me, dM=dm, dM_pred=pred, ratio=dm / pred if pred else None))
ok = all(r["ratio"] is not None and 0.7 <= r["ratio"] <= 1.3 for r in rows)
out = dict(rows=rows, all_ratio_in_0p7_1p3=bool(ok))
(C / "_band_ab/c2pin_hand0.json").write_text(json.dumps(out, indent=1))
for r in rows:
    print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
print("all in [0.7,1.3]:", ok)
