"""CFD ピンの軸 P 波増の切り分け ① (plan tooling-nozzle-cfd-pinned-initial-line §9 2026-10-05 事前登録、CFD 0 step):
pass 1 (pincal と同じ線、M_design 5.999584) と V0 の当てはめ壁の r″ 高周波残差 (局所 3 次 SG、窓 0.5 r_t) を帯別比較し、点上ゲート不合格点の x を出す。
usage: design/.venv-opt/bin/python cfdpin_noise_probe.py → _band_ab/cfdpin_noise_probe.json"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from cfd_initial_line import pinned_factory  # noqa: E402
from moc_wall_fit_ab import joint_fit  # noqa: E402
SRC = "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k"
W = {}
for nm, prob, pin in (("V0", "problem_d155_euler_c2final_n2400.yaml", False), ("pincal", "problem_d155_euler_c2final_n2400_pincal.yaml", True)):
    H = RA.HallThroat
    if pin:
        RA.HallThroat = pinned_factory(SRC, ["res_6000.h5"])
    try:
        d = RA.design_chain(RA.load_problem(C / prob))
    finally:
        RA.HallThroat = H
    s, _ = joint_fit(d["wall_inv"], d["R"], 1e-9); W[nm] = (d, s)
out = {}
xg = np.arange(1.5, 93.5, 1e-3); n = int(round(0.5 / 1e-3)) | 1
for nm, (d, s) in W.items():
    r2 = s(xg, 2); res = r2 - savgol_filter(r2, n, 3)
    tb = d["wall_inv"]; x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    dr = np.abs(s(x) - r); dth = np.abs(np.degrees(np.arctan(s(x, 1)) - th)); fail = (dr > 5e-6) | (dth > 0.005); fail[0] = False
    out[nm] = dict(hf={f"[{a},{b})": float(np.abs(res[(xg >= a) & (xg < b)]).max()) for a, b in ((2, 15), (15, 45), (45, 93))},
                   hf_rms={f"[{a},{b})": float(np.sqrt(np.mean(res[(xg >= a) & (xg < b)] ** 2))) for a, b in ((2, 15), (15, 45), (45, 93))},
                   fail_x=[float(v) for v in x[fail]], dth_max_xgt1=float(dth[x > 1].max()), dr_max_xgt1=float(dr[x > 1].max()),
                   moc_pt_spacing_15_45=float(np.median(np.diff(x[(x >= 15) & (x < 45)]))))
ratio = {k: out["pincal"]["hf"][k] / out["V0"]["hf"][k] for k in out["V0"]["hf"]}
out["ratio_pin_over_V0"] = ratio
out["line_noise_supported"] = bool(ratio["[15,45)"] >= 2.0 and any(v > 1.0 for v in out["pincal"]["fail_x"]))
out["line_noise_rejected"] = bool(ratio["[15,45)"] < 2.0 and all(v < 0.2 for v in out["pincal"]["fail_x"]))
(C / "_band_ab/cfdpin_noise_probe.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
