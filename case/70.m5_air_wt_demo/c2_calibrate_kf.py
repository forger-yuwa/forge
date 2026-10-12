"""case/70: C2 方式の k_f の較正 (手順書 §2 ④、case/45 `c2pin_solve.py` を物理スロート径で寸法を決める場合に合わせたもの)。

NS の run から E 法 (band_select="edge") で出口の δ_E(x_F) を測り、CONTUR 積分法の摩擦係数倍率 k_f を
δ_C(x_F; k_f) = δ_E(x_F) に合わせて解く (同じ r_t のまま。r_t は要求カードの生成器が k_f を入れた問題で物理スロート径から解き直す)。
結果は `c2_calibrate_kf.json` に書き、要求カードの `overrides: {deltastar_initializer.cf_scale: k_f}` に写して YAML を作り直す。

usage (case dir, design/.venv-opt):  python c2_calibrate_kf.py NS_RUN EULER_REF_RUN PROBLEM_YAML
"""
import json
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import _gam_or_gas, design_chain, load_problem  # noqa: E402
from forge_design.feedback.deltastar_integral import integral_bl  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next  # noqa: E402
from forge_design.metrics.deltastar import smooth_delta_quintic  # noqa: E402

ns, eu, prob = C / sys.argv[1], C / sys.argv[2], C / sys.argv[3]
summ = extract_and_merge(ns, eu, band_select="edge")
nx = read_delta_r_next(ns / "delta_r_next.csv")          # δ_E は未緩和の列 (第 2 列の delta_target は ω 緩和後; 取り違えない)
p = load_problem(prob)
d = design_chain(p)
rt = float(p.spec["r_throat"])
rt_ns = float(json.loads((ns / "prepare_info.json").read_text())["scale_m"])
if abs(rt / rt_ns - 1.0) > 1e-9:
    raise SystemExit(f"問題の r_t {rt} と NS の r_t {rt_ns} が違う (同じ r_t の問題で較正する)")
xF = float(d["wall_inv"][-1, 0])
dE = float(np.interp(xF, nx["x_rt"], nx["delta_E"]))
init = (p.raw or {}).get("deltastar_initializer") or {}


def dC(kf):
    res = integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt,
                      thermal_bc=p.wall_thermal_bc_integral, cf_scale=kf)
    f_s, _ = smooth_delta_quintic(res["x"], res["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)
    return float(f_s(np.r_[xF])[0])


k0 = float(init.get("cf_scale", 1.0))
k1 = k0 * 1.03
f0, f1 = dC(k0) - dE, dC(k1) - dE
hist = [(k0, f0 + dE), (k1, f1 + dE)]
for _ in range(15):
    k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
    f2 = dC(k2) - dE
    hist.append((k2, f2 + dE))
    k0, f0, k1, f1 = k1, f1, k2, f2
    if abs(f2) < 1e-7:
        break
out = {"ns_run": ns.name, "euler_ref": eu.name, "problem": prob.name, "r_t_m": rt, "x_F_rt": xF, "delta_E_xF_rt": dE,
       "delta_C_xF_kf1_rt": hist[0][1] if hist[0][0] == 1.0 else None, "k_f": float(k1), "k_f_resid_rt": float(f1),
       "k_f_hist": hist, "extract": {k: v for k, v in summ.items() if k != "massflow"},
       "note": "δ_C(x_F; k_f) = δ_E(x_F) を同じ r_t で解いた。要求カードの overrides deltastar_initializer.cf_scale に写す"}
(C / "c2_calibrate_kf.json").write_text(json.dumps(out, indent=1, default=str))
print(json.dumps({k: out[k] for k in ("delta_E_xF_rt", "k_f", "k_f_resid_rt", "x_F_rt")}, default=str))
