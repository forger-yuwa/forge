"""case/45 H (problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml) の参考の見積もり step 1: 問題 YAML の写し (初期線は Hall に置き換え:
CFD ピンの run は別のツリーにしかないので読まない) → design_chain → euler_reference_wall (H の上流を付けた設計壁)。
estimate_mesh と同じ式で y1,req(x) を作る (Tw 370 K 等温)。使い捨て。出力は c45/ の下 (step2_est.py --dir c45 で読める形)。"""
import sys, json, time, pickle
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path.home() / "forge-mknot/design"))
from forge_design.intake import wt_card as W
from forge_design.probdef import load_problem
from forge_design.evaluate.runner_axismach import design_chain, euler_reference_wall
from forge_design.evaluate.ic import _invert_area_ratio_gas, _isentropic_pressure_ratio_gas
from forge_design.gas import GasSemiPerfect
from forge_design.gas.composition import ResolvedSpeciesDB

OUT = Path.home() / "m55_mesh_est" / "c45"
p = load_problem(OUT / "problem_h_hall.yaml")
rt = float(p.spec["r_throat"])
t0 = time.time()
d = design_chain(p)
wall, rec = euler_reference_wall(p, d, rt)
print(f"chain {time.time() - t0:.1f}s  x_in {wall.x_in} x_e {wall.x_e}  design x_in {d['wall'].x_in}", flush=True)
pickle.dump(wall, open(OUT / "wall.pkl", "wb"))
import yaml
_py = yaml.safe_load(open(OUT / "problem_h_hall.yaml"))
Y = dict(_py["gas"]["species"])
json.dump(_py["mesh"]["wall_first_frac_table"], open(OUT / "h_frac_table.json", "w"))
Pt, Tt = float(p.spec["Pt"]), float(p.spec["Tt"])
Twall = float(p.spec["wall_thermal"]["Tw"])
gas = GasSemiPerfect(Y, Tt=Tt, db=ResolvedSpeciesDB.builtin())
# ---- estimate_mesh と同じ式 (wt_card.py の estimate_mesh、等温壁一様)
x = np.linspace(float(wall.x_in), float(wall.x_e), 8001)
rw = wall.r(x)
s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(rw)))]) * rt
AR = np.maximum(rw, 1.0) ** 2
Mx = _invert_area_ratio_gas(AR, x >= 0.0, gas)
T = np.where(Mx >= 1.0, gas.T_of_M(np.maximum(Mx, 1.0)), np.interp(np.minimum(Mx, 1.0), gas._Mu, gas._Tu))
P = Pt * _isentropic_pressure_ratio_gas(T, gas)
Rg = float(gas.R)
u = Mx * np.sqrt(gas.gamma(T) * Rg * T)
Tw = np.full_like(T, Twall)
Ts = T * (1.0 + 0.032 * Mx ** 2 + 0.58 * (Tw / T - 1.0))
rho_s, mu_s = P / (Rg * Ts), W._mu_sutherland(Ts)
rho_w, mu_w = P / (Rg * Tw), W._mu_sutherland(Tw)
with np.errstate(divide="ignore", invalid="ignore"):
    Re = rho_s * u * s / mu_s
    Cf = 0.0592 * Re ** -0.2
    tau = Cf * 0.5 * rho_s * u * u
    utau = np.sqrt(tau / rho_w)
    y1 = W.Y1PLUS_TARGET * mu_w / (rho_w * utau)
frac_req = y1 / (rw * rt)
np.savez(OUT / "sections.npz", x_rt=x, x_m=x * rt, Tw_K=Tw, y1_m=y1, y1_over_rw=frac_req)
m_dn, m_th = x > W.X_DOWNSTREAM_MIN, np.abs(x) <= W.X_THROAT_HALF
est = {"wall_first_frac": None, "wall_first_frac_throat": None, "wall_first_up_x0": None, "wall_first_up_x1": None,
       "frac_req_min_downstream": float(frac_req[m_dn].min()), "frac_req_min_throat": float(frac_req[m_th].min()),
       "note": "case/45 H の参考 (A は計算しない)"}
json.dump(est, open(OUT / "est.json", "w"), indent=1)
json.dump({"mesh": {"nj": 160}, "spec": {"r_throat": rt}}, open(OUT / "prob_parts.json", "w"), indent=1)
print(est, flush=True)
