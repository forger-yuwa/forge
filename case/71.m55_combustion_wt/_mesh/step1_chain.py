"""M5.5 (case/71) 見積もり step 1: カード (M_knot 2.80 固定、mesh.ni 外し) → 問題 → 設計チェーン 1 回 → estimate_mesh。
使い捨て (コードは変えない)。出力: est.json, sections.npz, wall_dense.npz, wall.pkl"""
import sys, json, time, pickle
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path.home() / "forge-mknot/design"))
from forge_design.intake import wt_card as W

OUT = Path.home() / "m55_mesh_est"
RT = 0.09990877093936801          # card_gen.log の r_t 反復 1 の固定点 (M_knot 2.80)
t0 = time.time()
card = W.load_card(OUT / "card_est.yaml")
prob = W.card_to_problem(card, RT)
print("geometry:", {k: prob["geometry"][k] for k in ("R", "r_inlet", "L_U", "L_pipe", "M_knot", "L_total_tol")}, flush=True)
print("mesh:", prob["mesh"], flush=True)
chain = W._design_chain_of(prob)
wall = chain["wall"]
print(f"chain {time.time() - t0:.1f}s  x_in {wall.x_in} x_e {wall.x_e} L_c {chain.get('L_c')}", flush=True)
try:
    pickle.dump(wall, open(OUT / "wall.pkl", "wb"))
    print("wall pickled", flush=True)
except Exception as ex:
    print("wall pickle failed:", ex, flush=True)
xd = np.linspace(float(wall.x_in), float(wall.x_e), 400001)
rd = np.asarray(wall.r(xd), dtype=float)
try:
    dr = np.asarray(wall.r(xd, 1), dtype=float)
except Exception:
    dr = np.gradient(rd, xd)
np.savez(OUT / "wall_dense.npz", x=xd, r=rd, dr=dr, x_in=float(wall.x_in), x_e=float(wall.x_e), L_U=float(wall.up.L_U))
t1 = time.time()
est = W.estimate_mesh(card, prob, chain)
print(f"estimate_mesh {time.time() - t1:.1f}s", flush=True)
sec = est.pop("sections")
np.savez(OUT / "sections.npz", **{k: np.asarray(v) for k, v in sec.items()})
json.dump(W._json_safe(est), open(OUT / "est.json", "w"), ensure_ascii=False, indent=1)
json.dump(W._json_safe({"mesh": prob["mesh"], "geometry": prob["geometry"], "spec": prob["spec"]}), open(OUT / "prob_parts.json", "w"),
          ensure_ascii=False, indent=1)
for k in ("wall_first_frac", "wall_first_frac_throat", "ni", "ni_required_for_ar", "blend_choice", "wall_first_up_x0", "wall_first_up_x1",
          "wall_first_frac_raw", "x_wall_first_frac_rt", "wall_first_frac_throat_raw", "x_wall_first_frac_throat_rt", "ar_max_wall_fast",
          "y1plus_est_max_with_blend", "x_y1plus_est_max_rt", "mesh_check", "nodes", "warnings"):
    print(k, ":", est.get(k), flush=True)
print(f"total {time.time() - t0:.1f}s", flush=True)
