import os, sys, copy, time
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "/home/sano/work/forge-mknot/design")
import yaml
from pathlib import Path
from forge_design.evaluate.mknot_select import evaluate_mknot
from forge_design.intake import wt_card as WC
ROOT = Path(sys.argv[1]).parent if len(sys.argv) > 1 else Path("/home/sano/work/forge-mknot")
CARD71 = Path(os.environ.get("CARD71", "/home/sano/work/forge-m55/case/71.m55_combustion_wt/card_m55.yaml"))
p = ROOT / "case/70.m5_air_wt_demo/_mknot_card_v4_3d/problem_c2b_auto_final.yaml"
probs = {"a": (yaml.safe_load(p.read_text()), p.parent, 2.05), "b": (WC.card_to_problem(WC.load_card(CARD71), 0.1, M_knot=2.80), CARD71.parent, 2.80)}
reps = int(os.environ.get("REPS", "2"))
for be in ("numpy", "numba"):
    os.environ["FORGE_MOC_KERNEL"] = be
    for w, (raw, base, MK) in probs.items():
        ts = []
        for _ in range(reps):
            t = time.perf_counter(); row = evaluate_mknot(copy.deepcopy(raw), MK, n_axis_inv=2400, base_dir=base); ts.append(time.perf_counter() - t)
        print(f"{be} {w}: " + " ".join(f"{t:.2f}" for t in ts) + f" s (ok {row.get('ok')}, θ_max {row.get('theta_max_deg')})", flush=True)
