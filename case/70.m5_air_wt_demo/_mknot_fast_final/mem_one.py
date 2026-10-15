import os, sys, copy
sys.path.insert(0, sys.argv[1])
import yaml
from pathlib import Path
from forge_design.evaluate.mknot_select import evaluate_mknot
from forge_design.intake import wt_card as WC
ROOT = Path(sys.argv[1]).parent
which, N = sys.argv[2], int(sys.argv[3])
if which == "a":
    p = ROOT / "case/70.m5_air_wt_demo/_mknot_card_v4_3d/problem_c2b_auto_final.yaml"; raw, base, MK = yaml.safe_load(p.read_text()), p.parent, 2.05
else:
    c = Path(os.environ["CARD71"]); raw, base, MK = WC.card_to_problem(WC.load_card(c), 0.1, M_knot=2.80), c.parent, 2.80
row = evaluate_mknot(copy.deepcopy(raw), MK, n_axis_inv=N, base_dir=base, L_total_tol=(0.1 if N == 600 else None))
print(which, N, "ok", row.get("ok"))
