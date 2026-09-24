#!/usr/bin/env python3
"""#10b: 各 run 単体の力係数時系列を作る。plan §6.2「#10b の再判定」手順 1。

**自前で係数を組まない** — `forge_design.metrics.sern_forces` の `force_history` /
`write_force_history_csv` が本番の runner (`runner_sern.py:885`) と同じ引数・同じ規約で作る
(`check_quasisteady --series-csv` 用の CSV をそのまま吐く)。

usage: v3sern_ct.py RUN --problem PROBLEM.yaml --out CSV
"""
import argparse, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.metrics.sern_forces import force_history, write_force_history_csv
from forge_design.probdef import load_problem

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--problem", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
rd = Path(a.run)
info = json.loads((rd/"prepare_info.json").read_text())
st = info["states"]; ex, en = st["exhaust"], st["ext"]; H = info["H_m"]
p = load_problem(a.problem)
xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
hist = force_history(rd, p_a=en["P"], F_ideal=info["F_ideal_N_per_m"], H=H,
                     x_ref=float(xr)*H, y_ref=float(yr)*H,
                     mdot_u_in=ex["ro"]*ex["u"]**2*H, p_in=ex["P"],
                     twall_on_fluid=(info.get("discretization","cell") == "cell"))
write_force_history_csv(a.out, hist)
print(f"# {a.out}: {len(hist)} 行" + (f"  終端 C_T {hist[-1].get('C_T'):.7f}  C_L {hist[-1].get('C_L'):.7f}  C_M {hist[-1].get('C_M'):.5f}" if hist else "  (空)"))
