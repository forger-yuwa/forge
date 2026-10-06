"""throat_moc_fix_probe.py の補足 (2026-10-07): axis_dx0 0.0075 は n_axis 2400 で壁フィルタ不合格 (半径が非単調) になるので、n_axis 4800 で 3 水準をそろえて回す。
usage: CASE_RUNS=<run_0062 のある case dir> python3 throat_moc_fix_probe_n4800.py → _band_ab/throat_moc_fix_probe.json の rows_n_axis_4800
"""
import json, os, sys, numpy as np
from pathlib import Path
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design")); sys.path.insert(0, str(C))
from forge_design.evaluate.runner_axismach import design_chain, load_problem
from forge_design.geometry.wall_axismach import joint_fit_wall
from throat_moc_restart_test import kernel
RUNS = Path(os.environ["CASE_RUNS"]); R = 2.0; RT_UM = 76.6539e3
def chain(k, dx0, ns, na):
    p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml"); p.geometry.update(n_start=ns, axis_dx0=dx0, n_axis_inv=na)
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    with kernel(k):
        return design_chain(p)
base = chain("prod", 0.03, 41, 2400)["wall"]; xe0 = float(chain("prod", 0.03, 41, 2400)["wall_inv"][-1, 0])
out = []
for k in ("prod", "K2c"):
    for dx0 in (0.03, 0.015, 0.0075):
        row = {"kernel": k, "dx0": dx0, "n_axis": 4800}
        for ns in (41, 321):
            try:
                d = chain(k, dx0, ns, 4800)
            except Exception as e:
                row[f"err{ns}"] = str(e)[:120]; continue
            tb = d["wall_inv"]; x, th = tb[:, 0], tb[:, 2]
            row[f"dth1_ns{ns}"] = float(np.degrees(th[1] - np.arctan(x[1] / R)))
            if ns == 41:
                s0, _ = joint_fit_wall(tb, R); xg = np.linspace(0, 0.3, 30001); row["bump_nomono"] = float(s0(xg, 2).max() - 0.5)
                W = d["wall"]; xx = np.linspace(-12, min(xe0, float(x[-1])), 200001)
                dr = (W.r(xx) - base.r(xx)) * RT_UM; row["max_dr_um"] = float(np.abs(dr).max()); row["x_max_dr"] = float(xx[np.argmax(np.abs(dr))])
                row["dr_exit_um"] = float(dr[-1]); m = (xx >= 40) & (xx <= 94)
                row["dth_TS_deg"] = float(np.abs(np.degrees(np.arctan(W.r(xx[m], 1)) - np.arctan(base.r(xx[m], 1)))).max())
        out.append(row); print(row, flush=True)
f = C / "_band_ab/throat_moc_fix_probe.json"; d = json.loads(f.read_text()); d["rows_n_axis_4800"] = out; f.write_text(json.dumps(d, indent=1, ensure_ascii=False))
