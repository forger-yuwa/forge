"""plan tooling-nozzle-throat-monotone-r2 の報告用: 旧壁 (problem_d155_ns_finemesh_recal_final.yaml) と単調壁 (…_mono.yaml) の
スロート近傍の r″ を、設計壁と物理壁 (prepare_ns と同じ経路: integral_delta_r → PhysicalNozzleWall) で比べる図 (CFD 0 step)。
usage: [CASE_RUNS=<run_0062 のある case dir>] python3 throat_mono_wall_fig.py [OUT_DIR (既定 _band_ab/throat_mono_compare)]
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, integral_delta_r, _gam_or_gas  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

_F = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
if os.path.exists(_F):
    font_manager.fontManager.addfont(_F); plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/throat_mono_compare"
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("CASE_RUNS", C))
KF = json.loads((C / "c2pin_solve_recal.json").read_text())["k_f"]


def walls(prob):
    p = load_problem(C / prob)
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    d = design_chain(p)
    _, drx, _ = integral_delta_r(p, d, {"model": "contur", "a_crocco": 1.0, "cf_scale": KF, "n_scale": 1.0})
    pw = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]),
                            _gam_or_gas(p), p.cp, offset="radial", delta_r_x=drx,
                            ramp=tuple(float(v) for v in p.geometry["pw_ramp"]), upstream="ramp")   # 2026-10-07 に既定が poly
    return d["wall"], pw, float(d["R"])


Wo, Po, R = walls("problem_d155_ns_finemesh_recal_final.yaml")
Wm, Pm, _ = walls("problem_d155_ns_finemesh_recal_final_mono.yaml")
x = np.linspace(-0.15, 0.3, 9001)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
PAIRS = (("設計壁 (MOC 点への当てはめ)", (Wo, Wm)), ("物理壁 (境界層ぶん外へ、生産と同じ経路)", (Po, Pm)))
for ax, (lab, (a, b)) in zip(axs, PAIRS):
    ax.plot(x, a.r(x, 2), color="#1d2b3a", lw=1.6, label="旧壁")
    ax.plot(x, b.r(x, 2), color="#0f766e", lw=1.6, ls="--", label="単調壁")
    ax.axhline(1 / R, color="#94a3b8", lw=.8, ls=":", label="1/R = 0.5")
    ax.axvline(0, color="#94a3b8", lw=.6)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(0.38, 0.52); ax.set_xlabel("x / r_t"); ax.set_ylabel("r″  [1/r_t]")
    ax.set_title(lab); ax.grid(alpha=.3); ax.legend(fontsize=9, loc="lower left")
fig.tight_layout(); fig.savefig(OUT / "fig_wall_r2_old_vs_mono.png", dpi=150); plt.close(fig)
xs = np.linspace(0, 0.3, 30001)
summ = {k: {"r2_max_0_0p3": float(w.r(xs, 2).max()), "x_r2_max": float(xs[np.argmax(w.r(xs, 2))])}
        for k, w in (("design_old", Wo), ("design_mono", Wm), ("physical_old", Po), ("physical_mono", Pm))}
xa = np.linspace(0, 95, 190001)
summ["design_max_abs_dr_um"] = float(np.abs(Wm.r(xa) - Wo.r(xa)).max() * 76.6539e3)
summ["physical_max_abs_dr_um"] = float(np.abs(Pm.r(xa) - Po.r(xa)).max() * 76.6539e3)
summ["physical_x_of_max_dr"] = float(xa[np.argmax(np.abs(Pm.r(xa) - Po.r(xa)))])
(OUT / "wall_r2_summary.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False))
print(json.dumps(summ, indent=1, ensure_ascii=False))
