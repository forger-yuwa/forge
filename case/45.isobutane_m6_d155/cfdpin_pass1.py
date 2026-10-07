"""CFD ピン pass 1 の形状判定 = plan tooling-nozzle-cfd-pinned-initial-line §6 V0b・V1 (事前登録, CFD 0 step)。
V0 の Euler 場 (run_0062 res_6000) から凍結した初期線・m*・(M_A, M′_A) で design_chain を回し、Hall 版 (現行) と比べる。
V0b law − 実測軸 (同じ V0 場の evenfit) ≤ 2e-3 on [x₀+0.3, x_K] / V1 c₀ ≤ 0.015°・a_θ ≤ 0.06°・V0 型当てはめの点上ゲート (i≥1) |Δr| ≤ 5e-6・|Δθ| ≤ 0.005°・r″ の山 (x∈[0,0.3]) < 0.534
usage: design/.venv-opt/bin/python cfdpin_pass1.py [V0 場の run] → _band_ab/cfdpin_pass1.json
"""
import json, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from cfd_initial_line import pinned_factory  # noqa: E402
from moc_wall_fit_ab import joint_fit  # noqa: E402

src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k")
PROB = C / "problem_d155_euler_c2final_n2400.yaml"
_cpl = MI.cplus_flux_wall; CAP = {}


def cap(levels, init_cum, mdot_star, gamma=1.4):
    CAP["mstar"] = float(mdot_star); return _cpl(levels, init_cum, mdot_star, gamma)


def run(pinned):
    MI.cplus_flux_wall = cap
    H = RA.HallThroat
    if pinned:
        RA.HallThroat = pinned_factory(src, ["res_6000.h5"])
    try:
        d = RA.design_chain(RA.load_problem(PROB)); d["mstar"] = CAP["mstar"]
        d["ht"] = RA.HallThroat(R=d["R"], gamma=d["gamma_hall"]) if pinned else None
    finally:
        RA.HallThroat = H; MI.cplus_flux_wall = _cpl
    return d


def wall_metrics(d):
    tb = d["wall_inv"]; x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    c0 = abs(np.degrees(np.arctan((r[1] - r[0]) / (x[1] - x[0])) - 0.5 * (th[0] + th[1])))
    mk = (x > 0) & (x < 0.08); a_th = float(np.degrees(np.polyfit(x[mk], np.tan(th[mk]), 1)[1]))
    s, _ = joint_fit(tb, d["R"], 1e-9)
    dr = np.abs(s(x) - r)[1:]; dth = np.abs(np.degrees(np.arctan(s(x, 1)) - th))[1:]
    xg = np.linspace(0, 0.3, 30001)
    return dict(c0_deg=float(c0), a_theta_deg=a_th, fit_dr_max=float(dr.max()), fit_dth_max_deg=float(dth.max()),
                fit_n_fail=int((dr > 5e-6).sum() + (dth > 0.005).sum()), fit_r2_max_0_0p3=float(s(xg, 2).max()), x_F=float(x[-1]), r_F=float(r[-1])), s


dH = run(False); dP = run(True)
mH, sH = wall_metrics(dH); mP, sP = wall_metrics(dP)
import json as _j
axs, _, _ = RA.axis_curve_node(src, float(_j.loads((src / "prepare_info.json").read_text())["scale_m"]))   # 実測軸 (同じ V0 場の evenfit、全長)


def law_gap(d):
    xs = np.linspace(d["x0"] + 0.3, d["x_K"], 400)
    return float(np.abs(np.array([float(d["law"](x)) for x in xs]) - axs(xs)).max())


gap, gapH = law_gap(dP), law_gap(dH)
dm = dP["mstar"] / dH["mstar"] - 1
xx = np.linspace(0, min(mH["x_F"], mP["x_F"]) - 0.01, 20001)
dshape = sP(xx) - sH(xx) * (1 + 0.5 * dm)       # m* による一様な膨らみ (r ∝ √m*) を除いた形状差
bands = {f"[{a},{b})": dict(maxabs=float(np.abs(dshape[(xx >= a) & (xx < b)]).max()), min=float(dshape[(xx >= a) & (xx < b)].min()), max=float(dshape[(xx >= a) & (xx < b)].max())) for a, b in ((0, 0.5), (0.5, 5), (5, 15), (15, 40), (40, 95))}
out = dict(source=str(src), hall=dict(mH, x0=dH["x0"], anchor=list(dH["anchor"]), mstar=dH["mstar"], law_gap_max=gapH),
           pinned=dict(mP, x0=dP["x0"], anchor=list(dP["anchor"]), mstar=dP["mstar"], law_gap_max=gap, x_E=dP["x_E"], x_K=dP["x_K"]),
           mstar_rel=float(dm), shape_change_bands=bands, rF_change=mP["r_F"] - mH["r_F"])
g = dict(V0b_law_gap=out["pinned"]["law_gap_max"] <= 2e-3, V1_c0=mP["c0_deg"] <= 0.015, V1_a_theta=abs(mP["a_theta_deg"]) <= 0.06,
         V1_point_gate=mP["fit_n_fail"] == 0, V1_r2_bump=mP["fit_r2_max_0_0p3"] < mH["fit_r2_max_0_0p3"])
out["gates"] = {k: bool(v) for k, v in g.items()}; out["all_pass"] = all(out["gates"].values())
(C / "_band_ab/cfdpin_pass1.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
