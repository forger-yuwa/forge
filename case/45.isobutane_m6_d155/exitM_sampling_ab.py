"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11e (codex diagnose 2026-10-05 exitM、事前登録): 出口コア M の標本位置 A/B と各量の時系列。
usage (AWS, case dir): python3 exitM_sampling_ab.py RUN REF_RUN → RUN/quantities_series.csv、_band_ab/exitM_sampling_ab_<RUN>.json
腕 A: 出口断面 (最終 i) の η∈[0.05,0.7] 節点の単純平均 (nozzle_report.metrics と同じ)。腕 B: 同じ M(η) を REF_RUN 最終場の出口断面の帯内 η 節点へ線形補間して単純平均。
δ_E/δ_C(x_F) は extract_and_merge (E 法、Euler 参照 run_0086) と CONTUR (c2pin_solve_fine.json の k_f) の比。"""
import csv, json, os, shutil, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field, eta_line, pspline  # noqa: E402
RUN, REF = C / sys.argv[1], C / sys.argv[2]
EU = C / "run_0086_euler_wallfit_pincal_r1_ext6k"


def core_eta(F):
    eta = F["R"][-1] / F["R"][-1, -1]; m = (eta >= 0.05) & (eta <= 0.7)
    return eta, m


Fr = load_field(REF, _res_files(REF)[-1]); eta_ref, m_ref = core_eta(Fr); eta_B = eta_ref[m_ref]


def one(res):
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    G = load_field(RUN, res); info = G["info"]; Md = float(info.get("Md", 6.0)); xE = float(info.get("x_E", 40.0)); xF = float(G["X"][-1, 0])
    xq = np.linspace(float(G["X"][0, 0]), xF, 2401); w_test = (xq >= xE + 2) & (xq <= xF - 1); w_ov = (xq >= xE - 15) & (xq <= xF)
    d = 100 * (eta_line(G, "M", 0.1, xq) / Md - 1); xx, v = xq[w_test], d[w_test]
    eta, m = core_eta(G); M = G["V"]["M"][-1]
    o = np.argsort(eta); A = float(M[m].mean()); B = float(np.interp(eta_B, eta[o], M[o]).mean())
    step = int(res.split("_")[1].split(".")[0])
    dd = Path(f"/tmp/exM_{RUN.name[:8]}_{step}"); shutil.rmtree(dd, ignore_errors=True); dd.mkdir()
    for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"): shutil.copy(RUN / f, dd / f)
    os.symlink((RUN / "nozzle.h5").resolve(), dd / "nozzle.h5"); os.symlink((RUN / res).resolve(), dd / res)
    extract_and_merge(dd, EU, band_select="edge"); nx = read_delta_r_next(dd / "delta_r_next.csv"); shutil.rmtree(dd)
    dE = float(np.interp(xF, nx["x_rt"], nx["delta_E"]))
    return dict(step=step, exitM_A=A, exitM_B=B, dM_BA=B - A, overshoot01=float(d[w_ov].max()),
                wave01=float(np.abs(v - pspline(xx, v)).max()), delta_E=dE, n_core_A=int(m.sum()), n_core_B=int(m_ref.sum()))


if __name__ == "__main__":
    from forge_design.feedback.deltastar_integral import integral_bl
    from forge_design.metrics.deltastar import smooth_delta_quintic
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    rows = list(ProcessPoolExecutor(3).map(one, _res_files(RUN)))
    p = load_problem(C / "problem_d155_ns_finemesh_pin_final.yaml"); dch = design_chain(p); xF = float(dch["wall_inv"][-1, 0])
    kf = json.load(open(C / "c2pin_solve_fine.json"))["k_f"]
    res = integral_bl(dch["wall"], dch["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), float(p.spec["r_throat"]),
                      thermal_bc=p.wall_thermal_bc_integral, cf_scale=kf)
    f_s, _ = smooth_delta_quintic(res["x"], res["delta_r"], knot_spacing=2.0, lam=1.0, positive=True); dC = float(f_s(np.r_[xF])[0])
    for r in rows: r["dE_over_dC"] = r["delta_E"] / dC
    with open(RUN / "quantities_series.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    t = rows[-5:]; dm = np.array([r["dM_BA"] for r in t])
    summ = dict(run=RUN.name, ref=REF.name, delta_C=dC, k_f=kf, tail5_steps=[r["step"] for r in t],
                tail5={k: [float(np.mean([r[k] for r in t])), float(np.ptp([r[k] for r in t]))] for k in ("exitM_A", "exitM_B", "dM_BA", "overshoot01", "wave01", "dE_over_dC")},
                all_tail_ge_0p000471=bool((dm >= 0.000471).all()), all_abs_lt_1e4=bool((np.abs(np.array([r["dM_BA"] for r in rows])) < 1e-4).all()))
    (C / "_band_ab" / f"exitM_sampling_ab_{RUN.name}.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))
