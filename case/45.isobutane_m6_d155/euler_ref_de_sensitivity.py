"""plan verification-case45-euler-total-enthalpy §6 E5: δ_E・C2 の参照場の感度 (2026-10-07 登録、plan 段レビュー M1)。

NS の 3 条件の dry の最終場を固定し、Euler の参照を本段 42000〜54000 の 13 枚の 1 枚ずつに差し替えて δ_E/δ_C(x_F) を求める。
経路は exitM_sampling_ab.py の one() と同じ (extract_and_merge(..., band_select="edge")、δ_C は c2pin_solve_recal.json の k_f と各条件の問題)。
変えるのは参照の Euler のスナップショットだけ (Euler の run の写しに prepare_info.json・nozzle.h5 と 1 枚の res だけを並べる)。

  前提: 54000 の 1 枚の値が、NS の連鎖の quantities_series.csv の step 80000 の dE_over_dC とビット一致 (経路の同一性)。
  合格: 各条件で 13 枚の幅 (最大 − 最小) ≤ 5e-4 (NS のゲート 1 ± 0.5 % の 1/10)。13 枚の平均 − 54000 の値も記録する。

usage (AWS の case dir で): python3 euler_ref_de_sensitivity.py → _band_ab/euler_ref_de_sensitivity.json
"""
import csv
import json
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))

CONDS = {"N0": ("run_0165_ns_n012_N0", "run_0164_euler_e4_recal_d1", "problem_d155_ns_n012_N0.yaml"),
         "N1": ("run_0166_ns_n012_N1", "run_0164_euler_e4_recal_d1", "problem_d155_ns_n012_N1.yaml"),
         "N2": ("run_0167_ns_n012_N2", "run_0174_euler_v5d_M_r1", "problem_d155_ns_n012_N2.yaml")}
NS_RES = "res_80000.h5"
EU_STEPS = list(range(42000, 54001, 1000))
SOLVE = "c2pin_solve_recal.json"
TOL = 5e-4


def delta_E(args):
    """exitM_sampling_ab.one() の δ_E の部分と同じ手順。参照は Euler の run の写し (res は 1 枚だけ)。"""
    ns_run, eu_run, step = args
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    from forge_design.report.nozzle_report import load_field
    RUN, EU = C / ns_run, C / eu_run
    G = load_field(RUN, NS_RES)
    xF = float(G["X"][-1, 0])
    with tempfile.TemporaryDirectory(prefix=f"e5_{ns_run[:8]}_{step}_") as td:
        td = Path(td)
        eu = td / "eu"; eu.mkdir()
        for f in ("prepare_info.json", "nozzle.h5"):
            os.symlink((EU / f).resolve(), eu / f)
        os.symlink((EU / f"res_{step}.h5").resolve(), eu / f"res_{step}.h5")
        dd = td / "ns"; dd.mkdir()
        for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
            shutil.copy(RUN / f, dd / f)
        os.symlink((RUN / "nozzle.h5").resolve(), dd / "nozzle.h5")
        os.symlink((RUN / NS_RES).resolve(), dd / NS_RES)
        extract_and_merge(dd, eu, band_select="edge")
        nx = read_delta_r_next(dd / "delta_r_next.csv")
    return ns_run, step, float(np.interp(xF, nx["x_rt"], nx["delta_E"]))


def delta_C(problem: str) -> float:
    """exitM_sampling_ab.py の __main__ と同じ δ_C(x_F)。"""
    from forge_design.feedback.deltastar_integral import integral_bl
    from forge_design.metrics.deltastar import smooth_delta_quintic
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas
    p = load_problem(C / problem); dch = design_chain(p); xF = float(dch["wall_inv"][-1, 0])
    kf = json.load(open(C / SOLVE))["k_f"]
    res = integral_bl(dch["wall"], dch["wall_inv"], _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), float(p.spec["r_throat"]),
                      thermal_bc=p.wall_thermal_bc_integral, cf_scale=kf)
    f_s, _ = smooth_delta_quintic(res["x"], res["delta_r"], knot_spacing=2.0, lam=1.0, positive=True)
    return float(f_s(np.r_[xF])[0])


def chain_value(ns_run: str):
    with open(C / ns_run / "quantities_series.csv") as f:
        rows = [r for r in csv.DictReader(f) if int(r["step"]) == 80000]
    return (float(rows[0]["dE_over_dC"]), float(rows[0]["delta_E"])) if rows else (None, None)


def main():
    jobs = [(ns, eu, s) for ns, eu, _ in CONDS.values() for s in EU_STEPS]
    got = {}
    with ProcessPoolExecutor(int(os.environ.get("E5_NPROC", "3"))) as ex:
        for ns, s, dE in ex.map(delta_E, jobs):
            got[(ns, s)] = dE
    out = {"item": "E5 δ_E・C2 の参照場の感度 (plan verification-case45-euler-total-enthalpy §6 E5)", "tol_width": TOL,
           "ns_res": NS_RES, "euler_steps": EU_STEPS, "solve_json": SOLVE, "conds": {}}
    ok_all = True
    for c, (ns, eu, prob) in CONDS.items():
        dC = delta_C(prob)
        r = np.array([got[(ns, s)] / dC for s in EU_STEPS])
        ch_ratio, ch_dE = chain_value(ns)
        last = float(r[-1])
        pre = bool(ch_ratio is not None and last == ch_ratio and got[(ns, EU_STEPS[-1])] == ch_dE)
        width = float(r.max() - r.min())
        ok = bool(pre and width <= TOL)
        ok_all &= ok
        out["conds"][c] = {"ns_run": ns, "euler_run": eu, "problem": prob, "delta_C": dC,
                           "precondition_bit_identical_to_chain": pre, "chain_dE_over_dC_80000": ch_ratio, "this_54000": last,
                           "dE_over_dC": {str(s): float(v) for s, v in zip(EU_STEPS, r)},
                           "width": width, "mean_minus_54000": float(r.mean() - last), "pass": ok,
                           "verdict": ("PASS" if ok else ("判定しない (経路の同一性の前提が不成立)" if not pre else "FAIL (幅 > 5e-4)"))}
        print(f"[E5] {c}: 前提 {pre}  幅 {width:.3e}  平均 − 54000 {float(r.mean() - last):+.3e}  → {out['conds'][c]['verdict']}", flush=True)
    out["verdict"] = "PASS" if ok_all else "NOT PASS (条件ごとの verdict を見る)"
    (C / "_band_ab" / "euler_ref_de_sensitivity.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(out["verdict"])


if __name__ == "__main__":
    main()
