"""単調拘束 (geometry.wall_fit_mono_r2 [0, 1.5]) が、MOC の軸処理を直した後 (analytic + converge) でも要るかの試算 (CFD 0 step)。
2026-10-07 ユーザ「順によろ」の 3 番目 (plan tooling-nozzle-throat-monotone-r2 の拘束は、legacy MOC の始点誤差で出た r″ の山を押さえるために入れた)。

  A: 生産の問題 problem_d155_ns_prod.yaml (拘束あり)
  B: A から wall_fit_mono_r2 だけを外したもの
両方で設計チェーン → 積分法の δ_r → 物理壁を作り、壁の差・r″ の形 (山・単調性の破れ)・形状ゲート (S6: [0, 1.5] で r″ の最大の増加 ≤ 0.002)・
出口半径を比べる。出力 _band_ab/mono_need_probe.json。

usage: design/.venv-opt/bin/python mono_need_probe.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
RUNS = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
PROB = "problem_d155_ns_prod.yaml"


def build(drop_mono: bool):
    from forge_design.evaluate.runner_axismach import build_physical_wall, design_chain, integral_delta_r, load_problem
    txt = (HERE / PROB).read_text().replace("initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k",
                                            f"initial_line_run: {RUNS / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    tmp = HERE / "_band_ab" / f"mono_need_{'B' if drop_mono else 'A'}.yaml"
    tmp.write_text(txt)
    p = load_problem(tmp)
    if drop_mono:
        p.geometry.pop("wall_fit_mono_r2")
    d = design_chain(p)
    _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    PW = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial")
    return p, d, PW


def r2_rise(w, lo=0.0, hi=1.5, n=30001):
    x = np.linspace(lo, hi, n)
    r2 = np.asarray(w.r(x, 2))
    run_min = np.minimum.accumulate(r2)
    return float((r2 - run_min).max()), float(x[np.argmax(r2 - run_min)]), float(r2.max()), float(x[np.argmax(r2)])


def main():
    pa, da, A = build(False)
    pb, db, B = build(True)
    rt = float(pa.spec["r_throat"])
    out = {"r_t_m": rt}
    for lab, lo, hi in (("throat", -0.1, 0.5), ("expansion", 0.5, 3.0), ("downstream", 3.0, float(min(A.x_e, B.x_e)))):
        x = np.linspace(lo, hi, 60001)
        res = {}
        for kind, a, b in (("physical", A, B), ("design", da["wall"], db["wall"])):
            d0 = np.asarray(b.r(x, 0)) - np.asarray(a.r(x, 0))
            th = np.degrees(np.arctan(np.asarray(b.r(x, 1))) - np.arctan(np.asarray(a.r(x, 1))))
            d2 = np.asarray(b.r(x, 2)) - np.asarray(a.r(x, 2))
            res[kind] = {"max_abs_dr_um": float(np.abs(d0).max() * rt * 1e6), "x_at": float(x[np.argmax(np.abs(d0))]),
                         "max_abs_dtheta_deg": float(np.abs(th).max()), "max_abs_dr2": float(np.abs(d2).max())}
        out[f"B_minus_A_{lab}"] = res
    for lab, (p, d, W) in (("A_mono", (pa, da, A)), ("B_free", (pb, db, B))):
        rise_p, xr_p, pk_p, xp_p = r2_rise(W)
        rise_d, xr_d, pk_d, xp_d = r2_rise(d["wall"])
        out[lab] = {"physical_r2_rise_0_1p5": rise_p, "x_rise": xr_p, "physical_r2_max": pk_p, "x_r2_max": xp_p,
                    "design_r2_rise_0_1p5": rise_d, "design_r2_max": pk_d, "design_x_r2_max": xp_d,
                    "S6_pass(<=0.002)": bool(rise_p <= 0.002), "validate": W.validate(),
                    "exit_radius_m": float(W.r(np.r_[W.x_e])[0] * float(p.spec["r_throat"])),
                    "throat": {"x": float(W.x_throat), "r": float(W.r_throat), "kappa": float(W.kappa_throat)},
                    "moc_gate_pass": (d.get("moc") or {}).get("gate", {}).get("pass")}
    out["exit_radius_diff_um"] = (out["B_free"]["exit_radius_m"] - out["A_mono"]["exit_radius_m"]) * 1e6
    (HERE / "_band_ab" / "mono_need_probe.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(out, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
