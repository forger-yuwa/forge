"""形状ゲート S (plan tooling-nozzle-throat-monotone-r2 §6 S1〜S8、§5.1 #4)。CFD 0 step。
対象: 最終問題 problem_d155_ns_finemesh_recal_final.yaml (run_0117/0118 の壁) の MOC 点群、λ 1e-9。
  現行壁 = キー無し (`geometry.wall_fit_mono_r2` 無し) / 単調壁 = `wall_fit_mono_r2: [0.0, 1.5]`。どちらも生産の設計チェーン (design_chain) で作る。
極値・増加量・∫(r‴)²・max|r⁗| は区間多項式で厳密に評価する (`r3_piecewise_exact`; 均等点ではない)。
S6 の物理壁は prepare_ns と同じ経路 (runner_axismach.integral_delta_r: integral_bl → smooth_delta_quintic → delta_r_from_table、
PhysicalNozzleWall、pw_ramp [−11, −6])。積分法の設定は run_0117 と同じ (contur, a_crocco 1, cf_scale = c2pin_solve_recal.json の k_f, n_scale 1)。
各項目に PASS / FAIL / RECORD と数値を出し、最後に VERDICT 行。閾値は plan §6 S1〜S8 の事前登録値 (ここで変えない)。
usage: [CASE_RUNS=<run_0062 / run_0117 のある case dir>] python3 throat_mono_shape_gate.py [OUT_JSON]
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate.runner_axismach import (_gam_or_gas, design_chain, integral_delta_r,  # noqa: E402
                                                   load_problem)
from forge_design.geometry.wall_axismach import PhysicalNozzleWall, joint_fit_wall, r3_piecewise_exact  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/throat_mono_shape_gate.json"
RUNS = Path(os.environ.get("CASE_RUNS", C))
PROB = C / "problem_d155_ns_finemesh_recal_final.yaml"
MONO = [0.0, 1.5]
BASE_COMMIT = "ebcbc632445d27997603b1e0e8d7437f5d16ed3c"   # 改修前の wall_axismach.py (S7 の基準)


def load(mono):
    p = load_problem(PROB)
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    if mono is not None:
        p.geometry["wall_fit_mono_r2"] = list(mono)
    return p


p0, pm = load(None), load(MONO)
d0, dm = design_chain(p0), design_chain(pm)
W0, Wm = d0["wall"], dm["wall"]
tb = dm["wall_inv"]
R = float(dm["R"])
x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
assert np.array_equal(tb, d0["wall_inv"]), "現行と単調で MOC 点群が違う (当てはめ以外が変わった)"
res = {}


def item(key, status, **vals):
    res[key] = {"status": status, **vals}
    print(f"{key:4s} {status:6s} " + "  ".join(f"{k}={v:.6g}" if isinstance(v, float) else f"{k}={v}" for k, v in vals.items()
                                            if not isinstance(v, (dict, list))))


# --- S1 単調性 ([0, 1.5]、区間多項式で厳密) ------------------------------------------------
m1 = r3_piecewise_exact(Wm._spl, None, *MONO)
m1_cur = r3_piecewise_exact(W0._spl, None, *MONO)
ok1 = m1["r3_max"] <= 1e-6 and m1["r2_max_increase"] <= 1e-7 and m1["r2_max"] <= 1.0 / R + 1e-9
item("S1", "PASS" if ok1 else "FAIL", r3_max=m1["r3_max"], x_r3_max=m1["x_r3_max"], r2_max_increase=m1["r2_max_increase"],
     r2_max=m1["r2_max"], r2_max_minus_1overR=m1["r2_max"] - 1.0 / R,
     r3_coef_max_unnormalized=dm["wall_fit"]["r3_coef_max_unnormalized"],
     limits={"r3_max": 1e-6, "r2_max_increase": 1e-7, "r2_max": "1/R + 1e-9"},
     kkt={k: dm["wall_fit"][k] for k in ("n_constraints", "n_active", "iters", "eq_resid", "ineq_max_normalized", "mu_min",
                                         "active_support")},
     current_wall={"r3_max": m1_cur["r3_max"], "r2_max_increase": m1_cur["r2_max_increase"], "r2_max": m1_cur["r2_max"],
                   "x_r2_max": m1_cur["x_r2_max"]})

# --- S2 接続 (x = 0 で上流 Hermite との跳び) -----------------------------------------------
z = np.array([0.0])
jump = [float(abs(Wm.up.r(z, k)[0] - Wm._spl(z, k)[0])) for k in (0, 1, 2)]
jump_cur = [float(abs(W0.up.r(z, k)[0] - W0._spl(z, k)[0])) for k in (0, 1, 2)]
item("S2", "PASS" if max(jump) <= 1e-8 else "FAIL", jump_r=jump[0], jump_r1=jump[1], jump_r2=jump[2], limit=1e-8,
     current_wall=dict(zip(("jump_r", "jump_r1", "jump_r2"), jump_cur)))

# --- S3 局所性 (x ≥ 0.3 で現行との差) ------------------------------------------------------
xs3 = np.unique(np.r_[np.linspace(0.3, Wm.x_e, 800001), Wm._spl.t[(Wm._spl.t >= 0.3)], W0._spl.t[(W0._spl.t >= 0.3)]])
dr3 = np.abs(Wm.r(xs3) - W0.r(xs3))
dth3 = np.abs(np.degrees(np.arctan(Wm.r(xs3, 1)) - np.arctan(W0.r(xs3, 1))))
item("S3", "PASS" if (dr3.max() <= 1e-7 and dth3.max() <= 1e-4) else "FAIL",
     dr_max_rt=float(dr3.max()), x_dr_max=float(xs3[np.argmax(dr3)]), dtheta_max_deg=float(dth3.max()),
     x_dtheta_max=float(xs3[np.argmax(dth3)]), limits={"dr_rt": 1e-7, "dtheta_deg": 1e-4}, n_eval=int(len(xs3)))

# --- S4 点上の忠実度 ----------------------------------------------------------------------
bound = float(np.degrees(th[1] - np.arctan(x[1] / R)))      # §4.1 の解析的下限 θ₁ − atan(x₁/R)


def fidelity(s):
    dth = np.degrees(np.arctan(s(x, 1)) - th)
    drr = s(x) - r
    return dict(dtheta_first_deg=float(dth[1]), dtheta_max_ge2_deg=float(np.abs(dth[2:]).max()),
                x_dtheta_max_ge2=float(x[2 + int(np.argmax(np.abs(dth[2:])))]),
                dr_max_ge2_rt=float(np.abs(drr[2:]).max()), dr_first_rt=float(drr[1]))


f4, f4c = fidelity(Wm._spl), fidelity(W0._spl)
gap = abs(abs(f4["dtheta_first_deg"]) - bound)
item("S4", "PASS" if gap <= 0.002 else "FAIL", x1=float(x[1]), bound_deg=bound, dtheta_first_deg=f4["dtheta_first_deg"],
     gap_to_bound_deg=gap, limit_deg=0.002,
     record_ge2={"mono": {k: v for k, v in f4.items() if k != "dtheta_first_deg"},
                 "current": f4c})

# --- S5 λ 感度 (S1 が成り立つか = 合否、差は記録) ------------------------------------------
g03 = np.linspace(0.0, 0.3, 30001)
xa = np.linspace(0.0, Wm.x_e, 400001)
s5 = {}
ok5 = True
for lam in (1e-10, 1e-8):
    dg = {}
    s_l, _ = joint_fit_wall(tb, R, lam=lam, mono_r2=MONO, diag=dg)
    ml = r3_piecewise_exact(s_l, None, *MONO)
    okl = ml["r3_max"] <= 1e-6 and ml["r2_max_increase"] <= 1e-7 and ml["r2_max"] <= 1.0 / R + 1e-9
    ok5 = ok5 and okl
    # 差は λ 1e-9 の単調壁に対して
    s5[f"{lam:g}"] = dict(S1_pass=bool(okl), r3_max=ml["r3_max"], r2_max_increase=ml["r2_max_increase"], r2_max=ml["r2_max"],
                          n_active=dg["n_active"], iters=dg["iters"],
                          dr_all_max_rt=float(np.abs(s_l(xa) - Wm._spl(xa)).max()),
                          dtheta_all_max_deg=float(np.abs(np.degrees(np.arctan(s_l(xa, 1)) - np.arctan(Wm._spl(xa, 1)))).max()),
                          d_r2_0_0p3=float(np.abs(s_l(g03, 2) - Wm._spl(g03, 2)).max()),
                          d_r3_0_0p3=float(np.abs(s_l(g03, 3) - Wm._spl(g03, 3)).max()),
                          d_r4_0_0p3=float(np.abs(s_l(g03, 4) - Wm._spl(g03, 4)).max()),
                          dtheta_first_deg=float(np.degrees(np.arctan(s_l(x[1], 1)) - th[1])))
item("S5", "PASS" if ok5 else "FAIL", note="合否は S1 が λ 1e-10・1e-8 で成り立つか。差 (λ 1e-9 の単調壁に対して) は記録", lam=s5)

# --- S6 物理壁 (生産経路) ----------------------------------------------------------------
kf = float(json.loads((C / "c2pin_solve_recal.json").read_text())["k_f"])
init_cfg = {"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}    # run_recal_cond.sh / run_0117 と同じ


def physical(p, d):
    res_init, drx, info = integral_delta_r(p, d, init_cfg)
    pw = p.geometry.get("pw_ramp")
    PW = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]),
                            _gam_or_gas(p), p.cp, dstar_x=None, offset="radial", delta_r_x=drx,
                            ramp=(None if pw is None else tuple(float(v) for v in pw)))
    return PW, drx, res_init, info


PWm, drm, rim, infm = physical(pm, dm)
PW0, dr0, ri0, inf0 = physical(p0, d0)
assert PWm.analytic and PW0.analytic, "物理壁が解析経路でない"


def phys_metrics(PW, W, drx, a, b):
    brk = np.r_[W._spl.t, drx.spline.t]
    m = r3_piecewise_exact(lambda xq, k: PW.r(xq, k), brk, a, b)
    mdr = r3_piecewise_exact(lambda xq, k: drx(xq, k), drx.spline.t, a, b)
    return m, mdr


m6, m6dr = phys_metrics(PWm, Wm, drm, *MONO)
m6c, m6cdr = phys_metrics(PW0, W0, dr0, *MONO)
e6, _ = phys_metrics(PWm, Wm, drm, 0.0, 3.0)
e6c, _ = phys_metrics(PW0, W0, dr0, 0.0, 3.0)
vm, v0 = PWm.validate(), PW0.validate()
xg6 = np.linspace(PWm.x_in, PWm.x_e, 800001)
d6 = np.abs(PWm.r(xg6) - PW0.r(xg6))
# 照合: 現行物理壁の δ_r が run_0117 の delta_r_initial.csv (AWS で prepare_ns が書いた値) と一致するか (経路の同一性の確認)
chk117 = None
f117 = RUNS / "run_0117_ns_recal_final_ext/delta_r_initial.csv"
if f117.exists():
    t117 = np.loadtxt(f117, delimiter=",", skiprows=1)
    chk117 = dict(n=int(len(t117)), same_x=bool(len(t117) == len(ri0["x"]) and np.allclose(t117[:, 0], ri0["x"], rtol=0, atol=1e-12)),
                  max_abs_delta_r_diff_rt=(float(np.abs(t117[:, 1] - ri0["delta_r"]).max()) if len(t117) == len(ri0["x"]) else None))
ok6 = m6["r2_max_increase"] <= 0.002 and vm == []
item("S6", "PASS" if ok6 else "FAIL", r2_max_increase_phys=m6["r2_max_increase"], limit=0.002,
     r2_max_phys=m6["r2_max"], x_r2_max_phys=m6["x_r2_max"], validate=vm,
     current_physical={"r2_max_increase": m6c["r2_max_increase"], "r2_max": m6c["r2_max"], "x_r2_max": m6c["x_r2_max"],
                       "validate": v0},
     diag_delta_r2_increase={"mono": m6dr["r2_max_increase"], "current": m6cdr["r2_max_increase"]},
     shape_diff_vs_current_physical={"dr_max_rt": float(d6.max()), "x_dr_max": float(xg6[np.argmax(d6)])},
     r2_interior_extrema_0_3={"mono_physical": e6["n_r2_extrema"], "current_physical": e6c["n_r2_extrema"]},
     throat={"mono": {"x_throat": PWm.x_throat, "r_throat": PWm.r_throat}, "current": {"x_throat": PW0.x_throat, "r_throat": PW0.r_throat}},
     delta_r_exit={"mono": infm["delta_r_exit"], "current": inf0["delta_r_exit"]},
     init_cfg=init_cfg, ramp_gate={"mono": PWm.ramp_gate, "current": PW0.ramp_gate}, check_vs_run_0117_delta_r=chk117)

# --- S7 既定のビット同一 ----------------------------------------------------------------
src = subprocess.run(["git", "-C", str(ROOT), "show", f"{BASE_COMMIT}:design/forge_design/geometry/wall_axismach.py"],
                     capture_output=True, text=True, check=True).stdout
ns = {"__name__": "wall_axismach_base"}
exec(compile(src.replace("from .wall_walldriven", "from forge_design.geometry.wall_walldriven"), "wall_axismach_base", "exec"), ns)
s_base, _ = ns["joint_fit_wall"](tb, R)
xs7 = np.linspace(0.0, W0.x_e, 800001)
dr7 = float(np.abs(W0._spl(xs7) - s_base(xs7)).max())
same_c = bool(np.array_equal(W0._spl.c, s_base.c) and np.array_equal(W0._spl.t, s_base.t))
item("S7", "PASS" if (dr7 == 0.0 and same_c) else "FAIL", dr_max=dr7, coef_knots_identical=same_c, base_commit=BASE_COMMIT)

# --- S8 滑らかさ ([0, 0.3]、厳密) ---------------------------------------------------------
m8, m8c = r3_piecewise_exact(Wm._spl, None, 0.0, 0.3), r3_piecewise_exact(W0._spl, None, 0.0, 0.3)
ok8 = m8["int_r3sq"] <= m8c["int_r3sq"] and m8["r4_absmax"] <= m8c["r4_absmax"]
item("S8", "PASS" if ok8 else "FAIL", int_r3sq_mono=m8["int_r3sq"], int_r3sq_current=m8c["int_r3sq"],
     r4_absmax_mono=m8["r4_absmax"], r4_absmax_current=m8c["r4_absmax"])

gated = [k for k in res if res[k]["status"] in ("PASS", "FAIL")]
verdict = "PASS" if all(res[k]["status"] == "PASS" for k in gated) else "FAIL"
fails = [k for k in gated if res[k]["status"] == "FAIL"]
out = {"plan": "plans/accepted/tooling-nozzle-throat-monotone-r2.md §6 S1〜S8", "problem": PROB.name, "mono_r2": MONO, "lam": 1e-9,
       "R": R, "r_t_mm": float(pm.spec["r_throat"]) * 1e3, "case_runs": str(RUNS), "items": res,
       "wall_fit_mono": {k: v for k, v in dm["wall_fit"].items() if k != "spline"},
       "VERDICT": verdict, "fail_items": fails}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
print(f"VERDICT: {verdict}" + (f" (FAIL: {', '.join(fails)})" if fails else "") + f"  -> {OUT}")
