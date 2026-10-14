"""積分法 CONTUR の contur_v2 (plan tooling-nozzle-isothermal-wall-chain §4.7・§5.1 #20) の試験。
  (a) 既定 (contur_v1) は変更前のモジュールとビット同一 (生産の問題、断熱と 300 K、生産の k_f)
  (b1) CPG 極限の式の一致: 同じ点 (T_e, M, T_w) で、エンタルピー形の積分量 (θ, δ*, θ_c, F_c) と断熱壁温が温度形と相対 1e-12 で一致
  (b2) CPG 極限の積分の一致: 同じ粘性 (今の Sutherland) の contur_v2 と contur_v1 の δ_r の差が、v1 自身の縁の格子の誤差
       (縁の格子 4000 と 64000 の差) 以下 (x ≥ 0.5 と試験部)。2026-10-08: 初版は「相対 1e-4」で不合格 (2.1e-4 / 5.2e-4) — 原因は
       縁の状態の補間と数値微分の離散化 (v1 自身の格子の誤差は 1.0〜2.3e-3、積分の許容差の影響は 2〜3e-4) で、式の差ではない
       ((b1) で一致)。1e-4 は積分器の誤差より小さく、試験の設計の誤り (plan §5.1 #20)
  (c) 生産の燃焼ガス: contur_v2 の断熱壁温が NS の断熱壁温 (run_0179 の窓の平均、ns_wall_T.csv) と試験部で ±0.5 %
  (d) 表の範囲外の温度 (T_w = 50 K) は例外 (端値で外挿しない)
  (e) 混合気の粘性が独立参照の値と一致 (transport_reference.py で 2026-10-08 に計算した値、相対 1e-9)
usage: design/.venv-opt/bin/python design/tests/run_contur_v2_tests.py [変更前の deltastar_integral.py]
       (変更前のモジュールが無ければ (a) は飛ばす)。生産の問題 (case/45 の _band_ab/delta_contur/prod_local.yaml) が無ければ (a)〜(d) を飛ばす。
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
CASE = ROOT / "case/45.isobutane_m6_d155"
PROB = CASE / "_band_ab/delta_contur/prod_local.yaml"
fails = []


def check(name, ok, detail=""):
    print(f"[{'OK ' if ok else 'NG '}] {name}  {detail}")
    if not ok:
        fails.append(name)


# (e) 混合気の粘性 (独立参照の値)
from forge_design.gas.transport import MixtureViscosity  # noqa: E402
REF = {200: 1.254353537e-05, 300: 1.780397924e-05, 700: 3.395440309e-05, 1000: 4.361272938e-05, 1470: 5.663785968e-05, 2000: 6.980256548e-05}
mv = MixtureViscosity({"CO2": 0.1677, "H2O": 0.0858, "O2": 0.0220, "N2": 0.7245},
                      {"CO2": "cea", "O2": "cea", "N2": "cea", "H2O": "custom:h2o_iapws_cea_v1"})
worst = max(abs(mv(T) / v - 1.0) for T, v in REF.items())
check("(e) 混合気の μ が独立参照と一致", worst < 1e-9, f"最大 {worst:.1e}")

if not PROB.is_file():
    print(f"生産の問題 {PROB} が無い — (a)〜(d) を飛ばす")
else:
    from forge_design.evaluate.runner_axismach import design_chain, load_problem, _gam_or_gas, contur_mu_fn
    from forge_design.feedback import deltastar_integral as DI
    from forge_design.feedback.deltastar import _sutherland
    p = load_problem(PROB); d = design_chain(p); rt = float(p.spec["r_throat"]); Tt = float(p.spec["Tt"]); Pt = float(p.spec["Pt"])
    kf = float(p.raw["deltastar_initializer"]["cf_scale"])
    walls = {"adiabatic": {"mode": "adiabatic"}, "Tw300": {"mode": "prescribed_temperature", "Tw": 300.0}}
    keys = ("x", "theta", "dstar_n", "H", "N", "delta", "Cf", "Tw", "Taw", "delta_r")
    # (a)
    old_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if old_path and old_path.is_file():
        spec = importlib.util.spec_from_file_location("forge_design.feedback.di_old", old_path)
        old = importlib.util.module_from_spec(spec); sys.modules[spec.name] = old; spec.loader.exec_module(old)
        for w, tbc in walls.items():
            a = old.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, Pt, Tt, rt, thermal_bc=tbc, cf_scale=kf)
            b = DI.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, Pt, Tt, rt, thermal_bc=tbc, cf_scale=kf)
            same = all(np.array_equal(np.asarray(a[k]), np.asarray(b[k])) for k in keys) and a["settings"] == b["settings"]
            check(f"(a) contur_v1 がビット同一 ({w})", same)
    else:
        print("変更前のモジュールが無い — (a) を飛ばす")
    # (b1) CPG 極限の式の一致 (同じ点)
    g = float(p.gamma)
    ec = DI.EdgeConditions(d["wall"], d["wall_inv"], g, p.cp, Pt, Tt, rt, version="contur_v2", mu_fn=_sutherland)
    wi = wt = 0.0; r = 0.72 ** (1 / 3)
    for M in (1.5, 3.0, 4.5, 6.0):
        Te = Tt / (1 + 0.5 * (g - 1) * M ** 2)
        taw1 = Te * (1 + r * 0.5 * (g - 1) * M ** 2)
        he = float(ec.H(Te)); taw2 = float(ec.Hinv(he + r * (float(ec.H(Tt)) - he)))
        wt = max(wt, abs(taw2 / taw1 - 1))
        for Tw in (300.0, 800.0, taw1):
            a = DI._profile_integrals(0.01, 7.0, Tw, taw1, Te, 0.5, 0.95, 1.0)
            b = DI._profile_integrals(0.01, 7.0, Tw, taw1, Te, 0.5, 0.95, 1.0, thermo=ec)
            wi = max(wi, max(abs(bb / aa - 1) for aa, bb in zip(a, b)))
    check("(b1) CPG 極限で式が一致 (同じ点)", wi < 1e-12 and wt < 1e-12, f"積分量 {wi:.1e}、T_aw {wt:.1e}")
    # (b2) CPG 極限の積分: v2 − v1 が v1 自身の縁の格子の誤差以下
    for w, tbc in walls.items():
        a4 = DI.integral_bl(d["wall"], d["wall_inv"], g, p.cp, Pt, Tt, rt, thermal_bc=tbc, cf_scale=kf)
        a64 = DI.integral_bl(d["wall"], d["wall_inv"], g, p.cp, Pt, Tt, rt, thermal_bc=tbc, cf_scale=kf, edge_n_grid=64000)
        b4 = DI.integral_bl(d["wall"], d["wall_inv"], g, p.cp, Pt, Tt, rt, thermal_bc=tbc, cf_scale=kf,
                            closure_version="contur_v2", mu_fn=_sutherland)
        res = []
        for m in (a4["x"] >= 0.5, (a4["x"] >= 40) & (a4["x"] <= 94)):
            dv = float(np.max(np.abs(b4["delta_r"][m] / a4["delta_r"][m] - 1)))
            dg = float(np.max(np.abs(a64["delta_r"][m] / a4["delta_r"][m] - 1)))
            res.append((dv, dg))
        ok = all(dv <= dg for dv, dg in res)
        check(f"(b2) CPG 極限で v2 − v1 ≤ v1 の格子の誤差 ({w})", ok,
              f"x≥0.5: {res[0][0]:.1e} ≤ {res[0][1]:.1e}; 試験部: {res[1][0]:.1e} ≤ {res[1][1]:.1e}")
    # (c) 生産の燃焼ガスの断熱壁温 vs NS
    tw_csv = CASE / "_band_ab/delta_contur/ns_wall_T.csv"
    if tw_csv.is_file():
        r2 = DI.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, Pt, Tt, rt, thermal_bc=walls["adiabatic"], cf_scale=kf,
                            closure_version="contur_v2", mu_fn=contur_mu_fn(p))
        ns = np.loadtxt(tw_csv, delimiter=",", skiprows=1)
        xq = np.linspace(40.0, 94.0, 109)
        rel = np.abs(np.interp(xq, r2["x"], r2["Taw"]) / np.interp(xq, ns[:, 0], ns[:, 1]) - 1.0)
        check("(c) contur_v2 の断熱壁温が NS と ±0.5 % (試験部)", float(rel.max()) <= 0.005, f"最大 {100 * rel.max():.3f} %")
    # (d) 範囲外
    try:
        DI.integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, Pt, Tt, rt, thermal_bc={"mode": "prescribed_temperature", "Tw": 50.0},
                       cf_scale=kf, closure_version="contur_v2", mu_fn=contur_mu_fn(p))
        check("(d) 範囲外の温度は例外", False, "例外が出なかった")
    except ValueError as e:
        check("(d) 範囲外の温度は例外", "範囲" in str(e), str(e)[:80])

print(f"FAIL 件数: {len(fails)}" + (f" ({', '.join(fails)})" if fails else ""))
sys.exit(1 if fails else 0)
