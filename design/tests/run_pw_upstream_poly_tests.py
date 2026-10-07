#!/usr/bin/env python3
"""物理壁の上流の多項式 (`geometry.pw_upstream: poly`) と寸法の逆算 (`solve_rt_throat`・CFD 前の `solve_rt`) の試験
(plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1・§4.2・§5 の 5)。

壁は case/45 の単調壁の生産問題 (`problem_d155_ns_finemesh_recal_final_mono.yaml`) の初期線だけ Hall に差し替えて作る (CFD ピンの
凍結源 run に依存しない)。δ_r は prepare_ns と同じ積分法の経路 (`integral_delta_r`、YAML の deltastar_initializer)。

1. キー: `_pw_upstream` (既定 poly・明示・不正値・poly と pw_ramp の併記・joint でない壁・ramp + single_bspline)、`_sizing_spec`、
   prepare_ns・Euler の prepare が run dir を作る前に止まる
2. poly の壁: [0, x_e] が ramp とビット同一 / Q の端条件と継ぎ目 / ゲート (|Q″ − H″| は密な点と一致) / 物理スロート (大域最小) が
   ramp の壁と 1e-6 r_t・1e-9 r_t 以内 / Q′ の実根は −L_U (重根) と物理スロートだけ / 既定ランプとの差の診断
3. 大域最小の探索器の負例: δ_r′(0) が正・0・負 (レビューの反例 x_t = +0.00155091)・下流に余分な極値・合成関数 (二つの最小・平らな底・
   区切り不足) — 位置は密な点の総当たりと比べる
4. 寸法の逆算: CFD 前の δ_r の供給が prepare_ns の経路と同一 / solve_rt (CFD 前) の自明な目標 / solve_rt_throat (NS 後、合成の
   δ_E の全分布) の往復 (書いた表から作り直した壁の最小半径と目標の差 ≤ 1e-9 m) / 反復の上限で不合格
5. 一般性: 標準の L_U 3.5・r_U 2.5 の縮流部 (ゲートを通るかを記録; 自動で ramp に切り替えない)

usage: design/.venv-opt/bin/python design/tests/run_pw_upstream_poly_tests.py   (所要 1〜2 分)
"""
import copy
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate.runner_axismach import (_gam_or_gas, _pw_upstream, _sizing_spec, build_physical_wall,  # noqa: E402
                                                   delta_r_from_table, design_chain, integral_delta_r, load_problem)
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.geometry.wall_axismach import (PhysicalNozzleWall, SingleBSplinePhysicalWall,  # noqa: E402
                                                 wall_global_min)

ROOT = Path(__file__).resolve().parents[2]
PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml"
FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc as e:  # noqa: PERF203
        return str(e)
    return None


if not PROB.exists():
    print(f"問題 YAML {PROB} が無い — 全体を飛ばす (合格扱いにしない)")
    sys.exit(2)

TMP = Path(tempfile.mkdtemp(prefix="pwpoly_"))


def hall_yaml(name: str, mut=None) -> Path:
    """生産問題の写し (初期線を Hall に、pw_ramp・pw_upstream を外す)。mut(g, spec) で geometry・spec を書き換える。"""
    import yaml
    raw = yaml.safe_load(PROB.read_text())
    g = raw["geometry"]
    g["initial_line"] = "hall"
    for k in ("initial_line_run", "initial_line_res", "pw_ramp", "pw_upstream"):
        g.pop(k, None)
    if mut is not None:
        mut(g, raw["spec"])
    out = TMP / f"{name}.yaml"
    out.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return out


# --- 1. キー -------------------------------------------------------------------------------------------------------
J = {"wall_repr": "joint"}
check("joint・キー無し → poly (既定)", _pw_upstream(J) == {"value": "poly", "source": "default", "requested": None})
check("joint・ramp を明示 → ramp", _pw_upstream({**J, "pw_upstream": "ramp", "pw_ramp": [-11, -6]})["value"] == "ramp")
check("joint・poly を明示 → poly (source explicit)", _pw_upstream({**J, "pw_upstream": "poly"}) == {"value": "poly", "source": "explicit", "requested": "poly"})
check("不正値 (null・大文字・空白・数値・真偽値・未知) は例外",
      all(raises(lambda v=v: _pw_upstream({**J, "pw_upstream": v})) is not None for v in (None, "Poly", " poly", "ramp ", 1, True, "hermite")))
for g, lab in (({**J, "pw_ramp": [-11, -6]}, "既定 poly と pw_ramp"), ({**J, "pw_upstream": "poly", "pw_ramp": [-11, -6]}, "明示 poly と pw_ramp"),
               ({**J, "pw_ramp": None}, "既定 poly と pw_ramp: null")):
    e = raises(lambda g=g: _pw_upstream(g))
    check(f"{lab} の併記は例外 ({(e or '')[:50]}…)", e is not None and "併記" in e)
e = raises(lambda: _pw_upstream({"wall_repr": "interp", "pw_upstream": "poly"}))
check(f"joint でない壁に poly を明示 → 例外 ({(e or '')[:40]}…)", e is not None)
check("joint でない壁・キー無し → 今の振る舞い (value None)", _pw_upstream({"wall_repr": "interp"})["value"] is None
      and _pw_upstream({})["value"] is None)
r_ = _pw_upstream({"wall_repr": "interp", "pw_upstream": "ramp"})
check(f"joint でない壁に ramp を明示 → 無効として記録 ({r_})", r_["value"] is None and "無効" in r_["source"])
e = raises(lambda: _pw_upstream({**J, "pw_upstream": "ramp", "pw_ramp": [-11, -6], "physical_wall_repr": "single_bspline"}))
check(f"ramp + single_bspline → 例外 ({(e or '')[:40]}…)", e is not None)
check("poly + single_bspline は受ける", _pw_upstream({**J, "physical_wall_repr": "single_bspline"})["value"] == "poly")
check("spec.sizing: 無し → None / exit・throat を受ける", _sizing_spec({}) is None
      and _sizing_spec({"sizing": {"method": "throat", "target_m": 0.0767}})["method"] == "throat")
check("spec.sizing の不正値 (method 不明・target_m が負・真偽値・文字列) は例外",
      all(raises(lambda v=v: _sizing_spec({"sizing": v})) is not None for v in
          ({"method": "area", "target_m": 1.0}, {"method": "exit", "target_m": -1.0}, {"method": "exit", "target_m": True},
           {"method": "exit", "target_m": "0.7"}, "exit", None)))
# prepare_ns / prepare は run dir を作る前に止まる
y_bad = hall_yaml("bad_poly_ramp", lambda g, s: g.__setitem__("pw_ramp", [-11.0, -6.0]))
rd = TMP / "never_ns"
e = raises(lambda: RA.prepare_ns(y_bad, rd))
check(f"prepare_ns: poly (既定) と pw_ramp の併記で run dir を作る前に例外 ({(e or '')[:30]}…)", e is not None and not rd.exists())
y_bad2 = hall_yaml("bad_value", lambda g, s: g.__setitem__("pw_upstream", "Poly"))
e = raises(lambda: RA.prepare_ns(y_bad2, rd))
check("prepare_ns: 不正値で run dir を作る前に例外", e is not None and not rd.exists())
y_bad3 = hall_yaml("bad_sizing", lambda g, s: s.__setitem__("sizing", {"method": "exit"}))
e = raises(lambda: RA.prepare_ns(y_bad3, rd))
check("prepare_ns: spec.sizing の不正値で run dir を作る前に例外", e is not None and not rd.exists())
y_poly = hall_yaml("euler_poly", lambda g, s: g.__setitem__("pw_upstream", "poly"))
e = raises(lambda: RA.prepare(y_poly, TMP / "never_euler"))
check(f"Euler の prepare: poly の明示は例外 (物理壁が無い) ({(e or '')[:30]}…)", e is not None and not (TMP / "never_euler").exists())

# --- 2. poly の壁 (case/45 の形) ---------------------------------------------------------------------------------
p = load_problem(hall_yaml("c45"))
d = design_chain(p)
res_init, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
w = d["wall"]
args = (w, d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
PP = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx)
PR = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, ramp=(-11.0, -6.0), upstream="ramp")
check("build_physical_wall (prepare_ns と同じ構築) は poly の壁を作る",
      build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=drx, offset="radial").pw_upstream == "poly")
xd = np.linspace(0.0, PP.x_e, 300001)
check("[0, x_e] が ramp の壁とビット同一 (r, r′, r″, r‴)", all(np.array_equal(PP.r(xd, n), PR.r(xd, n)) for n in range(4)))
L = float(w.up.L_U)
e0 = [float(w.r(np.r_[0.0], n)[0] + drx(np.r_[0.0], n)[0]) for n in range(3)]
check(f"Q の端条件: x = 0 で S + δ_r の (r, r′, r″) = {[round(v, 7) for v in e0]}", PP.upstream_poly["end_conditions"]["x_0"] == e0)
g = PP.upstream_gate
check(f"継ぎ目 (−L_U・0) の値・1 階・2 階の跳び ≤ 1e-8: {g['max_seam_jump']:.1e}", g["max_seam_jump"] <= 1e-8)
xs = np.r_[np.linspace(-L, 0.0, 200001)[:-1], -1e-13]                # x → 0⁻ (左の極限) を含める
d2_dense = float(np.max(np.abs(PP.r(xs, 2) - w.r(xs, 2))))
check(f"|Q″ − H″| の最大 (区間端と極値で厳密) {g['max_abs_d2_change_vs_H']:.9e} ≥ 密な点 {d2_dense:.9e} かつ差 ≤ 1e-9",
      g["max_abs_d2_change_vs_H"] >= d2_dense - 1e-15 and g["max_abs_d2_change_vs_H"] - d2_dense <= 1e-9)
check(f"ゲート合格 (正の半径 {g['positive_radius']}, 一意 {g['unique_min']}, 前 r′ ≤ 0 {g['monotone_before']}, 後 r′ ≥ 0 {g['monotone_after']}, "
      f"|Q″ − H″| {g['max_abs_d2_change_vs_H']:.3e} ≤ 5e-3)", g["pass"] and g["max_abs_d2_change_vs_H"] <= 5e-3)
dx, dr = PP.x_throat - PR.x_throat, PP.r_throat - PR.r_throat
check(f"物理スロート (大域最小) の ramp との差: 位置 {dx:.2e} ≤ 1e-6、半径 {dr:.2e} ≤ 1e-9 (x_t = {PP.x_throat:.7f})",
      abs(dx) <= 1e-6 and abs(dr) <= 1e-9)
check(f"物理スロートで r′ = 0 (|r′| {abs(float(PP.r(np.r_[PP.x_throat], 1)[0])):.1e} ≤ 1e-12)", abs(float(PP.r(np.r_[PP.x_throat], 1)[0])) <= 1e-12)
P = np.polynomial.Polynomial
q1 = P(PP._q_c).deriv()                                     # ξ 微分 (ξ = (x + L)/L)
rr = q1.roots()
rr = np.sort(rr[np.abs(rr.imag) <= 1e-6].real)
xr = -L + L * rr[(rr >= -1e-6) & (rr <= 1 + 1e-6)]
check(f"Q′ の定義域内の実根は −L_U (重根) と物理スロートだけ: {np.round(xr, 7).tolist()}",
      len(xr) >= 2 and np.all(np.abs(xr[:-1] + L) <= 1e-4) and abs(xr[-1] - PP.x_throat) <= 1e-9)
vr = g["vs_ramp_default"]
check(f"診断: 既定ランプ {vr['ramp_ref']} との差の最大 {vr['max_abs_diff_rt']:.3e} r_t (x = {vr['x_at_max_abs_diff']:.3f})",
      vr["max_abs_diff_rt"] > 0 and vr["ramp_ref_source"] == "default_pw_ramp")
check(f"validate() が空 ({PP.validate()})", PP.validate() == [])

# --- 3. 大域最小の探索器の負例 ---------------------------------------------------------------------------------------
xt_tab = res_init["x"]


def wall_with(fun):
    f = delta_r_from_table(xt_tab, fun(xt_tab))
    return f, (lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=f))


def brute_min(W, lo=-1.0, hi=1.0):
    """密な点の総当たり + 細かい格子での詰め (探索器と独立)。"""
    xx = np.linspace(lo, hi, 400001)
    i = int(np.argmin(W.r(xx)))
    for _ in range(4):
        h = xx[1] - xx[0]
        xx = np.linspace(xx[max(i - 2, 0)], xx[min(i + 2, len(xx) - 1)], 4001)
        i = int(np.argmin(W.r(xx)))
    return float(xx[i])


# 下流の伸び (x > 10 で 1e-6·(x − 10)³、x ≤ 10 で厳密に 0): 一定の δ_r だと出口の手前で設計壁の r′ (−2e-8) が出てゲートが別の理由で
# 止まるので、スロートから遠い所にだけ単調な伸びを足す (スロート近傍の δ_r は plan 段レビューの式のまま)
grow = lambda x: 1e-6 * np.maximum(x - 10.0, 0.0) ** 3  # noqa: E731
cases = {"δ_r′(0) > 0": lambda x: 0.0014 + 0.000776 * x * np.exp(-(x / 0.1) ** 2) + grow(x),
         "δ_r′(0) = 0": lambda x: 0.0014 + 0.0 * x + grow(x),
         "δ_r′(0) < 0 (レビューの反例)": lambda x: 0.0014 - 0.000776 * x * np.exp(-(x / 0.1) ** 2) + grow(x)}
xts = {}
for lab, fun in cases.items():
    f, mk = wall_with(fun)
    W = mk()
    xb = brute_min(W)
    xts[lab] = W.x_throat
    check(f"{lab}: δ_r′(0) = {float(f(np.r_[0.0], 1)[0]):+.2e}, x_t = {W.x_throat:+.8f} (総当たり {xb:+.8f}), 一意 {W.upstream_gate['unique_min']}, "
          f"ゲート {W.upstream_gate['pass']}", abs(W.x_throat - xb) <= 2e-6 and W.upstream_gate["pass"])
check("位置の符号: δ_r′(0) > 0 → x_t < 0 (Q の中)、= 0 → |x_t| ≤ 1e-9、< 0 → x_t > 0 (x = 0 より下流)",
      xts["δ_r′(0) > 0"] < 0 and abs(xts["δ_r′(0) = 0"]) <= 1e-9 and xts["δ_r′(0) < 0 (レビューの反例)"] > 0)
# レビューの値 +0.00155091 は CFD ピンの設計壁と別の x の表で求めたもの。ここ (Hall の設計壁・積分法の 1500 点の表) では幅 0.1 の
# ガウスが表の間隔 0.072 で鈍り δ_r′(0) = −7.15e-4 になるので位置は違う。合否は総当たりとの一致と符号で見る (U1 は生産の壁で記録)
print(f"info レビューの反例 (この壁・表): x_t = {xts['δ_r′(0) < 0 (レビューの反例)']:.8f} (レビューの値 +0.00155091 は別の壁・表)")
# 丸めの並び: δ_r′(0) が丸め程度に正 (根が x = 0 のすぐ上流 −1e-11 級) だと、根と継ぎ目 0 の r が丸めで並ぶ。全候補の最小で選ぶと
# 継ぎ目を拾い、根との小区間で r′ > 0 になって単調性を誤判定した (2026-10-07 U1 の負例で発見、極小の中から選ぶように直した)
a1 = float(drx(np.r_[0.0], 1)[0])
f_tie, mk_tie = wall_with(lambda x: drx(x) - a1 * x * np.exp(-x ** 2))
Wt = mk_tie()
check(f"丸めの並び (δ_r′(0) = {float(f_tie(np.r_[0.0], 1)[0]):.1e}): x_t = {Wt.x_throat:.2e} (≈ −δ_r′(0)/r″)、ゲート合格",
      Wt.upstream_gate["pass"] and abs(Wt.x_throat) <= 1e-9 and Wt._throat_diag["search"]["x_t_is_local_min"])
fB, mkB = wall_with(cases["δ_r′(0) < 0 (レビューの反例)"])
BB = SingleBSplinePhysicalWall(mkB())
check(f"1 本の B-spline も同じ探索器 (反例で x_t = {BB.x_throat:.8f}、元の壁との差 {BB.x_throat - xts['δ_r′(0) < 0 (レビューの反例)']:.1e})",
      abs(BB.x_throat - xts["δ_r′(0) < 0 (レビューの反例)"]) <= 1e-9 and BB._throat_diag["method"] == "single_bspline_global_min")
# 下流に余分な極値 (x ≈ 3 の δ_r の窪みで物理壁が局所的に縮む)
f_dip, mk_dip = wall_with(lambda x: 0.0014 + 0.000776 * x * np.exp(-(x / 0.1) ** 2) - 0.4 * np.exp(-((x - 3.0) / 0.3) ** 2) + grow(x))
e = raises(mk_dip)
check(f"下流に余分な極値 → ゲート不合格で例外 ({(e or '')[:30]}…「後で r′ ≥ 0 False」)", e is not None and "後で r′ ≥ 0 False" in e)
Wd = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f_dip, ramp=(-11.0, -6.0), upstream="ramp")
br = np.unique(np.r_[0.0, w._spl.t, f_dip.spline.t, f_dip.x_range, Wd.x_e]); br = br[(br >= 0) & (br <= Wd.x_e)]
sd = wall_global_min(Wd.r, br)
check(f"探索器: 余分な極値を単調性の違反として報告 (後の r′ の最小 {sd['min_r1_after']:.2e} at x = {sd['x_min_r1_after']:.3f})",
      (not sd["monotone_after"]) and 2.0 < sd["x_min_r1_after"] < 3.2 and not sd["pass"])
# 合成関数: 二つの同じ最小 / 平らな底 / 一つの最小 / 区切り不足
two = lambda x, n: {0: (x * x - 1.0) ** 2, 1: 4 * x * (x * x - 1.0), 2: 12 * x * x - 4.0}[n]  # noqa: E731
s2 = wall_global_min(two, np.linspace(-2.0, 2.0, 9))
check(f"合成 (x²−1)²: 二つの同じ最小 → 一意でない ({s2['n_at_min']} か所)", (not s2["unique"]) and s2["n_at_min"] == 2 and not s2["pass"])


def flat(x, n):
    x = np.asarray(x, dtype=float)
    a = np.abs(x) - 1.0
    v = np.where(a > 0, a ** 2, 0.0) if n == 0 else (np.where(a > 0, 2 * a * np.sign(x), 0.0) if n == 1 else np.where(a > 0, 2.0, 0.0))
    return v
sf = wall_global_min(flat, [-2.0, -1.0, 1.0, 2.0])
check(f"合成: 平らな底 (r′ ≡ 0 の区間が最小) → 一意でない ({sf['n_at_min']} か所)", (not sf["unique"]) and sf["n_flat_intervals"] == 1)
one = lambda x, n: {0: 1.0 + (x - 0.3) ** 2, 1: 2 * (x - 0.3), 2: 2.0 + 0 * x}[n]  # noqa: E731
s1 = wall_global_min(one, [-1.0, 0.0, 1.0])
check(f"合成: 一つの最小 → x_t = {s1['x_t']} (0.3)、一意・単調・合格", s1["pass"] and abs(s1["x_t"] - 0.3) <= 1e-14)
check("区切り不足 (区間内で多項式でない) → 例外", raises(lambda: wall_global_min(lambda x, n: np.cos(3 * np.asarray(x)) if n == 0
                                                       else -3 * np.sin(3 * np.asarray(x)) if n == 1 else -9 * np.cos(3 * np.asarray(x)),
                                                       [-1.0, 1.0])) is not None)

# --- 4. 寸法の逆算 --------------------------------------------------------------------------------------------------
yp = hall_yaml("sizing")
init = DL.resolve_integral_initializer(yp)
pS = load_problem(yp); dS = design_chain(pS)
S0 = float(pS.spec["r_throat"])
dA, _ = DL._sizing_delta_supplier(pS, init)(dS, S0)
_, dB, _ = integral_delta_r(pS, dS, pS.raw["deltastar_initializer"])
xq = np.linspace(-12.4, 95.0, 20001)
check(f"CFD 前の δ_r の供給 = prepare_ns の経路 (integral_delta_r、cf_scale {init.get('cf_scale')}) とビット同一",
      all(np.array_equal(dA(xq, n), dB(xq, n)) for n in range(3)))
W0 = build_physical_wall(pS, dS, S0, delta_r_x=dB, offset="radial")
R_exit0 = S0 * float(W0.r(np.r_[W0.x_e])[0])
rs = DL.solve_rt(yp, R_exit0)
check(f"solve_rt (CFD 前): 今の壁の出口半径 {R_exit0:.9f} m を目標にすると r_t = S0、残差 {rs['residual_m']:.1e}、δ の出典を記録",
      rs["r_t_m"] == S0 and abs(rs["residual_m"]) <= 1e-9 and "integral_delta_r" in rs["source"] and rs["pw_upstream"] == "poly"
      and rs["delta_r_source"]["settings"]["cf_scale"] == init.get("cf_scale"))
# NS 後: 合成の prev_run (δ_E の全分布、S_prev = 0.0768 m)
prev = TMP / "run_prev"
prev.mkdir()
S_prev = 0.0768
(prev / "prepare_info.json").write_text(json.dumps({"scale_m": S_prev}))
dE = dB(xt_tab) * 1.02 + 2e-5
np.savetxt(prev / "delta_r_next.csv", np.c_[xt_tab, dE, dE, dE, 0 * dE], delimiter=",", comments="",
           header=",".join(DL.DELTA_R_NEXT_COLUMNS) + " (test)")
R_t = 0.0768
out_csv = TMP / "delta_next_wall.csv"
rt_ = DL.solve_rt_throat(yp, R_t, prev_run=prev, delta_r_out=out_csv)
check(f"solve_rt_throat (NS 後): 収束 (r_t = {rt_['r_t_m']:.10f} m、残差 {rt_['residual_m']:.1e} m、{rt_['n_iter']} 回)、出口半径を記録",
      abs(rt_["residual_m"]) <= 1e-9 and rt_["exit_radius_m"] > 0 and rt_["sizing"] == "throat" and out_csv.exists())
# 往復: 書いた表を prepare_ns と同じ読み方で読み、解いた r_t の問題で壁を作り直す
yrt = hall_yaml("sizing_rt", lambda g, s: s.__setitem__("r_throat", rt_["r_t_m"]))
pR = load_problem(yrt); dR = design_chain(pR)
tbl = np.loadtxt(out_csv, delimiter=",", skiprows=1)
WR = build_physical_wall(pR, dR, float(pR.spec["r_throat"]), delta_r_x=delta_r_from_table(tbl[:, 0], tbl[:, 1]), offset="radial")
diff = float(pR.spec["r_throat"]) * WR.r_throat - R_t
check(f"solve_rt_throat (NS 後) の往復: 作り直した壁の最小半径と目標の差 {diff:.1e} m ≤ 1e-9 m", abs(diff) <= 1e-9)
check("solve_rt_throat (NS 後): 表は δ_E × (r_t/S_prev)^−0.2", np.allclose(tbl[:, 1], dE * (rt_["r_t_m"] / S_prev) ** -0.2, rtol=0, atol=1e-18))
e = raises(lambda: DL.solve_rt_throat(yp, 0.9 * R_t, prev_run=prev, delta_r_out=out_csv, max_iter=1), DL.SizingNotConverged)
check(f"反復の上限 (1 回) で不合格 (SizingNotConverged: {(e or '')[:40]}…)", e is not None)
check("NS 後で delta_r_out が無ければ例外 (寸法と壁を同じ関数で作るため)", raises(lambda: DL.solve_rt_throat(yp, R_t, prev_run=prev)) is not None)

# --- 5. 一般性: 標準の L_U 3.5・r_U 2.5 ------------------------------------------------------------------------------
p35 = load_problem(hall_yaml("lu35", lambda g, s: (g.__setitem__("L_U", 3.5), g.__setitem__("r_inlet", 2.5))))
d35 = design_chain(p35)
_, dr35, _ = integral_delta_r(p35, d35, p35.raw["deltastar_initializer"])
try:
    W35 = build_physical_wall(p35, d35, float(p35.spec["r_throat"]), delta_r_x=dr35, offset="radial")
    g35 = W35.upstream_gate
    print(f"info 標準 L_U 3.5・r_U 2.5: poly のゲート合格 (|Q″ − H″| {g35['max_abs_d2_change_vs_H']:.3e}, 物理スロート x {W35.x_throat:.6f})")
    check("標準 L_U 3.5・r_U 2.5: poly の壁 (ramp に切り替えていない)", W35.pw_upstream == "poly")
except ValueError as e:
    print(f"info 標準 L_U 3.5・r_U 2.5: poly のゲートで止まる ({str(e)[:160]}…)")
    check("標準 L_U 3.5・r_U 2.5: 止まるなら poly のゲートの例外 (自動で ramp に切り替えない)", "pw_upstream poly" in str(e))

print(f"FAIL 件数: {FAIL}")
sys.exit(1 if FAIL else 0)
