#!/usr/bin/env python3
"""joint 壁の r″ 単調拘束 (`joint_fit_wall(mono_r2=…)` / `JointFitCFDWall(mono_r2=…)` / `geometry.wall_fit_mono_r2`) の試験。
plan: plans/accepted/tooling-nozzle-throat-monotone-r2.md §5.1 #3 (§4.1・§4.2)。

1. 合成の MOC 風の点群 (run データ不要):
   - mono_r2=None は改修前の式 (下の `_ref_joint_fit_wall` = 改修前の `joint_fit_wall` の写し) と係数・ノットが完全一致
   - mono_r2=(0, 1.5) で [0, 1.5] の r‴ 最大 ≤ 1e-6・r″ の最大増加 ≤ 1e-7 (§6 S1 の形状用検査と同じ関数
     `r3_piecewise_exact`)、等式拘束の残差 ≤ 1e-9 (mono なし/あり)
   - 不正入力 (a ≥ b・x 範囲外・長さ ≠ 2・文字列・非有限・None 要素) は例外
   - `JointFitCFDWall` の `fit_diag` に mono_r2・KKT・spline が入り、spline から壁が復元できる
2. 最終問題 (problem_d155_ns_finemesh_recal_final.yaml、凍結源 run_0062 が要る):
   - キー無しの設計壁が改修前の式と max|Δr| = 0・係数完全一致 (§6 S7)
   - キーあり ([0, 1.5]) で S1 の形状用検査・等式残差 ≤ 1e-9、design_chain の wall_fit に記録が残る
   - wall_repr が joint 以外でキーがあれば例外

usage: python3 design/tests/run_joint_fit_mono_tests.py
  凍結源 run の場所: $CASE_RUNS (run_0062_euler_wallfit_fit_r1_ext6k のある case dir)。既定はリポジトリの case/45、
  無ければ /home/sano/work/forge/case/45.isobutane_m6_d155。どちらにも無ければ 2. を飛ばして exit 2 (合格扱いにしない)。
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.geometry.wall_axismach import (JointFitCFDWall, joint_fit_wall,  # noqa: E402
                                                 r3_piecewise_exact)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def rejects(name, fn, must=None):
    try:
        fn()
    except (ValueError, TypeError) as e:
        check(f"{name} を拒否 ({str(e)[:90]})", must is None or must in str(e))
        return
    check(f"{name} を拒否 (通ってしまった)", False)


def _ref_joint_fit_wall(wall_tbl, R, lam=1e-9, k=5, sig_r=1e-6, sig_th=1e-4, h0=0.0125, h1=0.5, x_g=6.0):
    """改修前 (commit 前の HEAD) の `joint_fit_wall` の写し。既定経路のビット同一の基準 (手で変えないこと)。"""
    tb = np.asarray(wall_tbl, dtype=float)
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (x_g - x0), 1.0)
        xs.append(xs[-1] + h0 + (h1 - h0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1])
    xi = xi[xi < xe - 0.5 * h1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]
    nc = len(t) - k - 1
    E = np.eye(nc)

    def D(xq, d):
        return np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1)
    w = np.gradient(x)
    w = w / w.mean()
    A = (B0.T * (w / sig_r ** 2)) @ B0 + (B1.T * (w / sig_th ** 2)) @ B1
    b = B0.T @ (w * r / sig_r ** 2) + B1.T @ (w * np.tan(th) / sig_th ** 2)
    if lam > 0:
        xq = np.linspace(x0, xe, 8000)
        B3 = D(xq, 3)
        A = A + lam * (B3.T * np.gradient(xq)) @ B3 / sig_r ** 2
    Cm = np.array([D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0],
                   D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]])
    dv = np.array([r[0], x0 / R, 1.0 / R, r[-1], np.tan(th[-1])])
    K = np.block([[A, Cm.T], [Cm, np.zeros((5, 5))]])
    c = np.linalg.solve(K, np.r_[b, dv])[:nc]
    return BSpline(t, c, k), nc


def s1_check(spl, a, b, R, tag):
    """§6 S1 の形状用検査 (形状ゲート throat_mono_shape_gate.py と同じ関数・同じ閾値)。"""
    m = r3_piecewise_exact(spl, None, a, b)
    check(f"{tag}: [{a}, {b}] の r‴ 最大 {m['r3_max']:.2e} ≤ 1e-6 (区間多項式で厳密)", m["r3_max"] <= 1e-6)
    check(f"{tag}: r″ の最大増加 {m['r2_max_increase']:.2e} ≤ 1e-7", m["r2_max_increase"] <= 1e-7)
    check(f"{tag}: r″ 最大 {m['r2_max']:.10f} ≤ 1/R + 1e-9", m["r2_max"] <= 1.0 / R + 1e-9)
    return m


def eq_resid(spl, tb, R):
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    v = [spl(x[0]) - r[0], spl(x[0], 1) - x[0] / R, spl(x[0], 2) - 1.0 / R, spl(x[-1]) - r[-1], spl(x[-1], 1) - np.tan(th[-1])]
    return float(np.max(np.abs(v)))


# --- 1. 合成の点群 ----------------------------------------------------------------
# スロート始点 (x=0, r=1)、r″(0)=1/R から始まり、直後 (x ≲ 0.06) に 1/R より強く曲がる (現行の当てはめで r″ の山ができる形)
R = 2.0
xs_ = [0.0]
while xs_[-1] < 12.0:
    xs_.append(xs_[-1] + 0.025 + 0.02 * xs_[-1])
xt = np.array(xs_)
r1_raw = xt / R + (0.002 * (2 * xt - xt ** 2 / 0.06)) * np.exp(-xt / 0.06)
# 下流は r′ が 0.35 付近で頭打ち (θ の最大) → 逆 MOC 壁らしく曲率が下がる
r1 = np.tanh(r1_raw / 0.35) * 0.35
r_ = 1.0 + np.r_[0.0, np.cumsum(0.5 * (r1[1:] + r1[:-1]) * np.diff(xt))]
tb_syn = np.c_[xt, r_, np.arctan(r1)]
for lam in (1e-10, 1e-9, 1e-8, 0.0):
    s_ref, n_ref = _ref_joint_fit_wall(tb_syn, R, lam=lam)
    dg = {}
    s_new, n_new = joint_fit_wall(tb_syn, R, lam=lam, diag=dg)
    check(f"合成 λ={lam:g}: mono なしは改修前と係数・ノットが完全一致", n_ref == n_new and np.array_equal(s_ref.c, s_new.c)
          and np.array_equal(s_ref.t, s_new.t))
    if lam > 0:   # λ = 0 (平滑化なし) は KKT 系の条件が悪く残差 ~1e-7 になる。生産は λ 1e-9 なので残差の合否は λ > 0 だけ
        er = eq_resid(s_new, tb_syn, R)   # BSpline の点評価で再計算 (diag は基底行列 @ c なので丸めが違う)
        check(f"合成 λ={lam:g}: mono なしの等式残差 {dg['eq_resid']:.1e} (再計算 {er:.1e}) ≤ 1e-9",
              dg["eq_resid"] <= 1e-9 and er <= 1e-9)
m0 = r3_piecewise_exact(_ref_joint_fit_wall(tb_syn, R)[0], None, 0.0, 1.5)
check(f"合成: 拘束なしでは r″ が [0, 1.5] で増える (最大増加 {m0['r2_max_increase']:.2e} > 1e-4; 試験の前提)", m0["r2_max_increase"] > 1e-4)
for lam in (1e-10, 1e-9, 1e-8):
    dg = {}
    sm, _ = joint_fit_wall(tb_syn, R, lam=lam, mono_r2=(0.0, 1.5), diag=dg)
    s1_check(sm, 0.0, 1.5, R, f"合成 λ={lam:g} mono")
    check(f"合成 λ={lam:g} mono: 等式残差 {dg['eq_resid']:.1e} ≤ 1e-9", dg["eq_resid"] <= 1e-9 and eq_resid(sm, tb_syn, R) <= 1e-9)
    check(f"合成 λ={lam:g} mono: KKT 記録 (有効 {dg['n_active']}・反復 {dg['iters']}・乗数最小 {dg['mu_min']})",
          dg["n_active"] >= 1 and dg["iters"] >= 2 and dg["mu_min"] is not None
          and dg["mu_min"] >= -1e-12 * max(1.0, abs(dg["mu_min"])) and dg["ineq_max_normalized"] <= 1e-10
          and len(dg["active_support"]) == dg["n_active"] and dg["mono_r2"] == [0.0, 1.5])
    check(f"合成 λ={lam:g} mono: 有効制約の台が (0, 1.5) にかかる",
          all(lo < 1.5 and hi > 0.0 for lo, hi in dg["active_support"]))
# 区間の外では拘束しない: [0.5, 1.5] だけ拘束したら [0, 0.5) の山は残りうる (拘束が区間に閉じていること)
dg = {}
sp_part, _ = joint_fit_wall(tb_syn, R, mono_r2=(0.5, 1.5), diag=dg)
mp = r3_piecewise_exact(sp_part, None, 0.5, 1.5)
check(f"合成 mono [0.5, 1.5]: その区間で r‴ 最大 {mp['r3_max']:.1e} ≤ 1e-6・増加 {mp['r2_max_increase']:.1e} ≤ 1e-7",
      mp["r3_max"] <= 1e-6 and mp["r2_max_increase"] <= 1e-7)
# r3_piecewise_exact (形状用検査) を密な点評価と照合: 厳密値は密な標本の上限・積分は台形則の極限と一致
s_chk, _ = _ref_joint_fit_wall(tb_syn, R)
mx = r3_piecewise_exact(s_chk, None, 0.0, 1.5)
xd = np.linspace(0.0, 1.5, 1500001)
r3d, r2d = s_chk(xd, 3), s_chk(xd, 2)
inc_d = float(np.max(r2d - np.minimum.accumulate(r2d)))
check(f"r3_piecewise_exact: r‴ 最大 {mx['r3_max']:.6g} ≥ 密標本 {r3d.max():.6g} (下回りは丸め 1e-10 まで、上回り ≤ 1e-6)", mx["r3_max"] >= r3d.max() - 1e-10
      and mx["r3_max"] - r3d.max() <= 1e-6)
check(f"r3_piecewise_exact: r″ 最大増加 {mx['r2_max_increase']:.6g} ≈ 密標本 {inc_d:.6g}", abs(mx["r2_max_increase"] - inc_d) <= 1e-9
      and mx["r2_max_increase"] >= inc_d - 1e-12)
check(f"r3_piecewise_exact: ∫(r‴)² {mx['int_r3sq']:.8g} ≈ 台形則 {np.trapz(r3d ** 2, xd) if hasattr(np, 'trapz') else np.trapezoid(r3d ** 2, xd):.8g}",
      abs(mx["int_r3sq"] - (np.trapz(r3d ** 2, xd) if hasattr(np, "trapz") else np.trapezoid(r3d ** 2, xd))) <= 1e-6 * max(1.0, mx["int_r3sq"]))
check(f"r3_piecewise_exact: max|r⁗| {mx['r4_absmax']:.6g} ≥ 密標本 {np.abs(s_chk(xd, 4)).max():.6g} (差 ≤ 1e-6 相対)",
      mx["r4_absmax"] >= np.abs(s_chk(xd, 4)).max() * (1 - 1e-9) and mx["r4_absmax"] <= np.abs(s_chk(xd, 4)).max() * (1 + 1e-6))
rejects("r3_piecewise_exact の区間 a ≥ b", lambda: r3_piecewise_exact(s_chk, None, 1.5, 1.5))
rejects("r3_piecewise_exact の区間 NaN", lambda: r3_piecewise_exact(s_chk, None, float("nan"), 1.5))
# 不正入力
for bad, why in (((1.5, 0.0), "a > b"), ((0.5, 0.5), "a = b"), ((-0.1, 1.5), "x 範囲の下"), ((0.0, 99.0), "x 範囲の上"),
                 ((0.0,), "長さ 1"), ((0.0, 1.0, 1.5), "長さ 3"), ("01", "文字列"), ((0.0, float("nan")), "NaN"),
                 ((0.0, float("inf")), "inf"), ((None, 1.5), "None 要素"), (1.5, "スカラー"), ((), "空")):
    rejects(f"mono_r2 = {bad!r} ({why})", lambda bad=bad: joint_fit_wall(tb_syn, R, mono_r2=bad))
rejects("JointFitCFDWall(mono_r2=(1.5, 0))", lambda: JointFitCFDWall(tb_syn, R, mono_r2=(1.5, 0.0)))
# 数値の書き方 (YAML から来うる型): list・numpy・整数・指数表記の文字列要素は数として解釈
for ok_in in ([0, 1.5], np.array([0.0, 1.5]), ("0", "1.5e0")):
    dg = {}
    s_ok, _ = joint_fit_wall(tb_syn, R, mono_r2=ok_in, diag=dg)
    check(f"mono_r2 = {ok_in!r} は [0.0, 1.5] として解釈", dg["mono_r2"] == [0.0, 1.5])
# YAML の書き方 (problem の geometry.wall_fit_mono_r2 は probdef.load_problem = yaml.safe_load で読まれる)
import yaml  # noqa: E402
for txt, want in (("k: [0.0, 1.5]", [0.0, 1.5]), ("k:\n  - 0\n  - 1.5", [0.0, 1.5]), ("k: [0, 1.5e0]", [0.0, 1.5]),
                  ("k: [0., 15e-1]", [0.0, 1.5]), ("k: ['0', '1.5']", [0.0, 1.5]), ("k: [0.0, 1.5]  # コメント", [0.0, 1.5])):
    v = yaml.safe_load(txt)["k"]
    dg = {}
    joint_fit_wall(tb_syn, R, mono_r2=v, diag=dg)
    check(f"YAML {txt!r} → {v!r} は {want} として解釈", dg["mono_r2"] == want)
for txt in ("k: '0, 1.5'", "k: [0.0]", "k: [1.5, 0.0]", "k: [.nan, 1.5]", "k: [-.inf, 1.5]", "k: {a: 0, b: 1.5}", "k: [0.0, 1.5, 2.0]", "k: 1.5", "k: [false, true]"):
    v = yaml.safe_load(txt)["k"]
    rejects(f"YAML {txt!r} ({v!r})", lambda v=v: joint_fit_wall(tb_syn, R, mono_r2=v))
# 重複キーは PyYAML では後勝ち (probdef.load_problem は yaml.safe_load の直読み; 本 plan の範囲外として既知)。ここでは事実だけ固定する
check("YAML の重複キーは後勝ちで読まれる (probdef の既知の制約)", yaml.safe_load("k: [0.0, 1.5]\nk: [0.0, 1.0]")["k"] == [0.0, 1.0])
# JointFitCFDWall の fit_diag と spline からの復元
wj = JointFitCFDWall(tb_syn, R, mono_r2=[0.0, 1.5])
fd = wj.fit_diag
s_rec = BSpline(np.array(fd["spline"]["t"]), np.array(fd["spline"]["c"]), fd["spline"]["k"])
xx = np.linspace(0.0, float(xt[-1]), 20001)
check("JointFitCFDWall: fit_diag の spline から壁が完全に復元できる", np.array_equal(s_rec(xx), wj.r(xx)))
check("JointFitCFDWall: fit_diag に mono_r2・n_active・eq_resid", fd["mono_r2"] == [0.0, 1.5] and fd["n_active"] >= 1
      and fd["eq_resid"] <= 1e-9)
check("JointFitCFDWall: fit_diag は JSON に書ける", isinstance(json.loads(json.dumps(fd)), dict))
check("JointFitCFDWall: validate が空 (mono)", wj.validate() == [])
w0 = JointFitCFDWall(tb_syn, R)
check("JointFitCFDWall(mono なし): fit_diag の mono_r2 は None、KKT の項目は無い",
      w0.fit_diag["mono_r2"] is None and "n_active" not in w0.fit_diag and "spline" in w0.fit_diag)

# --- 2. 最終問題 (run データ) ------------------------------------------------------
C45 = ROOT / "case/45.isobutane_m6_d155"
cands = [Path(os.environ["CASE_RUNS"])] if os.environ.get("CASE_RUNS") else [C45, Path("/home/sano/work/forge/case/45.isobutane_m6_d155")]
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402

p = load_problem(C45 / "problem_d155_ns_finemesh_recal_final.yaml")
src = next((c / p.geometry["initial_line_run"] for c in cands if (c / p.geometry["initial_line_run"] / p.geometry["initial_line_res"]).exists()), None)
if src is None:
    print(f"凍結源 run {p.geometry['initial_line_run']} が {[str(c) for c in cands]} に無い — 2. を飛ばす (CASE_RUNS で与える)")
    print("FAIL 件数:", FAIL)
    sys.exit(1 if FAIL else 2)
p.geometry["initial_line_run"] = str(src.resolve())
d = design_chain(p)
tb = d["wall_inv"]
Rf = float(d["R"])
W = d["wall"]
s_ref, _ = _ref_joint_fit_wall(tb, Rf)
xg = np.linspace(W.x_in, W.x_e, 400001)
dr0 = float(np.max(np.abs(W._spl(xg[xg >= 0]) - s_ref(xg[xg >= 0]))))
check(f"最終問題: キー無しの設計壁が改修前と max|Δr| = {dr0} (= 0)", dr0 == 0.0)
check("最終問題: キー無しの係数・ノットが改修前と完全一致", np.array_equal(W._spl.c, s_ref.c) and np.array_equal(W._spl.t, s_ref.t))
check(f"最終問題: キー無しの等式残差 {d['wall_fit']['eq_resid']:.1e} ≤ 1e-9", d["wall_fit"]["eq_resid"] <= 1e-9)
check("最終問題: キー無しの wall_fit.mono_r2 は None", d["wall_fit"]["mono_r2"] is None)
pm = load_problem(C45 / "problem_d155_ns_finemesh_recal_final.yaml")
pm.geometry["initial_line_run"] = str(src.resolve())
pm.geometry["wall_fit_mono_r2"] = [0.0, 1.5]
dm = design_chain(pm)
Wm = dm["wall"]
check("最終問題 mono: 点群は現行と同一 (当てはめだけが違う)", np.array_equal(dm["wall_inv"], tb))
s1_check(Wm._spl, 0.0, 1.5, Rf, "最終問題 mono")
fm = dm["wall_fit"]
check(f"最終問題 mono: 等式残差 {fm['eq_resid']:.1e} ≤ 1e-9", fm["eq_resid"] <= 1e-9)
check(f"最終問題 mono: wall_fit に mono_r2 {fm['mono_r2']}・有効 {fm['n_active']}・反復 {fm['iters']}・spline",
      fm["mono_r2"] == [0.0, 1.5] and fm["n_active"] >= 1 and len(fm["spline"]["c"]) == fm["n_cp"])
s_rec = BSpline(np.array(fm["spline"]["t"]), np.array(fm["spline"]["c"]), fm["spline"]["k"])
check("最終問題 mono: wall_fit の spline から設計壁が完全に復元できる", np.array_equal(s_rec(xg[xg >= 0]), Wm.r(xg[xg >= 0])))
check("最終問題 mono: x < 0 (上流 Hermite) は現行と同一", np.array_equal(Wm.r(xg[xg < 0]), W.r(xg[xg < 0])))
for repr_ in ("interp", "lsq"):
    pb = load_problem(C45 / "problem_d155_ns_finemesh_recal_final.yaml")
    pb.geometry.update({"initial_line_run": str(src.resolve()), "wall_repr": repr_, "wall_fit_mono_r2": [0.0, 1.5]})
    rejects(f"wall_repr {repr_} + wall_fit_mono_r2", lambda pb=pb: design_chain(pb), must="joint 専用")
pn = load_problem(C45 / "problem_d155_ns_finemesh_recal_final.yaml")
pn.geometry.update({"initial_line_run": str(src.resolve()), "wall_fit_mono_r2": None})
dn = design_chain(pn)
check("最終問題: wall_fit_mono_r2: null は拘束なし (キー無しと同一の壁)", np.array_equal(dn["wall"]._spl.c, W._spl.c))
pr = load_problem(C45 / "problem_d155_ns_finemesh_recal_final.yaml")
pr.geometry.update({"initial_line_run": str(src.resolve()), "wall_fit_mono_r2": [0.0, 1.0e3]})
rejects("最終問題: wall_fit_mono_r2 の上端が出口より下流", lambda: design_chain(pr), must="x 範囲")

print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
