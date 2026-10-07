#!/usr/bin/env python3
"""逆 MOC の軸上の解析極限 θ_r と予測修正の収束 (`geometry.moc_axis_limit` / `geometry.moc_corrector`) の試験。
plan: plans/accepted/discretization-moc-axis-limit-and-corrector.md §5.1 #3 (§4.0〜§4.3)。

1. θ_r の式: CPG の閉形式と一致 / 放射源流で 1/x / semi-perfect (case/45 のガス) で ν(M) 表の差分・
   連続の式 −½ d ln(ρu)/dx の数値微分と一致 (表の補間の違いによる差を記録)
2. 収束判定: converge の対は更新量 ≤ 1e-12 で止まり、最終残差 ≤ 1e-10。対ごとの結果が同じ段の他の対に依存しない。
   停止判定は有効な対だけ (NaN の対が上限まで回さない)。fixed2 は従来の回数
3. 負例 (§4.2 の 5 分類): 入力 NaN の対 = もともとの欠損、反復中の NaN = 反復の失敗、平行な特性線 = 幾何的棄却、
   上限到達 = 反復の失敗。ゲート (`moc_gate`) が反復の失敗・最終残差・壁の内側の幾何的棄却・接続不一致で不合格にする
4. 軸上の分岐: analytic は真の軸端点で θ_r (相手を借りない)、AXIS_LIMIT_FRAC の発火を数える、θ_r の無い軸上点は例外
5. 軸端の θ_r の組み立て (`axis_theta_r_init`): アンカー / 一般経路 (x_A ≠ x0) / 微分が無いときの ν スプライン
6. キーの検査: 不正値 (null・大文字・空白・数値・真偽値・空・リスト) は例外、YAML の書き方の揺れ (引用符・flow・merge
   キー・コロン前の空白・コメント) は同じ値、重複キーは後勝ち (probdef の既知の制約を固定)
7. §4.0 の判別 A/B (放射源流、粗い 2 解像度): 第 1 段の θ 誤差が B ≤ A/10、壁の誤差 B ≤ A、未収束 0
   (4 解像度の本番は `moc_axis_limit_radial.py`)
8. 既定のビット同一: 実装前のコミット (BASE) の package を `git archive` で取り出し、
   (a) 放射源流の網 (fill_levels) と (b) case/45 単調壁の問題の MOC 点群・当てはめ後の係数とノットが完全一致。
   (b) は凍結源 run_0062 が要る ($CASE_RUNS、既定 /home/sano/work/forge/case/45.isobutane_m6_d155)。

usage: python3 design/tests/run_moc_axis_limit_tests.py
  凍結源が無ければ 8 (b) と case/45 依存の項目を飛ばして exit 2 (合格扱いにしない)。
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from forge_design.evaluate.runner_axismach import _moc_keys, design_chain, load_problem  # noqa: E402
from forge_design.geometry import moc_kernel as MK  # noqa: E402
from forge_design.geometry.moc_inverse import (InverseMOC, axis_theta_r_init,  # noqa: E402
                                               inverse_design, moc_gate)
from forge_design.geometry.moc_kernel import (PAIR_CONVERGED, PAIR_DONE, PAIR_GEOM,  # noqa: E402
                                              PAIR_MAXITER, PAIR_MISSING, PAIR_NONFINITE,
                                              axis_theta_r, dnu_dM, interior_vec, pm_nu)
import moc_axis_limit_radial as RAD  # noqa: E402

FAIL = 0
SKIPPED = []
BASE = "170b0d75924fa695134be0840e32fb37e5754917"   # 実装前 (§5.1 #3 着手時の HEAD)。既定のビット同一の基準 (変えないこと)
C45 = ROOT / "case/45.isobutane_m6_d155"
PROB45 = "problem_d155_ns_finemesh_recal_final_mono.yaml"
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
HAVE45 = (C45 / PROB45).exists() and (RUNS / "run_0062_euler_wallfit_fit_r1_ext6k/res_6000.h5").exists()


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def rejects(name, fn, must=None):
    try:
        fn()
    except (ValueError, TypeError) as e:
        check(f"{name} を拒否 ({str(e)[:80]})", must is None or must in str(e))
        return
    check(f"{name} を拒否 (通ってしまった)", False)


def arr(*v):
    return [np.array(x, dtype=float) for x in v]


G = 1.4

# --- 1. θ_r の式 ---------------------------------------------------------------------------------
Ms = np.array([1.001, 1.05, 1.2, 1.5, 2.0, 3.0, 4.5, 6.0])
Mp = np.array([0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01])
for gam in (1.4, 1.3, 1.1):
    ref = 0.5 * (Ms ** 2 - 1.0) / (Ms * (1.0 + 0.5 * (gam - 1.0) * Ms ** 2)) * Mp
    e = float(np.max(np.abs(axis_theta_r(Ms, Mp, gam) / ref - 1.0)))
    check(f"θ_r: CPG γ {gam} の閉形式 ½(M²−1)M′/[M(1+(γ−1)M²/2)] と一致 (max rel {e:.1e})", e < 1e-13)
xs_r = np.linspace(1.2, 5.5, 37)
e = max(abs(float(axis_theta_r(RAD.src_M(x, 0.0), RAD.src_dMdx_axis(x), G)) * x - 1.0) for x in xs_r)
check(f"θ_r: 放射源流の軸上で 1/x (max rel {e:.1e})", e < 1e-10)

if (C45 / PROB45).exists():
    gas = load_problem(C45 / PROB45).gas_model
    Mq = np.linspace(1.2, 5.8, 24)
    # (i) dnu_dM (ν 表の 3 次スプライン微分) と ν 表の中心差分 (本体と同じ線形補間の ν)
    h = 2e-3
    fd_nu = (gas.nu(Mq + h) - gas.nu(Mq - h)) / (2 * h)
    e_nu = float(np.max(np.abs(dnu_dM(Mq, gas) / fd_nu - 1.0)))
    check(f"semi-perfect: ν_M (spline 微分) と ν 表の中心差分が一致 (max rel {e_nu:.2e} ≤ 1e-3)", e_nu < 1e-3)
    # (ii) θ_r と −½ d ln(ρu)/dx (ρu = mass_flux_density の表、連続の式) の数値微分 (M′ = 1 で比較)
    fd_F = (np.log(gas.mass_flux_density(Mq + h)) - np.log(gas.mass_flux_density(Mq - h))) / (2 * h)
    e_F = float(np.max(np.abs(axis_theta_r(Mq, 1.0, gas) / (-0.5 * fd_F) - 1.0)))
    check(f"semi-perfect: θ_r = −½ d ln(ρu)/dx (連続の式の数値微分) と一致 (max rel {e_F:.2e} ≤ 2e-3)", e_F < 2e-3)
    # case/45 のアンカー (M_A 1.21476, M′_A 0.63528) で spline 微分と線形補間の勾配の差 (諮問の独立計算 0.0156 %)
    MA = 1.2147584415478936
    i = int(np.searchsorted(gas._M, MA))
    lin = (gas._nu[i] - gas._nu[i - 1]) / (gas._M[i] - gas._M[i - 1])
    d_sl = float(dnu_dM(MA, gas) / lin - 1.0)
    print(f"info semi-perfect: case/45 のアンカーで ν_M の spline / 線形補間の勾配の差 {d_sl * 100:+.4f} % "
          f"(θ_r = {float(axis_theta_r(MA, 0.6352781033834302, gas)):.6f} rad/r_t)")
    check("semi-perfect: アンカーの spline/線形の勾配差が 0.05 % 以内 (熱力学近似誤差として記録)", abs(d_sl) < 5e-4)
    # 表の解像度依存 (plan §4.1: 反復の許容差とは別の熱力学近似誤差として記録)。n_tab 6000 と 24000 の θ_r の差
    from forge_design.gas import GasSemiPerfect  # noqa: E402
    _p45 = load_problem(C45 / PROB45)
    _Y, _, _ = _p45.gas_composition
    _th = {n: float(axis_theta_r(MA, 0.6352781033834302, GasSemiPerfect(_Y, Tt=float(_p45.spec["Tt"]), db=_p45.species_db,
                                                                         n_tab=n))) for n in (6000, 24000)}
    _d = _th[6000] / _th[24000] - 1.0
    print(f"info semi-perfect: θ_r の表の解像度依存 n_tab 6000 → 24000 で {_d:+.2e} (相対)")
    check("semi-perfect: θ_r の表の解像度依存が 1e-5 以内", abs(_d) < 1e-5)
else:
    SKIPPED.append("1 semi-perfect (case/45 の問題 YAML が無い)")

# --- 2. 収束判定 -----------------------------------------------------------------------------------
n_ax, n_st = 60, 13
init = RAD.build_init(n_ax, n_st)
thr, info = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                              M_line_axis=RAD.src_M(RAD.X1, 0.0))
A0 = np.array([[p.x, p.r, p.th, p.nu, p.M] for p in init[1:]])
B0 = np.array([[p.x, p.r, p.th, p.nu, p.M] for p in init[:-1]])
st = {}
q = interior_vec(*A0.T, *B0.T, G, 1.0, 2, axis_limit="analytic", thrA=thr[1:], thrB=thr[:-1],
                 corrector="converge", stats=st)
okc = st["status"] == PAIR_CONVERGED
check(f"converge: 全対が収束 ({int(okc.sum())}/{okc.size}、反復 max {int(st['n_iter'].max())})", bool(okc.all()))
check(f"converge: 最終残差 ≤ 1e-10 (幾何 {np.nanmax(st['resid_geom']):.1e}、適合式 {np.nanmax(st['resid_comp']):.1e})",
      float(np.nanmax(st["resid_geom"])) <= MK.RESID_TOL and float(np.nanmax(st["resid_comp"])) <= MK.RESID_TOL)
# 収束した値で単位過程をもう 1 回まわすと動かない (≤ 1e-12 級)
st1 = {}
q1 = interior_vec(*A0.T, *B0.T, G, 1.0, 2, axis_limit="analytic", thrA=thr[1:], thrB=thr[:-1],
                  corrector="converge", tol=1e-14, max_corr=80, stats=st1)
dq = float(np.max(np.abs(q1[2] - q[2])))
check(f"converge: 許容差 1e-14 で回し直した値との差 ≤ 1e-11 ({dq:.1e})", dq <= 1e-11)
# 対ごとの結果が同じ段の他の対に依存しない (1 対ずつ通した結果と一括の結果)
dmax = 0.0
for j in (0, 5, n_ax - 1, n_ax, n_ax + 3):
    qs = interior_vec(*A0[j:j + 1].T, *B0[j:j + 1].T, G, 1.0, 2, axis_limit="analytic",
                      thrA=thr[j + 1:j + 2], thrB=thr[j:j + 1], corrector="converge")
    dmax = max(dmax, float(np.max(np.abs(np.array(qs[:4])[:, 0] - np.array(q[:4])[:, j]))))
check(f"converge: 1 対ずつと一括の結果が一致 (CPG の Newton の丸めまで, max {dmax:.1e} ≤ 1e-14)", dmax <= 1e-14)
# fixed2 は従来どおり n_corr 回
stf = {}
interior_vec(*A0.T, *B0.T, G, 1.0, 2, stats=stf)
check("fixed2: 状態は DONE、反復回数 = n_corr (2)", bool(np.all(stf["status"] == PAIR_DONE)) and int(stf["n_iter"].max()) == 2)
# 停止判定は有効な対だけ: NaN の対を混ぜても他の対の反復回数は変わらず、NaN の対はもともとの欠損
An, Bn = A0.copy(), B0.copy()
An[3, 2] = np.nan
stn = {}
interior_vec(*An.T, *Bn.T, G, 1.0, 2, axis_limit="analytic", thrA=thr[1:], thrB=thr[:-1], corrector="converge", stats=stn)
check("converge: 入力 NaN の対は『もともとの欠損』で、他の対の反復回数は変わらない",
      int(stn["status"][3]) == PAIR_MISSING and np.array_equal(np.delete(stn["n_iter"], 3), np.delete(st["n_iter"], 3)))

# --- 3. 負例 (5 分類とゲート) ----------------------------------------------------------------------
# (a) 反復中の NaN (pm_mach_vec を差し替えて 2 回目の修正で 1 要素を NaN にする)
_orig = MK.pm_mach_vec
_calls = {"n": 0}


def _nan_injector(nu, g=1.4, *a, **k):
    _calls["n"] += 1
    out = np.array(_orig(nu, g, *a, **k), dtype=float)
    if _calls["n"] == 3 and out.size:
        out.flat[0] = np.nan
    return out


MK.pm_mach_vec = _nan_injector
try:
    stx = {}
    qx = interior_vec(*A0[20:24].T, *B0[20:24].T, G, 1.0, 2, axis_limit="analytic", thrA=thr[21:25],
                      thrB=thr[20:24], corrector="converge", stats=stx)
finally:
    MK.pm_mach_vec = _orig
check(f"負例: 反復中の NaN は『反復の失敗 (非有限)』で ok=False (状態 {stx['status'].tolist()})",
      int(stx["status"][0]) == PAIR_NONFINITE and not bool(qx[4][0]) and bool(np.all(stx["status"][1:] == PAIR_CONVERGED)))
# (b) 平行な特性線: θ_B − θ_A = μ_A + μ_B + 2μ_P (M 2 で 120°) → 予測子で tan(θ+μ) = tan(θ−μ)
nu2 = float(pm_nu(2.0, G))
Ap = arr([1.0], [0.5], [-np.pi / 3], [nu2], [2.0])
Bp = arr([1.2], [0.6], [np.pi / 3], [nu2], [2.0])
for corr in ("converge", "fixed2"):
    stp = {}
    qp = interior_vec(*Ap, *Bp, G, 1.0, 2, corrector=corr, stats=stp)
    check(f"負例 ({corr}): 平行な特性線は『幾何的棄却』で ok=False", int(stp["status"][0]) == PAIR_GEOM and not bool(qp[4][0]))
# (c) 上限到達: 収束に 2 回以上要る対 (軸の第 1 段) を上限 1 回で
stm = {}
qm = interior_vec(*A0[:3].T, *B0[:3].T, G, 1.0, 2, axis_limit="legacy", corrector="converge", max_corr=1, stats=stm)
check(f"負例: 上限到達は『反復の失敗 (上限)』で ok=False (状態 {stm['status'].tolist()})",
      bool(np.all(stm["status"] == PAIR_MAXITER)) and not bool(np.any(qm[4])))
# (d) 充填の集計とゲート: 上限 1 回の充填は反復の失敗で不合格、既定上限では合格
inv1 = InverseMOC(gamma=G, delta=1.0, corrector="converge", max_corr=1)
lev1 = inv1.fill_levels(init)
d1 = dict(inv1.last_diag)
g1 = moc_gate(d1, None, None)
check(f"ゲート: 反復の失敗があれば不合格 (上限到達 {d1['pairs']['iter_maxiter']} 対 → {g1['reasons']})",
      d1["pairs"]["iter_maxiter"] > 0 and g1["pass"] is False)
pc = d1["pairs"]
check("集計: 入力 = 欠損 + 幾何的棄却 + 反復の失敗 + 収束 (+ fixed2 の完了)",
      pc["input"] == sum(pc[k] for k in ("missing", "geom_parallel", "geom_below_axis", "iter_nonfinite",
                                          "iter_maxiter", "converged", "done_fixed")))
inv2 = InverseMOC(gamma=G, delta=1.0, axis_limit="analytic", corrector="converge")
lev2 = inv2.fill_levels(init, axis_thr=thr)
d2 = dict(inv2.last_diag)
g2 = moc_gate(d2, info, None)
check(f"ゲート: 既定上限 (50) の analytic + converge は合格 (失敗 0、残差 {g2['resid_max']:.1e})",
      g2["pass"] is True and d2["pairs"]["missing"] == 0)
# もともとの欠損 (初期前線の 1 点を NaN にすると、その点に依存する対が段ごとに死ぬ) は失敗に数えない
import copy  # noqa: E402
init_n = copy.deepcopy(init)
init_n[30].th = float("nan")
inv_n = InverseMOC(gamma=G, delta=1.0, axis_limit="analytic", corrector="converge")
inv_n.fill_levels(init_n, axis_thr=thr)
d_n = dict(inv_n.last_diag)
g_n = moc_gate(d_n, info, None)
check(f"ゲート: もともとの欠損 ({d_n['pairs']['missing']} 対) は反復の失敗に数えず合格",
      d_n["pairs"]["missing"] > 0 and d_n["pairs"]["iter_nonfinite"] + d_n["pairs"]["iter_maxiter"] == 0
      and g_n["pass"] is True)
# 最終残差のしきい値
d3 = dict(d2, resid=dict(d2["resid"], comp_max=2e-10))
check("ゲート: 最終残差 > 1e-10 で不合格", moc_gate(d3, info, None)["pass"] is False)
# 幾何的棄却の位置: 壁の外 (上) は許す、壁の内側は不合格
wall_t = np.array([[0.0, 1.0, 0, 0], [5.0, 1.0, 0, 0]])
d4 = dict(d2, pairs=dict(d2["pairs"], geom_parallel=1), _geom_all=[(3, 7, "parallel", 1.0, 1.2, 1.1, 1.3)])
d5 = dict(d4, _geom_all=[(3, 7, "parallel", 1.0, 0.4, 1.1, 0.5)])
d6 = dict(d4, _geom_all=[(3, 7, "parallel", 6.0, 0.4, 6.1, 0.5)])
check("ゲート: 幾何的棄却が壁の外 (両端が壁より上) なら許す", moc_gate(d4, info, wall_t)["pass"] is True)
check("ゲート: 幾何的棄却が壁の内側なら不合格", moc_gate(d5, info, wall_t)["pass"] is False)
check("ゲート: 幾何的棄却が壁の x 範囲の外 (網の端) なら許す", moc_gate(d6, info, wall_t)["pass"] is True)
check("ゲート: fixed2 は合否を出さない (pass = None)", moc_gate(dict(d2, corrector="fixed2"), info, None)["pass"] is None)

# --- 4. 軸上の分岐 ---------------------------------------------------------------------------------
# 軸節点 2 個の対: legacy は両端 0、analytic は両端 θ_r
sa, sb = {}, {}
interior_vec(*A0[:5].T, *B0[:5].T, G, 1.0, 2, stats=sa)
interior_vec(*A0[:5].T, *B0[:5].T, G, 1.0, 2, axis_limit="analytic", thrA=thr[1:6], thrB=thr[:5], stats=sb)
check("分岐: legacy の軸節点 2 個の対は両端とも 0 (axis_zero = 2/対)", int(sa["branch"]["axis_zero"].sum()) == 10)
check("分岐: analytic の軸節点 2 個の対は両端とも θ_r (axis_analytic = 2/対)", int(sb["branch"]["axis_analytic"].sum()) == 10)
# AXIS_LIMIT_FRAC: 軸外の A が r_B の 5 % 未満なら相手を借りる (発火を数える)。真の軸端点は先に θ_r
nuA = float(pm_nu(1.6, G))
Af = arr([1.00], [0.01], [0.002], [nuA], [1.6])
Bf = arr([1.05], [0.40], [0.08], [nuA + 0.01], [1.62])
sf = {}
interior_vec(*Af, *Bf, G, 1.0, 2, axis_limit="analytic", thrA=[np.nan], thrB=[np.nan], stats=sf)
check("分岐: 軸外の A (r 0.01 < 0.05·0.40) は AXIS_LIMIT_FRAC の発火として数える", bool(sf["branch"]["frac_A"][0]))
v_leg = MK._src_vec(np.array([0.0]), np.array([0.0]), np.array([0.4]), np.array([0.08]), np.array([0.3]), False)
v_ana = MK._src_vec(np.array([0.0]), np.array([0.0]), np.array([0.4]), np.array([0.08]), np.array([0.3]), True)
check(f"分岐: 真の軸端点は analytic で θ_r (0.3)、legacy は相手の sinθ/r ({float(v_leg[0]):.4f})",
      float(v_ana[0]) == 0.3 and abs(float(v_leg[0]) - np.sin(0.08) / 0.4) < 1e-15)
rejects("analytic で θ_r の無い軸上の既知点", lambda: interior_vec(*A0[:2].T, *B0[:2].T, G, 1.0, 2, axis_limit="analytic",
                                                         thrA=[np.nan, np.nan], thrB=[np.nan, np.nan]), "θ_r が無い")
rejects("analytic で axis_thr 無しの充填", lambda: InverseMOC(axis_limit="analytic").fill_levels(init), "axis_thr")
rejects("axis_thr の長さ違い", lambda: InverseMOC(axis_limit="analytic").fill_levels(init, axis_thr=thr[:-1]), "長さ")

# --- 5. 軸端の θ_r の組み立て ---------------------------------------------------------------------
check(f"θ_r 組立: 軸節点は軸則の解析微分、初期線の軸端は target の微分 ({info['source']})",
      info["source"] == {"axis_nodes": "law_derivative", "line_axis_end": "target_derivative@x_line"})
check(f"θ_r 組立: 放射源流で θ_r·x = 1 (max |·−1| {np.nanmax(np.abs(thr[:n_ax + 1] * np.r_[[p.x for p in init[:n_ax]], RAD.X1] - 1)):.1e})",
      float(np.nanmax(np.abs(thr[:n_ax + 1] * np.r_[[p.x for p in init[:n_ax]], RAD.X1] - 1))) < 1e-10)
check("θ_r 組立: 軸端点以外は NaN", bool(np.all(np.isnan(thr[n_ax + 1:]))))
M1, Mp1 = RAD.src_M(RAD.X1, 0.0), RAD.src_dMdx_axis(RAD.X1)
thr_a, info_a = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                                  axis_anchor=(RAD.X1, M1, Mp1), M_line_axis=M1)
check(f"θ_r 組立: アンカーが初期線の軸端にあればアンカー ({info_a['source']['line_axis_end']}、接続 ok {info_a['connection']['ok']})",
      info_a["source"]["line_axis_end"] == "anchor" and info_a["connection"]["ok"] and thr_a[n_ax] == thr[n_ax])
# 一般経路 x_A ≠ x0: 初期線の軸端を x_A と扱わない (アンカーを使わず target の微分)
thr_g, info_g = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                                  axis_anchor=(RAD.X1 + 0.3, 9.9, 9.9), M_line_axis=M1)
check(f"θ_r 組立: x_A ≠ x0 ではアンカーを使わない ({info_g['source']['line_axis_end']}、Δx "
      f"{info_g['connection']['dx_line_minus_xA']:+.2f})",
      info_g["source"]["line_axis_end"] == "target_derivative@x_line" and thr_g[n_ax] == thr[n_ax]
      and not info_g["connection"]["line_axis_is_x_A"])
# 接続不一致 (初期線の軸端の M が target とずれる) → ok False → ゲートの理由
_, info_b = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                              M_line_axis=M1 + 1e-6)
check(f"θ_r 組立: 軸端の M が 1e-6 ずれたら接続不一致 (dM {info_b['connection']['dM']:+.1e})",
      info_b["connection"]["ok"] is False and moc_gate(d2, info_b, None)["pass"] is False)
_, info_c = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                              axis_anchor=(RAD.X1, M1, Mp1 * (1 + 1e-6)), M_line_axis=M1)
check("θ_r 組立: アンカーの M′ が target の微分とずれたら接続不一致", info_c["connection"]["ok"] is False)
# 微分が無いときの代わり: ν の 3 次スプライン微分 (ν_M を掛けない)
thr_s, info_s = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=None, M_line_axis=M1)
e_s = float(np.max(np.abs(thr_s[:n_ax + 1] / thr[:n_ax + 1] - 1.0)))
check(f"θ_r 組立: 微分が無ければ ν スプライン ({info_s['source']})、解析値との差 {e_s:.1e} ≤ 1e-2",
      info_s["source"] == {"axis_nodes": "nu_spline", "line_axis_end": "nu_spline"} and e_s < 1e-2)
from scipy.interpolate import CubicSpline  # noqa: E402
_xx = np.r_[RAD.X1, [p.x for p in init[:n_ax]][::-1]]
_nn = np.r_[float(pm_nu(M1, G)), [p.nu for p in init[:n_ax]][::-1]]
_ref = 0.5 * np.sqrt(np.array([RAD.src_M(x, 0.0) for x in _xx]) ** 2 - 1) * CubicSpline(_xx, _nn).derivative()(_xx)
check("θ_r 組立: ν スプラインの経路は ½√(M²−1)·dν/dx (ν_M を掛けない)",
      float(np.max(np.abs(np.r_[thr_s[n_ax], thr_s[:n_ax][::-1]] - _ref))) < 1e-14)

# --- 6. キーの検査 ---------------------------------------------------------------------------------
check("キー: 無ければ既定 (legacy, fixed2)", _moc_keys({}) == ("legacy", "fixed2"))
check("キー: analytic + converge", _moc_keys({"moc_axis_limit": "analytic", "moc_corrector": "converge"})
      == ("analytic", "converge"))
for bad in (None, "Analytic", " analytic", "analytic ", "", True, 1, 0, 1e-12, ["analytic"], {"a": 1}, "legacy,converge",
            float("nan")):
    rejects(f"geometry.moc_axis_limit = {bad!r}", lambda b=bad: _moc_keys({"moc_axis_limit": b}))
for bad in (None, "converged", "fixed", "Fixed2", 2, -1, False):
    rejects(f"geometry.moc_corrector = {bad!r}", lambda b=bad: _moc_keys({"moc_corrector": b}))
_yaml_ok = {
    "block": "moc_axis_limit: analytic\nmoc_corrector: converge\n",
    "引用符 (値)": "moc_axis_limit: \"analytic\"\nmoc_corrector: 'converge'\n",
    "引用符 (キー)": "\"moc_axis_limit\": analytic\n'moc_corrector': converge\n",
    "コロン前の空白": "moc_axis_limit : analytic\nmoc_corrector   : converge\n",
    "flow": "{moc_axis_limit: analytic, moc_corrector: converge}\n",
    "merge キー": "_b: &b {moc_axis_limit: analytic, moc_corrector: converge}\n<<: *b\n",
    "行末コメント": "moc_axis_limit: analytic   # 解析極限\nmoc_corrector: converge # 収束\n",
}
for lab, txt in _yaml_ok.items():
    check(f"キー (YAML {lab}): analytic + converge と読む", _moc_keys(yaml.safe_load(txt)) == ("analytic", "converge"))
check("キー (YAML コメント行のキーは無いのと同じ = 既定)",
      _moc_keys(yaml.safe_load("# moc_axis_limit: analytic\nR: 2.0\n")) == ("legacy", "fixed2"))
for lab, txt in (("null", "moc_axis_limit: null\n"), ("~", "moc_axis_limit: ~\n"), ("空", "moc_axis_limit:\n"),
                 ("yes (YAML 1.1 の真偽値)", "moc_corrector: yes\n"), ("指数表記", "moc_corrector: 1e-12\n"),
                 ("リスト", "moc_axis_limit: [analytic]\n")):
    rejects(f"キー (YAML {lab})", lambda t=txt: _moc_keys(yaml.safe_load(t)))
# 重複キーは PyYAML では後勝ち (probdef.load_problem は yaml.safe_load の直読み — 範囲外の既知の制約。事実だけ固定)
check("キー (YAML 重複キー): 後勝ちで読まれる (probdef の既知の制約)",
      _moc_keys(yaml.safe_load("moc_axis_limit: analytic\nmoc_axis_limit: legacy\n")) == ("legacy", "fixed2"))
for bad in (0, -1, 1.5, True, "50"):
    rejects(f"InverseMOC(max_corr={bad!r})", lambda b=bad: InverseMOC(max_corr=b))
for bad in (0.0, -1e-12, float("nan"), float("inf"), "1e-12", None):
    rejects(f"interior_vec(converge, tol={bad!r})",
            lambda b=bad: interior_vec(*A0[:2].T, *B0[:2].T, G, 1.0, 2, corrector="converge", tol=b))
rejects("InverseMOC(axis_limit='bad')", lambda: InverseMOC(axis_limit="bad"))
rejects("InverseMOC(corrector='bad')", lambda: InverseMOC(corrector="bad"))
rejects("interior_vec(corrector='bad')", lambda: interior_vec(*A0[:2].T, *B0[:2].T, G, 1.0, 2, corrector="bad"))
rejects("inverse_design(axis_limit='bad')",
        lambda: inverse_design(None, lambda x: 2.0, 3.0, axis_limit="bad"))
if (C45 / PROB45).exists():
    _p = load_problem(C45 / PROB45)
    _p.geometry["moc_corrector"] = "Converge"
    _p.geometry["initial_line_run"] = "/nonexistent/run"          # キーの検査が先 (凍結源を読む前) であること
    rejects("design_chain: 不正な moc_corrector は凍結源を読む前に例外", lambda: design_chain(_p), "moc_corrector")

# --- 7. §4.0 の判別 A/B (粗い 2 解像度) -----------------------------------------------------------
lv = ((140, 25), (280, 49))
res = {arm: {f"{a}x{b}": RAD.run_arm(arm, a, b, with_restart=(a == 140)) for a, b in lv} for arm in RAD.ARMS}
jd = RAD.judge(res, lv)
r140 = jd["V0"]["first_level_theta_ratio_B_over_A"]
check(f"判別 A/B: 第 1 段の θ 誤差 B/A ≤ 1/10 ({', '.join(f'{k} {v:.1e}' for k, v in r140.items())})", jd["V0"]["c_first_level_le_0p1"])
check("判別 A/B: 壁の誤差 B ≤ A (各解像度)", jd["V0"]["c_wall_B_le_A"])
check(f"判別 A/B: 未収束の対 0 ({jd['V0']['unconverged_pairs']})", sum(jd["V0"]["unconverged_pairs"].values()) == 0)
check(f"判別 A/B: B の壁の観測次数 ≥ 1.7 (140→280: {jd['V0']['orders']['B'][0]:.3f})", jd["V0"]["orders"]["B"][0] >= 1.7)
rs = res["B"]["140x25"]["restart"]
check(f"再出発 (B): 列の 1 段目の作り直しの差 ≤ 1e-11 (dθ {rs['reproc_first_row']['dtheta_max']:.1e})、"
      f"充填し直した網 ≤ 1e-11 (dθ {rs['refill_column']['dtheta_max']:.1e})",
      rs["reproc_first_row"]["dtheta_max"] <= 1e-11 and rs["refill_column"]["dtheta_max"] <= 1e-11)
ra = res["A"]["140x25"]["restart"]
print(f"info 再出発 (A legacy+converge): 1 段目の作り直し dθ {ra['reproc_first_row']['dtheta_max']:.2e} "
      f"(軸端の代用が対の組み方で変わる分)")
# 流線経路 (fill_arrays) でも converge・analytic が通り、集計が残る
inv_s = InverseMOC(gamma=G, delta=1.0, axis_limit="analytic", corrector="converge")
pts_s = inv_s.fill(init, axis_thr=thr)
check(f"流線経路 (fill): analytic + converge で集計が残る (対 {inv_s.last_diag['pairs']['converged']}、失敗 0)",
      inv_s.last_diag["pairs"]["converged"] > 0
      and inv_s.last_diag["pairs"]["iter_nonfinite"] + inv_s.last_diag["pairs"]["iter_maxiter"] == 0)

# --- 8. 既定のビット同一 (実装前のコミットと比べる) ------------------------------------------------
_SNIP = r'''
import sys, numpy as np
sys.path.insert(0, sys.argv[1])          # 比べる package だけを先頭に (moc_axis_limit_radial は現行を足すので使わない)
out = {}
from forge_design.geometry import moc_inverse as MI
from forge_design.geometry.moc_kernel import _Pt, pm_nu
from forge_design.evaluate.ic import invert_area_ratio
assert MI.__file__.startswith(sys.argv[1]), MI.__file__
G, X1 = 1.4, 1.2
RT = X1 * np.tan(np.deg2rad(10.0))
def src_M(x, r):        # moc_axis_limit_radial.src_M と同じ
    R = max(np.hypot(x, r), 1.0001)
    return float(invert_area_ratio(np.array([R * R]), np.array([True]), G)[0])
init = ([_Pt(float(x), 0.0, 0.0, float(pm_nu(src_M(x, 0.0), G)), G) for x in np.linspace(5.5, X1, 141)[:-1]]
        + [_Pt(X1, float(r), float(np.arctan2(r, X1)), float(pm_nu(src_M(X1, r), G)), G) for r in np.linspace(0.0, RT, 25)])
out["lev"] = MI.InverseMOC(gamma=1.4, delta=1.0).fill_levels(init)
out["pts"] = MI.InverseMOC(gamma=1.4, delta=1.0).fill_arrays(init)
if sys.argv[4] == "1":
    from forge_design.evaluate.runner_axismach import design_chain, load_problem
    p = load_problem(sys.argv[5])
    p.geometry["initial_line_run"] = sys.argv[6]
    d = design_chain(p)
    out.update(wall_inv=d["wall_inv"], c=d["wall"]._spl.c, t=d["wall"]._spl.t)
np.savez(sys.argv[3], **out)
'''
py = sys.executable
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    try:
        subprocess.run(f"git -C {ROOT} archive {BASE} design/forge_design solver_density_cuda/data | tar -x -C {td}", shell=True, check=True,
                       capture_output=True)
        have_base = (td / "design/forge_design/geometry/moc_inverse.py").exists()
    except subprocess.CalledProcessError as e:
        have_base = False
        SKIPPED.append(f"8 既定のビット同一 (git archive に失敗: {e.stderr[:120]!r})")
    if have_base:
        (td / "snip.py").write_text(_SNIP)
        args = [str(PROB_PATH) for PROB_PATH in (C45 / PROB45,)] + [str(RUNS / "run_0062_euler_wallfit_fit_r1_ext6k")]
        res_files = {}
        for lab, pkg in (("base", td / "design"), ("new", ROOT / "design")):
            f = td / f"{lab}.npz"
            r = subprocess.run([py, str(td / "snip.py"), str(pkg), str(Path(__file__).resolve().parent), str(f),
                                "1" if HAVE45 else "0", *args], capture_output=True, text=True)
            if r.returncode != 0:
                check(f"ビット同一: {lab} の計算が失敗 ({r.stderr[-300:]})", False)
            res_files[lab] = f
        if all(f.exists() for f in res_files.values()):
            a, b = np.load(res_files["base"]), np.load(res_files["new"])
            same = lambda k: a[k].shape == b[k].shape and np.array_equal(a[k], b[k], equal_nan=True)  # noqa: E731
            check(f"ビット同一 (a): 放射源流の網 (fill_levels {a['lev'].shape}) が実装前と完全一致", same("lev"))
            check(f"ビット同一 (a): 放射源流の点群 (fill_arrays {a['pts'].shape}) が実装前と完全一致", same("pts"))
            if HAVE45:
                check(f"ビット同一 (b): case/45 単調壁の MOC 点群 ({a['wall_inv'].shape}) が実装前と完全一致", same("wall_inv"))
                check(f"ビット同一 (b): 当てはめ後の係数 ({a['c'].shape}) とノット ({a['t'].shape}) が完全一致",
                      same("c") and same("t"))
            else:
                SKIPPED.append(f"8 (b) case/45 のビット同一 (凍結源 {RUNS} が無い)")

# --- 9. 計算準備の入口でゲートを必須にする (2026-10-07 result 段レビュー M1) ---------------------------------------
import types as _types
import tempfile as _tempfile
from forge_design.evaluate import runner_axismach as _RA
_PC = _types.SimpleNamespace(geometry={"moc_corrector": "converge"})
_PF = _types.SimpleNamespace(geometry={})
_ok = {"moc": {"gate": {"applicable": True, "pass": True, "reasons": []}}}
_ng = {"moc": {"gate": {"applicable": True, "pass": False, "reasons": ["反復の失敗 1 対 (非有限 0・上限到達 1)"]}}}
try:
    _RA.require_moc_gate(_PC, _ok)
    check("入口のゲート: converge で合格なら通す", True)
except ValueError as e:
    check(f"入口のゲート: converge で合格なら通す ({e})", False)
rejects("入口のゲート: converge で不合格 (反復の失敗 1 対)", lambda: _RA.require_moc_gate(_PC, _ng), must="不合格")
rejects("入口のゲート: converge で診断が無い", lambda: _RA.require_moc_gate(_PC, {"moc": None}), must="診断が無い")
rejects("入口のゲート: converge で gate.applicable が偽", lambda: _RA.require_moc_gate(_PC, {"moc": {"gate": {"applicable": False, "pass": None}}}), must="診断が無い")
try:
    _RA.require_moc_gate(_PF, {"moc": {"gate": {"applicable": False, "pass": None, "reasons": []}}})
    check("入口のゲート: fixed2 (合否を出さない) は通す", True)
except ValueError as e:
    check(f"入口のゲート: fixed2 は通す ({e})", False)
# prepare (Euler) は design_chain の不合格で run dir を作らずに止まる (design_chain を差し替えて確かめる)
_PROB = C45 / "problem_d155_euler_t0cluster_u5em3.yaml"
if _PROB.is_file() and (C45 / "problem_d155_ns_n012_N2.yaml").is_file():
    _orig = _RA.design_chain
    try:
        _RA.design_chain = lambda p: _ng
        with _tempfile.TemporaryDirectory() as _td:
            _rd = Path(_td) / "run_x"
            try:
                _RA.prepare(_PROB, _rd, nsteps=10)
                check("prepare: MOC のゲート不合格で止まる (通ってしまった)", False)
            except ValueError as e:
                check(f"prepare: MOC のゲート不合格で止まり run dir を作らない ({str(e)[:60]})", not _rd.exists())
            _rn = Path(_td) / "run_y"
            try:
                _RA.prepare_ns(C45 / "problem_d155_ns_n012_N2.yaml", _rn, nsteps=10)
                check("prepare_ns: MOC のゲート不合格で止まる (通ってしまった)", False)
            except ValueError as e:
                check(f"prepare_ns: MOC のゲート不合格で止まり run dir を作らない ({str(e)[:60]})", not _rn.exists() and "MOC" in str(e))
    finally:
        _RA.design_chain = _orig
else:
    SKIPPED.append("9 prepare の入口のゲート (case/45 の converge の Euler 問題が無い)")

if not HAVE45:
    SKIPPED.append("case/45 依存の項目 (凍結源または問題 YAML が無い)")
for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
