#!/usr/bin/env python3
r"""放射源流 (厳密解) での軸上の極限 θ_r と予測修正の収束の判別 A/B (plan §4.0 / §6 V0・V2)。CFD 0 step。

計画: plans/active/discretization-moc-axis-limit-and-corrector.md §4.0・§6.0 V0・V2。
設定は `run_inverse_tests.py` §8 (c) と同じ放射源流 (γ 1.4、半角 10°、初期線 x = 1.2 の縦線、軸 5.5 → 1.2、
壁 = C⁺ 線上の流束閉包 `cplus_flux_wall`、壁の誤差は x = 1.4〜2.4 の 21 標本の最大相対誤差)。
軸上の厳密解は θ_r = 1/x (A/A* = x² から θ_r = −½ d ln(ρu)/dx = 1/x)。

腕 (違いは単位過程だけ):
  A  = legacy + converge    (軸端点の sinθ/r は相手の値で代用、修正子は収束まで)
  B  = analytic + converge  (軸端点はその点の θ_r、修正子は収束まで)
  V2 = legacy + fixed2      (現行: 予測 1 + 修正 2 回)
測る量 (解像度 n_axis × n_start ごと):
  - 第 1 段 (軸端点 2 個の対から作る点 L_1[i]) の θ・ν の厳密解からの誤差 (生成点の座標で比較)
  - 再出発の差: 網の列 (軸節点から上る C⁻) の対 (L_1[i], L_0[i]) を単位過程に通し直したときの L_1[i] との差
    (全列の最大) と、列 x ≈ 1.6 を初期線にして充填し直したときの網の節点の差 (同じ C⁺ を持つ節点)
  - 壁の誤差と収束次数、対の 5 分類・反復回数・最終残差・ゲート
判定 (§4.0 の事前登録をそのまま): B の第 1 段の θ の誤差が各解像度で A の 1/10 以下、壁の誤差も A 以下、
最細区間の壁の収束次数 ≥ 1.7 → 「軸端の源項が第 1 段の誤差の主因」を支持。V2 (§6 V2) は B の壁の誤差が
現行 (V2) 以下かつ最細区間の次数が現行以上。

usage: python3 design/tests/moc_axis_limit_radial.py [OUT_JSON] [--levels 140x25,280x49,...]
       既定の出力: case/45.isobutane_m6_d155/_band_ab/moc_axis_limit_v0_radial.json (plan の他の V の記録と同じ置き場)
       `run_moc_axis_limit_tests.py` は粗い 2 解像度で同じ関数を呼ぶ (判定の向きだけを試験に固定)。
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate.ic import invert_area_ratio  # noqa: E402
from forge_design.geometry import moc_kernel as MK  # noqa: E402
from forge_design.geometry.moc_inverse import (InverseMOC, _flux_along,  # noqa: E402
                                               axis_theta_r_init, cplus_flux_wall, moc_gate)
from forge_design.geometry.moc_kernel import _Pt, interior_vec, pm_nu  # noqa: E402

G = 1.4
THC = np.deg2rad(10.0)
X1 = 1.2
R_TOP = X1 * np.tan(THC)
X_AX0 = 5.5
XQ = np.linspace(1.4, 2.4, 21)                 # 壁の誤差の標本 (run_inverse_tests §8 (c) と同じ)
LEVELS = ((140, 25), (280, 49), (560, 97), (1120, 193))
ARMS = {"A": ("legacy", "converge"), "B": ("analytic", "converge"), "V2": ("legacy", "fixed2")}
X_RESTART = 1.6                                # 再出発に使う列の軸節点の x
OUT_DEFAULT = ROOT / "case/45.isobutane_m6_d155/_band_ab/moc_axis_limit_v0_radial.json"


def src_M(x, r):
    R = max(np.hypot(x, r), 1.0001)
    return float(invert_area_ratio(np.array([R * R]), np.array([True]), G)[0])


def src_dMdx_axis(x):
    """軸上の厳密な dM/dx (A/A* = x² と面積 – マッハ関係: d ln A/dM = (M²−1)/[M(1+(γ−1)M²/2)])。"""
    M = src_M(x, 0.0)
    return (2.0 / x) * M * (1.0 + 0.5 * (G - 1.0) * M * M) / (M * M - 1.0)


def build_init(n_ax: int, n_st: int):
    init = ([_Pt(float(x), 0.0, 0.0, float(pm_nu(src_M(x, 0.0), G)), G)
             for x in np.linspace(X_AX0, X1, n_ax + 1)[:-1]]
            + [_Pt(X1, float(r), float(np.arctan2(r, X1)), float(pm_nu(src_M(X1, r), G)), G)
               for r in np.linspace(0.0, R_TOP, n_st)])
    return init


def _cum0(init, n_ax):
    c = np.zeros(len(init))
    c[n_ax:] = _flux_along(init[n_ax:], G)
    return c


def first_level_error(lev, n_ax):
    """第 1 段の軸端点 2 個の対 (i = 0..n_ax−1: A = L_0[i+1]、B = L_0[i]、どちらも r = 0) から作った点の誤差。"""
    P = lev[1, :n_ax]
    ok = np.isfinite(P[:, 0])
    x, r, th, nu = P[ok, 0], P[ok, 1], P[ok, 2], P[ok, 3]
    th_ex = np.arctan2(r, x)
    nu_ex = np.array([float(pm_nu(src_M(a, b), G)) for a, b in zip(x, r)])
    eth, enu = th - th_ex, nu - nu_ex
    j = int(np.argmax(np.abs(eth)))
    return {"n": int(ok.sum()), "theta_abs_max": float(np.abs(eth).max()), "theta_rel_max": float(np.abs(eth / th_ex).max()),
            "nu_abs_max": float(np.abs(enu).max()), "at_max_theta": {"x": float(x[j]), "r": float(r[j]),
                                                                       "theta_rel": float(eth[j] / th_ex[j])}}


def column(lev, i):
    col = lev[:, i]
    ok = np.isfinite(col[:, 0]) & np.isfinite(col[:, 1])
    cut = int(np.argmin(ok)) if not ok.all() else len(ok)
    return np.array(col[:cut])


def restart_diff(inv, lev, init, thr, n_ax):
    """再出発の差 (V3 と同じ考え方の放射源流版)。

    (a) 全列 i で対 (A = L_1[i]、B = L_0[i]) を同じ単位過程に通し直し、L_1[i] との差 (θ・ν の最大)。
        列 = C⁻ なので、自己整合なら差は 0 (修正子の許容差まで)。
    (b) 軸節点 x ≈ X_RESTART の列 i を壁より下まで初期線にして充填し直し、列 < i の網の節点 (同じ C⁺ を持つもの)
        が元と一致するか (θ・ν・x・r の最大差)。"""
    Ls = []
    for i in range(1, n_ax):
        if np.isfinite(lev[1, i, 0]):
            Ls.append(i)
    Ls = np.array(Ls)
    A, B = lev[1, Ls], lev[0, Ls]
    kw = {}
    if inv.axis_limit == "analytic":
        kw = dict(axis_limit="analytic", thrA=np.full(len(Ls), np.nan), thrB=thr[Ls])
    if inv.corrector == "converge":
        kw.update(corrector="converge", tol=inv.tol, max_corr=inv.max_corr)
    q = interior_vec(*(A[:, c] for c in range(5)), *(B[:, c] for c in range(5)), G, 1.0, inv.n_corr, **kw)
    da = {"n_columns": int(len(Ls)), "dtheta_max": float(np.abs(q[2] - A[:, 2]).max()),
          "dnu_max": float(np.abs(q[3] - A[:, 3]).max()), "dx_max": float(np.abs(q[0] - A[:, 0]).max()),
          "dr_max": float(np.abs(q[1] - A[:, 1]).max())}
    # (b)
    i = int(np.argmin(np.abs(np.array([p.x for p in init[:n_ax]]) - X_RESTART)))
    col = column(lev, i)
    col = col[col[:, 1] <= col[:, 0] * np.tan(THC)]          # 壁 (レイ) より下
    new = init[:i] + [_Pt(float(a[0]), float(a[1]), float(a[2]), float(a[3]), G, float(a[4])) for a in col]
    inv2 = InverseMOC(gamma=G, delta=1.0, axis_limit=inv.axis_limit, corrector=inv.corrector)
    thr2 = None
    if inv.axis_limit == "analytic":
        thr2 = np.full(len(new), np.nan)
        thr2[:i + 1] = thr[:i + 1]                            # 列の軸端 = 元の軸節点 L_0[i] (同じ θ_r)
    lev2 = inv2.fill_levels(new, axis_thr=thr2)
    n2 = min(lev2.shape[0], lev.shape[0])
    m_max = i + len(col) - 1
    same_cplus = (np.arange(n2)[:, None] + np.arange(i)[None, :]) <= m_max
    db = {"column": i, "axis_x": float(col[0, 0]), "n_line": int(len(col))}
    for c, name in ((0, "dx_max"), (1, "dr_max"), (2, "dtheta_max"), (3, "dnu_max")):
        a, b = lev2[:n2, :i, c], lev[:n2, :i, c]
        both = np.isfinite(a) & np.isfinite(b) & same_cplus
        db[name] = float(np.abs(a[both] - b[both]).max()) if both.any() else float("nan")
        if c == 0:
            db["n_nodes_compared"] = int(both.sum())
    return {"reproc_first_row": da, "refill_column": db}


def run_arm(arm: str, n_ax: int, n_st: int, with_restart: bool = True) -> dict:
    axis_limit, corrector = ARMS[arm]
    init = build_init(n_ax, n_st)
    inv = InverseMOC(gamma=G, delta=1.0, axis_limit=axis_limit, corrector=corrector)
    thr, info = axis_theta_r_init(init, n_ax, lambda x: src_M(x, 0.0), G, target_dM=src_dMdx_axis,
                                  M_line_axis=src_M(X1, 0.0))
    t0 = time.time()
    lev = inv.fill_levels(init, axis_thr=thr if axis_limit == "analytic" else None)
    t_fill = time.time() - t0
    md = float(_flux_along(init[n_ax:], G)[-1])
    w = cplus_flux_wall(lev, _cum0(init, n_ax), md, G)
    ex = XQ * np.tan(THC)
    e_wall = float(np.max(np.abs(np.interp(XQ, w[:, 0], w[:, 1]) / ex - 1.0)))
    diag = dict(inv.last_diag)
    gate = moc_gate(diag, info if axis_limit == "analytic" else None, w)
    diag.pop("_geom_all", None)
    xs = np.array([p.x for p in init[:n_ax]] + [X1])
    thr_err = float(np.max(np.abs(np.r_[thr[:n_ax], thr[n_ax]] * xs - 1.0)))      # θ_r = 1/x との相対差
    out = {"arm": arm, "axis_limit": axis_limit, "corrector": corrector, "n_axis": n_ax, "n_start": n_st,
           "wall_err_max_rel": e_wall, "first_level": first_level_error(lev, n_ax),
           "theta_r_rel_err_vs_1_over_x": thr_err, "fill_seconds": t_fill,
           "pairs": diag["pairs"], "iters": {k: diag["iters"][k] for k in ("max", "mean", "n_pairs")},
           "resid": diag["resid"], "branch": diag["branch"], "gate": gate}
    if with_restart:
        out["restart"] = restart_diff(inv, lev, init, thr, n_ax)
    return out


def order(e_coarse, e_fine):
    return float(np.log2(e_coarse / e_fine))


def judge(res: dict, levels) -> dict:
    """§4.0 の判定 (事前登録の文言どおり) と §6 V2 の判定。res[arm][level_key]。"""
    keys = [f"{a}x{b}" for a, b in levels]
    A, B, V2 = res["A"], res["B"], res["V2"]
    ratio = {k: B[k]["first_level"]["theta_abs_max"] / A[k]["first_level"]["theta_abs_max"] for k in keys}
    c1 = all(v <= 0.1 for v in ratio.values())
    c2 = all(B[k]["wall_err_max_rel"] <= A[k]["wall_err_max_rel"] for k in keys)
    p = {arm: [order(res[arm][keys[j]]["wall_err_max_rel"], res[arm][keys[j + 1]]["wall_err_max_rel"])
               for j in range(len(keys) - 1)] for arm in ("A", "B", "V2")}
    c3 = p["B"][-1] >= 1.7
    unconv = {arm: sum(res[arm][k]["pairs"]["iter_nonfinite"] + res[arm][k]["pairs"]["iter_maxiter"] for k in keys)
              for arm in ("A", "B")}
    c0 = unconv["A"] == 0 and unconv["B"] == 0
    if not c0:
        v0 = "判定不能 (未収束の対あり — §4.0 の必須条件を満たさない)"
    elif c1 and c2 and c3:
        v0 = "支持 (軸端の源項が第 1 段の誤差の主因。実装を進める)"
    elif c1:
        v0 = "第 1 段だけ改善 (壁の条件を満たさない — 前提を棄却しユーザに報告)"
    else:
        v0 = "第 1 段も改善しない (適用箇所・微分・単位を点検)"
    v2_ok = all(B[k]["wall_err_max_rel"] <= V2[k]["wall_err_max_rel"] for k in keys) and p["B"][-1] >= p["V2"][-1]
    return {"V0": {"verdict": v0, "unconverged_pairs": unconv, "first_level_theta_ratio_B_over_A": ratio,
                   "c_first_level_le_0p1": c1, "c_wall_B_le_A": c2, "c_order_B_last_ge_1p7": c3,
                   "orders": p},
            "V2": {"verdict": "PASS" if v2_ok else "FAIL",
                   "wall_B_le_V2_all": all(B[k]["wall_err_max_rel"] <= V2[k]["wall_err_max_rel"] for k in keys),
                   "order_B_last": p["B"][-1], "order_V2_last": p["V2"][-1]}}


def main(argv):
    out = OUT_DEFAULT
    levels = LEVELS
    args = list(argv)
    if "--levels" in args:
        j = args.index("--levels")
        levels = tuple(tuple(int(v) for v in s.split("x")) for s in args[j + 1].split(","))
        del args[j:j + 2]
    if args:
        out = Path(args[0])
    t0 = time.time()
    res = {arm: {} for arm in ARMS}
    for (n_ax, n_st) in levels:
        for arm in ARMS:
            r = run_arm(arm, n_ax, n_st)
            res[arm][f"{n_ax}x{n_st}"] = r
            fl = r["first_level"]
            print(f"{arm:2s} {n_ax:5d}x{n_st:<4d} wall {r['wall_err_max_rel']:.4e}  第1段 θ {fl['theta_abs_max']:.3e} "
                  f"(rel {fl['theta_rel_max']:.3e})  ν {fl['nu_abs_max']:.3e}  反復 max {r['iters']['max']} "
                  f"mean {r['iters']['mean']:.2f}  失敗 {r['pairs']['iter_nonfinite'] + r['pairs']['iter_maxiter']}  "
                  f"再出発 (a) dθ {r['restart']['reproc_first_row']['dtheta_max']:.2e} "
                  f"(b) dθ {r['restart']['refill_column']['dtheta_max']:.2e}  ({r['fill_seconds']:.1f} s)", flush=True)
    jd = judge(res, levels)
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                           capture_output=True, text=True).stdout.strip()
    doc = {"plan": "plans/active/discretization-moc-axis-limit-and-corrector.md §4.0 (V0)・§6 V2",
           "commit": commit, "design_tree_dirty": dirty,
           "setup": {"gamma": G, "half_angle_deg": 10.0, "x_line": X1, "x_axis_start": X_AX0, "xq": XQ.tolist(),
                     "wall": "cplus_flux_wall", "theta_r_exact": "1/x",
                     "tol": MK.CORR_TOL, "max_corr": MK.CORR_MAX, "resid_tol": MK.RESID_TOL,
                     "x_restart_column": X_RESTART},
           "levels": [list(l) for l in levels], "results": res, "judge": jd, "elapsed_s": time.time() - t0}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1, ensure_ascii=False, default=float))
    print(f"V0: {jd['V0']['verdict']}  比 B/A (第 1 段 θ) "
          + " ".join(f"{k}:{v:.2e}" for k, v in jd["V0"]["first_level_theta_ratio_B_over_A"].items())
          + f"  次数 B {['%.3f' % v for v in jd['V0']['orders']['B']]} A {['%.3f' % v for v in jd['V0']['orders']['A']]}"
          + f" V2 {['%.3f' % v for v in jd['V0']['orders']['V2']]}")
    print(f"V2: {jd['V2']['verdict']}  -> {out} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(sys.argv[1:])
