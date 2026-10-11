#!/usr/bin/env python3
"""③ベル (D) 直接モードの C2 壁 `geometry/wall_bell_c2.py` の単体試験 (plan tooling-nozzle-bell-two-mode §5.1 #3)。

合格条件 (§5.1 #3): 拘束の残差 ≤ 1e-10、継ぎ目なしの C2 (スロートで収縮部と r・r'・r'' が一致)、同一入力でビット同一、
成立条件の負例試験 (r' < 0・θ' の符号変化 ≠ 1・|κ| > 1/R_t がそれぞれ INFEASIBLE になる)、収縮部の単調条件。
加えて: 拘束行列の階数 7・拘束の零空間で ∫(r''')² が正定値・KKT の解が最小 (零空間の摂動で J が増える)、
規約 (等間隔 16 区間の clamped 5 次・係数 21)、設計変数の範囲外・非有限は ValueError。

実行: design/.venv-opt の Python で `python design/tests/run_bell_c2_wall_tests.py`。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.wall_bell_c2 import (DEGREE, N_COEF, N_SPAN, RESID_TOL, BellC2Spec,  # noqa: E402
                                                BellC2Wall, cone15_length)

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc as e:
        return str(e) or type(e).__name__
    return None


LC = cone15_length(50.0)
LK = [f * LC for f in (0.6, 0.7, 0.8, 0.9, 1.0)]
GRID = [(L, te, tm, xf) for L in LK for tm in (25.0, 32.5, 40.0) for te in (5.0, 10.0, 15.0) for xf in (0.1, 0.2, 0.3)]

# ---------------------------------------------------------------- 寸法と規約
check("15° 円錐の長さ (ε 50) = 22.6575 r_t (plan §8-1 の 3.172 m @ r_t 0.14 m)", abs(LC - 22.657533530077966) < 1e-12,
      f"{LC:.10f} r_t = {LC * 0.14:.4f} m")
w = BellC2Wall(LK[1], 10.0, 32.5, 0.1)
check("ノット: 等間隔 16 区間の clamped 5 次 (端 6 重)",
      len(w.knots) == N_SPAN + 1 + 2 * DEGREE and np.all(w.knots[:DEGREE + 1] == 0.0)
      and np.all(w.knots[-DEGREE - 1:] == w.L) and np.allclose(np.diff(w.knots[DEGREE:-DEGREE]), w.L / N_SPAN, rtol=0, atol=1e-13))
check("係数 21 個", len(w.coef) == N_COEF == 21)
sp = BellC2Spec()
check("寸法 spec: ε 50・R_t 1・収縮比 3 (r_U = √3)・L_U 3・L_pipe 1",
      (sp.eps, sp.R_t, sp.contraction_ratio, sp.L_U, sp.L_pipe) == (50.0, 1.0, 3.0, 3.0, 1.0) and abs(sp.r_U - np.sqrt(3.0)) < 1e-15)

# ---------------------------------------------------------------- 拘束の残差・C2・階数・最小性 (V0 の格子 135 点すべて)
res_max = 0.0
jump_max = [0.0, 0.0, 0.0]
rank_ok = True
eig_min_rel = np.inf
for args in GRID:
    wg = BellC2Wall(*args)
    res_max = max(res_max, float(np.max(np.abs(wg.constraint_residuals()))))
    tj = wg.throat_jump()
    jump_max = [max(jump_max[i], abs(tj[i])) for i in range(3)]
    rank_ok &= wg.constraint_rank == 7
    eig_min_rel = min(eig_min_rel, wg.reduced_hessian_eig[0] / wg.reduced_hessian_eig[1])
check(f"拘束の残差 ≤ {RESID_TOL:g} (V0 の格子 135 点の最大)", res_max <= RESID_TOL, f"最大 {res_max:.2e}")
check("C2: スロートで収縮部と拡大部の r・r'・r'' の差 ≤ 1e-10 (135 点の最大)", max(jump_max) <= 1e-10,
      f"r {jump_max[0]:.1e}・r' {jump_max[1]:.1e}・r'' {jump_max[2]:.1e}")
check("拘束行列の階数 = 7 (135 点すべて)", rank_ok)
check("拘束の零空間で ∫(r''')² が正定値 (最小固有値 / 最大 > 1e-12、135 点の最小)", eig_min_rel > 1e-12, f"{eig_min_rel:.2e}")

# 最小性: 拘束を保つ摂動 (零空間) で J が増える
rng = np.random.default_rng(20261011)
_, _, Vt = np.linalg.svd(w.C)
N = Vt[7:].T
J0 = float(w.coef @ w.Q @ w.coef)
incr = []
for _ in range(20):
    z = rng.standard_normal(N.shape[1])
    c = w.coef + 1e-3 * (N @ z)
    incr.append(float(c @ w.Q @ c) - J0)
check("KKT の解が最小: 拘束を保つ摂動 20 方向で ∫(r''')² が増える", min(incr) > 0.0, f"最小の増分 {min(incr):.2e}")
check("拘束を保つ摂動で拘束の値は変わらない", float(np.max(np.abs(w.C @ (N @ rng.standard_normal(N.shape[1]))))) < 1e-12)

# 拡大部の内部は C4 (ノットで r'''' まで連続)
xk = w.knots[DEGREE + 1:-DEGREE - 1]
dj = max(abs(float(w.expansion(xx + 1e-12, nu=4)) - float(w.expansion(xx - 1e-12, nu=4))) for xx in xk)
check("拡大部の内部ノットで r'''' が連続 (単純ノットの 5 次 = C4)", dj < 1e-6, f"最大の差 {dj:.1e}")

# ---------------------------------------------------------------- 決定性 (ビット同一)
a = BellC2Wall(LK[2], 5.0, 32.5, 0.1)
b = BellC2Wall(LK[2], 5.0, 32.5, 0.1)
xs = np.linspace(-a.spec.L_U - a.spec.L_pipe, a.L, 5001)
check("同一入力でビット同一 (係数・ノット・標本の r・r'・r'')",
      a.coef.tobytes() == b.coef.tobytes() and a.knots.tobytes() == b.knots.tobytes()
      and all(np.asarray(a.r(xs, d)).tobytes() == np.asarray(b.r(xs, d)).tobytes() for d in (0, 1, 2)))
fa, fb = a.feasibility(), b.feasibility()
check("同一入力で成立判定も同一", fa == fb)

# ---------------------------------------------------------------- 成立条件: 正例と負例 (理由ごと)
pos = BellC2Wall(LK[1], 10.0, 32.5, 0.1).feasibility()
check("正例: L 70 %・θ_e 10°・θ_max 32.5°・x_m/L 0.1 は成立", pos["ok"] and not pos["violations"],
      f"r'最小 {pos['rp_min']:.2e}・符号変化 {pos['theta_sign_changes']}・|κ|R_t 最大 {pos['kappa_Rt_max']:.4f}")
check("正例: 符号変化の位置 = x_m (標本 1 点の幅で)", abs(pos["x_sign_change"][0] - 0.1 * LK[1]) <= LK[1] / 1999 * 1.01,
      f"{pos['x_sign_change']} vs x_m {0.1 * LK[1]:.4f}")
neg_rp = BellC2Wall(LK[4], 5.0, 25.0, 0.3).feasibility()
check("負例 r' < 0: L 100 %・θ_e 5°・θ_max 25°・x_m/L 0.3 (出口の手前で √ε を超えて戻る)",
      (not neg_rp["ok"]) and any(v.startswith("r' < 0") for v in neg_rp["violations"]), "; ".join(neg_rp["violations"]))
neg_sc = BellC2Wall(LK[0], 5.0, 25.0, 0.1).feasibility()
check("負例 θ' の符号変化 ≠ 1: L 60 %・θ_e 5°・θ_max 25°・x_m/L 0.1 (x_m の後で再び増える)",
      (not neg_sc["ok"]) and any("符号変化" in v for v in neg_sc["violations"])
      and not any(v.startswith("r' < 0") or "|κ|" in v for v in neg_sc["violations"]), "; ".join(neg_sc["violations"]))
neg_k = BellC2Wall(LK[0], 5.0, 40.0, 0.1).feasibility()
check("負例 |κ| > 1/R_t: L 60 %・θ_e 5°・θ_max 40°・x_m/L 0.1 (スロート直後で曲率が上限を超える)",
      (not neg_k["ok"]) and any("|κ|" in v for v in neg_k["violations"]), "; ".join(neg_k["violations"]))
# 標本の境界: 2000 点 (§5.1 #3)
check("成立判定の標本は 2000 点", pos["n_sample"] == 2000)

# ---------------------------------------------------------------- 収縮部
up = w.up
check("収縮部: U→T Hermite の単調条件 μ = L_U²/[R_t(r_U−1)] ≤ 20", up.mu <= 20.0 and not up.validate(),
      f"μ = {up.mu:.4f} (上限 L_U {np.sqrt(20 * (sp.r_U - 1.0)):.4f})")
xc = np.linspace(-sp.L_U, 0.0, 4001)
check("収縮部: 標本で r' ≤ 0 (単調収縮)・r ≥ 1", np.all(w.r(xc, 1) <= 1e-14) and np.all(w.r(xc) >= 1.0 - 1e-15))
xp = np.linspace(-sp.L_U - sp.L_pipe, -sp.L_U, 101)
check("入口直管: r = √3・r' = 0", np.all(w.r(xp) == sp.r_U) and np.all(w.r(xp[:-1], 1) == 0.0))
check("U 端 (x = −L_U) で直管と C2 (r = √3・r' = 0・r'' = 0)",
      abs(float(w.r(np.array([-sp.L_U]))[0]) - sp.r_U) < 1e-14 and abs(float(w.r(np.array([-sp.L_U]), 1)[0])) < 1e-14
      and abs(float(w.r(np.array([-sp.L_U]), 2)[0])) < 1e-12)
check("全体の最小半径がスロート (成立した壁)", float(np.min(w.r(np.linspace(w.x_in, w.L, 20001)))) >= 1.0 - 1e-12)

# ---------------------------------------------------------------- 入力の検査 (ValueError = 入力の誤り)
for name, fn in [("θ_max 19.9°", lambda: BellC2Wall(LK[0], 10.0, 19.9, 0.2)),
                 ("θ_max 45.1°", lambda: BellC2Wall(LK[0], 10.0, 45.1, 0.2)),
                 ("θ_e 1.9°", lambda: BellC2Wall(LK[0], 1.9, 30.0, 0.2)),
                 ("θ_e 20.1°", lambda: BellC2Wall(LK[0], 20.1, 30.0, 0.2)),
                 ("x_m/L 0.049", lambda: BellC2Wall(LK[0], 10.0, 30.0, 0.049)),
                 ("x_m/L 0.401", lambda: BellC2Wall(LK[0], 10.0, 30.0, 0.401)),
                 ("L = 0", lambda: BellC2Wall(0.0, 10.0, 30.0, 0.2)),
                 ("L = NaN", lambda: BellC2Wall(float("nan"), 10.0, 30.0, 0.2)),
                 ("θ_e = inf", lambda: BellC2Wall(LK[0], float("inf"), 30.0, 0.2)),
                 ("spec の ε ≤ 1", lambda: BellC2Wall(LK[0], 10.0, 30.0, 0.2, BellC2Spec(eps=1.0))),
                 ("r(x) の範囲外 (x > L)", lambda: w.r(np.array([w.L + 1e-9]))),
                 ("r(x) の範囲外 (x < 入口)", lambda: w.r(np.array([w.x_in - 1e-9])))]:
    msg = raises(fn)
    check(f"入力の検査: {name} は ValueError", msg is not None, msg or "例外なし")
check("範囲の端 (θ_max 20/45°・θ_e 2/20°・x_m/L 0.05/0.40) は受ける",
      raises(lambda: [BellC2Wall(LK[0], te, tm, xf) for te in (2.0, 20.0) for tm in (20.0, 45.0) for xf in (0.05, 0.40)]) is None)

# ---------------------------------------------------------------- V0 の格子の成立数 (記録)
nfeas = sum(BellC2Wall(*a_).feasibility()["ok"] for a_ in GRID)
print(f"info V0 の格子 (θ_max 3 × θ_e 3 × x_m/L 3 × L_k 5 = 135 点) の成立数: {nfeas}")

print("\nALL PASS" if FAIL == 0 else f"\n{FAIL} FAIL")
sys.exit(0 if FAIL == 0 else 1)
