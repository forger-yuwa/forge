"""ユーザ提案 (2026-10-07「スロートより上流に排除厚さを足しこむ処理をやめにしませんか」) の試算 (CFD 0 step): スロートより上流に δ_r を足さず、配管〜スロートを 5 次多項式 1 本にしたら、
今の物理壁 (run_0147) とどれだけ違うか。多項式の端条件: 配管端 x=−L_U で (r_U, 0, 0)、スロート x=0 で
下流の物理壁 (S + δ_r) の (r, r′, r″) に 2 階微分まで連続。
あわせて「スロートを動かさず下流で徐変」(δ_r を 0 から全量へ) の有効コア形状のずれ (s−1)·δ_r を見積もる。
usage: [DESIGN_DIR=<design/ の写し>] [CASE_RUNS=<run_0147・run_0062 のある case dir>] python3 upstream_poly_probe.py > _band_ab/upstream_poly_probe.txt
2026-10-07 の実測は HEAD 934d3086 の design/ の写しで回した (同じツリーで実装が進行中のため)。"""
import os
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("DESIGN_DIR", str(C.parents[1] / "design")))
from forge_design.evaluate.runner_axismach import design_chain, load_problem, delta_r_from_table, _gam_or_gas  # noqa: E402
from forge_design.geometry.wall_axismach import PhysicalNozzleWall  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
p = load_problem(C / "problem_d155_ns_finemesh_recal_final_mono.yaml")
p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
d = design_chain(p)
scale = float(p.spec["r_throat"]); um = scale * 1e6
tbl = np.loadtxt(RUNS / "run_0147_ns_mono_final/delta_r_initial.csv", delimiter=",", skiprows=1)
f = delta_r_from_table(tbl[:, 0], tbl[:, 1])
W = PhysicalNozzleWall(d["wall"], d["wall_inv"], scale, float(p.spec["Pt"]), float(p.spec["Tt"]),
                       _gam_or_gas(p), p.cp, offset="radial", delta_r_x=f, ramp=tuple(float(v) for v in p.geometry["pw_ramp"]), upstream="ramp")   # 今の物理壁 (ランプ); 2026-10-07 に既定が poly
dw = d["wall"]; LU = float(dw.up.L_U); rU = float(dw.up.r_U)
x0 = 0.0
e0 = [float(W.r(np.r_[x0 + 1e-12], n)[0]) for n in range(3)]   # 下流側の値 (x=0⁺)
# 5 次多項式 Q(x) = Σ a_k (x+L_U)^k、端条件 6 個
xa, xb = -LU, x0
A = np.zeros((6, 6)); b = np.zeros(6)
def row(xx, n):
    t = xx - xa
    return [0 if k < n else np.prod(range(k - n + 1, k + 1)) * t ** (k - n) for k in range(6)]
for i, (xx, n, v) in enumerate([(xa, 0, rU), (xa, 1, 0.0), (xa, 2, 0.0), (xb, 0, e0[0]), (xb, 1, e0[1]), (xb, 2, e0[2])]):
    A[i] = row(xx, n); b[i] = v
a = np.linalg.solve(A, b)
Q = lambda xx, n=0: sum(a[k] * (np.prod(range(k - n + 1, k + 1)) if n else 1) * (xx - xa) ** (k - n) for k in range(n, 6))  # noqa: E731
xs = np.linspace(-LU, 0, 24001)
dq = (Q(xs) - W.r(xs)) * um
print(f"L_U {LU}  r_U {rU}  throat end (r, r', r'') = {e0}")
print(f"Q − 今の物理壁 [µm]: min {dq.min():.1f} at x {xs[np.argmin(dq)]:.3f}, max {dq.max():.1f} at x {xs[np.argmax(dq)]:.3f}")
for xx in (-11, -9, -6, -3, -1, -0.5, -0.2, -0.05):
    print(f"  x {xx:6.2f}: Q−W {(Q(np.r_[xx]) - W.r(np.r_[xx]))[0] * um:9.2f} µm   δ_r(x) {f(np.r_[xx])[0] * um:8.1f} µm")
# 物理スロート (r′ = 0) の位置
from scipy.optimize import brentq
xt = brentq(lambda xx: Q(xx, 1), -0.5, 0.0)
print(f"新しい物理スロート x {xt * scale * 1e3:.4f} mm, r {Q(xt) :.7f} r_t (今 {W.x_throat * scale * 1e3:.4f} mm, {W.r_throat:.7f})")
print(f"Q の r′ 最大 {Q(xs, 1).max():.3e} (上流の単調性: r′ ≤ 0 が x<x_t で成り立つか {bool(np.all(Q(xs[xs < xt - 1e-6], 1) <= 0))})")
print(f"r″ の差 (Q″ − W″) 最大 |·| {np.abs(Q(xs, 2) - W.r(xs, 2)).max():.3e} 1/r_t")
# 徐変案: δ_r を x∈[0, L_b] で smoothstep で 0→1。有効コア (物理壁 − δ*) の設計壁からのずれ = (s−1)·δ_r
d0 = float(f(np.r_[0.0])[0])
print(f"δ_r(0) = {d0:.6f} r_t = {d0 * scale * 1e3:.4f} mm、有効スロート面積の変化 2δ_r/r = {2 * d0:.4%}")
for Lb in (1.0, 2.0, 5.0, 10.0, 20.0):
    xx = np.linspace(0, Lb, 4001); s = (xx / Lb) ** 3 * (10 - 15 * xx / Lb + 6 * (xx / Lb) ** 2)
    e = (s - 1) * f(xx)
    print(f"  L_b {Lb:5.1f} r_t: ずれの傾きの最大 {np.degrees(np.abs(np.gradient(e, xx)).max()):.4f}°")
