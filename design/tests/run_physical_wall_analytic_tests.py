#!/usr/bin/env python3
"""joint 壁の物理壁 (`PhysicalNozzleWall` の解析経路) の試験 = plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6b (事前登録)。

壁 = P1(b) の CFD ピン joint 壁 (case/45 n2400、run_0062 res_6000 から凍結した初期線、Md_moc_offset −4.16e-4)、
δ_r = case/45.isobutane_m6_d155/_band_ab/delta_r_c2final_run0051.csv (`delta_r_from_table`)。
物理壁 r_W = r_design + s·δ_r (s: x ≤ −11 で 0、[−11, −6] で 5 次 smoothstep、x ≥ −6 で 1) の合格条件:
  x ≤ −11 で設計壁と |Δr| ≤ 1e-12 / x ≥ −6 で (設計 + δ_r) と ≤ 1e-12 / [−11, −6] で r′ < 0 かつ |r″ − r″_design| ≤ 5e-3 /
  x=0 の r″ 跳び ≤ 1e-8、r″(0⁺) = 0.5 + δ_r″(0) ± 1e-3 / [0, 0.3] の r″ max ≤ 0.52 + max δ_r″ /
  r‴ max ([0, 0.3]) ≤ 設計壁の r‴ max + 0.05 / validate の非単調 0 件 (δ_r 表の範囲内)。
加えて δ_r の導関数 (差分との照合・表の範囲外 0) と入力の拒否。
2026-10-07 (plan tooling-nozzle-upstream-poly-and-throat-sizing §5 の 5): コードの既定が `pw_upstream: poly` になったので、ランプの
振る舞いを試す箇所は `upstream="ramp"` を明示する (poly は design/tests/run_pw_upstream_poly_tests.py)。

usage: design/.venv-opt/bin/python design/tests/run_physical_wall_analytic_tests.py [凍結源の run]
  凍結源の既定: $FORGE_CFDPIN_GRID_RUN か /home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k
  (git 管理外。無ければ exit 2 で止める — 合格扱いにしない)
"""
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.evaluate.runner_axismach import (_gam_or_gas, delta_r_from_table,  # noqa: E402
                                                   design_chain, load_problem)
from forge_design.geometry.wall_axismach import JointFitCFDWall, PhysicalNozzleWall  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


C45 = ROOT / "case/45.isobutane_m6_d155"
run = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "FORGE_CFDPIN_GRID_RUN", "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k"))
if not (run / "res_6000.h5").exists():
    print(f"凍結源の run が無い: {run} (引数か FORGE_CFDPIN_GRID_RUN で与える)")
    sys.exit(2)

# --- δ_r の導関数 ------------------------------------------------------------------
tb = np.loadtxt(C45 / "_band_ab/delta_r_c2final_run0051.csv", delimiter=",", skiprows=1)
f = delta_r_from_table(tb[:, 0], tb[:, 1])
lo_t, hi_t = float(tb[0, 0]), float(tb[-1, 0])
xq = np.linspace(-10.0, 90.0, 2001)
h = 1e-4
for k in (1, 2, 3):
    fd = (f(xq + h, k - 1) - f(xq - h, k - 1)) / (2 * h)
    err = float(np.max(np.abs(f(xq, k) - fd)))
    check(f"δ_r の {k} 階導関数 = 中心差分 (max 差 {err:.1e})", err <= 1e-5 * max(1.0, float(np.max(np.abs(fd)))))
xo = np.array([lo_t - 1.0, hi_t + 1.0])
check("δ_r の範囲外: 値は端値クリップ、導関数は 0",
      np.allclose(f(xo), [tb[0, 1], tb[-1, 1]], rtol=0, atol=1e-15) and all(np.all(f(xo, k) == 0.0) for k in (1, 2, 3)))

# --- P1(b) の壁と物理壁 ------------------------------------------------------------
p = load_problem(C45 / "problem_d155_euler_c2final_n2400.yaml")
p.geometry.update({"initial_line": "cfd", "initial_line_run": str(run), "initial_line_res": "res_6000.h5",
                   "Md_moc_offset": -4.16e-4, "wall_repr": "joint"})
d = design_chain(p)
w = d["wall"]
check("設計壁は JointFitCFDWall", isinstance(w, JointFitCFDWall))
args = (w, d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
PW = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, ramp=(-11.0, -6.0), upstream="ramp")   # case/45 の pw_ramp
check("joint 壁では解析経路が自動で選ばれる", PW.analytic)
check("upstream='ramp' の明示で旧来のランプ (pw_upstream ramp, source explicit)", PW.pw_upstream == "ramp" and PW.pw_upstream_source == "explicit")

xa = np.linspace(w.x_in, -11.0, 30001)
e1 = float(np.max(np.abs(PW.r(xa) - w.r(xa))))
check(f"x ≤ −11 で設計壁と |Δr| = {e1:.1e} ≤ 1e-12", e1 <= 1e-12)
xb = np.linspace(-6.0, w.x_e, 400001)
e2 = float(np.max(np.abs(PW.r(xb) - (w.r(xb) + f(xb)))))
check(f"x ≥ −6 で (設計 + δ_r) と |Δr| = {e2:.1e} ≤ 1e-12", e2 <= 1e-12)
xr = np.linspace(-11.0, -6.0, 30001)
r1max = float(np.max(PW.r(xr, 1)))
d2 = float(np.max(np.abs(PW.r(xr, 2) - w.r(xr, 2))))
check(f"[−11, −6] で r′ < 0 (max r′ {r1max:.3e})", r1max < 0.0)
check(f"[−11, −6] で |r″ − r″_design| = {d2:.2e} ≤ 5e-3", d2 <= 5e-3)
eps = 1e-9
jump = float(PW.r(np.r_[eps], 2)[0] - PW.r(np.r_[-eps], 2)[0])
r2p = float(PW.r(np.r_[0.0], 2)[0])
r2ref = 0.5 + float(f(0.0, 2))
check(f"x=0 の r″ の跳び {jump:.1e} ≤ 1e-8", abs(jump) <= 1e-8)
check(f"r″(0⁺) = {r2p:.6f} vs 0.5 + δ_r″(0) = {r2ref:.6f} (± 1e-3)", abs(r2p - r2ref) <= 1e-3)
xp = np.linspace(0.0, 0.3, 30001)
r2max = float(np.max(PW.r(xp, 2)))
lim2 = 0.52 + float(np.max(f(xp, 2)))
check(f"[0, 0.3] の r″ max {r2max:.5f} ≤ 0.52 + max δ_r″ = {lim2:.5f} (山 max − r″(0⁺) = {r2max - r2p:.4f})", r2max <= lim2)
r3 = float(np.max(np.abs(PW.r(xp, 3))))
r3d = float(np.max(np.abs(w.r(xp, 3))))
check(f"[0, 0.3] の r‴ max {r3:.4f} ≤ 設計壁 {r3d:.4f} + 0.05", r3 <= r3d + 0.05)
msgs = PW.validate()
check(f"validate: 非単調 0 件 (全メッセージ {msgs})", not any("非単調" in m for m in msgs))
xin = np.linspace(PW.x_throat, hi_t, 400001)
xout = np.linspace(hi_t, PW.x_e, 4001)
n_in = int(np.sum(np.diff(PW.r(xin)) < -1e-9))
n_out = int(np.sum(np.diff(PW.r(xout)) < -1e-9))
check(f"スロート下流の非単調: δ_r 表の範囲内 (≤ {hi_t:.3f}) {n_in} 件", n_in == 0)
print(f"info 表の範囲外 [{hi_t:.3f}, x_e={PW.x_e:.3f}] (δ_r は端値一定・導関数 0): 非単調 {n_out} 件")
print(f"info 物理スロート x_t={PW.x_throat:.6f} r_t={PW.r_throat:.6f} κ_t={PW.kappa_throat:.6f}")
check("物理スロートで r′ = 0", abs(float(PW.r(np.r_[PW.x_throat], 1)[0])) <= 1e-12)

# --- 入力の拒否 --------------------------------------------------------------------
for name, kw in (("δ_r 無し", {}), ("導関数を返せない δ_r", {"delta_r_x": lambda x: f(x)})):
    try:
        PhysicalNozzleWall(*args, offset="radial", **kw)
        check(f"解析経路で {name} を拒否", False)
    except ValueError:
        check(f"解析経路で {name} を拒否", True)
PWo = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, analytic=False)
check("analytic=False で従来経路 (κ_t 再推定) を選べる", (not PWo.analytic) and hasattr(PWo, "_herm_c") and PWo.pw_upstream is None)

# --- pw_ramp (plan §5.1 #9, codex result M5) -----------------------------------------
check(f"pw_ramp ゲート (case/45 [−11, −6]) を記録: {PW.ramp_gate}", PW.ramp_gate["pass"] and PW.ramp_gate["source"] == "pw_ramp")
from forge_design.geometry.wall_axismach import default_pw_ramp  # noqa: E402
PD = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, upstream="ramp")
x_on = default_pw_ramp(w)[0]
check(f"既定 pw_ramp = (r′ < −0.05 の最初の x {x_on:.4f}, −0.5·L_U = {-0.5 * w.up.L_U:g}) でゲート合格 ({PD.ramp_gate})",
      PD._ramp == (x_on, -0.5 * w.up.L_U) and PD.ramp_gate["pass"]
      and float(w.r(np.r_[x_on], 1)[0]) < -0.05 <= float(w.r(np.r_[x_on - 1e-4], 1)[0]))
for bad in ((-13.0, -6.0), (-6.0, -11.0), (-11.0,)):
    try:
        PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, ramp=bad, upstream="ramp")
        check(f"不正な pw_ramp {bad} を拒否", False)
    except ValueError:
        check(f"不正な pw_ramp {bad} を拒否", True)
# 標準 geometry (L_U 3.5, r_U 2.5): 既定ランプは [−3.x, −1.75] と短く δ_r·s″ が縮流部の曲率を壊す → ゲートで止まることが合格
p35 = load_problem(C45 / "problem_d155_euler_c2final_n2400.yaml")
p35.geometry.update({"wall_repr": "joint", "L_U": 3.5, "r_inlet": 2.5})
d35 = design_chain(p35)
args35 = (d35["wall"], d35["wall_inv"], float(p35.spec["r_throat"]), float(p35.spec["Pt"]), float(p35.spec["Tt"]),
          _gam_or_gas(p35), p35.cp)
try:
    P35 = PhysicalNozzleWall(*args35, offset="radial", delta_r_x=f, upstream="ramp")
    check(f"L_U 3.5: 既定 pw_ramp のゲートで止まる (通ってしまった: {P35.ramp_gate})", False)
except ValueError as e:
    check(f"L_U 3.5: 既定 pw_ramp {default_pw_ramp(d35['wall'])} のゲートで止まる ({str(e)[:120]}…)", "ゲート不合格" in str(e))

# --- 既定は poly (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1) -------------------------------
# この δ_r 表 (run0051) は出口 x_e の手前 (x = 95.10) で終わり、その先は δ_r 一定・導関数 0 なので、物理壁 = 設計壁の r′ が出口直前で
# わずかに負 (−1.3e-8) になる。poly のゲート (スロートの後で r′ ≥ 0、丸めの許容 1e-12) はこれを不合格にする (ランプの壁の validate は
# 4000 点の差分の −1e-9 の許容で通していた)。既定が poly であることは、このゲートで止まることで確かめる。
try:
    PhysicalNozzleWall(*args, offset="radial", delta_r_x=f)
    check("upstream 省略 (既定 poly): 表の外で r′ < 0 になる壁をゲートで止める", False)
except ValueError as e:
    check(f"upstream 省略 (既定 poly): 表の外で r′ < 0 になる壁をゲートで止める ({str(e)[:110]}…)",
          "pw_upstream poly" in str(e) and "後で r′ ≥ 0 False" in str(e))
tb2 = tb[tb[:, 0] <= 95.0]
x_ext = np.linspace(tb2[-1, 0], w.x_e, 6)[1:]
f2 = delta_r_from_table(np.r_[tb2[:, 0], x_ext], np.r_[tb2[:, 1], tb2[-1, 1] + 5.15e-3 * (x_ext - tb2[-1, 0])])   # 出口まで覆う表
PP = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f2)
check(f"upstream 省略は poly (既定): {PP.pw_upstream} / {PP.pw_upstream_source}, ゲート {PP.upstream_gate['pass']}",
      PP.pw_upstream == "poly" and PP.pw_upstream_source == "default" and PP.upstream_gate["pass"])
try:
    PhysicalNozzleWall(*args, offset="radial", delta_r_x=f, ramp=(-11.0, -6.0))
    check("既定 (poly) と pw_ramp の併記を拒否", False)
except ValueError as e:
    check(f"既定 (poly) と pw_ramp の併記を拒否 ({str(e)[:60]}…)", "併記" in str(e))

print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
