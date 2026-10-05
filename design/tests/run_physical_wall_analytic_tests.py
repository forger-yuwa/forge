#!/usr/bin/env python3
"""joint 壁の物理壁 (`PhysicalNozzleWall` の解析経路) の試験 = plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6b (事前登録)。

壁 = P1(b) の CFD ピン joint 壁 (case/45 n2400、run_0062 res_6000 から凍結した初期線、Md_moc_offset −4.16e-4)、
δ_r = case/45.isobutane_m6_d155/_band_ab/delta_r_c2final_run0051.csv (`delta_r_from_table`)。
物理壁 r_W = r_design + s·δ_r (s: x ≤ −11 で 0、[−11, −6] で 5 次 smoothstep、x ≥ −6 で 1) の合格条件:
  x ≤ −11 で設計壁と |Δr| ≤ 1e-12 / x ≥ −6 で (設計 + δ_r) と ≤ 1e-12 / [−11, −6] で r′ < 0 かつ |r″ − r″_design| ≤ 5e-3 /
  x=0 の r″ 跳び ≤ 1e-8、r″(0⁺) = 0.5 + δ_r″(0) ± 1e-3 / [0, 0.3] の r″ max ≤ 0.52 + max δ_r″ /
  r‴ max ([0, 0.3]) ≤ 設計壁の r‴ max + 0.05 / validate の非単調 0 件 (δ_r 表の範囲内)。
加えて δ_r の導関数 (差分との照合・表の範囲外 0) と入力の拒否。

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
PW = PhysicalNozzleWall(*args, offset="radial", delta_r_x=f)
check("joint 壁では解析経路が自動で選ばれる", PW.analytic)

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
check("analytic=False で従来経路 (κ_t 再推定) を選べる", (not PWo.analytic) and hasattr(PWo, "_herm_c"))

print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
