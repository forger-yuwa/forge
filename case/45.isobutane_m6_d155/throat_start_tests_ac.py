"""plan tooling-nozzle-throat-monotone-r2 §5.1 #10 の切り分け試験 (a)・(c) (CFD 0 step)。
(a) ピン初期線の壁足 θ = 0 強制を、追跡で得た壁足直近の値に替えて、MOC 壁点の始点付近の角度差 Δθ = θ − atan(x/R) が変わるか。
(c) 始点付近の Δθ(x) を x∈[0.01, 0.04] で a + b·x に当てはめ、b = 0 (曲率が 1/R) にする M′_A の倍率を逆算し、
    M′_A の取り方 (軸 evenfit への多項式窓フィットの半幅・次数) によるばらつきと比べる。
usage: python3 throat_start_tests_ac.py  (凍結元 run はローカル主ツリーを参照)
"""
import sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.feedback import cfd_initial_line as CI  # noqa: E402
M = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
R = 2.0
orig_build, orig_anchor = CI.CFDPinnedThroat._build_line, CI.CFDPinnedThroat.axis_anchor


def run(ns=161, foot="zero", f1=1.0):
    def build(self, ds):
        orig_build(self, ds)
        if foot == "traced":
            self._line_T = self._line_T.copy(); self._line_T[-1] = self._line_T[-2]
    def anchor(self, x):
        a, b, c = orig_anchor(self, x); return a, b * f1, c
    CI.CFDPinnedThroat._build_line, CI.CFDPinnedThroat.axis_anchor = build, anchor
    try:
        p = load_problem(C / "problem_d155_ns_finemesh_recal_final.yaml"); p.geometry["n_start"] = ns
        p.geometry["initial_line_run"] = str(M / p.geometry["initial_line_run"])
        tb = design_chain(p)["wall_inv"]
    finally:
        CI.CFDPinnedThroat._build_line, CI.CFDPinnedThroat.axis_anchor = orig_build, orig_anchor
    x, th = tb[:, 0], tb[:, 2]
    return x, np.degrees(th - np.arctan(x / R)), tb


def slope(x, d, lo=0.01, hi=0.04):
    m = (x >= lo) & (x <= hi); c = np.polyfit(x[m], d[m], 1); return c[0], c[1]   # [deg/r_t], 切片 [deg]


print("== (a) 壁足 θ: 強制 0 と追跡値")
for ns in (161, 321):
    for foot in ("zero", "traced"):
        x, d, tb = run(ns, foot)
        b, a = slope(x, d)
        print(f"  n_start {ns} foot={foot:6s}: 壁足 θ_in={np.degrees(tb[0, 2]):+.5f}°  x {np.round(x[1:5], 5)}  Δθ {np.round(d[1:5], 4)}  [0.01,0.04] 当てはめ a={a:+.4f}° b={b:+.3f}°/r_t")
print("== (c) b = 0 にする M′_A 倍率 (n_start 161)")
F = (0.99, 1.0, 1.01, 1.02)
B = []
for f in F:
    x, d, _ = run(161, "zero", f); b, a = slope(x, d); B.append(b)
    print(f"  M′_A ×{f}: a={a:+.4f}° b={b:+.3f}°/r_t")
k = np.polyfit(F, B, 1); f0 = -k[1] / k[0]
print(f"  → b = 0 の倍率 ≈ {f0:.5f} (M′_A を {100*(f0-1):+.2f} %)")
print("== (c) M′_A の取り方によるばらつき (run_0062 の軸 evenfit への多項式窓フィット)")
ht = CI.pinned_factory(M / "run_0062_euler_wallfit_fit_r1_ext6k", "res_6000.h5")(R, 1.2735422978978068)
xa = ht.x0_cfd; base = orig_anchor(ht, xa)[1]
vals = {}
for hw in (0.10, 0.15, 0.25, 0.35, 0.50):
    for deg in (2, 3, 4, 5):
        w = np.abs(ht._ax_x - xa) <= hw
        if w.sum() <= deg + 2: continue
        c = np.polyfit(ht._ax_x[w] - xa, ht._ax_M[w], deg); vals[(hw, deg)] = c[-2]
print("  基準 (半幅 0.25, 4 次) M′_A =", round(base, 6), " 窓内の軸点数 (半幅 0.25):", int((np.abs(ht._ax_x - xa) <= 0.25).sum()))
for (hw, deg), v in vals.items():
    print(f"  半幅 {hw:.2f} 次数 {deg}: M′_A = {v:.6f} ({100*(v/base-1):+.3f} %)")
v = np.array(list(vals.values())); print(f"  範囲 {100*(v.min()/base-1):+.3f} 〜 {100*(v.max()/base-1):+.3f} %")
