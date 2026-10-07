"""plan tooling-nozzle-throat-monotone-r2 §3 仮説 H: 始点付近の不整合が「入力の初期線」にあるか「MOC を通した後」に出るかを分ける (CFD 0 step)。
(1) 入力: 初期線 (ピン = run_0062 場からの追跡、Hall = 級数解) 上の θ を、壁足からの距離 1−r → 0 で見る (壁足は取り出しで θ=0 に強制)。
(2) 出力: MOC 壁点の角度差 Δθ_i = θ_i − atan(x_i/R) を n_start 41/81/161 で並べ、x → 0+ の傾向を見る。
usage: python3 throat_start_limit.py  (凍結元 run はローカル主ツリーを参照)
"""
import sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.feedback.cfd_initial_line import pinned_factory  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
M = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
pin = load_problem(C / "problem_d155_ns_finemesh_recal_final.yaml")
d = design_chain(dict.__class__ and pin.__class__ and (lambda q: (q.geometry.__setitem__("initial_line_run", str(M / q.geometry["initial_line_run"])), q)[1])(pin))
R, g = float(d["R"]), float(d["law"].gamma) if hasattr(d["law"], "gamma") else None
print("R", R)
ht = pinned_factory(M / "run_0062_euler_wallfit_fit_r1_ext6k", "res_6000.h5")(R, d.get("gamma_hall", 1.2735422978978068))
lr, lT, lx = ht._line_r, ht._line_T, ht._line_x
print("== (1a) ピン入力: 追跡した線の θ (取り出しで壁足 r=1 は θ=0 に強制)")
for s in (1e-6, 1e-4, 1e-3, 3e-3, 1e-2, 2.5e-2, 5e-2):
    i = np.argmin(np.abs(lr[:-1] - (1 - s)))
    print(f"  1−r={1-lr[i]:.1e}  x={lx[i]:.5f}  θ={np.degrees(lT[i]):+.5f}°")
print("  壁足 (強制):", np.degrees(lT[-1]), "°")
print("== (1b) 場の θ の壁際 (x=0 断面、η→1) と壁そのもの")
for e in (0.9, 0.97, 0.99, 0.999, 1.0):
    print(f"  η={e}: θ_field(x=0)={np.degrees(ht.theta(0.0, e * float(ht._rw(0.0)))):+.5f}°")
hall = HallThroat(R=R, gamma=1.2735422978978068)
xs, rr, MM, tt = hall.throat_characteristic(n=161)
print("== (1c) Hall 入力 (n=161): r→1 の θ")
for k in (-1, -2, -3, -5, -9, -17):
    print(f"  1−r={1-rr[k]:.4f}  x={xs[k]:.5f}  θ={np.degrees(tt[k]):+.5f}°")
print("== (2) MOC 出力の壁点の角度差 Δθ_i = θ_i − atan(x_i/R) [deg]")
for prob, lab in ((C / "problem_d155_ns_finemesh_recal_final.yaml", "pin"), (C / "problem_d155_euler_c2final_n2400.yaml", "hall")):
    for ns in (41, 81, 161, 321):
        p = load_problem(prob); p.geometry["n_start"] = ns
        if p.geometry.get("initial_line") == "cfd":
            p.geometry["initial_line_run"] = str(M / p.geometry["initial_line_run"])
        try:
            tb = design_chain(p)["wall_inv"]
        except Exception as e:
            print(f"  {lab} n_start {ns}: ERR {str(e)[:90]}"); continue
        x, th = tb[:, 0], tb[:, 2]; dth = np.degrees(th - np.arctan(x / R))
        print(f"  {lab} n_start {ns:3d}: x {np.round(x[:6], 5)}  Δθ {np.round(dth[:6], 4)}")
