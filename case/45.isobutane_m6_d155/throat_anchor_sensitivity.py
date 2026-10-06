"""plan tooling-nozzle-throat-monotone-r2 §3 仮説 H の候補 (A) の感度試験 (CFD 0 step):
軸アンカー (M_A, M′_A, M″_A) を少し変えて、MOC 壁点の始点付近の角度差 Δθ = θ − atan(x/R) が反応するかを見る。
逆設計では、始点付近の壁点から出る C⁻ は軸アンカー A の直後に着くので、始点付近の壁は初期線と A 直後の軸分布の両方で決まる (Goursat 問題) と見込む。
摂動は CFDPinnedThroat.axis_anchor の戻り値を差し替えるだけ (初期線・m* は不変)。
usage: python3 throat_anchor_sensitivity.py  (凍結元 run はローカル主ツリーを参照)
"""
import sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.feedback import cfd_initial_line as CI  # noqa: E402
M = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
orig = CI.CFDPinnedThroat.axis_anchor
R = 2.0


def run(dM=0.0, f1=1.0, f2=1.0, ns=161):
    def patched(self, x):
        a, b, c = orig(self, x)
        return a + dM, b * f1, c * f2
    CI.CFDPinnedThroat.axis_anchor = patched
    try:
        p = load_problem(C / "problem_d155_ns_finemesh_recal_final.yaml"); p.geometry["n_start"] = ns
        p.geometry["initial_line_run"] = str(M / p.geometry["initial_line_run"])
        d = design_chain(p)
    finally:
        CI.CFDPinnedThroat.axis_anchor = orig
    tb = d["wall_inv"]; x, th = tb[:, 0], tb[:, 2]
    return x, np.degrees(th - np.arctan(x / R)), d


x0, b0, d0 = run()
a = None
for k, v in d0.items():
    if "anchor" in k.lower() and not isinstance(v, (np.ndarray,)):
        print(" ", k, str(v)[:200])
print("基準 (n_start 161): x", np.round(x0[:6], 5), " Δθ", np.round(b0[:6], 4), " Δθ(x≈0.05)", round(float(np.interp(0.05, x0, b0)), 4))
for lab, kw in (("M′_A ×1.01", dict(f1=1.01)), ("M′_A ×0.99", dict(f1=0.99)),
                ("M″_A ×1.10", dict(f2=1.10)), ("M″_A ×0.90", dict(f2=0.90)),
                ("M_A +1e-4", dict(dM=1e-4)), ("M_A −1e-4", dict(dM=-1e-4))):
    try:
        x, b, _ = run(**kw)
        print(f"{lab:12s}: Δθ {np.round(b[:6], 4)}  変化 {np.round(b[:6] - b0[:6], 4)}  Δθ(x≈0.05) {float(np.interp(0.05, x, b)):+.4f}")
    except Exception as e:
        print(f"{lab:12s}: ERR {str(e)[:120]}")
