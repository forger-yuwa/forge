#!/usr/bin/env python3
r"""上流延長の感度の原因診断 (2026-09-30 disposition の判別 A/B、合否には使わない)。
同じ延長格子・流れ場・BC で、**追加区間 (x < −80R) の源項 S = Φ + u·∇p** だけを A: そのまま / B: 0 にして解き、
元の領域 (x ≥ −80R) の壁温の、延長なしの参照 (中間水準) に対する最大差 D_A・D_B を出す。
判定: D_B < 0.5 D_A なら「追加区間の源項が主因」を支持、それ以外は棄却。

    python3 ext_source_ab.py <run>
"""
import re, sys
from pathlib import Path
import numpy as np
import pipe_common as pc
import eval_conj as ev


def main():
    run = Path(sys.argv[1])
    pc.CASE = re.search(r"kind cht (A1|A2)", (run / "RUN_INPUTS.txt").read_text()).group(1)
    st = ev.last_step(run)
    xs_f, ys_f, G, _ = ev.load_fields_only(run, st)
    P0, xs0, ys0, jw0, _ = ev.build("A", xs_f, ys_f, G, 1, pc=pc); T0 = ev._solve_rowwise(P0, P0.T_in_profile)
    D = {}
    for tag in ("A", "B"):
        P, xs, ys, jw, _ = ev.build("A", xs_f, ys_f, G, 1, pc=pc, extend=True)
        if tag == "B":
            P.S[xs < -80 * pc.R - 1e-12, :] = 0.0
        T = ev._solve_rowwise(P, P.T_in_profile)
        i0 = int(np.argmin(np.abs(xs - xs0[0])))
        D[tag] = np.abs(T[i0:i0 + len(xs0), jw] - T0[:, jw0]).max()
    heat = (xs0 >= 0) & (xs0 <= pc.L_HEAT)
    rise = np.trapezoid(T0[heat, jw0] - pc.gc.T_IN, xs0[heat]) / pc.L_HEAT if hasattr(np, "trapezoid") else 0
    print(f"=== 延長感度の原因 A/B {run} ({pc.CASE}) step {st}")
    print(f"  D_A (源項あり) {D['A']:.4e} K ({D['A']/rise*100:.3f} % of 上昇)、D_B (追加区間の源項 0) {D['B']:.4e} K ({D['B']/rise*100:.3f} %)")
    v = "「追加区間の源項が主因」を支持" if D["B"] < 0.5 * D["A"] else "「追加区間の源項が主因」を棄却 (第 2 仮説を残す)"
    print(f"VERDICT: {v}  (D_B/D_A = {D['B']/D['A']:.3f})")


if __name__ == "__main__":
    main()
