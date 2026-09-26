#!/usr/bin/env python3
"""plan gradient-scalar-lsq-unification §5.1 #5h(b) (codex result-1 M3): FCT smoke の nSub 感度。

D_q = ‖q_lsq15 − q_gg15‖∞、E_q = ‖q_gg30 − q_gg15‖∞ + ‖q_lsq30 − q_lsq15‖∞ を、同じ物理時刻 (res_100・res_200) の
同じ節点で出す (ノルムは節点の最大絶対差に固定)。A: 全 q・全時刻で E_q < 0.1 D_q (保存・有界性ゲートは別に
check_passive_budget の 4 本の VERDICT を貼る)。B: 不成立。D_q = 0 は判定不能。

  python3 m3_compare.py GG15 LSQ15 GG30 LSQ30 [--out M3_fct.txt]
"""
import argparse
import os

import h5py
import numpy as np

Q = ("ro", "roUx", "roUy", "roe", "roY0", "roY1", "roY2", "roY3", "roY4", "roXi")


def load(run, st):
    with h5py.File(os.path.join(run, f"res_{st}.h5"), "r") as f:
        return {k: np.array(f["VALUE"][k]).astype(np.float64) for k in Q if k in f["VALUE"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gg15"); ap.add_argument("lsq15"); ap.add_argument("gg30"); ap.add_argument("lsq30")
    ap.add_argument("--out", default="M3_fct.txt")
    a = ap.parse_args()
    L = [f"# #5h(b) FCT smoke の nSub 15→30 感度\n- gg15 {a.gg15}\n- lsq15 {a.lsq15}\n- gg30 {a.gg30}\n- lsq30 {a.lsq30}",
         "| 時刻 | q | D_q = ‖lsq15−gg15‖∞ | E_q = ‖gg30−gg15‖∞+‖lsq30−lsq15‖∞ | E/D | 判定 (E < 0.1 D) |",
         "| --- | --- | --- | --- | --- | --- |"]
    ok, und = True, False
    for st in (100, 200):
        g15, l15, g30, l30 = (load(r, st) for r in (a.gg15, a.lsq15, a.gg30, a.lsq30))
        for q in Q:
            if q not in g15:
                continue
            D = np.abs(l15[q] - g15[q]).max()
            E = np.abs(g30[q] - g15[q]).max() + np.abs(l30[q] - l15[q]).max()
            if D == 0.0:
                L.append(f"| {st} | {q} | 0 | {E:.3e} | - | 判定不能 (D = 0) |"); und = True; continue
            good = E < 0.1 * D; ok &= good
            L.append(f"| {st} | {q} | {D:.3e} | {E:.3e} | {E / D:.3f} | {'ok' if good else 'NG'} |")
    v = "判定不能" if und else ("感度条件成立 (保存ゲートと合わせて A/B)" if ok else "B (感度条件不成立)")
    L.append(f"\nVERDICT #5h(b) 感度: {v}")
    txt = "\n".join(L) + "\n"
    open(a.out, "w").write(txt); print(txt)


if __name__ == "__main__":
    main()
