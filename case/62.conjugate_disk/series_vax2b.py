#!/usr/bin/env python3
r"""V-ax2b の系列判定 (plan §6 V-ax2b「合格」の格子部分)。

登録: **N_r=32 一様と非一様 (等比 1.1) でそれぞれ全界面節点 (端点を含む) の |T_w − T_w*| ≤ 0.5 % of 降下**、
かつ系列 **N_r = 8 / 16 / 32 一様で最大誤差が N_r 倍増ごとに 1/2.5 以下** (e16/e8 ≤ 0.4、e32/e16 ≤ 0.4)。

参考予測 (codex の独立計算、荷重移送だけの試験 — 連成計算の予測ではない): 0.9421 / 0.2983 / 0.08796 %。
`test_eval_vax2b.py` の合成 run でこの 3 値を再現済み (非一様は 0.1617 %)。

各 run の G-cons・連成保存・恒等式・準定常・G-if・流体残差は `eval_vax2b.py` と各ツールで別に判定する。

    python3 case/62.conjugate_disk/series_vax2b.py <run_r8u> <run_r16u> <run_r32u> --nonuniform <run_r32g> [--step N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "61.conjugate_annulus"))
import axcht  # noqa: E402
from eval_vax2b import NAME_CJ, PID_CJ, analytic  # noqa: E402


def err(run: Path, step, an):
    st = step if step is not None else axcht.last_step(run)
    d = axcht.wall_dump(run, NAME_CJ, PID_CJ, st)
    e = np.abs(d["Ts"] - an["Tw"]) / an["drop"] * 100
    return st, len(e), float(e.max()), float(d["xyz"][int(np.argmax(e)), 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("r8")
    ap.add_argument("r16")
    ap.add_argument("r32")
    ap.add_argument("--nonuniform", required=True, help="N_r=32 等比 1.1 の run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--tol", type=float, default=0.5, help="[%% of 降下] (登録値)")
    ap.add_argument("--ratio", type=float, default=0.4, help="倍増ごとの誤差比の上限 (= 1/2.5、登録値)")
    a = ap.parse_args()
    an = analytic()
    rows = {}
    for key, r, nr in (("8u", a.r8, 9), ("16u", a.r16, 17), ("32u", a.r32, 33), ("32g", a.nonuniform, 33)):
        st, n, e, rpos = err(Path(r), a.step, an)
        if n != nr:
            raise SystemExit(f"{r}: 界面節点 {n} (N_r={nr - 1} の run ではない)")
        rows[key] = e
        print(f"N_r={key:4s} {r}: step {st}, max |T_w−T_w*| = {e:.5f} % of 降下 (r = {rpos*1e3:.3f} mm)")
    ok = []
    ok.append(rows["32u"] <= a.tol)
    ok.append(rows["32g"] <= a.tol)
    q1, q2 = rows["16u"] / rows["8u"], rows["32u"] / rows["16u"]
    ok.append(q1 <= a.ratio)
    ok.append(q2 <= a.ratio)
    print(f"\n  N_r=32 一様   {rows['32u']:.5f} % <= {a.tol}   {'PASS' if ok[0] else '**FAIL**'}")
    print(f"  N_r=32 非一様 {rows['32g']:.5f} % <= {a.tol}   {'PASS' if ok[1] else '**FAIL**'}")
    print(f"  e16/e8  = {q1:.4f} <= {a.ratio}   {'PASS' if ok[2] else '**FAIL**'}")
    print(f"  e32/e16 = {q2:.4f} <= {a.ratio}   {'PASS' if ok[3] else '**FAIL**'}")
    v = "PASS" if all(ok) else "FAIL"
    print(f"\nVERDICT: {v}")
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
