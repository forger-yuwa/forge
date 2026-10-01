#!/usr/bin/env python3
r"""case/65 cfl 判別 (plan `boundary-cht-conjugate-flat-plate.md` §4.6、事前登録どおり): limiter 0 のもとで cfl_pseudo 2 が使えるか。

    python3 ab_judge_cfl2lim0.py <run_0020 (limiter 0・cfl 2、6000 step)>

区間は 2000 step × 3 (最初の区間は切替の過渡)。「cfl 2 採用」はすべて満たすとき:
1. NaN/Inf なし (残差・系列)
2. `rms_ro/rms_roUx/rms_roUy/rms_roe` の区間平均が 3 区間で厳密減少
3. 最終区間の 4 列が run_0019 の最終区間 (累計 [44000,48000)) 以下
4. 最終区間の規格化標準偏差 (節点最大): 上流 slip P ≤ 3.378e-5、前縁帯の界面 q ≤ 1.139e-2、下流 slip P ≤ 1.354e-4
   (limiter 2・cfl 2 の対照 run_0015 最終区間の 1/10)、かつ上流 slip P (最大変動の節点) の最大ピークのパワー比 < 50 %
1 つでも外れたら本番は cfl 0.5。最終行に `DECISION: cfl_pseudo=2.0` または `DECISION: cfl_pseudo=0.5` を出す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plate_common as pc  # noqa: E402
from ab_judge_cfl import COLS, Q_INF, fmt  # noqa: E402
from ab_judge_lim import resid  # noqa: E402
from ab_series import q_scale, refuse  # noqa: E402

LI = 2000
REF19 = {"rms_ro": 1.185e-11, "rms_roUx": 8.703e-09, "rms_roUy": 1.901e-09, "rms_roe": 3.405e-06}
LIM = {"上流 slip の P": 3.378e-5, "前縁帯の界面 q": 1.139e-2, "下流 slip の P": 1.354e-4}


def main():
    run = Path(sys.argv[1])
    print(f"=== cfl 判別 (limiter 0、cfl_pseudo 2)  {run.name}  2000 step × 3")
    R, _ = resid(run, 0, 3 * LI, LI)
    d = np.load(run / "ab_series.npz", allow_pickle=True)
    st = d["step"]
    for nm in d.files:
        if d[nm].dtype.kind == "f" and not np.isfinite(d[nm]).all():
            refuse(f"系列 {nm} に非有限値")
    if len(st) < 3 * LI or (np.diff(st[:3 * LI]) != 1).any():
        refuse(f"系列が {len(st)} 点・不連続")
    g = d["group"]; wx = d["wall_x"] / pc.L; qs = q_scale()
    tail = slice(2 * LI, 3 * LI)
    up = d["P"][tail][:, g == "up"] / Q_INF
    amp = {"上流 slip の P": np.std(up, axis=0).max(),
           "前縁帯の界面 q": np.std(d["q"][tail][:, wx <= 0.02 + 1e-9] / qs, axis=0).max(),
           "下流 slip の P": np.std(d["P"][tail][:, g == "down"] / Q_INF, axis=0).max()}
    s = up[:, np.argmax(np.std(up, axis=0))]; s = s - s.mean()
    pw = np.abs(np.fft.rfft(s)) ** 2
    pfrac = float(pw[1:].max() / pw[1:].sum()) if pw[1:].sum() > 0 else 0.0
    f = np.fft.rfftfreq(len(s)); per = 1 / f[1:][np.argmax(pw[1:])]
    ok = True
    c2 = all(R[c][0] > R[c][1] > R[c][2] for c in COLS)
    ok &= c2
    print("  残差の区間平均:")
    for c in COLS:
        print(f"    {c:9s} {fmt(R[c])}  (run_0019 最終 {REF19[c]:.3e})")
    print(f"  {'PASS' if c2 else 'FAIL'}  条件 2: 4 列が 3 区間で厳密減少")
    c3 = all(R[c][-1] <= REF19[c] for c in COLS)
    ok &= c3
    print(f"  {'PASS' if c3 else 'FAIL'}  条件 3: 最終区間 ≤ run_0019 最終区間 (" + ", ".join(f"{c} {R[c][-1]/REF19[c]:.3f}" for c in COLS) + ")")
    c4 = all(amp[k] <= LIM[k] for k in LIM) and pfrac < 0.5
    ok &= c4
    print(f"  {'PASS' if c4 else 'FAIL'}  条件 4: " + ", ".join(f"{k} {amp[k]:.3e} (≤ {LIM[k]:.3e})" for k in LIM)
          + f"、上流 P の最大ピーク 周期 {per:.1f} step・パワー比 {pfrac:.3f} (< 0.5)")
    print("  PASS  条件 1: NaN/Inf なし")
    print(f"VERDICT: {'cfl 2 採用' if ok else 'cfl 0.5 採用 (cfl 2 の条件を満たさない)'}")
    print(f"DECISION: cfl_pseudo={'2.0' if ok else '0.5'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
