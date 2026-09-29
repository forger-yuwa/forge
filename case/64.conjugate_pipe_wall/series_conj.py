#!/usr/bin/env python3
r"""case/64 (A)・case/65 (C) の準定常の系列 (plan `boundary-cht-conjugate-benchmarks.md` §4.7)。forge の出力だけを使う。

全スナップショットについて、窓内の全壁節点の規格化した界面温度・界面熱流束と積分量を系列 CSV にし、
`check_quasisteady.py --series-csv --tail 0.5 --drift D --osc D` を登録の閾値で回す。

- 規格化尺度は**最終スナップショットの forge の値で固定** (系列の途中で尺度を変えない)。
  A: 温度 (T_i − T_in)/上昇 (加熱区間の長さ平均)、q_i/q_o (q_o = 加熱区間の外面から入る熱流束の平均、固体ダンプの q_hole から)。
  C: θ_i = (T_i − T_∞)/(T_h − T_∞)、q_i/(窓内平均)。
- 窓: A は**全長の壁節点**、C は x/L ∈ [0.2, 0.9]。熱流束も**窓内の全節点** (2026-09-30 codex diagnose M2: 初版は |q| < 5 % の節点を外していた)。
- 判定: 各系列を登録尺度で規格化したうえで、末尾半分の (i) 線形トレンドの幅 |傾き × 区間| と (ii) 変動幅 max − min を、
  **系列の平均でなく登録尺度 (=1) に対して** D と比べる (平均 0 近傍の節点で相対変動が発散しないため。check_quasisteady の
  drift/fluct の分母を登録尺度に置き換えたもの)。閾値 D = 比較許容の 1/5: A 温度 0.002・q 0.004、C θ 0.002・q 0.006。

    python3 series_conj.py A <run>
    python3 series_conj.py C <run>
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "65.conjugate_flat_plate"))
TOOL = HERE.parents[1] / "solver_density_cuda" / "tools" / "check_quasisteady.py"
trap = getattr(np, "trapezoid", None) or np.trapz


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=["A", "C"]); ap.add_argument("run")
    a = ap.parse_args()
    run = Path(a.run)
    if a.case == "A":
        import pipe_common as pc
        wall = "wall_3"; solid_pid = 3
    else:
        import plate_common as pc
        wall = "plate_5"; solid_pid = 5
    gc = pc.gc
    steps = sorted(int(m.group(1)) for p in run.glob(f"res_{wall}_*.h5") for m in [re.match(rf"res_{wall}_(\d+)\.h5$", p.name)] if m)
    if len(steps) < 8:
        print(f"REFUSED: スナップショット {len(steps)} 枚 (< 8)\nVERDICT: REFUSED"); return 2
    def wall_at(st):
        with h5py.File(run / f"res_{wall}_{st}.h5", "r") as w:
            c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
            q = -np.asarray(w["VALUE/iface_q_eff"][:], float); T = np.asarray(w["VALUE/iface_Tw_bc"][:], float)
        o = np.argsort(c[:, 0]); return c[o, 0], T[o], q[o]
    x, Tl, ql = wall_at(steps[-1])
    if a.case == "A":
        win = np.ones_like(x, bool); heat = (x >= -1e-12) & (x <= pc.L_HEAT + 1e-9)
        tscale = trap(Tl[heat] - gc.T_IN, x[heat]) / pc.L_HEAT
        with h5py.File(run / f"res_solid_{solid_pid}_{steps[-1]}.h5", "r") as s:
            qhole = float(np.asarray(s["VALUE/q_hole"][:], float).sum())
        qscale = abs(qhole) / (pc.R_O * pc.L_HEAT)            # per rad の熱量 / 外面の長さ (per rad)
        tnorm = lambda T: (T - gc.T_IN) / tscale
        Dt, Dq = 0.002, 0.004
    else:
        win = (x >= 0.2 * pc.L - 1e-12) & (x <= 0.9 * pc.L + 1e-12)
        qscale = trap(ql[win], x[win]) / (x[win][-1] - x[win][0])
        tnorm = lambda T: (T - gc.T_IN) / pc.DT_H
        Dt, Dq = 0.002, 0.006
    qsel = win.copy()
    rows = []
    for st in steps:
        xx, T, q = wall_at(st)
        if not np.allclose(xx, x): print("REFUSED: 壁節点が時刻で違う\nVERDICT: REFUSED"); return 2
        integ = (pc.R * trap(q[x <= 0], x[x <= 0]) if a.case == "A" else trap(q[win], x[win]) / (x[win][-1] - x[win][0])) / (qscale * (pc.R_O * pc.L_HEAT if a.case == "A" else 1.0))
        rows.append(np.concatenate([[st], tnorm(T[win]), q[qsel] / qscale, [integ]]))
    R = np.array(rows)
    nT, nQ = int(win.sum()), int(qsel.sum())
    colsT = [f"T_n{i}" for i in range(nT)]; colsQ = [f"q_n{i}" for i in range(nQ)]
    out = run / "series_conj.csv"
    np.savetxt(out, R, delimiter=",", comments="", header=",".join(["step"] + colsT + colsQ + ["integral"]), fmt=["%d"] + ["%.12e"] * (nT + nQ + 1))
    res = []
    tail = R[len(R) // 2:]
    st_t = tail[:, 0]
    def judge(block, D, nm):
        span = st_t.max() - st_t.min()
        slope = np.polyfit(st_t, block, 1)[0] if len(st_t) > 1 else np.zeros(block.shape[1])
        drift = np.abs(slope * span); fl = np.ptp(block, axis=0)
        worst = float(max(drift.max(), fl.max())); nbad = int(((drift > D) | (fl > D)).sum())
        return (nm, D, f"最大 {worst:.3e} (drift {drift.max():.3e}、変動 {fl.max():.3e})", nbad, block.shape[1], nbad == 0)
    res.append(judge(tail[:, 1:1 + nT], Dt, "温度"))
    res.append(judge(tail[:, 1 + nT:], Dq, "熱流束・積分量"))
    print(f"=== 準定常 {a.case} {run}  スナップショット {len(steps)} 枚 (step {steps[0]}–{steps[-1]})、窓の温度節点 {nT}、熱流束節点 {nQ}")
    ok = True
    for nm, D, ov, nbad, n, good in res:
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {nm}: {ov}  (閾値超え {nbad} / {n} 系列、D {D})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
