#!/usr/bin/env python3
r"""case/64 (A)・case/65 (C) の準定常の系列 (plan `boundary-cht-conjugate-benchmarks.md` §4.7)。forge の出力だけを使う。

全スナップショットについて、窓内の全壁節点の規格化した界面温度・界面熱流束と積分量を系列 CSV にし、
`check_quasisteady.py --series-csv --tail 0.5 --drift D --osc D` を登録の閾値で回す。

- 規格化尺度は**最終スナップショットの forge の値で固定** (系列の途中で尺度を変えない)。
  A: 温度 (T_i − T_in)/上昇 (加熱区間の長さ平均)、q_i/q_o (q_o = 加熱区間の外面から入る熱流束の平均、固体ダンプの q_hole から)。
  C: θ_i = (T_i − T_∞)/(T_h − T_∞)、q_i/(窓内平均)。
- 窓: A は加熱区間 + 予熱域 (−40R ≤ x ≤ L_h)、C は x/L ∈ [0.2, 0.9]。熱流束は |q_i| ≥ 0.05 × 尺度の節点だけ (0 近傍は相対変動が意味を失う)。
- 閾値 D = 比較許容の 1/5: A 温度 0.002・q 0.004、C θ 0.002・q 0.006。**絶対温度そのままの系列には掛けない**。

    python3 series_conj.py A <run>
    python3 series_conj.py C <run>
"""
from __future__ import annotations

import argparse
import re
import subprocess
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
        win = (x >= -40 * pc.R - 1e-12) & (x <= pc.L_HEAT + 1e-9); heat = (x >= -1e-12) & (x <= pc.L_HEAT + 1e-9)
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
    qsel = win & (np.abs(ql) >= 0.05 * abs(qscale))
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
    for cols, D, nm in ((colsT, Dt, "温度"), (colsQ + ["integral"], Dq, "熱流束・積分量")):
        p = subprocess.run([sys.executable, str(TOOL), "--series-csv", str(out), "--series-cols", ",".join(cols), "--tail", "0.5",
                            "--drift", str(D), "--osc", str(D)], capture_output=True, text=True)
        o = p.stdout + p.stderr
        (run / f"series_conj_QS_{'T' if nm == '温度' else 'q'}.txt").write_text(o)
        ov = [l for l in o.splitlines() if "OVERALL" in l]
        nbad = sum(1 for l in o.splitlines() if re.search(r"DRIFTING|OSCILLATING|TRANSIENT|NONFINITE", l))
        res.append((nm, D, ov[-1].strip() if ov else "(OVERALL 無し)", nbad, len(cols)))
    print(f"=== 準定常 {a.case} {run}  スナップショット {len(steps)} 枚 (step {steps[0]}–{steps[-1]})、窓の温度節点 {nT}、熱流束節点 {nQ}")
    ok = True
    for nm, D, ov, nbad, n in res:
        good = "ALL STEADY" in ov
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {nm}: {ov}  (非 STEADY {nbad} / {n} 系列、drift/osc {D})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
