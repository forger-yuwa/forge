#!/usr/bin/env python3
r"""case/64 (A)・case/65 (C) の準定常の系列 (plan `boundary-cht-conjugate-benchmarks.md` §4.7)。forge の出力だけを使う。

全スナップショットについて、窓内の全壁節点の規格化した界面温度・界面熱流束と積分量を系列 CSV にし、
`check_quasisteady.py --series-csv --tail 0.5 --drift D --osc D` を登録の閾値で回す。

- 規格化尺度は**最終スナップショットの forge の値で固定** (系列の途中で尺度を変えない)。
  A: 温度 (T_i − T_in)/上昇 (加熱区間の長さ平均)、q_i/q_o (q_o = 加熱区間の外面から入る熱流束の平均、固体ダンプの q_hole から)。
  C: θ_i = (T_i − T_∞)/(T_h − T_∞)、q_i/(窓内平均)。
- 窓: A は**全長の壁節点**、C は x/L ∈ [0.2, 0.9]。熱流束も**窓内の全節点** (2026-09-30 codex diagnose M2: 初版は |q| < 5 % の節点を外していた)。
- 判定: 登録尺度で規格化した系列を正式ツール `check_quasisteady.py --abs-scale 1` に渡す (drift/fluct の分母を系列の平均でなく登録尺度にする。
  平均 0 近傍の節点で相対変動が発散しないため)。系列・座標・尺度の非有限・0、壁節点の欠落・重複は REFUSED (2026-09-30 result レビュー M3)。
  固体の効果 (厚さ方向の温度差・軸方向熱量) の時系列も含める (M5)。閾値 D = 比較許容の 1/5: A 温度 0.002・q 0.004、C θ 0.002・q 0.006。

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
    import eval_conj as ev
    if a.case == "A":
        import pipe_common as pc
        wall = "wall_3"; pid = 3
    else:
        import plate_common as pc
        wall = "plate_5"; pid = 5
    pc.CASE = re.search(r"kind cht (A1|A2|C1|C2)", (run / "RUN_INPUTS.txt").read_text()).group(1)
    gc = pc.gc; k_s = pc.KS_RATIO[pc.CASE] * gc.K_F
    steps = sorted(int(m.group(1)) for p in run.glob(f"res_{wall}_*.h5") for m in [re.match(rf"res_{wall}_(\d+)\.h5$", p.name)] if m)
    if len(steps) < 8:
        ev.refuse(f"スナップショット {len(steps)} 枚 (< 8)")
    # 期待節点集合 (メッシュから) と照合 (2026-09-30 result レビュー M3)
    xs_f, ys_f, _, _ = ev.load_fields_only(run, steps[-1])
    exp = ev.expected_wall_x(a.case, xs_f, ys_f, pc)
    def wall_at(st):
        with h5py.File(run / f"res_{wall}_{st}.h5", "r") as w:
            c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
            q = -np.asarray(w["VALUE/iface_q_eff"][:], float); T = np.asarray(w["VALUE/iface_Tw_bc"][:], float)
        o = np.argsort(c[:, 0]); x = c[o, 0]
        if len(x) != len(np.unique(np.round(x, 10))): ev.refuse(f"step {st}: 壁ダンプの節点に重複")
        if len(x) != len(exp) or not np.allclose(x, exp, rtol=0, atol=1e-9): ev.refuse(f"step {st}: 壁ダンプの節点が期待集合と一致しない")
        if not (np.isfinite(T).all() and np.isfinite(q).all()): ev.refuse(f"step {st}: 壁ダンプに非有限値")
        return x, T[o], q[o]
    x, Tl, ql = wall_at(steps[-1])
    xs_s, ys_s, Ts_l, Qtot_l = ev.solid_from_dump(run, steps[-1], pid, k_s)
    if a.case == "A":
        win = np.ones_like(x, bool); heat = (x >= -1e-12) & (x <= pc.L_HEAT + 1e-9)
        tscale = trap(Tl[heat] - gc.T_IN, x[heat]) / pc.L_HEAT
        qscale = Qtot_l / (pc.R_O * pc.L_HEAT)
        tnorm = lambda T: (T - gc.T_IN) / tscale
        Dt, Dq = 0.002, 0.004
    else:
        win = (x >= 0.2 * pc.L - 1e-12) & (x <= 0.9 * pc.L + 1e-12)
        qscale = trap(ql[win], x[win]) / (x[win][-1] - x[win][0])
        tnorm = lambda T: (T - gc.T_IN) / pc.DT_H
        Dt, Dq = 0.002, 0.006
    for nm, v in (("温度の尺度", tscale if a.case == "A" else pc.DT_H), ("熱流束の尺度", qscale)):
        if not (np.isfinite(v) and abs(v) > 0): ev.refuse(f"{nm}が 0 または非有限 ({v})")
    rows = []
    for st in steps:
        xx, T, q = wall_at(st)
        xs_s2, ys_s2, Ts, Qtot = ev.solid_from_dump(run, st, pid, k_s)
        dTs, Qax = ev.solid_indicators(a.case, pc, xs_s2, ys_s2, Ts, k_s, Qtot)
        integ = (pc.R * trap(q[x <= 0], x[x <= 0]) / Qtot) if a.case == "A" else trap(q[win], x[win]) / (x[win][-1] - x[win][0]) / qscale
        rows.append(np.concatenate([[st], tnorm(T[win]), q[win] / qscale, [integ, dTs / (tscale if a.case == "A" else pc.DT_H), Qax]]))
    R = np.array(rows)
    if not np.isfinite(R).all(): ev.refuse("系列に非有限値")
    nT = int(win.sum())
    colsT = [f"T_n{i}" for i in range(nT)]; colsQ = [f"q_n{i}" for i in range(nT)] + ["integral", "dTs_norm", "Qax_frac"]
    out = run / "series_conj_allnodes.csv"
    np.savetxt(out, R, delimiter=",", comments="", header=",".join(["step"] + colsT + colsQ), fmt=["%d"] + ["%.12e"] * (2 * nT + 3))
    print(f"=== 準定常 {a.case} ({pc.CASE}) {run}  スナップショット {len(steps)} 枚 (step {steps[0]}–{steps[-1]})、"
          f"評価窓内の全節点 {nT} (期待集合と一致)、系列 {out.name}")
    ok = True
    for cols, D, nm, tag in ((colsT, Dt, "温度", "T"), (colsQ, Dq, "熱流束・積分量・固体の効果", "q")):
        p = subprocess.run([sys.executable, str(TOOL), "--series-csv", str(out), "--series-cols", ",".join(cols), "--tail", "0.5",
                            "--drift", str(D), "--osc", str(D), "--abs-scale", "1"], capture_output=True, text=True)
        o = p.stdout + p.stderr
        (run / f"series_conj_QS_{tag}.txt").write_text(o)
        ov = [l for l in o.splitlines() if "OVERALL" in l or "->" in l]
        nbad = sum(1 for l in o.splitlines() if re.search(r"DRIFTING|OSCILLATING|TRANSIENT|NONFINITE", l))
        good = p.returncode == 0 and nbad == 0 and any("STEADY" in l for l in ov)
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {nm}: check_quasisteady --abs-scale 1 --drift/--osc {D} → 非 STEADY {nbad} / {len(cols)} 系列 (rc {p.returncode})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
