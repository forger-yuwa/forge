"""plan tooling-nozzle-throat-monotone-r2 §5.1 #10 の事前試算 (CFD 0 step、ユーザ依頼 2026-10-07「その数分で終わる作業やって」)。
MOC の軸側の誤差を直したら何が変わるかを、単調壁の生産問題 (problem_d155_ns_finemesh_recal_final_mono.yaml) で測る。
条件: 単位過程 {prod (生産), K2c (試験用: 修正子を収束 + 軸上の sinθ/r の代用をやめる)} × axis_dx0 {0.03 (生産), 0.015, 0.0075}。
  n_axis_inv は 2400 のまま。各条件を n_start 41 (生産) と 321 (始点付近の観察用) で回す。
測る量:
  - 始点付近の角度差 Δθ = θ − atan(x/R) の第 1 点 (n_start 41 で x ≈ 0.025、321 で x ≈ 0.003)
  - 単調拘束なしの joint 当てはめの r″ の山 (x ∈ [0, 0.3] の最大 − 1/R) と、単調拘束ありの第 1 点の流れ角ずれ (代償)
  - 設計壁 (単調拘束あり) の生産条件からの変化: 半径の最大差 [µm] とその位置、出口の半径・壁角の差
  - MOC 点群の出口半径 r_e/r_t (面積比、r_t の解き直しが要るかの目安)
試験用の単位過程は throat_moc_restart_test.kernel (with を抜けると生産に戻る)。生産コードは変えない。
usage: [CASE_RUNS=<run_0062 のある case dir>] python3 throat_moc_fix_probe.py → _band_ab/throat_moc_fix_probe.json
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design")); sys.path.insert(0, str(C))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.geometry.wall_axismach import joint_fit_wall  # noqa: E402
from throat_moc_restart_test import kernel  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", C))
PROB = "problem_d155_ns_finemesh_recal_final_mono.yaml"
RT_UM = 76.6539e3
R = 2.0


def chain(kname, dx0, n_start):
    p = load_problem(C / PROB)
    p.geometry.update(n_start=int(n_start), axis_dx0=float(dx0))
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    t0 = time.time()
    with kernel(kname):
        d = design_chain(p)
    return d, time.time() - t0


rows, walls = [], {}
for kname in ("prod", "K2c"):
    for dx0 in (0.03, 0.015, 0.0075):
        for ns in (41, 321):
            try:
                d, dt = chain(kname, dx0, ns)
            except Exception as e:  # noqa: BLE001
                rows.append(dict(kernel=kname, dx0=dx0, n_start=ns, error=f"{type(e).__name__}: {str(e)[:200]}")); print(rows[-1], flush=True)
                continue
            tb = d["wall_inv"]; x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
            row = dict(kernel=kname, dx0=dx0, n_start=ns, seconds=round(dt, 1),
                       x1=float(x[1]), dth1_deg=float(np.degrees(th[1] - np.arctan(x[1] / R))),
                       re_moc_rt=float(r[-1]), xe_moc_rt=float(x[-1]), th_exit_deg=float(np.degrees(th[-1])))
            if ns == 41:
                s0, _ = joint_fit_wall(tb, R)                     # 単調拘束なし (旧壁の当てはめ)
                xg = np.linspace(0, 0.3, 30001)
                row["bump_unconstrained"] = float(s0(xg, 2).max() - 1 / R)
                W = d["wall"]                                     # 単調拘束あり (生産)
                row["mono_first_point_dev_deg"] = float(np.degrees(np.arctan(W._spl(np.r_[x[1]], 1)[0]) - th[1]))
                walls[(kname, dx0)] = (W, x[-1])
            rows.append(row); print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)

base, xe0 = walls[("prod", 0.03)]
xs = np.linspace(-12.0, xe0, 200001)
r0 = base.r(xs)
for row in rows:
    if row.get("n_start") != 41 or "error" in row:
        continue
    W, xe = walls[(row["kernel"], row["dx0"])]
    xx = xs[xs <= min(xe, xe0)]
    dr = (W.r(xx) - base.r(xx)) * RT_UM
    dth = np.degrees(np.arctan(W.r(xx, 1)) - np.arctan(base.r(xx, 1)))
    m = (xx >= 40) & (xx <= 94)
    row.update(wall_max_abs_dr_um=float(np.abs(dr).max()), x_of_max_dr=float(xx[np.argmax(np.abs(dr))]),
               wall_dr_exit_um=float(dr[-1]), wall_max_abs_dtheta_test_section_deg=float(np.abs(dth[m]).max()),
               exit_x_shift_rt=float(xe - xe0))
(C / "_band_ab").mkdir(exist_ok=True)
(C / "_band_ab/throat_moc_fix_probe.json").write_text(json.dumps({"problem": PROB, "rows": rows}, indent=1, ensure_ascii=False))
print("\n%-5s %-7s %-4s %10s %12s %12s %12s %12s %12s" % ("kern", "dx0", "ns", "dth1[deg]", "bump(no-mono)", "mono_dev1", "max|dr|[um]", "dr_exit[um]", "dth_TS[deg]"))
for r_ in rows:
    if "error" in r_:
        print(r_["kernel"], r_["dx0"], r_["n_start"], r_["error"]); continue
    print("%-5s %-7g %-4d %10.4f %12s %12s %12s %12s %12s" % (r_["kernel"], r_["dx0"], r_["n_start"], r_["dth1_deg"],
          f"{r_['bump_unconstrained']:.4f}" if "bump_unconstrained" in r_ else "", f"{r_['mono_first_point_dev_deg']:.4f}" if "mono_first_point_dev_deg" in r_ else "",
          f"{r_['wall_max_abs_dr_um']:.2f}" if "wall_max_abs_dr_um" in r_ else "", f"{r_['wall_dr_exit_um']:.2f}" if "wall_dr_exit_um" in r_ else "",
          f"{r_['wall_max_abs_dtheta_test_section_deg']:.4f}" if "wall_max_abs_dtheta_test_section_deg" in r_ else ""))
