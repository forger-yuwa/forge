#!/usr/bin/env python3
"""V0-b の補足: (D) の成立する x_m の集合の細かさ (走査 0.025 r_t で取りこぼす細い成立区間があるか)。

V0-b (`run_v0b.py`) の走査範囲の確認 (0.01 r_t 刻み) で、`design_bell_c2` の走査 (0.025 r_t 刻み、141 点) では不成立の
3 点が 0.01 刻みでは成立し、成立 31 点のうち 3 点は成立区間の幅が 0.01 r_t 未満だった。ここでは 60 点すべてで x_m を
1.0〜3.0 r_t の 0.001 r_t 刻みで振り (V0-b の 0.01 刻みの確認で成立する x_m は全体で 1.08〜2.75 r_t)、成立区間の数・幅と、
`design_bell_c2` の判定 (OK / INFEASIBLE) との対応を記録する。判定ロジックは変えない (記録だけ)。

出力: `v0b_sliver_check.csv`。実行: design/.venv-opt の Python で `python case/68.bell_two_mode/v0_moc_gate/v0b_sliver_check.py`。
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.geometry.wall_bell_c2 import BellC2Spec, BellC2Wall, cone15_length, design_bell_c2  # noqa: E402

EPS, R_T = 50.0, 1.0
L_FRACS = (0.6, 0.7, 0.8, 0.9, 1.0)
THETA_MAX = (20.0, 25.0, 32.5, 40.0)
THETA_E = (5.0, 10.0, 15.0)
X_LO, X_HI, DX = 1.0, 3.0, 0.001


def main() -> int:
    spec = BellC2Spec(eps=EPS, R_t=R_T)
    Lc = cone15_length(EPS)
    xs = X_LO + DX * np.arange(int(round((X_HI - X_LO) / DX)) + 1)
    rows = []
    for fr in L_FRACS:
        L = fr * Lc
        for tm in THETA_MAX:
            for te in THETA_E:
                pid = f"D_L{int(round(fr * 100))}_tm{tm:g}_te{te:g}"
                d = design_bell_c2(L, te, tm, spec)
                ok = np.array([BellC2Wall(L, te, tm, spec=spec, x_m=float(x)).feasibility()["ok"] for x in xs])
                runs, i = [], 0
                while i < len(xs):
                    if ok[i]:
                        j = i
                        while j + 1 < len(xs) and ok[j + 1]:
                            j += 1
                        runs.append((float(xs[i]), float(xs[j]), (j - i + 1) * DX))
                        i = j + 1
                    else:
                        i += 1
                w_chosen = None
                if d.ok:
                    for (a, b, w) in runs:
                        if a - DX <= d.x_m <= b + DX:
                            w_chosen = w
                rows.append({"point": pid, "L_frac": fr, "theta_max_deg": tm, "theta_e_deg": te,
                             "design_status": d.status, "x_m": d.x_m, "n_fine_intervals": len(runs),
                             "fine_measure_rt": float(ok.sum() * DX),
                             "fine_width_min_rt": min((w for _, _, w in runs), default=None),
                             "fine_width_max_rt": max((w for _, _, w in runs), default=None),
                             "chosen_interval_width_rt": w_chosen,
                             "fine_intervals": ";".join(f"{a:.3f}-{b:.3f}" for a, b, _ in runs[:10])})
                r = rows[-1]
                print(f"{pid:20s} {d.status:10s} 細かい成立区間 {r['n_fine_intervals']} 個・幅の計 {r['fine_measure_rt']:.4f}・"
                      f"選んだ区間の幅 {w_chosen}  {r['fine_intervals']}", flush=True)
    with open(HERE / "v0b_sliver_check.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    miss = [r["point"] for r in rows if r["design_status"] != "OK" and r["n_fine_intervals"] > 0]
    narrow = [r["point"] for r in rows if r["design_status"] == "OK" and (r["chosen_interval_width_rt"] or 0) < 0.025]
    print(f"\n0.025 刻みの走査で INFEASIBLE だが 0.001 刻みで成立区間がある点: {len(miss)} {miss}")
    print(f"OK の点のうち選んだ成立区間の幅 < 0.025 r_t (走査の刻み): {len(narrow)} {narrow}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
