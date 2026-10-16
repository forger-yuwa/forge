#!/usr/bin/env python3
r"""case/71 run_0002_ns の最終の場 (既定 res_183000.h5) の非有限・物理性の検査 (plan tooling-nozzle-card-local-mesh §5.1 #3b)。

- VALUE の全ての浮動小数点の量の非有限の数
- ρ・P・T > 0 (最小値と 0 以下の数)、k・ω ≥ 0
- 静圧 ≤ 全圧: 全温・全圧は `solver_density_cuda/tools/total_quantities.py` の `total_state` (VALUE/h0 から逆算) で作る
  (スクリプトで T + u²/2c_p を組まない; AGENTS.md「出力と後処理の原則」)。P > P0 の数と P/P0 の最大、P0 の最大と入口の Pt、
  T0 の範囲と入口の Tt (冷却壁なので T0 > Tt は境界層の全温の行き過ぎか数値の誤差; 数だけ記録する)。

使い方 (AWS A): python3 final_field_check.py --tools ~/forge-cardmesh/solver_density_cuda/tools --run <run> [--step 183000] --out-json <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--step", type=int, default=183000)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.tools).resolve()))
    from total_quantities import total_state

    run = Path(a.run).resolve()
    res = run / f"res_{a.step}.h5"
    with h5py.File(res, "r") as h:
        V = {k: np.array(d[:]) for k, d in h["VALUE"].items()}
    nonfinite = {k: int(np.count_nonzero(~np.isfinite(v))) for k, v in V.items() if v.dtype.kind == "f"}
    pos = {}
    for k in ("ro", "P", "T", "k", "omega", "roK", "roOmega"):
        v = V[k]
        pos[k] = {"min": float(v.min()), "max": float(v.max()), "n_le_0": int(np.count_nonzero(v <= 0)),
                  "n_lt_0": int(np.count_nonzero(v < 0))}
    st = total_state(str(run), str(res))
    T0, P0 = np.asarray(st["T0"], float), np.asarray(st["P0"], float)
    P, T = V["P"].astype(float), V["T"].astype(float)
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    Pt = float(bc["inlet"]["floats"]["Pt"]); Tt = float(bc["inlet"]["floats"]["Tt"])
    ratio = P / P0
    out = {
        "run": str(run), "res": res.name, "n_nodes": int(P.size),
        "nonfinite_total": int(sum(nonfinite.values())), "nonfinite_by_var": nonfinite,
        "positivity": pos,
        "total_quantities": {"tool": "solver_density_cuda/tools/total_quantities.py total_state", "method": st["method"],
                             "h0_includes_k": st["includes_k"], "species_source": st.get("species_source"), "warn": st.get("warn"),
                             "T0": {"min": float(T0.min()), "median": float(np.median(T0)), "max": float(T0.max())},
                             "P0": {"min": float(P0.min()), "median": float(np.median(P0)), "max": float(P0.max())},
                             "nonfinite_T0": int(np.count_nonzero(~np.isfinite(T0))), "nonfinite_P0": int(np.count_nonzero(~np.isfinite(P0)))},
        "static_le_total": {"n_P_gt_P0": int(np.count_nonzero(P > P0)), "max_P_over_P0": float(ratio.max()),
                            "n_P_gt_P0_rel1e-9": int(np.count_nonzero(P > P0 * (1 + 1e-9))),
                            "n_T_gt_T0": int(np.count_nonzero(T > T0))},
        "inlet": {"Pt": Pt, "Tt": Tt, "P0_max_over_Pt": float(P0.max() / Pt), "n_P0_gt_Pt": int(np.count_nonzero(P0 > Pt)),
                  "n_P0_gt_Pt_rel1e-6": int(np.count_nonzero(P0 > Pt * (1 + 1e-6))),
                  "T0_max_minus_Tt": float(T0.max() - Tt), "n_T0_gt_Tt_plus_1K": int(np.count_nonzero(T0 > Tt + 1.0))},
    }
    # 位置 (node の構造格子: index = i·nj + j、j = 0 が軸・j = nj − 1 が壁; prepare_info.json の mesh.ni・nj)
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.array(h["/MESH/COORD"][:]).reshape(-1, 3)
    umag = np.hypot(V["Ux"], V["Uy"])

    def where(mask):
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            return {"n": 0}
        i, j = idx // nj, idx % nj
        jw = nj - 1 - j
        return {"n": int(idx.size), "x_m": [float(xy[idx, 0].min()), float(xy[idx, 0].max())],
                "i": [int(i.min()), int(i.max())], "j_from_wall": {"min": int(jw.min()), "max": int(jw.max()),
                                                                   "n_wall_row": int(np.count_nonzero(jw == 0))},
                "n_axis_row": int(np.count_nonzero(j == 0)), "umag": [float(umag[idx].min()), float(umag[idx].max())],
                "T": [float(T[idx].min()), float(T[idx].max())]}
    out["where"] = {"note": "index = i·nj + j (j = 0 軸、j = nj − 1 壁)", "ni": ni, "nj": nj,
                    "P_gt_P0": where(P > P0), "P0_gt_Pt_rel1e-6": where(P0 > Pt * (1 + 1e-6)),
                    "T0_gt_Tt_plus_1K": where(T0 > Tt + 1.0)}
    q = P0 / Pt
    out["inlet"]["P0_over_Pt_quantiles"] = {p: float(np.percentile(q, p)) for p in (50, 90, 99, 99.9, 100)}
    out["verdict"] = ("OK (非有限 0、ρ・P・T > 0、静圧 ≤ 全圧)"
                      if out["nonfinite_total"] == 0 and all(pos[k]["n_le_0"] == 0 for k in ("ro", "P", "T"))
                      and out["static_le_total"]["n_P_gt_P0"] == 0 else "NG (欄を見る)")
    Path(a.out_json).write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
