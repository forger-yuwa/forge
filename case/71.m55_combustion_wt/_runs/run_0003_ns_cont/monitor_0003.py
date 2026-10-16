#!/usr/bin/env python3
r"""case/71 run_0003_ns_cont の監視の量の時系列 (forge 0 step)。plan tooling-nozzle-card-local-mesh §5.1 #3d の「監視」:
S の h0・P0、壁の積分熱流、入口・出口のエネルギー流束、入口の列の静温の超過。

- 全温・全圧は `solver_density_cuda/tools/total_quantities.py` の `total_state` (VALUE/h0 から逆算)。
- S = run_0002_ns の res_183000 で P0 > Pt (1 + 1e−6) の節点の集合 (#3c の S と同じ定義、節点の index で固定)。
  各時点で S の平均の ΔP0/Pt・h0 − h0_in、その時点の超過の点数と最大。
- エネルギー流束: 入口・出口の列の折れ線を通る ∫ ρ h0 (U·n) 2πr ds (`stop_rule.mass_flows` と同じ式に h0 を掛けたもの) [W]。
- 壁の積分熱流: res_wall_3_<n>.h5 の qwall を壁の節点の x の順の折れ線で 2πr ds に台形で積む [W] (符号はソルバの規約のまま)。
- 入口の列 (i = 0) の静温 T > Tt の点数と最大の T − Tt。
- 収支: (E_in − E_out + Q_wall)/E_in と (E_in − E_out − Q_wall)/E_in の両方を出す (qwall の符号の規約をここで決めない)。

使い方 (AWS A): python3 monitor_0003.py --tools ~/forge-cardmesh/solver_density_cuda/tools --design ~/forge-cardmesh/design \
    --run <run_0003_ns_cont> --ref-run <run_0002_ns> --steps 2000,20000,...,200000 --out-csv <csv> --out-json <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

THR = 1e-6
H0_IN = 1147573.5007131016   # h_mix(Tt, 入口の Y) (#3c p0_ab.json・inlet_face_check.json と同じ値)


def col_flux(X, R, W, UX, UY, S: float) -> tuple:
    """入口・出口の列の折れ線を通る ∫ W (U·n) 2πr ds (`stop_rule.mass_flows` の式で ρ を W に替えたもの)。"""
    out = []
    for i in (0, X.shape[0] - 1):
        xx = X[i] * S; rr = R[i] * S
        a = W[i] * UX[i] * rr; b = W[i] * UY[i] * rr
        out.append(2 * np.pi * float(np.sum(0.5 * (a[1:] + a[:-1]) * np.diff(rr) - 0.5 * (b[1:] + b[:-1]) * np.diff(xx))))
    return out[0], out[1]


def wall_heat(path: Path) -> dict:
    with h5py.File(path, "r") as h:
        xyz = np.array(h["MESH/COORD"][:]).reshape(-1, 3).astype(float)
        q = np.array(h["VALUE/qwall"][:]).astype(float)
    o = np.argsort(xyz[:, 0], kind="stable")
    x, r, qq = xyz[o, 0], xyz[o, 1], q[o]
    ds = np.hypot(np.diff(x), np.diff(r))
    f = qq * r
    Q = 2 * np.pi * float(np.sum(0.5 * (f[1:] + f[:-1]) * ds))
    return {"Q_wall": Q, "qwall_min": float(qq.min()), "qwall_max": float(qq.max()), "n_wall_nodes": int(qq.size)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--ref-run", required=True)
    ap.add_argument("--ref-step", type=int, default=183000)
    ap.add_argument("--steps", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.tools).resolve()))
    sys.path.insert(0, str(Path(a.design).resolve()))
    from total_quantities import total_state
    from forge_design.evaluate.stop_rule import RunGeometry, read_values, nonfinite_count

    run, ref = Path(a.run).resolve(), Path(a.ref_run).resolve()
    geo = RunGeometry(run)
    ni, nj = geo.ni, geo.nj
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    Pt = float(bc["inlet"]["floats"]["Pt"]); Tt = float(bc["inlet"]["floats"]["Tt"])

    st_ref = total_state(str(ref), str(ref / f"res_{a.ref_step}.h5"))
    S = np.asarray(st_ref["P0"], float) > Pt * (1 + THR)
    with h5py.File(ref / "nozzle.h5", "r") as h:
        same_mesh = bool(np.array_equal(np.array(h["/MESH/COORD"][:]), np.array(h5py.File(run / "nozzle.h5", "r")["/MESH/COORD"][:])))

    rows = []
    jj = np.tile(np.arange(nj), ni); ii = np.repeat(np.arange(ni), nj)
    for n in [int(s) for s in a.steps.split(",")]:
        res = run / f"res_{n}.h5"
        V = read_values(res)
        st = total_state(str(run), str(res))
        if int(st["includes_k"]) != 0:
            raise SystemExit(f"{res}: h0_includes_k={st['includes_k']} (この台本は k を含まない h0 だけを扱う)")
        P0 = np.asarray(st["P0"], float); T0 = np.asarray(st["T0"], float)
        h0 = V["h0"].astype(float); T = V["T"].astype(float); P = V["P"].astype(float)
        sh = (ni, nj)
        Ein, Eout = col_flux(geo.X, geo.R, (V["ro"].astype(float) * h0).reshape(sh), V["Ux"].astype(float).reshape(sh),
                             V["Uy"].astype(float).reshape(sh), geo.S)
        Min, Mout = col_flux(geo.X, geo.R, V["ro"].astype(float).reshape(sh), V["Ux"].astype(float).reshape(sh),
                             V["Uy"].astype(float).reshape(sh), geo.S)
        wh = wall_heat(run / f"res_wall_3_{n}.h5")
        ex = P0 > Pt * (1 + THR)
        inl = ii == 0
        row = {"step": n, "nonfinite": nonfinite_count(V),
               "n_P_gt_P0": int(np.count_nonzero(P > P0)),
               "n_P0_exceed": int(np.count_nonzero(ex)), "max_dP0_over_Pt": float(P0.max() / Pt - 1),
               "n_P0_exceed_axis_row": int(np.count_nonzero(ex & (jj == 0))),
               "S_mean_dP0_over_Pt": float(np.mean(P0[S] / Pt - 1)), "S_mean_h0_minus_h0in": float(np.mean(h0[S] - H0_IN)),
               "n_T0_gt_Tt_plus_1K": int(np.count_nonzero(T0 > Tt + 1.0)),
               "inlet_n_T_gt_Tt": int(np.count_nonzero(T[inl] > Tt)), "inlet_max_T_minus_Tt": float(np.max(T[inl] - Tt)),
               "mdot_in": Min, "mdot_out": Mout, "E_in": Ein, "E_out": Eout, "Q_wall": wh["Q_wall"],
               "bal_plus": (Ein - Eout + wh["Q_wall"]) / Ein, "bal_minus": (Ein - Eout - wh["Q_wall"]) / Ein,
               "qwall_min": wh["qwall_min"], "qwall_max": wh["qwall_max"]}
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)

    keys = list(rows[0].keys())
    with open(a.out_csv, "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(repr(r[k]) for k in keys) + "\n")
    Path(a.out_json).write_text(json.dumps({"run": str(run), "ref": f"{ref}/res_{a.ref_step}.h5", "same_mesh_coords": same_mesh,
                                            "S_n": int(S.sum()), "Pt": Pt, "Tt": Tt, "h0_in": H0_IN, "thr": THR,
                                            "rows": rows}, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
