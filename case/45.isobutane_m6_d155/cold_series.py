"""冷却壁の腕の状態の水準の時系列 (plan time_integration-implicit-thermal-jacobian §6.0 の run_0224・切り戻しの判定、
plan tooling-nozzle-isothermal-wall-chain §5.1 #27)。各スナップショット res_<n>.h5 から:
  - 欠損 = 2π Σ res_ro [kg/s] (extraFields に res_ro がある run だけ。run_0223 の判定と同じ定義)
  - θ_r・δ_loc (x = 40・70・94、cold_xcheck.reduce_fields と同じ式・同じ帯の外縁)、Q_w [W]、入口・出口の ṁ (列の折れ線)
状態の水準 (事前登録): 欠損 ≤ 0.1 kg/s かつ θ_r のドリフト ≤ 0.05 %/2 万 step。ドリフトは末尾 2 万 step の両端の差を 2 万 step あたりに換算する。
usage (AWS の case dir): python3 cold_series.py <run> [<run> ...] [--every N] [--from S] → _band_ab/cold_pair/series_<run>.json
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402

DRIFT_WIN = 20000


def snapshot(run: Path, st: int, yb_x, yb):
    xy, ro, ux, uy, T, tke = XC.load_forge(run, st)
    ni, nj, S = XC.mesh_info(len(ro), run)
    o = XC.reduce_fields(xy, ro, ux, uy, T, tke, False, ni, nj, S, yb_x, yb)
    rec = {"step": st, "nonfinite": int(sum(np.count_nonzero(~np.isfinite(a)) for a in (ro, ux, uy, T, tke))),
           "Q_w": float(o["Q_w"]), "mdot_in": float(o["mdot_in"]), "mdot_out": float(o["mdot_out"])}
    for xs in XC.XS:
        i = int(np.argmin(np.abs(o["x"] - xs)))
        rec[f"theta_r_{int(xs)}"] = float(o["theta_r"][i]); rec[f"delta_loc_{int(xs)}"] = float(o["delta_loc"][i])
    with h5py.File(run / f"res_{st}.h5", "r") as h:
        rec["deficit"] = float(np.sum(h["VALUE/res_ro"][:]) * 2 * math.pi) if "VALUE/res_ro" in h else None
    return rec


def drift(rows, key):
    """末尾 DRIFT_WIN step の両端の相対差を DRIFT_WIN あたりに換算 [%]。窓が足りなければ None。"""
    st = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if st[-1] - st[0] < DRIFT_WIN:
        return None
    j = int(np.searchsorted(st, st[-1] - DRIFT_WIN))
    return float(100.0 * (v[-1] / v[j] - 1.0) * DRIFT_WIN / (st[-1] - st[j]))


def main(run: Path, every: int, start: int):
    yb_x, yb = XC.common_yb()
    steps = sorted(int(m.group(1)) for p in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", p.name)))
    steps = [s for s in steps if s >= start and (every <= 0 or s % every == 0)]
    rows = []
    for st in steps:
        rows.append(snapshot(run, st, yb_x, yb))
        r = rows[-1]
        print(f"[series] {run.name} {st}: 欠損 {r['deficit'] if r['deficit'] is None else round(r['deficit'], 4)} kg/s、"
              f"θ_r(40/70/94) {r['theta_r_40']:.5e} {r['theta_r_70']:.5e} {r['theta_r_94']:.5e}、Q_w {r['Q_w'] / 1e6:.4f} MW、"
              f"ṁ 入口 {r['mdot_in']:.4f} 出口 {r['mdot_out']:.4f}、非有限 {r['nonfinite']}", flush=True)
    out = {"run": run.name, "rows": rows,
           "drift_pct_per_20k": {k: drift(rows, k) for k in ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")} if rows else {},
           "level_rule": "欠損 ≤ 0.1 kg/s かつ θ_r のドリフト ≤ 0.05 %/2 万 step (plan time_integration-implicit-thermal-jacobian §6.0)"}
    p = XC.OUTD / f"series_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(f"[series] → {p}  ドリフト [%/2 万 step]: {out['drift_pct_per_20k']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+"); ap.add_argument("--every", type=int, default=0); ap.add_argument("--from", dest="start", type=int, default=0)
    a = ap.parse_args()
    for r in a.runs:
        main(HERE / r, a.every, a.start)
