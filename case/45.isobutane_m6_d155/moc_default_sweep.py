"""MOC のキーの既定の切り替え (2026-10-07: legacy + fixed2 → analytic + converge、plans/accepted/discretization-moc-axis-limit-and-corrector.md
§9) の波及の確認 (CFD 0 step)。MOC のキーを書いていない axis-Mach の問題 YAML (case/41〜45) を全部、設計チェーンだけで 2 通りに流す。

  A: moc_axis_limit: legacy・moc_corrector: fixed2 を明示 (切り替え前の設計と同じ)
  B: キー無し (新しい既定 analytic + converge)

記録: B の MOC のゲート (applicable・pass・reasons)、A と B の設計壁の差 (共通の x の範囲で max|Δr| [r_t]・max|Δθ| [°]・x_F の差)、例外。
CFD でピン止めした初期線の run (initial_line_run) は手元の主ツリー (/home/sano/work/forge/case/<case>/) から読む。無ければ skip。

usage: design/.venv-opt/bin/python moc_default_sweep.py [--np 4] [--only case/41.wind_tunnel_design] → _band_ab/moc_default_sweep.json
"""
import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAIN = Path("/home/sano/work/forge")
sys.path.insert(0, str(REPO / "design"))


def targets(only=None):
    out = []
    for p in sorted(REPO.glob("case/*/problem*.yaml")):
        if only and not str(p.relative_to(REPO)).startswith(only):
            continue
        t = p.read_text()
        if "type: wind_tunnel_axisym_axismach" not in t or "moc_axis_limit" in t or "moc_corrector" in t:
            continue
        out.append(p)
    return out


def _load(path: Path):
    from forge_design.evaluate.runner_axismach import load_problem
    p = load_problem(path)
    g = p.geometry
    run = g.get("initial_line_run")
    if run and not Path(run).is_absolute():
        cand = [path.parent / run, MAIN / path.parent.relative_to(REPO) / run]
        hit = next((c for c in cand if c.exists()), None)
        if hit is None:
            raise FileNotFoundError(f"initial_line_run {run} が手元に無い")
        g["initial_line_run"] = str(hit)
    return p


def one(path_str: str) -> dict:
    from forge_design.evaluate.runner_axismach import design_chain
    path = Path(path_str)
    rec = {"problem": str(path.relative_to(REPO))}
    t0 = time.time()
    try:
        pa = _load(path)
        pa.geometry["moc_axis_limit"], pa.geometry["moc_corrector"] = "legacy", "fixed2"
        da = design_chain(pa)
        pb = _load(path)
        db = design_chain(pb)
        g = (db.get("moc") or {}).get("gate") or {}
        rec["gate"] = {"applicable": g.get("applicable"), "pass": g.get("pass"), "reasons": g.get("reasons")}
        wa, wb = np.asarray(da["wall_inv"], float), np.asarray(db["wall_inv"], float)
        lo, hi = max(wa[0, 0], wb[0, 0]), min(wa[-1, 0], wb[-1, 0])
        x = np.linspace(lo, hi, 4001)
        ra, rb = np.interp(x, wa[:, 0], wa[:, 1]), np.interp(x, wb[:, 0], wb[:, 1])
        ta, tb = np.interp(x, wa[:, 0], wa[:, 2]), np.interp(x, wb[:, 0], wb[:, 2])
        rec["wall_diff"] = {"max_abs_dr_rt": float(np.abs(rb - ra).max()), "x_at_max": float(x[np.argmax(np.abs(rb - ra))]),
                            "max_abs_dtheta_deg": float(np.degrees(np.abs(tb - ta).max())),
                            "x_F": [float(wa[-1, 0]), float(wb[-1, 0])], "r_F": [float(wa[-1, 1]), float(wb[-1, 1])],
                            "dr_exit_rel": float(wb[-1, 1] / wa[-1, 1] - 1.0)}
        rec["status"] = "ok"
    except FileNotFoundError as e:
        rec["status"] = "skip"
        rec["reason"] = str(e)
    except Exception as e:  # noqa: BLE001 — 例外は記録して続ける
        rec["status"] = "error"
        rec["reason"] = f"{type(e).__name__}: {e}"
        rec["trace_tail"] = traceback.format_exc()[-800:]
    rec["sec"] = round(time.time() - t0, 1)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    ts = targets(a.only)
    print(f"[sweep] {len(ts)} 問題", flush=True)
    rows = []
    with ProcessPoolExecutor(a.np) as ex:
        for r in ex.map(one, [str(t) for t in ts]):
            rows.append(r)
            g = r.get("gate") or {}
            wd = r.get("wall_diff") or {}
            print(f"[sweep] {r['status']:5s} {r['problem']}  gate {g.get('pass')}  max|Δr| {wd.get('max_abs_dr_rt')}  {r.get('reason', '')[:80]}  ({r['sec']} s)", flush=True)
    summ = {"n": len(rows), "by_status": {s: sum(r["status"] == s for r in rows) for s in ("ok", "skip", "error")},
            "gate_fail": [r["problem"] for r in rows if r["status"] == "ok" and (r["gate"]["pass"] is not True)],
            "max_abs_dr_rt_over_all": max((r["wall_diff"]["max_abs_dr_rt"] for r in rows if r["status"] == "ok"), default=None)}
    out = HERE / "_band_ab" / "moc_default_sweep.json"
    out.write_text(json.dumps({"summary": summ, "rows": rows}, indent=1, ensure_ascii=False))
    print(json.dumps(summ, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
