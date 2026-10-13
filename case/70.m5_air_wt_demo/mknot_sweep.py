"""case/70: 長さを固定 (from_length) して M_knot だけを振り、設計の量を比べる (設計計算だけ、CFD なし)。

見る量 (case/42.isobutane_wt/compare_axislaw_D.py と同じ診断):
  - 壁での余裕 min(μ_w − θ_w) [°] (`wall_margin_diagnostics`)
  - 網の余裕: 特性線の網の内部の三角形の最小 sin 角・向きの反転の数 (`characteristic_topology_diagnostics`)
  - 最大壁角 θ_max [°]、設計壁の曲率の最大 |κ| と曲率の変化率の最大 |dκ/ds|・∫(dκ/ds)² ds (スロート窓の外 x > 0.25 r_t)
  - L_c・x_K・出口 x_F (長さは from_length で固定)
usage (case dir、design/.venv-opt):  python mknot_sweep.py PROBLEM.yaml OUT.json [--n-axis N] M_K [M_K ...]
  --n-axis N: geometry.n_axis_inv を N に上書き (解像度の A/B; plan tooling-nozzle-mknot-auto の諮問 2026-10-13)
"""
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))


def one(args):
    prob_path, MK, n_over = args
    import forge_design.evaluate.runner_axismach as RA
    from forge_design.geometry.moc_diagnostics import (characteristic_topology_diagnostics, wall_curvature_diagnostics,
                                                        wall_margin_diagnostics)
    from forge_design.probdef import load_problem
    raw = yaml.safe_load(Path(prob_path).read_text())
    raw.setdefault("geometry", {})["M_knot"] = float(MK)
    if n_over:
        raw["geometry"]["n_axis_inv"] = int(n_over)
    raw.pop("_card", None)
    cap = {}
    orig = RA.inverse_design

    def wrapped(*a, **k):
        r = orig(*a, **k)
        cap["res"], cap["n_axis"] = r, k.get("n_axis")
        return r

    RA.inverse_design = wrapped
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", dir=C, delete=False) as f:
            yaml.safe_dump(raw, f, allow_unicode=True, sort_keys=False)
            tmp = f.name
        p = load_problem(tmp)
        d = RA.design_chain(p)
        w = d["wall_inv"]
        x_A = float(d["x_A"])
        n_axis = int(cap.get("n_axis") or p.geometry.get("n_axis_inv", 500))
        topo = characteristic_topology_diagnostics(cap["res"]["levels"], n_axis - 1, n_crossing_sample=20,
                                                   seam_i_pad=90, seam_k_max=110, seam_x_max=x_A + 1.0, wall=w[:, :2])
        mg = wall_margin_diagnostics(w)
        wc = wall_curvature_diagnostics(d["wall"], 0.0, float(d["wall"].x_e), x_split=0.25)
        return {"M_K": float(MK), "ok": True, "L_c": float(d["L_c"]), "x_K": float(d["x_K"]), "x_F": float(d["qa"]["x_F"]),
                "theta_max_deg": float(d["qa"]["theta_max_deg"]), "margin_deg": mg["min_mu_minus_theta_deg"], "x_at_min_margin": mg["x_at_min_margin"],
                "min_sin_angle": float(topo.get("min_sin_angle_interior", float("nan"))),
                "n_flip_interior": int(topo.get("n_orientation_flip_interior", -1)),
                "n_crossings_sampled": int(topo.get("n_crossings_sampled", -1)), "pileup_pairs_interior": int(topo.get("pileup_pairs_interior", -1)),
                "n_axis": n_axis, "qa_violations": list(d["qa"].get("violations") or []),
                "max_kappa_down": float(wc.get("downstream_max_abs_kappa", float("nan"))),
                "max_dkappa_ds_down": float(wc.get("downstream_max_abs_dkappa_ds", float("nan"))),
                "J_dkappa_down": float(wc.get("downstream_J_dkappa_ds2", float("nan"))),
                "x_at_max_dkappa_down": float(wc.get("downstream_x_at_max_dkappa_ds", float("nan")))}
    except Exception as e:  # noqa: BLE001 — 失敗も結果として残す
        return {"M_K": float(MK), "ok": False, "error": f"{type(e).__name__}: {str(e)[:300]}"}
    finally:
        RA.inverse_design = orig
        try:
            Path(tmp).unlink()
        except Exception:  # noqa: BLE001
            pass


if __name__ == "__main__":
    argv = list(sys.argv[1:])
    n_over = None
    if "--n-axis" in argv:
        k = argv.index("--n-axis"); n_over = int(argv[k + 1]); del argv[k:k + 2]
    prob, out = argv[0], Path(argv[1])
    mks = [float(v) for v in argv[2:]]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(one, [(prob, m, n_over) for m in mks]))
    out.write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(json.dumps(r))
