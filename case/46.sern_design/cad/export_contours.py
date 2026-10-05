#!/usr/bin/env python3
"""接続模型 (hex_junction_model.py) に渡す MOC 輪郭を CSV に書き出す (plan tooling-sern-mesh-blocking §5.1 B4-1 / §4.14-3)。

生産設計 (既定 `problem_3d_prod_m6on_wallres_lswx08.yaml`) の逆設計輪郭を、runner_sern3d.prepare と同じ書式
(列 `x_m,y_m`、単位 m) で OUT_DIR/ramp_contour.csv と OUT_DIR/cowl_contour.csv に書く。
**ランプは `mesh.ramp_fillet` を適用した後の 2D 輪郭** (legacy `mesh_sern3d.generate_sern_mesh3d` と同じ式:
接点 x1 = −t, x2 = t cosθ, t = R tan(θ/2)、丸め区間は円弧 y = 1 + R − sqrt(R² − (x + t)²) と折れ線の大きい方)。
円弧区間は密な点列 (既定 400 点) にし、その後ろに ramp_xy の x > x2 の点を続ける。カウルは cowl_xy をそのまま書く。

usage (system python3。gmsh venv には yaml が無い):
  cd design && python3 ../case/46.sern_design/cad/export_contours.py OUT_DIR [--problem YAML] [--n-arc 400]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402


def fillet_ramp(rx, ry, R_f, n_arc):
    """ランプ角部の丸めを適用した輪郭 (/H)。mesh_sern3d.generate_sern_mesh3d の y_top と同じ式"""
    if R_f <= 0.0:
        return np.stack([rx, ry], 1), dict(R_f=0.0)
    th = float(np.arctan2(ry[1] - ry[0], rx[1] - rx[0])); t_f = R_f * np.tan(0.5 * th)
    xf1, xf2 = -t_f, t_f * np.cos(th)
    xa = np.linspace(xf1, xf2, n_arc)
    lin = np.where(xa < 0.0, 1.0, np.interp(xa, rx, ry))
    arc = 1.0 + R_f - np.sqrt(np.clip(R_f * R_f - (xa + t_f) ** 2, 0.0, None))
    ya = np.maximum(lin, arc)
    keep = rx > xf2 + 1e-12
    xy = np.vstack([np.stack([xa, ya], 1), np.stack([rx[keep], ry[keep]], 1)])
    return xy, dict(R_f=R_f, theta_r0_deg=float(np.degrees(th)), t_f=float(t_f), x_f1=float(xf1), x_f2=float(xf2), n_arc=int(n_arc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--problem", default=str(ROOT / "case/46.sern_design/problem_3d_prod_m6on_wallres_lswx08.yaml"))
    ap.add_argument("--n-arc", type=int, default=400)
    a = ap.parse_args()
    p = load_problem(a.problem); d0 = R2.design_snapshot(p)
    _, d, _, _ = R2.design_from_problem(p, design=d0)
    H = float(p.spec["H_m"]); R_f = float(p.raw.get("mesh", {}).get("ramp_fillet", 0.0))
    ramp, meta = fillet_ramp(d.ramp_xy[:, 0], d.ramp_xy[:, 1], R_f, a.n_arc)
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    np.savetxt(out / "ramp_contour.csv", ramp * H, delimiter=",", header="x_m,y_m", comments="", fmt="%.17g")
    np.savetxt(out / "cowl_contour.csv", d.cowl_xy * H, delimiter=",", header="x_m,y_m", comments="", fmt="%.17g")
    meta.update(problem=str(a.problem), H_m=H, L_cowl_H=float(d.cowl_xy[-1, 0]), y_te_H=float(d.cowl_xy[-1, 1]),
                L_ramp_H=float(d.L_ramp), n_ramp=int(len(ramp)), n_cowl=int(len(d.cowl_xy)))
    json.dump(meta, open(out / "contours_meta.json", "w"), indent=1, ensure_ascii=False)
    print(json.dumps(meta, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
