"""R4 (codex M6): 固定形状 (smoke 設計, 加速点, 3D Euler node, ext_top あり) で遠方境界・出口距離・格子の感度を測る。
各変種の YAML を run_dir に書き、runner_sern3d を順に回して metrics.json のノズル力係数を表にする。
判定: base に対する |ΔC_T| < 0.002 (plan §6 の許容; 最適化の差より小さい) を領域独立の目安にする。

  python3 case/46.sern_design/r4_domain_study.py <out_root>   (例 run_0102_r4_domain)
"""
import json, re, subprocess, sys
from pathlib import Path
C = Path(__file__).resolve().parent; D = C.parent.parent / "design"
PY = "/home/sano/work/forge/design/.venv-opt/bin/python"
out_root = Path(sys.argv[1] if len(sys.argv) > 1 else C / "run_0102_r4_domain").resolve(); out_root.mkdir(parents=True, exist_ok=True)
base = (C / "problem_r2_euler_node_3d_accel_vehicle.yaml").read_text()
# ext_top (2D の mesh: キー) を 3D にも付ける
base = base.replace("  first_wall_frac: 0.004\n", "  first_wall_frac: 0.004\n  ext_top: 1\n  top_depth: 2.0\n  nj_ext_top: 33\n  vehicle_taper: 0.35\n  vehicle_wedge_deg: 3.0\n  vehicle_clearance: 0.02\n  first_top_frac: 0.02\n", 1)
assert "ext_top: 1" in base
def variant(name, **kw):
    s = base
    for k, v in kw.items():
        s2 = re.sub(rf"^(\s*){k}: [\d.]+.*$", lambda m: f"{m.group(1)}{k}: {v}", s, count=1, flags=re.M)
        assert s2 != s, (name, k); s = s2
    return s
VARIANTS = [("base", {}), ("Zext3", {"Z_ext": 3.0, "nz_out": 21}), ("xout4", {"x_out_extra": 4.0, "ni_plume": 150}),
            ("bot15", {"bot_depth": 1.5}), ("top4", {"top_depth": 4.0, "nj_ext_top": 45}),
            ("fine", {"ni_noz": 75, "ni_plume": 138, "nj_top": 61, "nj_bot": 39, "nz_in": 31, "nz_out": 19, "nj_ext_top": 41})]
rows = []
for name, kw in VARIANTS:
    rd = out_root / name
    if not (rd / "metrics.json").exists():
        yml = out_root / f"problem_{name}.yaml"; yml.write_text(variant(name, **kw))
        if rd.exists():
            subprocess.run(["rm", "-rf", str(rd)])
        log = out_root / f"{name}.log"
        with open(log, "w") as fh:
            subprocess.run([PY, "-m", "forge_design.evaluate.runner_sern3d", str(yml), str(rd), "--op", "accel", "--soft-steps", "1500"],
                           cwd=D, env={**__import__("os").environ, "PYTHONPATH": str(D)}, stdout=fh, stderr=subprocess.STDOUT)
    if (rd / "metrics.json").exists():
        m = json.loads((rd / "metrics.json").read_text()); g = m["gates"]; info = json.loads((rd / "prepare_info.json").read_text())
        rows.append({"name": name, "cells": info["mesh"]["cells"], "C_T": m["C_T"], "C_L": m["C_L"], "C_M": m["C_M"], "C_T_vehicle": m.get("C_T_vehicle"),
                     "C_T_vehicle_top": m.get("C_T_vehicle_top"), "gate": g["verdict"], "residual": g["residual"]["verdict"],
                     "steady": all(v["verdict"] == "STEADY" for v in g["steadiness"]["series"].values()), **kw})
    else:
        rows.append({"name": name, "gate": "NO RESULT", **kw})
    (out_root / "domain_study.json").write_text(json.dumps(rows, indent=1))
b = next((r for r in rows if r["name"] == "base" and "C_T" in r), None)
md = ["# R4 領域独立性 (3D Euler, 加速点, ext_top)", "", "| variant | cells | C_T | ΔC_T | C_L | ΔC_L | C_M | ΔC_M | C_T_veh | C_T_veh_top | gate | residual | steady |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
for r in rows:
    if "C_T" in r and b:
        md.append(f"| {r['name']} | {r['cells']} | {r['C_T']:.5f} | {r['C_T']-b['C_T']:+.5f} | {r['C_L']:.4f} | {r['C_L']-b['C_L']:+.4f} | {r['C_M']:.3f} | {r['C_M']-b['C_M']:+.3f} | {r['C_T_vehicle']:.4f} | {r['C_T_vehicle_top']:.4f} | {r['gate']} | {r['residual']} | {r['steady']} |")
    else:
        md.append(f"| {r['name']} | | | | | | | | | | {r['gate']} | | |")
(out_root / "domain_study.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
