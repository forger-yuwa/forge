"""V9 (f): cplus × characteristic 以外の分岐 (streamline・flux の terminal_exit、cplus × lip) の出口・壁が変更の前後でビット一致か。"""
import sys, json, hashlib, tempfile
from pathlib import Path
import numpy as np, yaml
root = sys.argv[1]; sys.path.insert(0, root)
from forge_design.evaluate.runner_axismach import design_chain
from forge_design.probdef import load_problem
out = {}
for wm, em in (("streamline", "characteristic"), ("flux", "characteristic"), ("cplus", "lip"), ("cplus", "characteristic")):
    B = {"name": "inv", "type": "wind_tunnel_axisym_axismach", "gas": {"gamma": 1.4, "cp": 1004.5},
         "spec": {"Pt": 1.0e6, "Tt": 300.0, "M_design": 4.0, "r_throat": 0.05},
         "dv": {"L_c": {"value": 6.0}},
         "geometry": {"R": 2.0, "start_line": "throat_char", "wall_mode": wm, "exit_mode": em, "n_axis_inv": 300, "n_start": 41,
                      "Lc_mode": "explicit"},
         "mesh": {}, "evaluate": {}}
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "p.yaml"; p.write_text(yaml.safe_dump(B, sort_keys=False))
        try:
            d = design_chain(load_problem(p))
            out[f"{wm}/{em}"] = {"x_F": repr(float(d["exit"]["x_F"])), "wall_inv_sha": hashlib.sha256(np.ascontiguousarray(d["wall_inv"]).tobytes()).hexdigest()[:16]}
        except Exception as e:
            out[f"{wm}/{em}"] = {"err": repr(e)[:200]}
print(json.dumps(out, ensure_ascii=False))
