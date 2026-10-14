"""case/70: edge の帯位置の切り替わりの断面で、梯子の R(y) = δ_r(1.25 y)/δ_r(y) − 1 の形を見る
(plan tooling-nozzle-cooled-wall-ns-settling §5.1 #3、forge 0 step、読むだけ; 今の抽出器をそのまま呼ぶ)。
usage (case dir、design/.venv-opt):  python ladder_probe.py NS_RUN STEP EULER_RUN OUT.json
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402

ns, step, eu, out = (C / sys.argv[1]).resolve(), int(sys.argv[2]), (C / sys.argv[3]).resolve(), Path(sys.argv[4])
orig = DS._load_structured


def load(run_dir):
    run_dir = Path(run_dir)
    if run_dir.resolve() != ns:
        return orig(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = run_dir / f"res_{step}.h5"
    with h5py.File(res) as f:
        q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
    return dict(info=info, S=S, res=res.name, x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S, q=q.reshape(ni, nj))


DS._load_structured = load
d = DS.deltastar_from_core_matched_euler(ns, eu, band_select="edge", return_ladders=True)
e = d["edge_detail"]
rows = []
for k in range(len(e["x"])):
    y = np.asarray(e["ladders_y"][k], float); dd = np.asarray(e["ladders_d"][k], float)
    rows.append({"x": float(e["x"][k]), "y_star": float(e["y_star"][k]), "y_b": float(e["y_b"][k]), "y_ref": float(e["y_ref"][k]),
                 "kappa_med": float(e["kappa_med"][k]), "ladder_y": y.tolist(), "ladder_d": dd.tolist()})
out.write_text(json.dumps({"ns": str(ns), "step": step, "euler": str(eu), "rows": rows}))
print(out, len(rows))
