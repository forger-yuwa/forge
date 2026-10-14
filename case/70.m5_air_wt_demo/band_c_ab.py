"""case/70: 帯位置の係数 c (y_b = c × 縁 y*) を変えたときの、帯の切り替わりへの強さと δ_E の偏りの測定
(plan tooling-nozzle-cooled-wall-ns-settling §5.1 #3、forge 0 step、読むだけ)。

run_0010 の 50k と 55k の場 (帯位置の切り替わりで x 2〜12 の δ_E が 5e-4 r_t 動いた 2 枚) を、Euler run_0012 に対して
同じ修正済みの抽出 (main の実装) に c = C で通し、使用域 [0, x_e] の 50k→55k の差と、c = 1.25 (既定) からの偏りを出す。

usage (case dir、design/.venv-opt):  python band_c_ab.py C STEP OUT_DIR
"""
import json
import sys
from pathlib import Path

import numpy as np

C_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(C_DIR.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge  # noqa: E402

c, step, out = float(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
ns, eu = (C_DIR / "run_0010_ns_c2b_final").resolve(), (C_DIR / "run_0012_euler_c2b").resolve()
orig = DS._load_structured


def load(run_dir):
    run_dir = Path(run_dir)
    if run_dir.resolve() != ns:
        return orig(run_dir)
    import h5py
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
od = out / f"c{c:g}" / str(step)
extract_and_merge(ns, eu, band_select="edge", edge_c=c, out_dir=od)
print(od)
