"""edge の縁 y* の推定に「持続」(交点の後 K 段続けて |R| < ε) を要求したときの、帯位置の切り替わりと δ_E への効き
(plan tooling-nozzle-cooled-wall-ns-settling §5.1 #3、forge 0 step、読むだけ; 測定用の差し替えで、実装ではない)。

K = 1 は今の実装 (最初に |R| < ε を満たす段) と同じ。NS_RUN の res_<STEP>.h5 (STEP = 0 なら最新) を Euler に対して抽出し、
OUT_DIR/K<K>/<tag>/ に delta_r_next.csv・delta_r_edge.csv を書く。

usage (design/.venv-opt):  python band_persist_ab.py NS_RUN EULER_RUN STEP K OUT_DIR TAG
"""
import inspect
import json
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge  # noqa: E402

ns, eu, step, K, out, tag = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), int(sys.argv[3]), int(sys.argv[4]), Path(sys.argv[5]), sys.argv[6]
if step > 0:
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
if K > 1:
    src = inspect.getsource(DS.edge_band_positions)
    old = "        j = int(np.argmax(Rabs < 0.0)) if np.any(Rabs < 0.0) else -1\n"
    assert src.count(old) == 1, "edge_band_positions の交点の行が見つからない"
    new = ("        _ok = Rabs < 0.0\n"
           f"        _K = {K}\n"
           "        _c = [jj for jj in range(len(_ok) - _K + 1) if bool(np.all(_ok[jj:jj + _K]))]\n"
           "        j = _c[0] if _c else -1\n")
    g = dict(DS.__dict__)
    exec(compile(src.replace(old, new), "patched", "exec"), g)
    DS.edge_band_positions = g["edge_band_positions"]
od = out / f"K{K}" / tag
extract_and_merge(ns, eu, band_select="edge", out_dir=od)
print(od)
