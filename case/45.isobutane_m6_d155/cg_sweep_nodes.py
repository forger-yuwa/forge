"""plan tooling-nozzle-core-grid §4.42 M-d: FORGE_LINE_DUMP_NODES に渡す節点 (1 ライン 1 節点) を選ぶ。
縮流部 x ∈ [−8, −2] r_t で、壁の x が 0.05 r_t ごとの目標に最も近い列と、島の x ∈ [−6, −3.25] の全列。
各列の壁から 80 節点目 (島の帯の中、ラインに必ず乗る) を渡す。節点の並びは cg_mode_shape.py と同じ (列 i の節点は i·nj + jj、壁は jj = nj − 1)。
usage: python3 cg_sweep_nodes.py <run> → 標準出力にカンマ区切りの節点番号"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run = sys.argv[1]
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
src = {"Gref": "run_0619_cg_gref_b0_ctl", "Grefc": "run_0625_cg_grefc"}
seed = next(HERE / s / "res_60000.h5" for s in src.values() if json.loads((HERE / s / "prepare_info.json").read_text())["mesh"]["ni"] == ni)
with h5py.File(seed, "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
nj = len(xy) // ni
if ni * nj != len(xy):
    raise SystemExit("節点数が ni の倍数でない")
xw = xy[np.arange(ni) * nj + nj - 1, 0] / S
cols = set()
for t in np.arange(-7.95, -2.0, 0.05):
    cols.add(int(np.argmin(np.abs(xw - t))))
cols |= set(np.where((xw >= -6.0) & (xw <= -3.25))[0].tolist())
cols = sorted(c for c in cols if -8.0 <= xw[c] <= -2.0)
print(",".join(str(c * nj + nj - 1 - 80) for c in cols))
print(f"{len(cols)} 列 (x {xw[cols[0]]:.3f}〜{xw[cols[-1]]:.3f})", file=sys.stderr)
