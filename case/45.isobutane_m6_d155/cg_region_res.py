"""plan tooling-nozzle-core-grid §4.36: 1 出力の、縮流部 x∈[−8, −2] r_t と全域の、単位体積あたりの残差の体積重み付き RMS。
  q = sqrt(Σ V·(R/V)² / Σ V)  (R = res_ro・res_roe ほか、V = volume)。同じ物理領域・同じ定義で格子の違う run を比べる。
usage: python3 cg_region_res.py <run> <step> → _band_ab/core_grid/region_res_<run>.jsonl に 1 行足す"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, step = sys.argv[1], int(sys.argv[2])
S = json.loads((HERE / run / "prepare_info.json").read_text())["scale_m"]
with h5py.File(HERE / run / f"res_{step}.h5", "r") as h:
    x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
    V = np.asarray(h["VALUE/volume"][:], float)
    R = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega") if "VALUE/" + k in h}
rec = {"step": step, "nonfinite": int(sum(np.count_nonzero(~np.isfinite(v)) for v in R.values()))}
mc = (x >= -8) & (x <= -2)   # §4.80: 縮流部の外 (outside) を足した (既存のキーの値は不変)
for lab, m in (("contraction", mc), ("all", np.ones_like(x, bool)), ("outside", ~mc)):
    for k, r in R.items():
        rec[f"{lab}_{k}"] = float(np.sqrt(np.sum(V[m] * (r[m] / V[m]) ** 2) / np.sum(V[m])))
with open(HERE / "_band_ab" / "core_grid" / f"region_res_{run}.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
