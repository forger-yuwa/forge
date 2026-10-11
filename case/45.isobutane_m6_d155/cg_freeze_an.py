"""plan tooling-nozzle-core-grid §4.38 F1: SST の凍結の実効性を、続いた 2 つの出力の差で測る (境界の節点を除く)。
usage: python3 cg_freeze_an.py <run> <step_prev> <step> → _band_ab/core_grid/freeze_check_<run>.jsonl に 1 行
  roK・roOmega の変化 max|Δ|/max|Q| と RMS(Δ)/RMS(Q)、vis_turb の変化 RMS(Δμt)/RMS(μt)、dt_local の比の中央値"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
ni = int(json.loads((HERE / run / "prepare_info.json").read_text())["mesh"]["ni"])
def L(k):
    with h5py.File(HERE / run / f"res_{k}.h5", "r") as h:
        return {n: np.asarray(h["VALUE/" + n][:], float) for n in ("roK", "roOmega", "vis_turb", "dt_local") if "VALUE/" + n in h}
A, B = L(a), L(b); n = len(A["roK"]); nj = n // ni
jj = np.tile(np.arange(nj), ni); ii = np.repeat(np.arange(ni), nj)
inner = (jj > 0) & (jj < nj - 1) & (ii > 0) & (ii < ni - 1)
rec = {"steps": [a, b]}
for k in ("roK", "roOmega", "vis_turb"):
    if k in A:
        d = (B[k] - A[k])[inner]; q = A[k][inner]
        rec[k] = {"max_rel": float(np.abs(d).max() / max(np.abs(q).max(), 1e-300)), "rms_rel": float(np.sqrt(np.mean(d ** 2)) / max(np.sqrt(np.mean(q ** 2)), 1e-300))}
if "dt_local" in A:
    rec["dt_ratio_median"] = float(np.median(B["dt_local"][inner] / np.maximum(A["dt_local"][inner], 1e-300)))
with open(HERE / "_band_ab" / "core_grid" / f"freeze_check_{run}.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
