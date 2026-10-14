"""plan tooling-nozzle-core-grid §4.59: 1 出力の物理性 (ρ・P・T の最小、非有限の数)。usage: python3 cg_phys.py <run> <step> → _band_ab/core_grid/phys_<run>.jsonl に 1 行足す"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, step = sys.argv[1], int(sys.argv[2])
with h5py.File(HERE / run / f"res_{step}.h5", "r") as h:
    v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("ro", "P", "T") if "VALUE/" + k in h}
rec = {"step": step, "nonfinite": int(sum(np.count_nonzero(~np.isfinite(a)) for a in v.values()))}
rec.update({f"min_{k}": float(np.nanmin(a)) for k, a in v.items()})
with open(HERE / "_band_ab" / "core_grid" / f"phys_{run}.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
