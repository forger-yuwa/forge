# 診断用: run_0005_ns の最後の場 (res_200000) から関門と同じ固定の条件で δ_E を抽出し、入れた δ (δ_in) との差を測る。
# run dir には書かない (出力は argv[2])。関門の判定ではない (関門は止まったまま)。
import json, sys
from pathlib import Path
import numpy as np
from forge_design.feedback import deltastar_loop as DL
from forge_design.evaluate import stop_rule as SR
run = Path(sys.argv[1]); od = Path(sys.argv[2]); od.mkdir(parents=True, exist_ok=True)
info = json.loads((run / "prepare_info.json").read_text())
eref = SR.euler_ref_spec(info.get("euler_ref")) or info.get("euler_ref")
print("euler_ref", eref)
st = DL.handoff_extract_settings()
summ = DL.extract_and_merge(run, eref, out_dir=od, ns_res="res_200000.h5", **st)
d = DL.read_delta_r_next(od / "delta_r_next.csv")
rt = float(info.get("r_throat") or info.get("physical_wall", {}).get("r_throat") or 0.0)
x = d["x_rt"]; dE = d["delta_E"]; di = d["delta_in_prev"]; dt = d["delta_target"]
fin = np.isfinite(dE)
print(json.dumps({k: summ[k] for k in ("n_stations", "n_ok", "n_hard_ok", "n_held", "use_over_in_median", "use_over_in_p10_p90",
                                        "max_abs_change", "delta_r_throat_use", "delta_r_throat_in")}, indent=1, default=str))
print("r_t (prepare_info)", rt)
np.save(od / "cmp.npy", np.c_[x, di, dE, dt])
