"""plan tooling-nozzle-core-grid §4.104 P-hs (事前登録、codex 諮問 2026-10-14 hi-settle): 1 出力の停止規則の量。
設計チェーンの停止規則 `design/forge_design/evaluate/stop_rule.py` (feature/nozzle-cooled-wall-settling 719dfbf4 を `~/forge-coregrid/_ext/f5design_719dfbf4` に固定して読む)
の `snapshot` (入口・出口の質量流量と収支 r、x/x_E 0.4・0.7・0.95 の θ_r、出口コア M、無次元の流量) に、
参考として M6 の見張りの断面 x = 40・70・94 r_t の θ_r (同じ帯・同じ式) を足す。帯 y_b(x) は共通の起点で固定したもの (`cg_hs_band.py`)。
usage: python3 cg_hs_stoprow.py <run> <step> → _band_ab/core_grid/stoprows_<run>.jsonl に 1 行足す"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
F5 = Path.home() / "forge-coregrid" / "_ext" / "f5design_719dfbf4" / "design"
sys.path.insert(0, str(F5))
from forge_design.evaluate import stop_rule as SR  # noqa: E402
run, step = sys.argv[1], int(sys.argv[2])
band = json.loads((HERE / "_band_ab" / "core_grid" / "hs_band.json").read_text())
geo = SR.RunGeometry(HERE / run)
V = SR.read_values(HERE / run / f"res_{step}.h5")
row = SR.snapshot(geo, step, band, values=V)
if not row["nonfinite_all"]:
    sh = (geo.ni, geo.nj)
    for x0, th in SR.theta_r_sections(geo.X, geo.R, V["ro"].astype(float).reshape(sh), V["Ux"].astype(float).reshape(sh),
                                      np.asarray(band["x"]), np.asarray(band["y_b"]), [40.0, 70.0, 94.0]):
        row[f"m9_theta_r_{round(x0)}"] = th
    row["deficit_kg_s"] = row["mdot_in"] - row["mdot_out"]
with open(HERE / "_band_ab" / "core_grid" / f"stoprows_{run}.jsonl", "a") as fh:
    fh.write(json.dumps(row) + "\n")
