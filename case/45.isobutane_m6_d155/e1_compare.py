"""E1 の事前の確認 (plan time_integration-line-viscous-jacobian §6.8 改訂): 1 step 目のライン行列の書き出しを比べる。
  (1) マスク 7 と 5: D・rhs (sweep 0)・拘束の行が一致し、Kprev/Knext の差が行 4 (熱伝導) だけ
  (2) マスク 0 と 値 0 + キー 7: D が相対 1e-12 以内、K が一致
usage: python3 e1_compare.py <dump_m7> <dump_m5> <dump_m0> <dump_v0k7> → 標準出力と <dump_m7>/e1_compare.json
"""
import json
import sys
from pathlib import Path

import numpy as np


def load(d):
    meta = {}
    for line in (Path(d) / "meta.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        name, r, c = line.split(); meta[name] = (int(r), int(c))
    return {k: np.fromfile(Path(d) / f"{k}.f64").reshape(v) for k, v in meta.items()}


def rel(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


m7, m5, m0, v7 = (load(x) for x in sys.argv[1:5])
out = {}
assert np.array_equal(m7["node_line"], m5["node_line"]) and np.array_equal(m7["node_line"], v7["node_line"])
out["m7_vs_m5_D"] = rel(m7["D"], m5["D"])
out["m7_vs_m5_rhs_s0"] = rel(m7["rhs_s0"], m5["rhs_s0"])
for K in ("Kprev", "Knext"):
    d = (m7[K] - m5[K]).reshape(-1, 5, 5)
    out[f"m7_vs_m5_{K}_rows0to3_maxabs"] = float(np.max(np.abs(d[:, :4, :])))
    out[f"m7_vs_m5_{K}_row4_maxabs"] = float(np.max(np.abs(d[:, 4, :])))
    out[f"m0_vs_v0k7_{K}"] = rel(m0[K], v7[K])
out["m0_vs_v0k7_D"] = rel(m0["D"], v7["D"])
out["m0_vs_v0k7_rhs_s0"] = rel(m0["rhs_s0"], v7["rhs_s0"])
ok = (out["m7_vs_m5_D"] == 0.0 and out["m7_vs_m5_rhs_s0"] <= 1e-10 and out["m7_vs_m5_Kprev_rows0to3_maxabs"] == 0.0
      and out["m7_vs_m5_Knext_rows0to3_maxabs"] == 0.0 and out["m0_vs_v0k7_D"] <= 1e-12
      and out["m0_vs_v0k7_Kprev"] <= 1e-12 and out["m0_vs_v0k7_Knext"] <= 1e-12)
out["VERDICT"] = "比較は有効" if ok else "比較は無効 (事前の確認に外れた)"
print(json.dumps(out, indent=1, ensure_ascii=False))
(Path(sys.argv[1]) / "e1_compare.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
