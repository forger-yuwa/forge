"""plan tooling-nozzle-core-grid §4.59 の資格 (A・B の最初の外反復のライン行列の書き出し): D・dt_local・体積が同一、初期の右辺 rhs_s0 の差が丸めの水準 (成分別 ≤ 1e−8)、
K (Kprev・Knext) の差がエネルギー行 (行 4) だけにあり 0 でない。usage: python3 cg_lhs_qual.py <A> <B> → _band_ab/core_grid/lhs_qual.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; A, B = sys.argv[1:3]
def dg(run):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    return lambda n: np.fromfile(HERE / run / "linedump" / f"{n}.f64").reshape(*meta[n])
ga, gb = dg(A), dg(B)
q = {"nodes_equal": bool(np.array_equal(ga("node_line"), gb("node_line"))),
     "D_equal": bool(np.array_equal(ga("D"), gb("D"))), "dt_vol_equal": bool(np.array_equal(ga("dt_vol"), gb("dt_vol")))}
ra, rb = ga("rhs_s0"), gb("rhs_s0")
q["rhs_s0_rel"] = [float(np.linalg.norm(ra[:, j] - rb[:, j]) / max(np.linalg.norm(ra[:, j]), 1e-300)) for j in (0, 1, 2, 4)]
for k in ("Kprev", "Knext"):
    d = (gb(k) - ga(k)).reshape(-1, 5, 5)
    q[f"{k}_diff_rows"] = [float(np.max(np.abs(d[:, i, :]))) for i in range(5)]
q["K_diff_only_row4"] = all(max(q[f"{k}_diff_rows"][:4]) == 0.0 for k in ("Kprev", "Knext")) and max(q["Kprev_diff_rows"][4], q["Knext_diff_rows"][4]) > 0
q["ok"] = q["nodes_equal"] and q["D_equal"] and q["dt_vol_equal"] and all(v <= 1e-8 for v in q["rhs_s0_rel"]) and q["K_diff_only_row4"]
(HERE / "_band_ab" / "core_grid" / "lhs_qual.json").write_text(json.dumps(q, indent=1)); print(json.dumps(q))
