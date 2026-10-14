"""§6.10 の起動の条件 (plan time_integration-line-viscous-jacobian): 1 step 目の書き出しでマスク 7・15・5 を比べる。外れたら終了コード 1 (比較の run を起動しない)。
  (i) D が 7・15 で完全一致 (ii) 拘束のフラグが一致 (iii) rhs (sweep 0) の差 ≤ 1e-10
  (iv) K: 行 0〜3 と行 4 の列 1〜4 が 7・15 で完全一致、行 4 の列 0 は 15 と 5 で完全一致、7 と 15 では違う
加えて A (マスク 7) の 1 step 目 (sweep 0) について、熱伝導の K の作用を列に分けた量を出す (K7 − K5 の行 4 = 熱伝導の K、ΔQ は dqnew_s0 / relax)。
usage: python3 e1b_gate.py <dump_m7> <dump_m15> <dump_m5>
"""
import json
import sys
from pathlib import Path

import numpy as np


def load(d):
    meta, relax = {}, None
    for line in (Path(d) / "meta.txt").read_text().splitlines():
        if line.startswith("#"):
            relax = float(line.split("implicitRelax")[1].split()[0]); continue
        n, r, c = line.split(); meta[n] = (int(r), int(c))
    a = {k: np.fromfile(Path(d) / f"{k}.f64").reshape(v) for k, v in meta.items()}
    a["_relax"] = relax
    return a


m7, m15, m5 = (load(x) for x in sys.argv[1:4])
o = {}
o["i_D_identical"] = bool(np.array_equal(m7["D"], m15["D"]))
o["ii_flags_identical"] = bool(np.array_equal(m7["flags_wall_iso_axis"], m15["flags_wall_iso_axis"]) and np.array_equal(m7["node_line"], m15["node_line"]))
o["iii_rhs_s0_rel"] = float(np.max(np.abs(m7["rhs_s0"] - m15["rhs_s0"])) / max(np.max(np.abs(m15["rhs_s0"])), 1e-300))
ok_iv = True
for K in ("Kprev", "Knext"):
    a7, a15, a5 = (x[K].reshape(-1, 5, 5) for x in (m7, m15, m5))
    same_rest = np.array_equal(a7[:, :4, :], a15[:, :4, :]) and np.array_equal(a7[:, 4, 1:], a15[:, 4, 1:])
    col0_15_eq_5 = np.array_equal(a15[:, 4, 0], a5[:, 4, 0])
    col0_7_ne_15 = bool(np.any(a7[:, 4, 0] != a15[:, 4, 0]))
    o[f"iv_{K}"] = {"rows0to3_and_row4_cols1to4_identical_7_15": bool(same_rest), "row4_col0_identical_15_5": bool(col0_15_eq_5), "row4_col0_differs_7_15": col0_7_ne_15}
    ok_iv = ok_iv and same_rest and col0_15_eq_5 and col0_7_ne_15
ok = o["i_D_identical"] and o["ii_flags_identical"] and o["iii_rhs_s0_rel"] <= 1e-10 and ok_iv
# A の熱伝導の K の作用を列に分ける (sweep 0、ライン内の隣の補正で)。K_heat = (K7 − K5) の行 4
nl = m7["node_line"]; lines = nl[:, 1].astype(int)
dq = m7["dqnew_s0"] / m7["_relax"]
parts = {"rho": [], "mom": [], "E": [], "sum": []}
for K, sh in (("Kprev", -1), ("Knext", +1)):
    Kh = (m7[K] - m5[K]).reshape(-1, 5, 5)[:, 4, :]
    for i in range(len(nl)):
        j = i + sh
        if j < 0 or j >= len(nl) or lines[j] != lines[i]:
            continue
        v = Kh[i] * dq[j]
        parts["rho"].append(v[0]); parts["mom"].append(v[1:4].sum()); parts["E"].append(v[4]); parts["sum"].append(v.sum())
o["A_heatK_action_sweep0_maxabs"] = {k: float(np.max(np.abs(v))) for k, v in parts.items()}
o["A_heatK_action_sweep0_rms"] = {k: float(np.sqrt(np.mean(np.square(v)))) for k, v in parts.items()}
o["GATE"] = "通過" if ok else "不通過 (比較の run を起動しない)"
print(json.dumps(o, indent=1, ensure_ascii=False))
(Path(sys.argv[1]) / "e1b_gate.json").write_text(json.dumps(o, indent=1, ensure_ascii=False))
sys.exit(0 if ok else 1)
