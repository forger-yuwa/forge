"""水準だけ (|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) ドリフト ≤ 0.05 %/2 万 step) に初めて入った出力を既存の系列で数える (plan time_integration-line-implicit-speed §6.15、事後の集計)。usage (AWS の case dir): python3 levelreach.py"""
import json, numpy as np
from pathlib import Path
D = Path("_band_ab/cold_pair")
def series(r):
    p = Path(r) / "m9_watch.json"
    rows = json.loads(p.read_text())["rows"] if p.exists() else json.loads((D / f"series_{r}.json").read_text())["rows"]
    return sorted([x for x in rows if x["step"] > 0], key=lambda x: x["step"])
def drift(rows, key):
    s_ = np.array([r["abs"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return 100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j])
K = ("theta_r_40", "theta_r_70", "theta_r_94")
def reach(parts, per_run):
    allrows = []
    for r, base in parts:
        rs = series(r); print(f"  {r}: step {rs[0]['step']}..{rs[-1]['step']} (通算 {base + rs[0]['step']}..{base + rs[-1]['step']}、{len(rs)} 点)")
        cur = []
        for x in rs:
            y = dict(x); y["abs"] = base + x["step"]; cur.append(y); allrows.append(y)
            win = cur if per_run else allrows
            dr = [drift(win, k) for k in K]
            if y.get("deficit") is not None and None not in dr and abs(y["deficit"]) <= 0.1 and all(abs(d) <= 0.05 for d in dr):
                return y["abs"], [round(d, 3) for d in dr], round(y["deficit"], 4)
    return None
L0 = [("run_0223_ns_coldmesh_tw300_linedir_tj5_cap50", 0), ("run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext", 15000), ("run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2", 75000)]
L5 = [("run_0353_m9_L5", 0)]
for per_run in (True, False):
    print("ドリフトの窓:", "run の中だけ (見張りと同じ)" if per_run else "run をまたいで連結")
    for name, parts in (("L0 ライン (値 0・上限 50)", L0), ("L5 ライン (粘性・上限なし)", L5)):
        print(name); print("  →", reach(parts, per_run))
