"""§6.9 の判定 (plan time_integration-line-implicit-speed、2026-10-10): 3 腕の到達の step・総時間・E2・分解能。
P・L0 は既存の系列 (事後の集計)、L5 は run_0353_m9_L5 (ライン) と run_0354_m9_L5cut (point) の見張りの系列。
usage: python3 m9_judge.py → _band_ab/cold_pair/m9_judge.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
U = json.loads((D / "m9_unit.json").read_text())
def load(chain):
    rows = []
    for r, base in chain:
        for x in json.loads((D / f"series_{r}.json").read_text())["rows"]:
            if x["step"] == 0 and rows: continue
            y = dict(x); y["abs"] = base + x["step"]; rows.append(y)
    return rows
def drift(rows, key):
    st = np.array([r["abs"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if st[-1] - st[0] < 20000: return None
    j = int(np.searchsorted(st, st[-1] - 20000)); return 100 * (v[-1] / v[j] - 1) * 20000 / (st[-1] - st[j])
def reach(rows):
    for n in range(len(rows)):
        pre = rows[: n + 1]; r = pre[-1]
        dr = [drift(pre, k) for k in ("theta_r_40", "theta_r_70", "theta_r_94")]
        if r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr):
            return r
    return None
out = {}
P = reach(load([("run_0217_ns_coldmesh_tw300_cfl4_ext3", 440000), ("run_0263_ns_coldmesh_tw300_cfl4_ext4", 640000)]))
out["P"] = {"source": "既存の系列 (事後の集計)", "reach": P["abs"], "res_steps": 10000,
            "phases": [{"mode": "P", "steps": P["abs"]}], "restarts": 4}
L0l = reach(load([("run_0223_ns_coldmesh_tw300_linedir_tj5_cap50", 0), ("run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext", 15000), ("run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2", 75000)]))
L0c = reach(load([("run_0262_ns_coldmesh_tw300_cutback_point", 0)]))
out["L0"] = {"source": "既存の系列 (事後の集計)", "line_reach": L0l["abs"], "switch_at": 135000, "cut_reach_point_steps": L0c["abs"], "res_steps": 5000,
             "phases": [{"mode": "L0", "steps": 135000}, {"mode": "P", "steps": L0c["abs"]}], "restarts": 4, "final": L0c}
st5 = json.loads((HERE / "run_0353_m9_L5" / "m9_watch.json").read_text())
out["L5"] = {"source": "新しい run (事前登録)", "line_status": st5["status"], "line_reach": st5.get("reach_step"), "res_steps": 5000}
if st5["status"] == "REACHED":
    stc = json.loads((HERE / "run_0354_m9_L5cut" / "m9_watch.json").read_text())
    out["L5"].update(cut_status=stc["status"], cut_reach_point_steps=stc.get("reach_step"))
    if stc["status"] == "REACHED":
        c = [r for r in json.loads((D / "series_run_0354_m9_L5cut.json").read_text())["rows"] if r["step"] == stc["reach_step"]][0]
        out["L5"].update(phases=[{"mode": "L5", "steps": st5["reach_step"]}, {"mode": "P", "steps": stc["reach_step"]}], restarts=2, final=c)
for k, a in out.items():
    if "phases" not in a: a["verdict"] = "未到達・失敗"; continue
    a["total_s"] = sum(ph["steps"] * U[ph["mode"]]["median_ms"] / 1000 for ph in a["phases"]) + a["restarts"] * U["P"]["median_startup_s"]
    a["uncert_s"] = sum(a["res_steps"] * U[ph["mode"]]["median_ms"] / 1000 for ph in a["phases"])
    if k != "P":
        f = a["final"]
        a["E2"] = {key: (f[key] / P[key] - 1) * 100 for key in ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")}
        a["E2_pass"] = all(abs(v) <= 0.1 for v in a["E2"].values())
def cmp(x, y):
    a, b = out[x], out[y]
    if "total_s" not in a or "total_s" not in b: return "判別不能 (未到達)"
    if a.get("E2_pass") is False or b.get("E2_pass") is False: return "比べない (別の解)"
    d = a["total_s"] - b["total_s"]; u = a["uncert_s"] + b["uncert_s"]
    return f"{x} が速い" if d < -u else f"{x} が遅い" if d > u else "判別不能 (差が分解能以内)"
out["compare"] = {"L5_vs_L0": cmp("L5", "L0"), "L0_vs_P": cmp("L0", "P"), "L5_vs_P": cmp("L5", "P")}
(D / "m9_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float)); print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("final",)} if isinstance(v, dict) else v for k, v in out.items()}, indent=1, ensure_ascii=False, default=float))
