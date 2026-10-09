"""§6.11 の判定 v2 (plan time_integration-line-implicit-speed、2026-10-10、codex plan-2 の反映)。
終わり E = point の段の出力で (水準: |欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step、窓はその段の中、step 0 は使わない)
           かつ θ_r(40/70/94)・Q_w が参照 R (m9_ref.json) から相対 0.1 % 以内、を初めて満たした出力。
総時間 (推定、共通の環境の単価で換算) = Σ 段 (step 数 × 1 step の時間) + 出力の回数 × 出力 1 回の費用 + 起動の回数 × 起動と終了の時間 (各段の構成の値)。
分解能 = Σ 段 (出力の間隔 × 1 step の時間) + Σ 段 (step 数 × 単価の 3 本の幅)。
P・L0 は既存の系列の事後の集計 (古いバイナリの step 数を新しいバイナリでも同じと見なす、という仮定)。L5 は run_0353_m9_L5・run_0354_m9_L5cut の見張りの状態。
品質: 新しい段の check_convergence (DIVERGED は比較不可、NOT CONVERGED は「水準への到達の比較」に限る) と、到達の窓 [N − 2 万, N] の θ_r・Q_w の
check_quasisteady (系列 CSV を窓で切り出して --tail 1.0 --drift 0.0005 --min-snaps 5) を記録する。
比較の前提 (codex plan-3 M3): 見張りの到達の step と系列からの再計算が一致、新しい段の check_convergence の VERDICT がある (DIVERGED でない)、準定常の判定が正常に終わる
(OVERALL の行がある)。どれかが欠けたら比較不可。NOT CONVERGED・NOT ALL STEADY は「水準と E2 への到達時間」の比較に限って扱う。
usage (AWS の case dir): python3 m9_judge.py → _band_ab/cold_pair/m9_judge.json"""
import csv, json, math, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
U = json.loads((D / "m9_unit.json").read_text()); R = json.loads((HERE / "m9_ref.json").read_text())
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
def series(r):
    p = HERE / r / "m9_watch.json"
    rows = json.loads(p.read_text())["rows"] if p.exists() else json.loads((D / f"series_{r}.json").read_text())["rows"]
    return [x for x in rows if x["step"] > 0]
def chain(parts):
    out = []
    for r, base in parts:
        for x in series(r): y = dict(x); y["abs"] = base + x["step"]; out.append(y)
    return out
def drift(rows, key):
    s_ = np.array([r["abs"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return 100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j])
def reach(rows):
    for n in range(len(rows)):
        pre = rows[: n + 1]; r = pre[-1]
        if r.get("nonfinite") or r.get("nonfinite_all"): return {"diverged_at": r["abs"]}
        dr = [drift(pre, k) for k in KEYS[:3]]
        e2 = {k: 100 * (r[k] / R[k] - 1) for k in KEYS}
        if r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr) and all(abs(v) <= 0.1 for v in e2.values()):
            return {"abs": r["abs"], "drift": dr, "e2": e2, "window": [x for x in pre if x["abs"] >= r["abs"] - 20000]}
    return None
def qs(window, tag):
    p = D / f"m9_qs_{tag}.csv"
    with open(p, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step"] + list(KEYS))
        for x in window: w.writerow([x["abs"]] + [x[k] for k in KEYS])
    o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p), "--series-cols", ",".join(KEYS),
                        "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
    txt = (o.stdout + o.stderr).strip()
    ov = [l for l in txt.splitlines() if "OVERALL" in l]
    return {"rc": o.returncode, "overall": ov[-1] if ov else None, "ok": bool(ov), "tail": txt.splitlines()[-6:]}
def conv(run):
    v = HERE / run / "CONVERGENCE_VERDICT.txt"
    t = v.read_text(errors="replace") if v.exists() else ""
    return "DIVERGED" if "DIVERGED" in t else "NOT CONVERGED" if "NOT CONVERGED" in t else "PASS" if "-> PASS" in t else "記録なし"
def total(phases, restarts):
    t = sum(ph["steps"] * U[ph["mode"]]["median_ms"] / 1000 + ph["outputs"] * U[ph["mode"]]["median_output_ms"] / 1000 for ph in phases)
    t += sum(U[m]["median_startup_s"] for m in restarts)
    u = sum(ph["interval"] * U[ph["mode"]]["median_ms"] / 1000 + ph["steps"] * U[ph["mode"]]["spread_ms"] / 1000 for ph in phases)
    return t, u
out = {"ref": R}
# P (事後の集計): run_0183 res_100000 → 0191 (4 万) → 0198 → 0208 → 0217 (各 20 万、1 万ごと)。系列は 0217 (通算 44 万〜) と 0263 だけ。
rP = reach(chain([("run_0217_ns_coldmesh_tw300_cfl4_ext3", 440000), ("run_0263_ns_coldmesh_tw300_cfl4_ext4", 640000)]))
if rP and "abs" in rP:
    nP = rP["abs"]; tP, uP = total([{"mode": "P", "steps": nP, "outputs": nP // 10000, "interval": 10000}], ["P"] * (4 if nP <= 640000 else 5))
    out["P"] = {"source": "既存の系列 (事後の集計)", "reach": nP, "e2": rP["e2"], "total_s": tP, "uncert_s": uP, "quasisteady": qs(rP["window"], "P")}
# L0 (事後の集計): ライン 13.5 万 step (0223・0224・0252) → point の run_0262。
rC = reach(chain([("run_0262_ns_coldmesh_tw300_cutback_point", 0)]))
if rC and "abs" in rC:
    tL, uL = total([{"mode": "L0", "steps": 135000, "outputs": 6 + 12 + 12, "interval": 5000}, {"mode": "P", "steps": rC["abs"], "outputs": rC["abs"] // 5000, "interval": 5000}], ["L0", "L0", "L0", "P"])
    out["L0"] = {"source": "既存の系列 (事後の集計)", "line_steps": 135000, "point_steps": rC["abs"], "e2": rC["e2"], "total_s": tL, "uncert_s": uL, "quasisteady": qs(rC["window"], "L0")}
else:
    out["L0"] = {"source": "既存の系列 (事後の集計)", "status": "point の段で未到達"}
# L5 (事前登録)
s5 = json.loads((HERE / "run_0353_m9_L5" / "m9_watch.json").read_text())
out["L5"] = {"source": "新しい run (事前登録)", "line": {"status": s5["status"], "reach_step": s5.get("reach_step"), "convergence": conv("run_0353_m9_L5")}}
if s5["status"] == "REACHED":
    c5 = json.loads((HERE / "run_0354_m9_L5cut" / "m9_watch.json").read_text())
    out["L5"]["point"] = {"status": c5["status"], "reach_step": c5.get("reach_step"), "convergence": conv("run_0354_m9_L5cut")}
    if c5["status"] == "REACHED":
        n1, n2 = s5["reach_step"], c5["reach_step"]
        t5, u5 = total([{"mode": "L5", "steps": n1, "outputs": n1 // 5000, "interval": 5000}, {"mode": "P", "steps": n2, "outputs": n2 // 5000, "interval": 5000}], ["L5", "P"])
        rows = chain([("run_0354_m9_L5cut", 0)]); r5 = reach(rows)
        out["L5"]["recheck"] = {"reach_from_series": r5.get("abs") if r5 else None, "watch_reach": n2, "match": bool(r5 and r5.get("abs") == n2)}
        if out["L5"]["recheck"]["match"]:
            out["L5"].update(total_s=t5, uncert_s=u5, e2=r5["e2"], quasisteady=qs(r5["window"], "L5"))
def evidence(k):
    """比較の前提がそろっているか (codex plan-3 M3)。欠けていればその理由を返す。"""
    a = out.get(k, {})
    if "total_s" not in a: return "未到達・失敗・到達の再計算の不一致"
    if a.get("quasisteady") is None or not a["quasisteady"]["ok"]: return "準定常の判定が正常に終わらない"
    if k == "L5":
        for ph in ("line", "point"):
            c = a.get(ph, {}).get("convergence")
            if c in (None, "記録なし"): return f"{ph} の段の check_convergence が無い"
            if c == "DIVERGED": return f"{ph} の段が DIVERGED"
    return None
def cmp(x, y):
    a, b = out.get(x, {}), out.get(y, {})
    for k in (x, y):
        e = evidence(k)
        if e: return f"比較不可 ({k}: {e})"
    d = a["total_s"] - b["total_s"]; u = a["uncert_s"] + b["uncert_s"]
    note = "" if all("ALL STEADY" in (out[k]["quasisteady"]["overall"] or "") and "NOT ALL" not in (out[k]["quasisteady"]["overall"] or "") for k in (x, y)) else " [到達の窓の準定常は未確認 = 水準と E2 への到達時間の比較に限る]"
    return (f"{x} が速い ({a['total_s']/3600:.2f} h vs {b['total_s']/3600:.2f} h)" if d < -u else f"{x} が遅い ({a['total_s']/3600:.2f} h vs {b['total_s']/3600:.2f} h)" if d > u else "判別不能 (差が分解能以内)") + note
out["compare"] = {"L5_vs_L0": cmp("L5", "L0"), "L0_vs_P": cmp("L0", "P"), "L5_vs_P": cmp("L5", "P")}
(D / "m9_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk not in ("window",)} if isinstance(v, dict) else v) for k, v in out.items()}, indent=1, ensure_ascii=False, default=float))
