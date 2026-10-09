"""§6.13 の判定 (plan time_integration-line-implicit-speed、2026-10-10): implicitRelax 0.85・1.0 × 値 0・粘性入りの 4 腕の総時間を、緩和 0.7 の腕 (§6.11 の L0・L5、m9_judge.json) と P と比べる。
各腕: ライン段 (緩和 r) の到達 n1 (水準) と point 段の到達 n2 (水準 + E2) を見張りの状態から取り、系列から再計算して一致を確かめる。
総時間 (推定) = n1 × 1 step (L0 か L5 の単価、緩和は 1 step の費用を変えない = 仮定) + n2 × 1 step (P) + 出力の回数 × 出力 1 回 + 起動 2 回。分解能・比較の前提は §6.11 と同じ。
usage: python3 m9r_judge.py → _band_ab/cold_pair/m9r_judge.json"""
import csv, json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
U = json.loads((D / "m9_unit.json").read_text()); R = json.loads((HERE / "m9_ref.json").read_text())
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
def rows_of(r):
    return [x for x in json.loads((HERE / r / "m9_watch.json").read_text())["rows"] if x["step"] > 0]
def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return 100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j])
def first(rows, e2):
    for n in range(len(rows)):
        pre = rows[: n + 1]; r = pre[-1]
        if r.get("nonfinite") or r.get("nonfinite_all"): return None
        dr = [drift(pre, k) for k in KEYS[:3]]
        ok = r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr)
        if ok and e2: ok = all(abs(100 * (r[k] / float(R[k]) - 1)) <= 0.1 for k in KEYS)
        if ok: return {"step": r["step"], "window": [x for x in pre if x["step"] >= r["step"] - 20000]}
    return None
def qs(window, tag):
    p = D / f"m9r_qs_{tag}.csv"
    with open(p, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step"] + list(KEYS))
        for x in window: w.writerow([x["step"]] + [x[k] for k in KEYS])
    o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p), "--series-cols", ",".join(KEYS),
                        "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
    ov = [l for l in (o.stdout + o.stderr).splitlines() if "OVERALL" in l]
    return {"overall": ov[-1] if ov else None, "ok": bool(ov)}
def conv(run):
    v = HERE / run / "CONVERGENCE_VERDICT.txt"; t = v.read_text(errors="replace") if v.exists() else ""
    return "DIVERGED" if "DIVERGED" in t else "NOT CONVERGED" if "NOT CONVERGED" in t else "PASS" if "-> PASS" in t else "記録なし"
ARMS = {"L0_r085": ("run_0364_m9r_L0_r085", "run_0365_m9r_L0_r085cut", "L0", 0.85), "L0_r100": ("run_0366_m9r_L0_r100", "run_0367_m9r_L0_r100cut", "L0", 1.0),
        "L5_r085": ("run_0368_m9r_L5_r085", "run_0369_m9r_L5_r085cut", "L5", 0.85), "L5_r100": ("run_0370_m9r_L5_r100", "run_0371_m9r_L5_r100cut", "L5", 1.0)}
out = {}
for k, (rl, rp, m, rx) in ARMS.items():
    a = {"line_run": rl, "point_run": rp, "mode": m, "relax": rx}
    try:
        sl = json.loads((HERE / rl / "m9_watch.json").read_text()); a["line"] = {"status": sl["status"], "reach": sl.get("reach_step"), "fail_step": sl.get("fail_step"), "convergence": conv(rl)}
        if sl["status"] == "REACHED":
            sp = json.loads((HERE / rp / "m9_watch.json").read_text()); a["point"] = {"status": sp["status"], "reach": sp.get("reach_step"), "convergence": conv(rp)}
            if sp["status"] == "REACHED":
                fl, fp = first(rows_of(rl), False), first(rows_of(rp), True)
                n1, n2 = sl["reach_step"], sp["reach_step"]
                a["recheck"] = {"line": fl and fl["step"], "point": fp and fp["step"], "match": bool(fl and fp and fl["step"] == n1 and fp["step"] == n2)}
                if a["recheck"]["match"]:
                    t = n1 * U[m]["median_ms"] / 1000 + n2 * U["P"]["median_ms"] / 1000 + (n1 // 5000) * U[m]["median_output_ms"] / 1000 + (n2 // 5000) * U["P"]["median_output_ms"] / 1000 \
                        + U[m]["median_startup_s"] + U["P"]["median_startup_s"]
                    u = 5000 * (U[m]["median_ms"] + U["P"]["median_ms"]) / 1000 + (n1 * U[m]["spread_ms"] + n2 * U["P"]["spread_ms"]) / 1000
                    a.update(total_s=t, uncert_s=u, quasisteady=qs(fp["window"], k))
    except FileNotFoundError as e:
        a["status"] = f"状態ファイルなし ({e.filename})"
    out[k] = a
base = json.loads((D / "m9_judge.json").read_text()) if (D / "m9_judge.json").exists() else {}
def cmp(a, b, name_b):
    if "total_s" not in a or not isinstance(b, dict) or "total_s" not in b: return f"判別不能 (どちらかが未到達・失敗・未集計: {name_b})"
    if not a["quasisteady"]["ok"] or (b.get("quasisteady") and not b["quasisteady"].get("ok", True)): return "比較不可 (準定常の判定が正常に終わらない)"
    if any(x.get("convergence") == "DIVERGED" for x in (a.get("line", {}), a.get("point", {}))): return "比較不可 (DIVERGED)"
    d = a["total_s"] - b["total_s"]; u = a["uncert_s"] + b["uncert_s"]
    note = "" if ("NOT ALL" not in (a["quasisteady"]["overall"] or "") and all(x.get("convergence") == "PASS" for x in (a["line"], a["point"]))) else " [水準と E2 への到達時間の比較に限る]"
    return (f"速い ({a['total_s']/3600:.2f} h vs {b['total_s']/3600:.2f} h)" if d < -u else f"遅い ({a['total_s']/3600:.2f} h vs {b['total_s']/3600:.2f} h)" if d > u else "判別不能 (差が分解能以内)") + note
out["compare"] = {}
for k, a in out.items():
    if k == "compare": continue
    out["compare"][f"{k}_vs_relax0.7"] = cmp(a, base.get(a["mode"]), f"緩和 0.7 の {a['mode']}")
    out["compare"][f"{k}_vs_P"] = cmp(a, base.get("P"), "P")
(D / "m9r_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
print(json.dumps(out, indent=1, ensure_ascii=False, default=float))
