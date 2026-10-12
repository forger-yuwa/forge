"""plan architecture-float-state-double-geometry §6.56 (T2、事前登録、codex diagnose 2026-10-12 `2026-10-12-pref-result-path-diagnose.md`) の判定。B で回す (A の 2 本は B へ写してから)。
腕: FP64 = run_0752_t2_d1 (A)・run_0754_t2_d2 (B)、float = run_0753_t2_f1 (A)・run_0755_t2_f2 (B)。いずれも §6.51 の上限 1 の run の res_30000 から再開。
停止規則 (§6.11 と同じ、m9_watch.py): 2,500 step ごと、|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step を 2 出力連続で満たした 2 つ目で止める。
  上限は共通の起点から通算 20 万 step (この run では 170,000 step)。
比較が成り立つ条件 (§6.11 と同じ、満たさなければ判定不能): 到達を保存した系列から計算し直し、出力の間隔 2500・末尾 2 万 step に 9 点以上・必要な残差の列・非有限なし。
判定: 4 本とも到達し、FP64 と float の全組 (2 × 2) で θ_r(40/70/94) の相対差 ≤ 0.05 %、符号つきの Q_w ≤ 0.1 %、末尾 2 万 step の各 rms_* の中央値 ≤ FP64 の 2 倍、
  到達時の 2π·Σ|res_ro| ≤ FP64 の 2 倍 (分母は FP64 の到達時の値) → 結果 A (停止規則での比較にフラグなし)。未到達またはフラグあり → 結果 B。
記録: 到達の step、同じ精度の 2 本の差、check_convergence、check_quasisteady (末尾 2 万 step、θ_r は drift/osc 1e-4、Q_w は 2e-4)。→ _band_ab/cold_pair/t2_judge.json"""
import csv, json, math, subprocess, sys, h5py, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
ARMS = {"fp64": ("run_0752_t2_d1", "run_0754_t2_d2"), "float": ("run_0753_t2_f1", "run_0755_t2_f2")}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
NEED = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega")
OUT = {"arms": {}, "flags": [], "undecidable": [], "pairs": {}}
def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return float(100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j]))
def reach(rows):
    run_ok = 0
    for n in range(len(rows)):
        r = rows[n]
        if r.get("nonfinite") or r.get("nonfinite_all"): return ("DIV", r["step"])
        dr = [drift(rows[: n + 1], k) for k in KEYS]
        win = sum(1 for x in rows[: n + 1] if r["step"] - 20000 <= x["step"] <= r["step"])
        ok = r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(d) <= 0.05 for d in dr) and win >= 9
        run_ok = run_ok + 1 if ok else 0
        if run_ok >= 2: return ("REACH", r["step"])
    return None
A = {}
for prec, runs in ARMS.items():
    for name in runs:
        run = HERE / name
        if not (run / "m9_watch.json").exists(): OUT["undecidable"].append(f"{name}: m9_watch.json が無い"); continue
        st = json.loads((run / "m9_watch.json").read_text())
        rows = sorted([r for r in st["rows"] if r["step"] > 0], key=lambda r: r["step"])
        steps = [r["step"] for r in rows]
        a = {"prec": prec, "status": st["status"], "reach_step_watch": st.get("reach_step"), "n_rows": len(rows)}
        if steps != [2500 * (i + 1) for i in range(len(steps))]: OUT["undecidable"].append(f"{name}: 出力の間隔が 2500 でない")
        d = reach(rows); a["reach_recomputed"] = d
        if d is not None and d[0] == "DIV": OUT["undecidable"].append(f"{name}: 非有限 (step {d[1]})")
        if d is None or d[0] != "REACH":
            a["reached"] = False
            if st["status"] not in ("CENSORED",): OUT["undecidable"].append(f"{name}: 到達せず、見張りの状態が {st['status']} (上限まで回り切っていない)")
        else:
            a["reached"] = True; n = d[1]
            if st.get("reach_step") != n: OUT["undecidable"].append(f"{name}: 見張りの到達 {st.get('reach_step')} と計算し直した到達 {n} が違う")
            R = next(r for r in rows if r["step"] == n); a["at_reach"] = {k: R[k] for k in KEYS + ("Q_w", "deficit")}
            a["window_points"] = sum(1 for x in rows if n - 20000 <= x["step"] <= n)
            a["drift_at_reach"] = [drift([x for x in rows if x["step"] <= n], k) for k in KEYS]
            hr = [x for x in csv.DictReader(open(run / "residual_history.csv")) if x["phase"] == "outer_end"]
            cols = [c for c in hr[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
            miss = [c for c in NEED if c not in cols] + ([] if any(c.startswith("rms_roY") for c in cols) else ["rms_roY*"])
            if miss: OUT["undecidable"].append(f"{name}: 残差の列が無い {miss}")
            tail = [x for x in hr if n - 20000 <= int(x["step"]) + 1 <= n]
            vals = {c: [float(x[c]) for x in tail] for c in cols}
            if any(not math.isfinite(v) for c in cols for v in vals[c]): OUT["undecidable"].append(f"{name}: 末尾の残差に非有限")
            a["rms_median"] = {c: float(np.median(vals[c])) if vals[c] else None for c in cols}
            p = run / f"res_{n}.h5"
            if not p.exists(): OUT["undecidable"].append(f"{name}: 到達の出力 res_{n} が無い")
            else:
                with h5py.File(p, "r") as h:
                    rr = np.asarray(h["VALUE/res_ro"][:], np.float64) if "VALUE/res_ro" in h else None
                if rr is None or not np.all(np.isfinite(rr)): OUT["undecidable"].append(f"{name}: res_{n} の res_ro が無いか非有限")
                else: a["sum_abs_res_ro_2pi"] = float(2 * math.pi * np.sum(np.abs(rr)))
        OUT["arms"][name] = a; A[name] = a
        print(f"== {name} ({prec}): 見張り {a['status']} 到達 {a['reach_step_watch']}、計算し直し {a['reach_recomputed']}、出力 {a['n_rows']} 点")
if OUT["undecidable"]:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"])
elif not all(A[n]["reached"] for runs in ARMS.values() for n in runs):
    verdict = "結果 B (未到達: " + ", ".join(n for runs in ARMS.values() for n in runs if not A[n]["reached"]) + ")"
else:
    for fn in ARMS["float"]:
        for dn in ARMS["fp64"]:
            f, r = A[fn], A[dn]; key = f"{fn} 対 {dn}"; pr = {}
            for k in KEYS + ("Q_w",):
                dv = 100 * (f["at_reach"][k] - r["at_reach"][k]) / abs(r["at_reach"][k]); lim = 0.1 if k == "Q_w" else 0.05
                pr[k] = dv
                if abs(dv) > lim: OUT["flags"].append(f"{key}: {k} {dv:+.4f} %")
            for c, vr in r["rms_median"].items():
                vf = f["rms_median"].get(c)
                if vr is None or vf is None: OUT["flags"].append(f"{key}: {c} の中央値が無い"); continue
                if (vr == 0 and vf != 0) or (vr > 0 and vf > 2 * vr): OUT["flags"].append(f"{key}: {c} 中央値 float {vf:.3e} > 2 × FP64 {vr:.3e}")
            sf, sr = f["sum_abs_res_ro_2pi"], r["sum_abs_res_ro_2pi"]; pr["sum_abs_res_ro_ratio"] = sf / sr
            if sf > 2 * sr: OUT["flags"].append(f"{key}: 2π·Σ|res_ro| {sf:.3g} > 2 × {sr:.3g}")
            OUT["pairs"][key] = pr
            print(f"  {key}: " + "  ".join(f"{k} {v:+.4f} %" for k, v in pr.items() if k != "sum_abs_res_ro_ratio") + f"  Σ|res_ro| 比 {pr['sum_abs_res_ro_ratio']:.2f}")
    verdict = "結果 A (停止規則での比較にフラグなし)" if not OUT["flags"] else "結果 B (フラグあり: " + "; ".join(OUT["flags"][:6]) + ")"
# 記録: 同じ精度の 2 本の差
for prec, (n1, n2) in ARMS.items():
    if A.get(n1, {}).get("reached") and A.get(n2, {}).get("reached"):
        OUT["arms"][f"{prec} の 2 本の差 [%]"] = {k: 100 * (A[n2]["at_reach"][k] - A[n1]["at_reach"][k]) / abs(A[n1]["at_reach"][k]) for k in KEYS + ("Q_w",)}
print(f"== VERDICT §6.56: {verdict}")
OUT["verdict"] = verdict
print("== 記録: 判定ツール")
for prec, runs in ARMS.items():
    for name in runs:
        run = HERE / name
        if name not in A: continue
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["arms"][name]["check_convergence"] = ov[-1] if ov else "(行なし)"
        print(f"  {name} check_convergence: {OUT['arms'][name]['check_convergence']}")
        a = A[name]
        if a.get("reached"):
            n = a["reach_recomputed"][1]; st = json.loads((run / "m9_watch.json").read_text())
            rows = [x for x in st["rows"] if n - 20000 <= x["step"] <= n]
            sc = run / "t2_series_tail.csv"
            with open(sc, "w") as fh:
                fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
                for x in sorted(rows, key=lambda x: x["step"]): fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
            for cols, dr in ((",".join(KEYS), "0.0001"), ("Q_w", "0.0002")):
                p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", cols,
                                    "--tail", "1.0", "--drift", dr, "--osc", dr], capture_output=True, text=True)
                ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l or "->" in l]
                OUT["arms"][name][f"check_quasisteady {cols}"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
                print(f"  {name} check_quasisteady ({cols}): {OUT['arms'][name][f'check_quasisteady {cols}']}")
o = HERE / "_band_ab/cold_pair/t2_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
