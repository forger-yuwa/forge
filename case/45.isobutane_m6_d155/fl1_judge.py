"""plan architecture-float-state-double-geometry §6.20 (事前登録: hoop の修正 (キー 1) の格子での V4 のやり直し、§5.1 #15 (1)) の判定。fl1.sh と hp7.sh の後に AWS の case/45 で回す。
判定の式と閾値は §6.11 (v4_judge.py) と同じ。腕だけを float = run_0484_fl1_f32_k1、FP64 = run_0483_hp7_k1 に替えた。
- 比較が成り立つ条件: 到達を保存した系列 (m9_watch.json) から計算し直して確かめる。出力の間隔 2500・末尾 2 万 step に 9 点以上・
  必要な残差の列・非有限なし・発散なし。満たさなければ判定不能。
- 判定: float が上限までに水準に入らなければ FAIL。各腕の到達時の値で、θ_r(40/70/94) ≤ 0.05 %、符号つきの Q_w ≤ 0.1 %、
  末尾 2 万 step の各 rms_* の中央値 ≤ FP64 の 2 倍 (FP64 が 0 なら float も 0)、到達時の 2π·Σ|res_ro| ≤ FP64 の 2 倍。分母は FP64 の到達時の値の絶対値。
- 記録: 到達の step、check_convergence と check_quasisteady (系列、θ_r は drift/osc 1e-4、Q_w は 2e-4、--tail 1.0)、commit_loss.csv の統計。
結果は標準出力と _band_ab/cold_pair/fl1_judge.json。末尾に V4 (v4_judge.json) の差と並べる。"""
import csv, json, math, subprocess, sys, h5py, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
ARMS = {"float": HERE / "run_0484_fl1_f32_k1", "fp64": HERE / "run_0483_hp7_k1"}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
NEED = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega")
OUT = {"arms": {}, "flags": [], "undecidable": []}
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
for name, run in ARMS.items():
    st = json.loads((run / "m9_watch.json").read_text())
    rows = sorted([r for r in st["rows"] if r["step"] > 0], key=lambda r: r["step"])
    steps = [r["step"] for r in rows]
    a = {"status": st["status"], "reach_step_watch": st.get("reach_step"), "n_rows": len(rows)}
    if steps != [2500 * (i + 1) for i in range(len(steps))]: OUT["undecidable"].append(f"{name}: 出力の間隔が 2500 でない")
    d = reach(rows); a["reach_recomputed"] = d
    if d is None or d[0] != "REACH":
        a["reached"] = False
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
        tail = [x for x in hr if n - 20000 <= int(x["step"]) + 1 <= n]   # residual_history の step は res の番号より 1 小さい
        vals = {c: [float(x[c]) for x in tail] for c in cols}
        if any(not math.isfinite(v) for c in cols for v in vals[c]): OUT["undecidable"].append(f"{name}: 末尾の残差に非有限")
        a["rms_median"] = {c: float(np.median(vals[c])) if vals[c] else None for c in cols}
        with h5py.File(run / f"res_{n}.h5", "r") as h:
            rr = np.asarray(h["VALUE/res_ro"][:], np.float64) if "VALUE/res_ro" in h else None
        if rr is None or not np.all(np.isfinite(rr)): OUT["undecidable"].append(f"{name}: res_{n} の res_ro が無いか非有限")
        else: a["sum_abs_res_ro_2pi"] = float(2 * math.pi * np.sum(np.abs(rr)))
    OUT["arms"][name] = a; A[name] = a
    print(f"== {name}: 見張り {a['status']} 到達 {a['reach_step_watch']}、計算し直し {a['reach_recomputed']}、出力 {a['n_rows']} 点")
f, r = A["float"], A["fp64"]
if not r.get("reached"): OUT["undecidable"].append("FP64 の参照が水準に入っていない")
if OUT["undecidable"]:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"])
elif not f.get("reached"):
    verdict = "FAIL (float が上限までに水準に入らない)"
else:
    for k in KEYS + ("Q_w",):
        dv = 100 * (f["at_reach"][k] - r["at_reach"][k]) / abs(r["at_reach"][k]); lim = 0.1 if k == "Q_w" else 0.05
        print(f"  {k:11s} float {f['at_reach'][k]:.6g}  FP64 {r['at_reach'][k]:.6g}  差 {dv:+.4f} %  (上限 {lim} %)")
        if abs(dv) > lim: OUT["flags"].append(f"{k} {dv:+.4f} %")
    for c, vr in r["rms_median"].items():
        vf = f["rms_median"].get(c)
        if vr is None or vf is None: OUT["flags"].append(f"{c} の中央値が無い"); continue
        if (vr == 0 and vf != 0) or (vr > 0 and vf > 2 * vr): OUT["flags"].append(f"{c} 中央値 float {vf:.3e} > 2 × FP64 {vr:.3e}")
        print(f"  {c:12s} 末尾 2 万 step の中央値 float {vf:.3e}  FP64 {vr:.3e}  比 {vf / vr if vr else float('nan'):.2f}")
    sf, sr = f.get("sum_abs_res_ro_2pi"), r.get("sum_abs_res_ro_2pi")
    print(f"  2π·Σ|res_ro| float {sf:.4g}  FP64 {sr:.4g}  比 {sf / sr:.2f}")
    if sf > 2 * sr: OUT["flags"].append(f"2π·Σ|res_ro| {sf:.3g} > 2 × {sr:.3g}")
    verdict = "登録した水準へ到達し、停止時の比較でフラグなし" if not OUT["flags"] else "到達、フラグあり: " + "; ".join(OUT["flags"])
print(f"== VERDICT §6.20: {verdict}")
print(f"   到達の step: float {f.get('reach_recomputed')}  FP64 {r.get('reach_recomputed')}")
OUT["verdict"] = verdict
print("== 記録: 判定ツール")
for name, run in ARMS.items():
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    print(f"  {name} check_convergence: {ov[-1] if ov else '(行なし)'}")
    a = A[name]
    if a.get("reached"):
        n = a["reach_recomputed"][1]; st = json.loads((run / "m9_watch.json").read_text())
        rows = [x for x in st["rows"] if n - 20000 <= x["step"] <= n]
        sc = run / "fl1_series_tail.csv"
        with open(sc, "w") as fh:
            fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
            for x in sorted(rows, key=lambda x: x["step"]): fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
        for cols, dr in ((",".join(KEYS), "0.0001"), ("Q_w", "0.0002")):
            p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", cols,
                                "--tail", "1.0", "--drift", dr, "--osc", dr], capture_output=True, text=True)
            ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
            print(f"  {name} check_quasisteady ({cols}): {ov[-1] if ov else '(行なし) rc=' + str(p.returncode)}")
    cl = run / "commit_loss.csv"
    if cl.exists():
        rows = list(csv.DictReader(open(cl)))
        last = max(int(x["step"]) for x in rows) if rows else None
        print(f"  {name} commit_loss.csv の最後の step {last}:")
        for x in rows:
            if int(x["step"]) == last and x["region"] == "interior":
                print(f"     {x['qty']:8s} n_req {x['n_req']}  n_lost {x['n_lost']}  S_act/S_req {float(x['S_act']) / max(float(x['S_req']), 1e-300):.4f}  S_lostfrac {x['S_lostfrac']}")
o = HERE / "_band_ab/cold_pair/fl1_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
print("== 記録: V4 (キー 0 の格子、§6.16) との並び (float − FP64)/FP64 [%]")
try:
    v4 = json.loads((HERE / "_band_ab/cold_pair/v4_judge.json").read_text())
    for k in KEYS + ("Q_w",):
        old = 100 * (v4["arms"]["float"]["at_reach"][k] - v4["arms"]["fp64"]["at_reach"][k]) / abs(v4["arms"]["fp64"]["at_reach"][k])
        new = 100 * (A["float"]["at_reach"][k] - A["fp64"]["at_reach"][k]) / abs(A["fp64"]["at_reach"][k]) if A["float"].get("reached") and A["fp64"].get("reached") else float("nan")
        print(f"  {k:11s} V4 {old:+.3f} %  今回 {new:+.3f} %")
except Exception as e:
    print("  (並べられない)", e)
