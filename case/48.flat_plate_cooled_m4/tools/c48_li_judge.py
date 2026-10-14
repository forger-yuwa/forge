#!/usr/bin/env python3
"""plan architecture-float-state-double-geometry §6.37 (事前登録、codex diagnose 2026-10-11) の判定: case/48 に case/45 の陰解法の設定一式を入れた float と FP64 の
**有限期間の再現判定** (定常解の比較ではない)。c48_li.sh の後に AWS A の case/48 で回す。結果は c48_li_judge.json。
腕: float = run_0061_li_f32a・run_0062_li_f32b、FP64 = run_0063_li_f64a・run_0064_li_f64b。系列は各 run の prec_series.csv (step 0〜48,000、共通の double の座標) と
prec_series_fixje.csv (縁の点 je を共通の起点 (step 0) の値に固定した抽出)。
前提 (欠ければ判定不能。結果を見て変えない):
  1. 4 本とも RUN_RC 0・48,000 step、残差の全行・全列が有限、系列 25 点が有限、check_convergence が実行できる (VERDICT は記録)。
  2. 線: 4 本の forge_run.log の "[lineImplicit] lines=… covered CVs=…" が同じで、被覆率 ≥ 95 % (評価断面の境界層が線に載ることの代用。直接の確認は無い)。
  3. 初期化後: 4 本の step 0 の θ(0.3/0.6/0.9) の精度間の差が |·| ≤ 0.0005 (相対)。res_0 の保存量 7 つの精度間の最大相対差を記録する。
判定 (主判定は θ(0.3)・θ(0.6)・θ(0.9)、窓 W1 = 24,000〜36,000・W2 = 36,000〜48,000):
  Δ_r(s) = (θ_r(s) − θ_r(0)) / θ_FP64(0)。e = (float の 2 本・窓の全点の Δ の平均) − (FP64 の同じ平均) (符号つき。case/45 と同じ向きは正)。
  E = (float の 2 本の窓の全点の、その平均からの最大偏差) + (FP64 の同じもの)。E は観測の範囲で信頼区間ではない。
  結果 A (再現): 同じ θ の断面で両窓とも e − E ≥ 0.005。その断面で窓の中に je が step 0 から変わった run があれば、je を固定した抽出でも A が保たれること。
  結果 B (再現せず): 3 断面・両窓とも |e| + E ≤ 0.0005。
  それ以外は判別不能。逆符号の大差・q_w や Cf・δ* だけの大差は記録し、A にも B にも読み替えない。
記録: check_quasisteady (窓 24,000〜48,000、θ・δ* は drift/osc 1e-4、Cf・q_w は 2e-4)、§6.21 の旧判定器の結果 (c48_li_judge_old.json)、方向別 dt のログ。"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parents[1]
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
F32 = ("run_0061_li_f32a", "run_0062_li_f32b"); F64 = ("run_0063_li_f64a", "run_0064_li_f64b")
XS = (0.3, 0.6, 0.9); TH = [f"theta_{x}" for x in XS]
REC = [f"{q}_{x}" for x in XS for q in ("dstar", "Cf", "qw")] + ["CD", "HF"]
WIN = {"W1": (24000, 36000), "W2": (36000, 48000)}
Q7 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega")
OUT = {"undecidable": [], "pre": {}, "table": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def rd(run, name):
    p = HERE / run / name
    if not p.exists(): return None
    rows = list(csv.DictReader(open(p)))
    return {int(r["step"]): {k: float(v) for k, v in r.items() if k != "step"} for r in rows}
S, SF, lines = {}, {}, {}
for r in F32 + F64:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}")
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 48000: und(f"{r}: 48,000 step に届かない")
    S[r] = rd(r, "prec_series.csv"); SF[r] = rd(r, "prec_series_fixje.csv")
    for nm, s in (("prec_series.csv", S[r]), ("prec_series_fixje.csv", SF[r])):
        if s is None or len(s) != 25 or any(not math.isfinite(v) for row in s.values() for v in row.values()): und(f"{r}: {nm} が 25 点そろって有限でない")
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(d)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    if not ov: und(f"{r}: check_convergence が実行できない")
    OUT["pre"][f"{r} check_convergence"] = ov[-1] if ov else None
    log = (d / "forge_run.log").read_text(errors="replace")
    m = re.search(r"\[lineImplicit\] lines=(\d+)\s+covered CVs=(\d+)/(\d+) \(([\d.]+)%\)\s+maxLen=(\d+)", log)
    lines[r] = m.groups() if m else None
    OUT["records"][f"{r} 方向別 dt のログ"] = [l[:200] for l in log.splitlines() if "lineDt" in l or "Directional" in l][:6]
# 2. 線
OUT["pre"]["線 (lines, covered, total, %, maxLen)"] = lines
if any(v is None for v in lines.values()): und("線の記録 ([lineImplicit] lines=…) が無い run がある")
elif len(set(lines.values())) != 1: und(f"線の接続が 4 本で同じでない: {lines}")
elif float(next(iter(lines.values()))[3]) < 95.0: und(f"線の被覆率が 95 % 未満: {next(iter(lines.values()))}")
# 3. 初期化後
if not OUT["undecidable"]:
    t0 = {k: (np.mean([S[r][0][k] for r in F32]), np.mean([S[r][0][k] for r in F64])) for k in TH}
    OUT["pre"]["step 0 の θ の精度間の相対差"] = {k: (a - b) / b for k, (a, b) in t0.items()}
    if any(abs((a - b) / b) > 0.0005 for a, b in t0.values()): und(f"初期化後の θ の精度間の差が 0.0005 を超える: {OUT['pre']['step 0 の θ の精度間の相対差']}")
    with h5py.File(HERE / F32[0] / "res_0.h5", "r") as ha, h5py.File(HERE / F64[0] / "res_0.h5", "r") as hb:
        OUT["records"]["res_0 の保存量の精度間の最大相対差"] = {q: float(np.max(np.abs(np.asarray(ha["VALUE/" + q][:], float) - np.asarray(hb["VALUE/" + q][:], float)) / np.maximum(np.abs(np.asarray(hb["VALUE/" + q][:], float)), 1e-30))) for q in Q7 if q in ha["VALUE"]}
verdict = None
def ev(ser, k):
    base = np.mean([ser[r][0][k] for r in F64]); out = {}
    for w, (w0, w1) in WIN.items():
        pts = {g: np.array([(ser[r][s][k] - ser[r][0][k]) / base for r in runs for s in ser[r] if w0 <= s <= w1]) for g, runs in (("f32", F32), ("f64", F64))}
        mu = {g: v.mean() for g, v in pts.items()}
        E = sum(np.max(np.abs(v - mu[g])) for g, v in pts.items())
        out[w] = dict(e=float(mu["f32"] - mu["f64"]), E=float(E), mu32=float(mu["f32"]), mu64=float(mu["f64"]))
    return out
if not OUT["undecidable"]:
    ra, rb = {}, []
    print("量          窓  e (float − FP64)   E          e − E      |e| + E    ")
    for k in TH + REC:
        t = ev(S, k); OUT["table"][k] = t
        for w, v in t.items():
            print(f"{k:11s} {w}  {v['e']:+.6f}   {v['E']:.6f}   {v['e'] - v['E']:+.6f}   {abs(v['e']) + v['E']:.6f}{'' if k in TH else '  (記録)'}")
        if k in TH:
            ra[k] = all(v["e"] - v["E"] >= 0.005 for v in t.values()); rb.append(all(abs(v["e"]) + v["E"] <= 0.0005 for v in t.values()))
    # je の切り替わり
    for k in TH:
        x = k.split("_")[1]; jc = f"je_{x}"
        moved = [r for r in F32 + F64 if any(S[r][s][jc] != S[r][0][jc] for s in S[r] if s >= 24000)]
        OUT["records"][f"{k} の窓で je が step 0 から変わった run"] = moved
        if ra.get(k) and moved:
            tf = ev(SF, k); OUT["table"][k + " (je 固定)"] = tf
            ok = all(v["e"] - v["E"] >= 0.005 for v in tf.values())
            print(f"  {k}: je が変わった run {moved} → je 固定の抽出で A {'保たれる' if ok else '保たれない'}")
            ra[k] = ok
    big_rev = [k for k in TH + REC for v in OUT["table"][k].values() if -v["e"] - v["E"] >= 0.005]
    OUT["records"]["逆符号の大差 (−e − E ≥ 0.005)"] = sorted(set(big_rev))
    OUT["records"]["θ 以外の大差 (|e| − E ≥ 0.005、記録)"] = sorted({k for k in REC for v in OUT["table"][k].values() if abs(v["e"]) - v["E"] >= 0.005})
    verdict = (f"結果 A (再現): {[k for k, v in ra.items() if v]} で両窓とも e − E ≥ 0.005" if any(ra.values())
               else "結果 B (再現せず): 3 断面・両窓とも |e| + E ≤ 0.0005" if all(rb) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.37: {verdict}")
# 記録: check_quasisteady (窓)、旧判定器
for r in F32 + F64:
    d = HERE / r; sc = d / "li_window_series.csv"
    if S.get(r):
        with open(sc, "w", newline="") as f:
            w = csv.writer(f); keys = TH + REC; w.writerow(["step"] + keys)
            for s in sorted(x for x in S[r] if x >= 24000): w.writerow([s] + [S[r][s][k] for k in keys])
        q = [subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(d), "--series-csv", str(sc), "--series-cols", ",".join(cs),
                             "--drift", tol, "--osc", tol, "--tail", "1.0"], capture_output=True, text=True).stdout for cs, tol in
             ((TH + [f"dstar_{x}" for x in XS], "0.0001"), ([f"{q}_{x}" for x in XS for q in ("Cf", "qw")] + ["CD", "HF"], "0.0002"))]
        OUT["records"][f"{r} check_quasisteady"] = [l for t in q for l in t.splitlines() if "VERDICT" in l or "STEADY" in l or "DRIFT" in l or "OSC" in l][:20]
p = subprocess.run([sys.executable, str(HERE / "tools/c48_prec_judge.py"), "--arms", f"f32={','.join(F32)};f64={','.join(F64)}", "--out", "c48_li_judge_old.json"], capture_output=True, text=True, cwd=HERE)
OUT["records"]["§6.21 の旧判定器 (記録)"] = [l for l in p.stdout.splitlines() if "VERDICT" in l][-1:] or p.stderr[-300:]
(HERE / "c48_li_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", HERE / "c48_li_judge.json")
