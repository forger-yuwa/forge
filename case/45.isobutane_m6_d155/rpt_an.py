"""plan architecture-float-state-double-geometry §6.48 (事前登録、codex diagnose 2026-10-11) の判定: point の陰解法に固定した float と FP64 の比較 (案 P)。B で回す。
D1・D2 = run_0713・0714 (FP64)、F1・F2 = run_0715・0716 (float)。使い方: rpt_an.py --stage 1 (30,000 step、窓 20,000〜25,000・25,000〜30,000) / --stage 2 (120,000 step、
窓 100,000〜110,000・110,000〜120,000)。→ _band_ab/cold_pair/rpt_judge_s{1,2}.json
Δ = 100 × (窓の平均 − 自分の step 0) / |旧 A (run_0487) の step 0|、D = (+1.892, +1.647, +1.431) %。N = max(|ΔD1 − ΔD2|, |ΔF1 − ΔF2|)。
結果 A: その段の両窓・3 断面で min_{i,j} sign(D)(ΔF_i − ΔD_j) − N ≥ 0.5|D|。
結果 B (段 2 だけ): 両窓・3 断面で max_{i,j} |ΔF_i − ΔD_j| + N ≤ 0.1|D|、かつ 4 本とも check_convergence が PASS、かつ θ_r の 100,000〜120,000 の check_quasisteady (1e-4) が STEADY。
それ以外は判別不能 (段 1 で A でなければ判定は段 2 に持ち越し)。
前提: 各精度の中で res_0 の 9 保存量がビット一致 (キーの完全一致)。設定の差は旧 run (0487・0488) からの 4 キーの削除だけ。全期間の残差・使う出力が有限。"""
import argparse, csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--stage", type=int, choices=(1, 2), required=True); A = ap.parse_args()
RUNS = {"D1": "run_0713_pt_d1", "D2": "run_0714_pt_d2", "F1": "run_0715_pt_f1", "F2": "run_0716_pt_f2"}
BASE = {"D1": "run_0487_tr_fp64", "D2": "run_0487_tr_fp64", "F1": "run_0488_tr_f32", "F2": "run_0488_tr_f32"}
WIN = {1: {"W1": (20000, 25000), "W2": (25000, 30000)}, 2: {"W3": (100000, 110000), "W4": (110000, 120000)}}[A.stage]
STEPS = [0] + sorted({s for w0, w1 in WIN.values() for s in range(w0, w1 + 1, 500)})
LAST = max(w1 for _, w1 in WIN.values())
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
REMOVED = ("lineDtDirectionalCap: 50.0, ", "implicitThermalJacobian: 5, ", "lineImplicit: 1, ", "lineDtDirectional: 1, ")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"stage": A.stage, "undecidable": [], "checks": {}, "records": {}, "delta": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def q9(path):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"] if k in Q9}
for k, r in RUNS.items():
    d = HERE / r
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < LAST: und(f"{r}: {LAST} step に届かない ({int(rows[-1]['step']) + 1})")
    for s in STEPS:
        p = d / f"res_{s}.h5"
        if not p.exists(): und(f"{r}: res_{s} が無い"); break
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][q].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][q][:])) for q in h["VALUE"]): und(f"{r}: res_{s} に非有限"); break
    a_, b_ = (d / "solverConfig.yaml").read_text().splitlines(), (HERE / BASE[k] / "solverConfig.yaml").read_text().splitlines()
    diff = [(x, y) for x, y in zip(a_, b_) if x != y]
    # 許す差: deltaT の行の 4 キーの削除と、nStepOuter (旧 30,000 → 120,000、登録の時点で書き漏らした必然の差、run の前に追記)
    rest = [(x, y) for x, y in diff if not (x.strip() == "last: {nStepOuter: 120000}" and y.strip() == "last: {nStepOuter: 30000}")]
    ok = len(a_) == len(b_) and len(rest) == 1
    if ok:
        y = rest[0][1]
        for t in REMOVED: y = y.replace(t, "")
        ok = (y == rest[0][0])
    if not ok: und(f"{r}: {BASE[k]} との設定の差が 4 キーの削除だけでない: {diff[:2]}")
    OUT["checks"][f"{r} 設定の差"] = diff
if not OUT["undecidable"] or A.stage == 1:
    for g in (("D1", "D2"), ("F1", "F2")):
        a, b = q9(HERE / RUNS[g[0]] / "res_0.h5"), q9(HERE / RUNS[g[1]] / "res_0.h5")
        if set(a) != set(Q9) or set(b) != set(Q9): und(f"{g}: res_0 の保存量のキーが 9 つそろわない")
        elif any(not np.array_equal(a[q], b[q]) for q in Q9): und(f"{g}: res_0 の保存量が精度の中でビット一致しない")
    a, b = q9(HERE / RUNS["D1"] / "res_0.h5"), q9(HERE / RUNS["F1"] / "res_0.h5")
    OUT["records"]["初期化の差 (D1 と F1 の res_0、最大相対差)"] = {q: float(np.max(np.abs(a[q].astype(float) - b[q].astype(float)) / np.maximum(np.abs(a[q].astype(float)), 1e-30))) for q in Q9 if q in a and q in b}
tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
A0 = {int(k): v for k, v in tj["series"]["A"].items()}[0]
yb_x, yb = XC.common_yb()
verdict = None
if not OUT["undecidable"]:
    ser = {}
    for k, r in RUNS.items():
        ser[k] = {}
        for s in STEPS:
            x = CS.snapshot(HERE / r, s, yb_x, yb)
            if x.get("nonfinite") or any(not np.isfinite(x[q]) for q in KEYS + (QW,)): und(f"{r}: step {s} の系列に非有限"); break
            ser[k][s] = {q: float(x[q]) for q in KEYS + (QW,)}
    OUT["series"] = ser
if not OUT["undecidable"]:
    def dl(k, q, w):
        w0, w1 = WIN[w]; return 100 * (np.mean([ser[k][s][q] for s in ser[k] if w0 <= s <= w1]) - ser[k][0][q]) / abs(A0[q])
    ra, rb = [], []
    print(f"段 {A.stage}  窓  断面         ΔD1      ΔD2      ΔF1      ΔF2      N       min sgn(F−D)−N  max|F−D|+N  0.5|D|  0.1|D|")
    for w in WIN:
        for q in KEYS + (QW,):
            v = {k: dl(k, q, w) for k in RUNS}; sg, Dk = np.sign(D[q]), abs(D[q])
            N = max(abs(v["D1"] - v["D2"]), abs(v["F1"] - v["F2"]))
            pairs = [(f, d) for f in ("F1", "F2") for d in ("D1", "D2")]
            lo = min(sg * (v[f] - v[d]) for f, d in pairs) - N; hi = max(abs(v[f] - v[d]) for f, d in pairs) + N
            A_, B_ = lo >= 0.5 * Dk, hi <= 0.1 * Dk
            OUT["delta"][f"{w}|{q}"] = dict(**v, N=N, lo=lo, hi=hi, resultA=bool(A_), resultB_part=bool(B_))
            print(f"       {w}  {q:11s} {v['D1']:+7.3f}  {v['D2']:+7.3f}  {v['F1']:+7.3f}  {v['F2']:+7.3f}  {N:.3f}  {lo:+.3f}          {hi:.3f}       {0.5 * Dk:.3f}   {0.1 * Dk:.3f}  {'A' if A_ else ''}{'b' if B_ else ''}{' (副)' if q == QW else ''}")
            if q != QW: ra.append(A_); rb.append(B_)
    for k, r in RUNS.items():
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
        sc = HERE / r / f"rpt_series_s{A.stage}.csv"
        with open(sc, "w", newline="") as fh:
            wcsv = csv.writer(fh); wcsv.writerow(["step"] + list(KEYS) + [QW])
            for s in sorted(t for t in ser[k] if t > 0): wcsv.writerow([s] + [ser[k][s][q] for q in KEYS + (QW,)])
        qs = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r), "--series-csv", str(sc), "--series-cols", ",".join(KEYS),
                             "--drift", "0.0001", "--osc", "0.0001", "--tail", "1.0"], capture_output=True, text=True)
        OUT["records"][f"{r} check_quasisteady"] = [l.strip() for l in qs.stdout.splitlines() if "csv" in l or any(t in l for t in KEYS)][:4]
    if all(ra): verdict = "結果 A: point に固定しても大差が残る (第 1 仮説を支持、「B0 の構成が大差に必須」を退ける。特定の演算が原因とは言わない)"
    elif A.stage == 2 and all(rb):
        conv_ok = all("PASS" in str(OUT["records"].get(f"{r} check_convergence", "")) for r in RUNS.values())
        qs_ok = all(any("-> STEADY" in l for l in OUT["records"].get(f"{r} check_quasisteady", [])) for r in RUNS.values())
        OUT["records"]["結果 B の付帯条件 (収束・準定常)"] = [conv_ok, qs_ok]
        verdict = ("結果 B: point では大差が残らない (第 2 仮説を支持。B0 のどの設定が効いたかは未確定)" if conv_ok and qs_ok
                   else "判別不能 (差は小さいが、収束・準定常の付帯条件を満たさない)")
    else: verdict = "判別不能" + (" (段 2 へ持ち越し)" if A.stage == 1 else "")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.48 段 {A.stage}: {verdict}")
(HERE / f"_band_ab/cold_pair/rpt_judge_s{A.stage}.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print(f"→ _band_ab/cold_pair/rpt_judge_s{A.stage}.json")
