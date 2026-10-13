"""plan tooling-nozzle-core-grid §4.96 P-ts の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 thermofloat-stop): 熱物性の経路が B0 の止まる先を変えるか。
4 本 (T1a・T0a・T0b・T1b、各 60,000 step、1000 ごと) の 8 量の系列 (gref_series_<run>_<n>.csv) から、主指標 θ_r[65,75]、補助 θ_r[84,94]・Q_w。
百分率の分母は point の止まった値 (grefx_refpoint_judge.json の ref_point の μ)。窓 W1 = (20k, 40k]・W2 = (40k, 60k] (各 20 出力)。
D_w = (T0 の 2 本の窓平均の平均 − T1 の 2 本の平均) の百分率。F = 両経路・両窓の同じ経路の複製の差の最大。
u (run・窓) = max(|窓の前半 10k の平均 − 後半 10k の平均|、窓の末尾が単調なら classify_series の漸近値と窓平均の差)。U = max(u_T0) + max(u_T1)。H = F + U。
資格: 入力の sha が 4 本で同じ・restart ビット一致・実効の thermoFloat・非有限なし・ρ・P・T > 0・完走・各 run の各窓で主指標が STEADY (単調なら漸近値が出せる)・主指標で H ≤ 0.025 %。
結果 A: 主指標で両窓の D_w の符号が同じかつ |D_w| − H ≥ 0.1 %。結果 B: 主指標で両窓 |D_w| + H < 0.1 %。それ以外・資格不足は判定不能。補助量にも同じ判定 (記述)。
usage: python3 cg_tfs_judge.py → _band_ab/core_grid/tfs_judge.json"""
import csv, hashlib, json, re, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
sys.path.insert(0, str(HERE / ".." / ".." / "solver_density_cuda" / "tools"))
from check_quasisteady import classify_series, _monotone_limit
RUNS = {"T1": ["run_0905_tfs_t1a", "run_0908_tfs_t1b"], "T0": ["run_0906_tfs_t0a", "run_0907_tfs_t0b"]}
TF = {"T1": "1", "T0": "0"}; W = ((20000, 40000), (40000, 60000)); Q = ("theta_r_w65_75", "theta_r_w84_94", "Q_w")
ref = json.loads((OUTD / "grefx_refpoint_judge.json").read_text())["stats"]["ref_point"]["q"]
MU = {q: float(ref[q]["mu"]) for q in Q}
def series(run):
    rows = {}
    for n in range(1000, 60001, 1000):
        p = OUTD / f"gref_series_{run}_{n}.csv"
        if p.exists():
            r = list(csv.DictReader(open(p)))
            if r: rows[n] = r[-1]
    return rows
def logt(run): return "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists())
def cfg(run, key):
    m = re.search(r"'%s' in '[^']*': ([^\s]+)" % key, (HERE / run / "forge_run.log").read_text(errors="replace")); return m.group(1) if m else None
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
S = {r: series(r) for path in RUNS for r in RUNS[path]}
allruns = [r for path in RUNS for r in RUNS[path]]
qual = {"same_input": len({sha(HERE / r / "nozzle.h5") for r in allruns}) == 1,
        "restart_bitwise": "VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / RUNS["T1"][0] / "restart.log").read_text(),
        "thermoFloat_effective": all((cfg(r, "thermoFloat") or "").startswith(TF[p]) for p in RUNS for r in RUNS[p]),
        "completed": all((HERE / r / "RUN_RC").exists() and (HERE / r / "RUN_RC").read_text().strip() == "0" and len(S[r]) == 60 for r in allruns),
        "finite_physical": True}
for r in allruns:
    for n, row in S[r].items():
        if row.get("nonfinite", "0") not in ("0", "0.0"): qual["finite_physical"] = False
    p = OUTD / f"phys_{r}.jsonl"
    if p.exists():
        for l in open(p):
            x = json.loads(l)
            if x.get("nonfinite", 0) or not (x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0): qual["finite_physical"] = False
def win(run, q, lo, hi):
    pts = [(n, float(S[run][n][q])) for n in sorted(S[run]) if lo < n <= hi]
    return np.array([a for a, _ in pts], float), np.array([b for _, b in pts], float)
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.96", "mu_point": MU, "q": {}}
steady_ok = True
for q in Q:
    rec = {"windows": {}}
    for wi, (lo, hi) in enumerate(W):
        mean = {}; u = {}; verdict = {}
        for path in RUNS:
            for r in RUNS[path]:
                st, v = win(r, q, lo, hi)
                if len(v) != 20: mean[r] = None; u[r] = None; verdict[r] = "points"; continue
                mean[r] = float(v.mean())
                half = (st <= lo + 10000)
                u1 = abs(v[half].mean() - v[~half].mean())
                vd, det, _ = classify_series(st, v, 0.4, 0.05, 0.10, 4)
                mono = _monotone_limit(st, v)
                u2 = 0.0; unstable = False
                if mono is not None:
                    lim, _, rel = mono
                    d = np.diff(v[-max(5, int(np.ceil(0.5 * len(v)))):])
                    with np.errstate(divide="ignore", invalid="ignore"): rr = d[1:] / d[:-1]
                    rr = rr[np.isfinite(rr) & (rr > 0) & (rr < 1.0)]
                    unstable = rr.size == 0          # 増分が減衰しない = 漸近値を出せない
                    u2 = abs(lim - mean[r])
                u[r] = 100 * max(u1, u2) / MU[q]; verdict[r] = vd + (" (漸近値を出せない)" if unstable else "")
                if q == Q[0] and (vd != "STEADY" or unstable): steady_ok = False
        if any(v is None for v in mean.values()):
            rec["windows"][f"W{wi+1}"] = {"missing": True}; continue
        m0 = np.mean([mean[r] for r in RUNS["T0"]]); m1 = np.mean([mean[r] for r in RUNS["T1"]])
        D = 100 * (m0 - m1) / MU[q]
        F = 100 * max(abs(mean[RUNS[p][0]] - mean[RUNS[p][1]]) for p in RUNS) / MU[q]
        rec["windows"][f"W{wi+1}"] = {"D": D, "F": F, "u": u, "verdict": verdict, "mean": mean,
                                       "T0_minus_point_pct": 100 * (m0 - MU[q]) / MU[q], "T1_minus_point_pct": 100 * (m1 - MU[q]) / MU[q]}
    out["q"][q] = rec
# F は両窓の最大、U は経路ごとの u の最大の和
for q, rec in out["q"].items():
    ws = [w for w in rec["windows"].values() if not w.get("missing")]
    if len(ws) == 2:
        F = max(w["F"] for w in ws); U = sum(max(w["u"][r] for w in ws for r in RUNS[p]) for p in RUNS)
        rec["F"], rec["U"], rec["H"] = F, U, F + U
        Ds = [w["D"] for w in ws]
        rec["verdict"] = ("A" if (np.sign(Ds[0]) == np.sign(Ds[1]) and all(abs(d) - (F + U) >= 0.1 for d in Ds)) else
                          "B" if all(abs(d) + (F + U) < 0.1 for d in Ds) else "判定不能")
qual["primary_steady"] = steady_ok
qual["primary_H_le_0.025pct"] = "H" in out["q"][Q[0]] and out["q"][Q[0]]["H"] <= 0.025
pv = out["q"][Q[0]].get("verdict")
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif pv == "A":
    res = "A: この起点・B0 では、熱物性の経路が θ_r[65,75] の準定常の平均を実行間差以上に (≥ 0.1 %) 変えた"
elif pv == "B":
    res = "B: この条件・期間では、θ_r[65,75] に 0.1 % 以上の経路依存を検出しなかった" + ("" if all(out["q"][q].get("verdict") == "B" for q in Q[1:]) else " (補助量の一部は B でない: 「経路は効かない」と総括しない)")
else:
    res = "判定不能 (中間・符号が窓で変わる)"
out["qualification"] = qual; out["result"] = res
(OUTD / "tfs_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for q, rec in out["q"].items():
    for wn, w in rec["windows"].items():
        if w.get("missing"): print(q, wn, "点が足りない"); continue
        print(f"{q:15s} {wn} D {w['D']:+.4f} %  F {w['F']:.4f}  T0−point {w['T0_minus_point_pct']:+.3f} %  T1−point {w['T1_minus_point_pct']:+.3f} %  u " + " ".join(f"{k[-3:]} {v:.4f}" for k, v in w["u"].items()) + "  " + " ".join(f"{k[-3:]} {v}" for k, v in w["verdict"].items()))
    print(f"{q:15s} F {rec.get('F', float('nan')):.4f}  U {rec.get('U', float('nan')):.4f}  H {rec.get('H', float('nan')):.4f}  → {rec.get('verdict')}")
print("資格:", qual); print("結果:", res)
