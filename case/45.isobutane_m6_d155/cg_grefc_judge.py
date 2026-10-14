"""plan tooling-nozzle-core-grid §4.36 の判定 (対照 = Gref の B0 の継続、B = Grefc)。
usage: python3 cg_grefc_judge.py <対照の run> <Grefc の run> → _band_ab/core_grid/grefc_judge.json
主判定 (残差の減衰): 縮流部 x∈[−8, −2] r_t の単位体積あたりの残差の体積重み付き RMS (cg_region_res.py)。窓 W1 = 20k〜40k・W2 = 40k〜60k の出力の平均。
  対照の停滞 = ρ・ρE とも W2/W1 ≥ 1/2。
  結果 A: 対照が停滞し、両窓・ρ と ρE の両方で Grefc/対照 ≤ 1/3 → 第 1 仮説を支持 (「平均流の停滞の緩和にスロートの粗化は必須」を退ける)。
         SST (roK・roOmega) が Grefc で悪化していれば「平均流のみ改善、連成系は未解決」と付ける。
  結果 B: 対照が停滞し、両窓・両成分で Grefc/対照 ≥ 2/3 → 第 1 仮説を退ける。
  それ以外 (対照も下がる・割れる・移送の過渡が支配) → 判定不能。延長しない。
副判定 (参照値への到達、記録): Grefc の 8 量の 40k〜60k の窓平均と point の止まった Gref (GrefP3 の 180k〜200k) の差 E、Grefc の U・d、check_quasisteady の STEADY。
  到達 = 8 量とも STEADY かつ |E| + U_c + U_ref + p_c + p_ref ≤ 0.1 % (相対) / 0.005 %pt (軸・η 0.1 の M)。一度横切っただけでは到達と呼ばない。"""
import csv, glob, json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"; TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
A, B = sys.argv[1:3]
REQ8 = ("theta_r_w65_75", "theta_r_w84_94", "delta_loc_w65_75", "delta_loc_w84_94", "Q_w", "exit_M_line", "eta0.0_mean_pct", "eta0.1_mean_pct")
P = {"theta": 0.010, "delta": 0.023, "Q_w": 0.009, "exit_M_line": 0.0006, "eta_pctpt": 0.0014}     # Gref と同じ (Grefc の試験部の列は Gref とビット一致、Q_w の x の差 +0.0001 %)
def kind(k):
    return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"
def region(run):
    rows = [json.loads(l) for l in open(OUTD / f"region_res_{run}.jsonl")]
    return sorted({r["step"]: r for r in rows}.values(), key=lambda r: r["step"])
W = ((20000, 40000), (40000, 60000))
RA, RB = region(A), region(B); out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.36", "A": A, "B": B, "region": {}}
def wmean(R, key, lo, hi):
    v = [r[key] for r in R if lo <= r["step"] <= hi]; return (float(np.mean(v)), len(v)) if v else (None, 0)
nonfinite = any(r["nonfinite"] for r in RA + RB)
stall = True; ratios = {}
for comp in ("res_ro", "res_roe", "res_roUx", "res_roUy", "res_roK", "res_roOmega"):
    key = f"contraction_{comp}"; ra = [wmean(RA, key, *w) for w in W]; rb = [wmean(RB, key, *w) for w in W]
    rr = [(b[0] / a[0]) if a[0] and b[0] is not None else None for a, b in zip(ra, rb)]
    out["region"][comp] = {"A": [x[0] for x in ra], "B": [x[0] for x in rb], "n": [x[1] for x in ra], "B_over_A": rr,
                           "A_W2_over_W1": (ra[1][0] / ra[0][0]) if ra[0][0] else None,
                           "all_B_over_A": [(wmean(RB, f"all_{comp}", *w)[0] / wmean(RA, f"all_{comp}", *w)[0]) for w in W]}
    ratios[comp] = rr
    if comp in ("res_ro", "res_roe"):
        stall &= (out["region"][comp]["A_W2_over_W1"] or 0) >= 0.5
mean_ok = all(r is not None and r <= 1 / 3 for c in ("res_ro", "res_roe") for r in ratios[c])
mean_bad = all(r is not None and r >= 2 / 3 for c in ("res_ro", "res_roe") for r in ratios[c])
sst_worse = any(r is not None and r > 1 for c in ("res_roK", "res_roOmega") for r in ratios[c])
out["result"] = ("発散 (非有限)" if nonfinite else "判定不能 (対照が停滞していない)" if not stall else
                 ("A: 縮流部の局所の粗化だけで平均流の停滞が大きく緩和 (第 1 仮説を支持、スロートの粗化は必須でない)" + ("。ただし SST は Grefc で悪化 (平均流のみ改善、連成系は未解決)" if sst_worse else "")) if mean_ok else
                 "B: 縮流部の局所の粗化では 6 万 step 内に大きく改善しない (第 1 仮説を退ける)" if mean_bad else "判定不能 (中間・割れる)")
# 副判定: 参照値への到達
def merged(run):
    rows = []
    for f in sorted(glob.glob(str(OUTD / f"gref_series_{run}_*.csv")), key=lambda p: int(p.rsplit("_", 1)[1][:-4])):
        rows += list(csv.DictReader(open(f)))
    with open(OUTD / f"gref_series_{run}.csv", "w") as fh:
        fh.write(",".join(rows[0].keys()) + "\n"); [fh.write(",".join(r.values()) + "\n") for r in rows]
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
def load(nm):
    rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv"))); return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
def tm(D, k, lo, hi):
    m = (D["step"] >= lo) & (D["step"] <= hi); return float(np.trapezoid(D[k][m], D["step"][m]) / (hi - lo))
SB = merged(B); merged(A); P3 = load("GrefP3"); GA = load("GrefA")
reach = True; out["reach"] = {}
tmp = OUTD / "_qs_grefc.csv"; m = (SB["step"] >= 20000) & (SB["step"] <= 60000)
with open(tmp, "w") as fh:
    fh.write("step," + ",".join(REQ8) + "\n")
    for i in np.where(m)[0]:
        fh.write(f"{SB['step'][i]:.0f}," + ",".join(f"{(1 + SB[k][i] / 100) if k.startswith('eta') else SB[k][i]:.12g}" for k in REQ8) + "\n")
for k in REQ8:
    D0 = abs(tm(GA, k, 40000, 60000)); thr = 0.005 if kind(k) == "eta_pctpt" else 0.1
    sc = (lambda x: x) if kind(k) == "eta_pctpt" else (lambda x, D0=D0: 100 * x / D0)
    mu_b = tm(SB, k, 40000, 60000); mu_r = tm(P3, k, 180000, 200000)
    mr = (P3["step"] >= 160000) & (P3["step"] <= 200000)
    U_b = sc(float(np.abs(SB[k][m] - mu_b).max())); U_r = sc(float(np.abs(P3[k][mr] - mu_r).max())); d_b = sc(abs(tm(SB, k, 20000, 40000) - mu_b))
    rel = (thr / 100) / abs(1 + np.mean(SB[k][m]) / 100) if kind(k) == "eta_pctpt" else (thr / 100) * D0 / abs(np.mean(SB[k][m]))
    pr = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k, "--tail", "1.0",
                         "--drift", f"{rel / 8:.6g}", "--osc", f"{rel / 2:.6g}"], capture_output=True, text=True)
    E = sc(mu_b - mu_r); tot = abs(E) + U_b + U_r + 2 * P[kind(k)]
    out["reach"][k] = {"E": E, "U_c": U_b, "U_ref": U_r, "d_c": d_b, "steady": pr.returncode == 0, "budget": tot, "thr": thr,
                       "B0_ref_minus_point": sc(tm(GA, k, 40000, 60000) - mu_r)}
    reach &= pr.returncode == 0 and tot <= thr
tmp.unlink()
out["reach_result"] = "到達 (8 量 STEADY・予算内)" if reach else "到達とは言えない"
(OUTD / "grefc_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for c, v in out["region"].items():
    print(f"縮流部 {c:11s} 対照 {v['A']} Grefc {v['B']} Grefc/対照 {v['B_over_A']} 対照 W2/W1 {v['A_W2_over_W1']} (全域の比 {v['all_B_over_A']})")
for k, v in out["reach"].items():
    print(f"  {k:16s} E {v['E']:+.4f} U_c {v['U_c']:.4f} U_ref {v['U_ref']:.4f} d_c {v['d_c']:.4f} STEADY {v['steady']} 予算 {v['budget']:.4f} (閾値 {v['thr']}) (B0 の Gref − point {v['B0_ref_minus_point']:+.3f})")
print("結果:", out["result"]); print("副判定:", out["reach_result"])
