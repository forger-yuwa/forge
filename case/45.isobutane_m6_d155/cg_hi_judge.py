"""plan tooling-nozzle-core-grid §4.98 P-hi の判定 (事前登録、ユーザ提案・codex 諮問 2026-10-14 half-inlet): 入口半径を半分にした形 H 対 今の形 C。
本段 (C = run_0912_hi_c_main、H = run_0915_hi_h_main、1000 ごと) の縮流部 x/r_t ∈ [−8, −2] の q = √(ΣV(R/V)²/ΣV) (ρ・ρE)。
窓 W1 = (20k, 40k]・W2 = (40k, 60k] (各 20 出力) の平均 q̄、q̄_H/q̄_C と各腕の d = q̄(W2)/q̄(W1)。
資格: (1) C の d ≥ 0.5 (ρ・ρE、対照が停滞を再現)。(2) 6 本とも完走・非有限なし・ρ・P・T > 0。
  (3) 本段の最後の場で、主流 (壁からの距離 ≥ 局所半径の 20 %) の M の中央値の比 H/C が 入口 x ∈ [−12.5, −12] で ≥ 3、x ∈ [−6, −3.25] で ≥ 1.2。
  (4) H の縮流部の ρ の q の系列と θ_r[65,75]・Q_w の系列が W1 ∪ W2 で classify_series (既定) の DRIFTING・TRANSIENT-UNSETTLED でない。
結果 A: 両窓・ρ・ρE で q̄_H/q̄_C ≤ 1/3、H の d ≤ 1.1、かつ両窓で縮流部の ρUx・ρUy・ρk・ρω と縮流部の外の ρ・ρUx・ρUy・ρE の q̄_H/q̄_C ≤ 1.5。
結果 B: 両窓・ρ・ρE で q̄_H/q̄_C ≥ 2/3、かつ H の d ≥ 0.8。それ以外・資格不足は判定不能。
usage: python3 cg_hi_judge.py → _band_ab/core_grid/hi_judge.json"""
import csv, json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
sys.path.insert(0, str(HERE / ".." / ".." / "solver_density_cuda" / "tools"))
from check_quasisteady import classify_series
C, H = "run_0912_hi_c_main", "run_0915_hi_h_main"
ALL = ["run_0910_hi_c_s1", "run_0911_hi_c_s2", C, "run_0913_hi_h_s1", "run_0914_hi_h_s2", H]
W = ((20000, 40000), (40000, 60000))
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
R = {r: jl(f"region_res_{r}.jsonl") for r in ALL}; PH = {r: jl(f"phys_{r}.jsonl") for r in ALL}
def wm(run, key, lo, hi):
    v = [x[key] for x in R[run] if lo < x["step"] <= hi and key in x]; return float(np.mean(v)) if len(v) == 20 else None
qual = {}
qual["completed_finite_physical"] = all((HERE / r / "RUN_RC").exists() and (HERE / r / "RUN_RC").read_text().strip() == "0" for r in ALL) \
    and not any(x["nonfinite"] for r in ALL for x in R[r] + PH[r]) and all(x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for r in ALL for x in PH[r])
d = {run: {c: (wm(run, f"contraction_{c}", *W[1]) / wm(run, f"contraction_{c}", *W[0])) if wm(run, f"contraction_{c}", *W[0]) else None for c in ("res_ro", "res_roe")} for run in (C, H)}
qual["control_stall"] = all(v is not None and v >= 0.5 for v in d[C].values())
# (3) マッハ数
def mach(run):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / "res_60000.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("Ux", "Uy", "sonic", "wall_dist")}
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    nj = xy.shape[0] // ni; X = (xy[:, 0] / S).reshape(ni, nj)[:, ::-1]; Rr = (xy[:, 1] / S).reshape(ni, nj)[:, ::-1]
    M = (np.hypot(v["Ux"], v["Uy"]) / v["sonic"]).reshape(ni, nj)[:, ::-1]; wd = (v["wall_dist"] / S).reshape(ni, nj)[:, ::-1]
    rw = Rr[:, 0][:, None]; core = wd >= 0.2 * rw
    out = {}
    for lab, (a, b) in {"inlet": (-12.5, -12.0), "stallband_x": (-6.0, -3.25), "x-6": (-6.25, -5.75), "x-4": (-4.25, -3.75), "x-3.25": (-3.5, -3.0)}.items():
        m = (X >= a) & (X <= b) & core; out[lab] = float(np.median(M[m]))
    return out
MC, MH = mach(C), mach(H)
mr = {k: MH[k] / MC[k] for k in MC}
qual["mach_increase"] = mr["inlet"] >= 3 and mr["stallband_x"] >= 1.2
# (4) 過渡
def series(run, q):
    st, v = [], []
    for n in range(21000, 60001, 1000):
        p = OUTD / f"gref_series_{run}_{n}.csv"
        if p.exists():
            r = list(csv.DictReader(open(p)))
            if r: st.append(n); v.append(float(r[-1][q]))
    return st, v
tr = {}
st = [x["step"] for x in R[H] if 20000 < x["step"] <= 60000]; vv = [x["contraction_res_ro"] for x in R[H] if 20000 < x["step"] <= 60000]
tr["contraction_res_ro"] = classify_series(st, vv, 0.4, 0.05, 0.10, 4)[0]
for q in ("theta_r_w65_75", "Q_w"):
    s_, v_ = series(H, q); tr[q] = classify_series(s_, v_, 0.4, 0.05, 0.10, 4)[0] if len(v_) == 40 else "points"
qual["H_not_transient"] = all(v not in ("DRIFTING", "TRANSIENT-UNSETTLED", "NONFINITE", "points") for v in tr.values())
# 比
main = ("res_ro", "res_roe")
side = [("contraction", c) for c in ("res_roUx", "res_roUy", "res_roK", "res_roOmega")] + [("outside", c) for c in ("res_ro", "res_roUx", "res_roUy", "res_roe")]
ratio = {}
for wi, (lo, hi) in enumerate(W):
    rr = {}
    for rg, c in [("contraction", m) for m in main] + side + [("throatwall", m) for m in main]:
        a, b = wm(H, f"{rg}_{c}", lo, hi), wm(C, f"{rg}_{c}", lo, hi)
        rr[f"{rg}_{c}"] = (a / b) if (a is not None and b) else None
    ratio[f"W{wi+1}"] = rr
missing = any(ratio[w][f"contraction_{c}"] is None for w in ratio for c in main)
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif missing:
    res = "判定不能 (証拠不足)"
else:
    mr_ = [ratio[w][f"contraction_{c}"] for w in ratio for c in main]
    sr = [ratio[w][f"{rg}_{c}"] for w in ratio for rg, c in side]
    if all(v <= 1 / 3 for v in mr_) and all(v is not None and v <= 1.1 for v in d[H].values()) and all(v is not None and v <= 1.5 for v in sr):
        res = "A: この形の変更 (入口半径の半減) で、縮流部の B0 の停滞が大きく緩んだ"
    elif all(v >= 2 / 3 for v in mr_) and all(v is not None and v >= 0.8 for v in d[H].values()):
        res = "B: 入口半径の半減で停滞を大きく緩められる、を退ける (低マッハは無関係とは言わない)"
    else:
        res = "判定不能 (中間・成分で割れる・悪化の移動)"
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.98", "C": C, "H": H, "qualification": qual, "d": d, "ratio_H_over_C": ratio,
       "mach_C": MC, "mach_H": MH, "mach_ratio": mr, "H_transient_verdicts": tr, "result": res,
       "series_contraction_ro": {run: [(x["step"], x["contraction_res_ro"]) for x in R[run] if x["step"] % 5000 == 0] for run in (C, H)}}
(OUTD / "hi_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
fmt = lambda v: "-" if v is None else f"{v:.3g}"
for w, rr in ratio.items(): print(w, "  ".join(f"{k} {fmt(v)}" for k, v in rr.items()))
print("d:", {r: {c: fmt(v) for c, v in dd.items()} for r, dd in d.items()})
print("M (主流の中央値) C:", {k: round(v, 4) for k, v in MC.items()}); print("M H:", {k: round(v, 4) for k, v in MH.items()}); print("比:", {k: round(v, 2) for k, v in mr.items()})
print("H の過渡:", tr); print("縮流部 ρ の q:", out["series_contraction_ro"])
print("資格:", qual); print("結果:", res)
