"""plan tooling-nozzle-core-grid §4.88 P-k の判定 (事前登録、codex 諮問 2026-10-13 pdt-interp): Γ = I (precondEps 1.0) で前処理のカーネル (K3) と標準のカーネル (K0)。
主指標: 縮流部 (x/r_t −8〜−2) の ρ・ρE の q。窓 W1 = 10k〜14k・W2 = 16k〜20k (各 5 出力) の平均 q̄、d = q̄(W2)/q̄(W1)。
資格: (1) restart ビット一致・設定の差は lowMachPrecond の 1 行・K3 のログにスイッチの行。(2) step 0 の全残差列の相対差 ≤ 1e−6。
  (3) 島の step 1000 の dt_local(K3)/dt_local(K0) の中央値 0.98〜1.02。(4) precondEps の実効値 1・K3 の lowMachPrecond 3。(5) 非有限なし・ρ・P・T > 0・完走。
結果 A: 両窓で縮流部の ρ・ρE の q̄_K3/q̄_K0 ≥ 2。
結果 B: 両窓で縮流部の ρ・ρE の q̄_K3/q̄_K0 と d_K3/d_K0 がともに 0.9〜1.1、かつ throatwall・縮流部の外の全保存量の q̄_K3/q̄_K0 ≤ 1.5。
それ以外・資格不足は判定不能。
usage: python3 cg_pk_judge.py <K0 の run> <K3 の run> → _band_ab/core_grid/pk_judge.json"""
import csv, json, re, subprocess, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
K0, K3 = sys.argv[1:3]; W = ((10000, 14000), (16000, 20000))
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
R = {K0: jl(f"region_res_{K0}.jsonl"), K3: jl(f"region_res_{K3}.jsonl")}; PH = {K0: jl(f"phys_{K0}.jsonl"), K3: jl(f"phys_{K3}.jsonl")}
def wm(run, key, lo, hi):
    v = [r[key] for r in R[run] if lo <= r["step"] <= hi and key in r]; return float(np.mean(v)) if len(v) == 5 else None
def logtext(run): return "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists())
def cfg_val(run, key):
    m = re.search(r"'%s' in '[^']*': ([^\s]+)" % key, (HERE / run / "forge_run.log").read_text(errors="replace")); return m.group(1) if m else None
bit = all("VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / r / "restart.log").read_text() for r in (K0, K3))
dlines = [l for l in subprocess.run(["diff", str(HERE / K0 / "solverConfig.yaml"), str(HERE / K3 / "solverConfig.yaml")], capture_output=True, text=True).stdout.splitlines() if l.startswith(">")]
q1 = bit and len(dlines) == 1 and "lowMachPrecond: 3, " in dlines[0] and logtext(K3).count("FORGE_DIAG_PRECOND_NO_DT_SCALE=1:") >= 1
def step0(run):
    rows = list(csv.reader(open(HERE / run / "residual_history.csv"))); h = rows[0]
    r0 = next(r for r in rows[1:] if r[0] == "0" and r[1] == "-1")
    return {h[i]: float(r0[i]) for i in range(len(h)) if h[i].startswith("rms_")}
a0, b0 = step0(K0), step0(K3)
rel0 = {k: (abs(a0[k] - b0[k]) / max(abs(a0[k]), abs(b0[k])) if max(abs(a0[k]), abs(b0[k])) > 0 else 0.0) for k in a0}
q2 = all(v <= 1e-6 for v in rel0.values())
def load(run, n, keys):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}; x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
    nj = x.size // ni
    return {k: a.reshape(ni, nj)[:, ::-1] for k, a in v.items()}, x.reshape(ni, nj)[:, ::-1], ni, nj
g0, X, ni, nj = load(K0, 1000, ("dt_local",)); g3, _, _, _ = load(K3, 1000, ("dt_local",))
isl = np.zeros((ni, nj), bool); cols = (X[:, 0] >= -6.0) & (X[:, 0] <= -3.25); isl[np.ix_(cols, np.arange(60, 101))] = True
dtr = float(np.median(g3["dt_local"][isl] / g0["dt_local"][isl])); q3 = 0.98 <= dtr <= 1.02
eff = {"precondEps_K0": cfg_val(K0, "precondEps"), "precondEps_K3": cfg_val(K3, "precondEps"), "lowMachPrecond_K3": cfg_val(K3, "lowMachPrecond")}
q4 = eff["lowMachPrecond_K3"] == "3" and all(v is not None and abs(float(v) - 1.0) < 1e-12 for v in (eff["precondEps_K0"], eff["precondEps_K3"]))
q5 = not any(r["nonfinite"] for run in (K0, K3) for r in R[run] + PH[run]) and all(x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for run in (K0, K3) for x in PH[run]) \
    and all((HERE / r / "RUN_RC").exists() and (HERE / r / "RUN_RC").read_text().strip() == "0" for r in (K0, K3))
qual = {"same_start_config_switch": bool(q1), "same_rhs_step0": bool(q2), "same_dt": bool(q3), "effective_values": bool(q4), "finite_physical_completed": bool(q5)}
comps = ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega")
tab = {}
for rg in ("contraction", "throatwall", "outside"):
    for c in comps:
        k = f"{rg}_{c}"; w = {run: [wm(run, k, *W[0]), wm(run, k, *W[1])] for run in (K0, K3)}
        lv = [(w[K3][i] / w[K0][i]) if (w[K3][i] is not None and w[K0][i]) else None for i in range(2)]
        d = {run: (w[run][1] / w[run][0]) if (w[run][0] and w[run][1] is not None) else None for run in (K0, K3)}
        tab[k] = {"level_ratio_K3_over_K0": lv, "d_K0": d[K0], "d_K3": d[K3], "d_ratio": (d[K3] / d[K0]) if (d[K0] and d[K3] is not None) else None}
main = ["contraction_res_ro", "contraction_res_roe"]
missing = any(v is None for k in tab for v in tab[k]["level_ratio_K3_over_K0"] + [tab[k]["d_ratio"]])
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif missing:
    res = "判定不能 (証拠不足: 窓の点が足りない・ゼロ分母)"
else:
    A = all(v >= 2 for k in main for v in tab[k]["level_ratio_K3_over_K0"])
    B = all(0.9 <= v <= 1.1 for k in main for v in tab[k]["level_ratio_K3_over_K0"] + [tab[k]["d_ratio"]]) and \
        all(v <= 1.5 for k in tab if not k.startswith("contraction") for v in tab[k]["level_ratio_K3_over_K0"])
    res = ("A: Γ = I でも前処理のカーネルは縮流部の残差を 2 倍以上にする (カーネルの違いを優先の候補に)" if A else
           "B: Γ = I では前処理のカーネルと標準のカーネルが同程度 (カーネルの違い説を退け、β < 1 の時間項を優先)" if B else
           "判定不能 (中間・成分で割れる・悪化の移動)")
series = {run: [{"step": r["step"], "contraction_ro": r.get("contraction_res_ro"), "throatwall_ro": r.get("throatwall_res_ro")} for r in R[run]] for run in (K0, K3)}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.88", "K0": K0, "K3": K3, "qualification": qual, "effective": eff, "dt_ratio_island_step1000": dtr,
       "step0_rel_diff_max": max(rel0.values()), "table": tab, "series": series, "result": res}
(OUTD / "pk_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
fmt = lambda v: "-" if v is None else f"{v:.3g}"
for k, t in tab.items():
    print(f"{k:24s} 水準比 K3/K0 {[fmt(v) for v in t['level_ratio_K3_over_K0']]} d_K0 {fmt(t['d_K0'])} d_K3 {fmt(t['d_K3'])}")
for run, s in series.items():
    print(run, "縮流部 ρ:", " ".join(f"{x['step']}:{fmt(x['contraction_ro'])}" for x in s if x["step"] % 2000 == 0))
print("dt 比:", dtr, " 実効値:", eff); print("資格:", qual); print("結果:", res)
