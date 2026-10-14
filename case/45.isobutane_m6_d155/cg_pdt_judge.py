"""plan tooling-nozzle-core-grid §4.86 P-dt の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 precond-worsens): lowMachPrecond 3 で、刻みの拡大あり (S+) と なし (S0)。
領域: 縮流部 (x/r_t −8〜−2) と throatwall (x/r_t −2〜0・壁から 0〜10 節点目)。窓 W1 = 10k〜14k・W2 = 16k〜20k (各 5 出力、非重複) の平均 q̄、d = q̄(W2)/q̄(W1)。
傾き: 各窓 5 点の log q の最小二乗。
資格: (1) 2 本の restart_field がビット一致、設定ファイルが同一、S0 のログにスイッチ有効の行があり S+ に無い。
  (2) step 0 の全残差列が両腕で相対 1e−6 以内。
  (3) 島 (x/r_t −6〜−3.25・壁から 60〜100 節点目) の step 1000 で dt(S+)/dt(S0) を S+ の状態からの期待の倍率で割った値の中央値が 0.85〜1.15、
      かつ |c_S+/c_S0 − 1| の中央値 ≤ 0.01。
  (4) S+ が §4.83 の P3 を再現: 10000・15000・20000 の縮流部の ρ の q の S+/P3 が 0.67〜1.5、かつ S+ の縮流部の ρ・ρE の d > 1。
  (5) 全出力で非有限なし・ρ・P・T > 0。
結果 A: 両領域の ρ・ρE で S0 の d ≤ 0.8 かつ両窓の傾きが負、かつ縮流部・throatwall・縮流部の外の全保存量の残差で S0 の d ≤ 1.1。
結果 B: 両領域の ρ・ρE で S0 の d ≥ 1.1 かつ d_S0/d_S+ ≥ 0.8。それ以外・資格不足は判定不能。
usage: python3 cg_pdt_judge.py <S+ の run> <S0 の run> [P3 の run (既定 run_0689_pp_p3)] → _band_ab/core_grid/pdt_judge.json"""
import csv, json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
SP, S0 = sys.argv[1:3]; P3 = sys.argv[3] if len(sys.argv) > 3 else "run_0689_pp_p3"; EPS = 0.15
W = ((10000, 14000), (16000, 20000))
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
R = {SP: jl(f"region_res_{SP}.jsonl"), S0: jl(f"region_res_{S0}.jsonl"), P3: jl(f"region_res_{P3}.jsonl")}
PH = {SP: jl(f"phys_{SP}.jsonl"), S0: jl(f"phys_{S0}.jsonl")}
def pts(run, key, lo, hi): return [(r["step"], r[key]) for r in R[run] if lo <= r["step"] <= hi and key in r]
def wm(run, key, lo, hi):
    v = pts(run, key, lo, hi); return float(np.mean([b for _, b in v])) if len(v) == 5 else None
def slope(run, key, lo, hi):
    v = pts(run, key, lo, hi)
    if len(v) != 5 or any(b <= 0 for _, b in v): return None
    return float(np.polyfit([a for a, _ in v], np.log([b for _, b in v]), 1)[0])
def dval(run, key):
    a, b = wm(run, key, *W[0]), wm(run, key, *W[1]); return (b / a) if (a and b is not None) else None
def logtext(run): return "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists())
# 資格
bit = all("VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / r / "restart.log").read_text() for r in (SP, S0))
same_cfg = (HERE / SP / "solverConfig.yaml").read_text() == (HERE / S0 / "solverConfig.yaml").read_text()
sw = {r: logtext(r).count("FORGE_DIAG_PRECOND_NO_DT_SCALE=1:") for r in (SP, S0)}
q1 = bit and same_cfg and sw[S0] >= 1 and sw[SP] == 0
def step0(run):
    rows = list(csv.reader(open(HERE / run / "residual_history.csv"))); h = rows[0]
    r0 = next(r for r in rows[1:] if r[0] == "0" and r[1] == "-1")
    return {h[i]: float(r0[i]) for i in range(len(h)) if h[i].startswith("rms_")}
a0, b0 = step0(SP), step0(S0)
rel0 = {k: (abs(a0[k] - b0[k]) / max(abs(a0[k]), abs(b0[k])) if max(abs(a0[k]), abs(b0[k])) > 0 else 0.0) for k in a0}
q2 = all(v <= 1e-6 for v in rel0.values())
def load(run, n, keys):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
        x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
    nj = x.size // ni
    return {k: a.reshape(ni, nj)[:, ::-1] for k, a in v.items()}, x.reshape(ni, nj)[:, ::-1], ni, nj
gp, X, ni, nj = load(SP, 1000, ("dt_local", "Ux", "Uy", "sonic")); g0, _, _, _ = load(S0, 1000, ("dt_local", "sonic"))
isl = np.zeros((ni, nj), bool); cols = (X[:, 0] >= -6.0) & (X[:, 0] <= -3.25); isl[np.ix_(cols, np.arange(60, 101))] = True
c = np.maximum(gp["sonic"], 1e-8); vm = np.hypot(gp["Ux"], gp["Uy"]); Ur = np.minimum(c, np.maximum(vm, EPS * c)); beta = (Ur / c) ** 2
cp = 0.5 * np.sqrt((1 - beta) ** 2 * vm ** 2 + 4 * Ur ** 2); expect = np.maximum((vm + c) / np.maximum(0.5 * (1 + beta) * vm + cp, 1e-30), 1.0)
meas = gp["dt_local"] / g0["dt_local"]; mo = float(np.median(meas[isl] / expect[isl])); close = float(np.median(np.abs(gp["sonic"][isl] / g0["sonic"][isl] - 1)))
dtq = {"island_median_measured_ratio": float(np.median(meas[isl])), "island_median_expected_ratio": float(np.median(expect[isl])),
       "island_median_measured_over_expected": mo, "island_median_abs_sonic_rel_diff": close}
q3 = 0.85 <= mo <= 1.15 and close <= 0.01
rep = {}
for n in (10000, 15000, 20000):
    a = [r["contraction_res_ro"] for r in R[SP] if r["step"] == n]; b = [r["contraction_res_ro"] for r in R[P3] if r["step"] == n]
    rep[n] = (a[0] / b[0]) if (a and b and b[0] > 0) else None
dSP = {k: dval(SP, f"contraction_{k}") for k in ("res_ro", "res_roe")}
q4 = all(v is not None and 0.67 <= v <= 1.5 for v in rep.values()) and all(v is not None and v > 1 for v in dSP.values())
q5 = not any(r["nonfinite"] for run in (SP, S0) for r in R[run] + PH[run]) and all(x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for run in (SP, S0) for x in PH[run]) \
    and all((HERE / r / "RUN_RC").exists() and (HERE / r / "RUN_RC").read_text().strip() == "0" for r in (SP, S0))
qual = {"same_start_config_switch": bool(q1), "same_rhs_step0": bool(q2), "dt_scaling_toggled": bool(q3), "Splus_reproduces_P3": bool(q4), "finite_physical_completed": bool(q5)}
# 判定
regs = ("contraction", "throatwall"); comps = ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega")
tab = {f"{rg}_{c}": {"d_S0": dval(S0, f"{rg}_{c}"), "d_Splus": dval(SP, f"{rg}_{c}"), "slope_S0": [slope(S0, f"{rg}_{c}", *w) for w in W]}
       for rg in regs + ("outside",) for c in comps}
main = [f"{rg}_{c}" for rg in regs for c in ("res_ro", "res_roe")]
missing = any(tab[k]["d_S0"] is None or tab[k]["d_Splus"] is None for k in tab)
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif missing:
    res = "判定不能 (証拠不足: 窓の点が足りない)"
else:
    A = all(tab[k]["d_S0"] <= 0.8 and all(s is not None and s < 0 for s in tab[k]["slope_S0"]) for k in main) and all(tab[k]["d_S0"] <= 1.1 for k in tab)
    B = all(tab[k]["d_S0"] >= 1.1 and tab[k]["d_S0"] / tab[k]["d_Splus"] >= 0.8 for k in main)
    res = ("A: 刻みの拡大を止めると近壁と縮流部の残差の増加が消える (刻みの拡大への依存を支持)" if A else
           "B: 刻みの拡大を止めても増加が残る (「拡大を外せば解消する」を退け、β < 1 の境界と Γ を優先の候補に)" if B else
           "判定不能 (中間・成分で割れる・悪化の移動)")
series = {run: [{"step": r["step"], "contraction_ro": r.get("contraction_res_ro"), "throatwall_ro": r.get("throatwall_res_ro")} for r in R[run]] for run in (SP, S0)}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.86", "Splus": SP, "S0": S0, "P3": P3, "qualification": qual, "switch_lines": sw,
       "step0_rel_diff_max": max(rel0.values()), "dt": dtq, "Splus_over_P3": rep, "table": tab, "series": series, "result": res}
(OUTD / "pdt_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
fmt = lambda v: "-" if v is None else f"{v:.3g}"
for k, t in tab.items():
    print(f"{k:24s} d_S0 {fmt(t['d_S0'])} d_S+ {fmt(t['d_Splus'])} 傾き S0 {[fmt(s) for s in t['slope_S0']]}")
for run, s in series.items():
    print(run, "縮流部 ρ:", " ".join(f"{x['step']}:{fmt(x['contraction_ro'])}" for x in s if x["step"] % 2000 == 0))
    print(run, "throatwall ρ:", " ".join(f"{x['step']}:{fmt(x['throatwall_ro'])}" for x in s if x["step"] % 2000 == 0))
print("dt:", dtq); print("S+/P3:", rep); print("資格:", qual); print("結果:", res)
