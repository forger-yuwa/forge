"""plan tooling-nozzle-core-grid §4.90 P-g の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 gamma-term)。
凍結の状態 Q0 (run_0693 の res_20000) で、GI (Γ = I) と GB (局所 β の Γ) の nStepInner 5・400・800 の 1 回目の更新 δ = Q1 − Q0 と、
その応答 R(Q0 + δ) (各 run の res_2 の残差。出力の残差は更新前、保存量は更新後)、線形性の確認 R(Q0 + δ★/2) (半分の run の res_1)。
圧力残差 R_p/c² (r = ∂p/∂Q、r と c² は Q0 で固定)、単位体積あたりの体積重み付き RMS。
L(δ) = 1 − ‖R_p(Q0 + δ)‖/‖R_p(Q0)‖、射影 ⟨R_p(Q0 + δ), R_p(Q0)⟩/‖R_p(Q0)‖²。u₅ = max(ρ, ρE) RMS_V(δ₅ − δ★)/RMS_V(δ★)。δ★ = 800 sweep。
領域: 縮流部 (主)・島・壁際の内部 (x −2〜0・壁から 1〜10)・壁の節点 (x −2〜0・壁から 0)。
資格: (1) restart・設定・スイッチの行・半分の書き込み・SST 凍結 (ρk・ρω が Q0 から不変)。(2) 6 本の step 0 の残差が相対 1e−6 以内。
  (3) 縮流部の ρ・ρE の RMS(δ₈₀₀ − δ₄₀₀)/RMS(δ₈₀₀) ≤ 1e−8 (各腕)。(4) 縮流部の線形性の誤差 ≤ 0.05 (各腕)。(5) 非有限なし・RUN_RC 0。
結果 A: 両腕で u₅ ≤ 0.01 かつ 縮流部で L_GB(δ★) ≤ ½ L_GI(δ★)。
結果 B: GB で u₅ ≥ 0.3 かつ |L_GB(δ★) − L_GI(δ★)| ≤ 0.1 |L_GI(δ★)| かつ (L_GI(δ★) − L_GB(δ★)) ≤ 0.2 (L_GI(δ₅) − L_GB(δ₅))。
それ以外・信号不足 (|L_GI(δ★)| < 0.001)・資格不足は判定不能。
usage: python3 cg_pg_judge.py → _band_ab/core_grid/pg_judge.json"""
import csv, json
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
BASE = "run_0693_pdt_off"; SEED = HERE / BASE / "res_20000.h5"
RUNS = {"GI": {5: "run_0696_pg_gi5", 400: "run_0697_pg_gi400", 800: "run_0698_pg_gi800", "h": "run_0702_pg_gih"},
        "GB": {5: "run_0699_pg_gb5", 400: "run_0700_pg_gb400", 800: "run_0701_pg_gb800", "h": "run_0703_pg_gbh"}}
CONS = ("ro", "roUx", "roUy", "roUz", "roe"); RES = ("res_ro", "res_roUx", "res_roUy", "res_roe")
info = json.loads((HERE / BASE / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def rd(path, keys):
    with h5py.File(path, "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
with h5py.File(SEED, "r") as h:
    x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
Q0 = rd(SEED, CONS + ("roK", "roOmega", "P", "sonic", "h0", "Ux", "Uy", "volume"))
nj = x.size // ni; jw = nj - 1 - np.arange(x.size) % nj; V = Q0["volume"]
REG = {"contraction": (x >= -8) & (x <= -2), "island": (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100),
       "wall_interior": (x >= -2) & (x <= 0) & (jw >= 1) & (jw <= 10), "wall_node": (x >= -2) & (x <= 0) & (jw == 0)}
c2 = Q0["sonic"] ** 2; gam = c2 * Q0["ro"] / Q0["P"]; kap = gam - 1; ek = 0.5 * (Q0["Ux"] ** 2 + Q0["Uy"] ** 2); chi = c2 - kap * (Q0["h0"] - ek)
def Rp(R):   # 圧力残差 / c²、単位体積あたり
    return ((chi + kap * ek) * R["res_ro"] - kap * Q0["Ux"] * R["res_roUx"] - kap * Q0["Uy"] * R["res_roUy"] + kap * R["res_roe"]) / c2 / V
def nrm(z, m): return float(np.sqrt(np.sum(V[m] * z[m] ** 2) / np.sum(V[m])))
def ip(a, b, m): return float(np.sum(V[m] * a[m] * b[m]) / np.sum(V[m]))
D, R1, R2, NF = {}, {}, {}, 0
for arm, rr in RUNS.items():
    for k, run in rr.items():
        if k == "h":
            R1[(arm, k)] = rd(HERE / run / "res_1.h5", RES); continue
        q1 = rd(HERE / run / "res_1.h5", CONS + ("roK", "roOmega") + RES)
        D[(arm, k)] = {c: q1[c] - Q0[c] for c in CONS}; D[(arm, k)]["frozen_max"] = max(float(np.max(np.abs(q1[c] - Q0[c]))) for c in ("roK", "roOmega"))
        R1[(arm, k)] = {c: q1[c] for c in RES}; R2[(arm, k)] = rd(HERE / run / "res_2.h5", RES)
        NF += sum(int(np.count_nonzero(~np.isfinite(v))) for v in list(q1.values()) + list(R2[(arm, k)].values()))
# 資格
logs = {run: "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists()) for rr in RUNS.values() for run in rr.values()}
q1 = all("VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / run / "restart.log").read_text() for rr in RUNS.values() for run in rr.values()) \
    and all(l.count("FORGE_DIAG_PRECOND_NO_DT_SCALE=1:") >= 1 for l in logs.values()) \
    and all("Q0 + δ★/2 を書いた" in (HERE / rr["h"] / "half.log").read_text() for rr in RUNS.values()) \
    and all(D[k]["frozen_max"] == 0.0 for k in D)
def step0(run):
    rows = list(csv.reader(open(HERE / run / "residual_history.csv"))); h = rows[0]; r0 = next(r for r in rows[1:] if r[0] == "0" and r[1] == "-1")
    return {h[i]: float(r0[i]) for i in range(len(h)) if h[i].startswith("rms_")}
s0 = [step0(RUNS[a][k]) for a in RUNS for k in (5, 400, 800)]
rel0 = max((abs(s[c] - s0[0][c]) / max(abs(s0[0][c]), 1e-300)) for s in s0 for c in s0[0] if s0[0][c] != 0)
q2 = rel0 <= 1e-6
cauchy = {arm: {rg: {c: nrm(D[(arm, 800)][c] - D[(arm, 400)][c], m) / max(nrm(D[(arm, 800)][c], m), 1e-300) for c in ("ro", "roe")} for rg, m in REG.items()} for arm in RUNS}
q3 = all(cauchy[a]["contraction"][c] <= 1e-8 for a in RUNS for c in ("ro", "roe"))
b0 = {arm: Rp(R1[(arm, 5)]) for arm in RUNS}   # R_p(Q0) (各腕の 5 sweep の run の res_1)
lin = {}
for arm in RUNS:
    full = Rp(R2[(arm, 800)]) - b0[arm]; half = Rp(R1[(arm, "h")]) - b0[arm]
    lin[arm] = {rg: nrm(half - 0.5 * full, m) / max(nrm(0.5 * full, m), 1e-300) for rg, m in REG.items()}
q4 = all(lin[a]["contraction"] <= 0.05 for a in RUNS)
q5 = NF == 0 and all((HERE / run / "RUN_RC").read_text().strip() == "0" for rr in RUNS.values() for run in rr.values())
qual = {"setup_switch_half_freeze": bool(q1), "same_rhs_step0": bool(q2), "inner_converged_800": bool(q3), "linear_response": bool(q4), "finite_completed": bool(q5)}
# 量
tab = {}
for arm in RUNS:
    for k in (5, 400, 800):
        a = Rp(R2[(arm, k)]); b = b0[arm]
        tab[f"{arm}{k}"] = {rg: {"L": 1 - nrm(a, m) / max(nrm(b, m), 1e-300), "proj": ip(a, b, m) / max(ip(b, b, m), 1e-300)} for rg, m in REG.items()}
u5 = {arm: {rg: max(nrm(D[(arm, 5)][c] - D[(arm, 800)][c], m) / max(nrm(D[(arm, 800)][c], m), 1e-300) for c in ("ro", "roe")) for rg, m in REG.items()} for arm in RUNS}
LI5, LI8 = tab["GI5"]["contraction"]["L"], tab["GI800"]["contraction"]["L"]; LB5, LB8 = tab["GB5"]["contraction"]["L"], tab["GB800"]["contraction"]["L"]
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif abs(LI8) < 0.001:
    res = "判定不能 (信号不足: |L_GI(δ★)| < 0.001)"
elif u5["GI"]["contraction"] <= 0.01 and u5["GB"]["contraction"] <= 0.01 and LB8 <= 0.5 * LI8:
    res = "A: 両腕とも 5 sweep でほぼ解けているのに、Γ (β < 1) の十分解いた補正は圧力残差の低減が Γ = I の半分以下 (第 1 仮説を支持)"
elif u5["GB"]["contraction"] >= 0.3 and abs(LB8 - LI8) <= 0.1 * abs(LI8) and (LI8 - LB8) <= 0.2 * (LI5 - LB5):
    res = "B: Γ (β < 1) は 5 sweep の解き残しが大きく、十分解くと Γ = I と同程度の低減に戻る (第 2 仮説を支持)"
else:
    res = "判定不能 (中間・領域で割れる)"
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.90", "qualification": qual, "step0_rel_max": rel0, "cauchy_800_vs_400": cauchy, "linearity": lin,
       "u5": u5, "L_and_projection": tab, "result": res}
(OUTD / "pg_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for k, t in tab.items():
    print(f"{k:6s} " + "  ".join(f"{rg} L {v['L']:+.4f} 射影 {v['proj']:.4f}" for rg, v in t.items()))
print("u5:", {a: {rg: round(v, 4) for rg, v in d.items()} for a, d in u5.items()})
print("δ★ の収束 (縮流部):", {a: {c: f"{cauchy[a]['contraction'][c]:.2e}" for c in ("ro", "roe")} for a in RUNS})
print("線形性 (縮流部):", {a: round(lin[a]["contraction"], 4) for a in RUNS}, " step0 の差:", rel0)
print("資格:", qual); print("結果:", res)
