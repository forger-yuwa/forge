"""plan architecture-float-state-double-geometry §6.55 (事前登録、codex diagnose `2026-10-12-quant-result-next-diagnose.md`) の判定:
G0c・B0 (上限 50)・継続場 `_ic/cont_q32.h5` から、`space.pRef` 0 / 5,499,000 を両精度で比べる。B で回す。
腕: C = FP64、Cp = FP64 + pRef、F = float、Fp = float + pRef (バイナリは §6.54 と同じ fgeom11)、各 2 本・3,000 step・250 ごと。
評価域・窓は §6.54 と同じ (x/r_t ∈ [−9, −6]・壁から 40〜90 節点、W1 = 2,000・2,250・2,500、W2 = 2,500・2,750・3,000)。
E0_i = (T_Fi − T̄_C) の評価域・窓の RMS、Ep_i = (T_Fpi − T̄_Cp) の RMS、N = max(RMS(C1 − C2), RMS(Cp1 − Cp2))、G = max_ij RMS(Cp_i − C_j) [K]。
前提 (欠ければ判別不能): 両窓で min E0 ≥ 10 N かつ G ≤ 0.05 × min E0。閉性 (`rpc_gate.json`) g ≤ 1e-12。バイナリ・入力の sha、設定 (C・F は run_0717 から 2 行、
  Cp・Fp はさらに pRef だけ)、完走、必須の残差の列がそろい全行が有限、保存した全場の全量が有限、同じ精度の中で pRef の有無によらず res_0 の 9 保存量と座標がビット一致。
結果 A: 両窓で max Ep + N + G ≤ 0.5 × min E0 → 第 1 仮説 (pRef だけで登録窓の温度差を少なくとも半減) を支持し、「pRef 単独では半減できない」を退ける。
結果 B: 両窓で min Ep − N − G > 0.5 × max E0 → 第 1 仮説を退ける (圧力に関わる誤差全般の除外の根拠にはしない)。
それ以外は判別不能。記録: 最大 |ΔT|、rms_ro の窓平均、check_convergence、評価域の平均温度の check_quasisteady。→ _band_ab/cold_pair/rpc_judge.json"""
import csv, glob, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
ARMS = {"C": ("run_0743_pc_c1", "run_0747_pc_c2"), "Cp": ("run_0744_pc_cp1", "run_0748_pc_cp2"),
        "F": ("run_0745_pc_f1", "run_0749_pc_f2"), "Fp": ("run_0746_pc_fp1", "run_0750_pc_fp2")}
ALLR = [r for g in ARMS.values() for r in g]
WIN = {"W1": (2000, 2250, 2500), "W2": (2500, 2750, 3000)}
SAVED = list(range(0, 3001, 250))
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
NEED = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
PREF = "5499000.0"; NJ = 121
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "E": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
log = (HERE / "rpc.log").read_text(errors="replace") if (HERE / "rpc.log").exists() else ""
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 継続場 ([0-9a-f]{16}) pRef (\S+)", log)
OUT["checks"]["バイナリ・入力・pRef"] = hd[-1] if hd else None
if not hd or hd[-1] != ("1f9ff32232da8823", "74954e7be5ef38a8", "3baafa3d23d31e5e", PREF): und(f"バイナリ・入力・pRef が登録と違う ({hd[-1] if hd else None})")
gp = HERE / "_band_ab/cold_pair/rpc_gate.json"
gate = json.loads(gp.read_text()) if gp.exists() else None
OUT["checks"]["閉性 (rpc_gate.json)"] = gate.get("closure") if gate else None
if not gate or gate.get("fail") or gate.get("closure", {}).get("g_max", 1.0) > 1e-12: und("閉性の関門が通っていない、または記録が無い")
base = (HERE / "run_0717_c50_d1/solverConfig.yaml").read_text().replace("nStepOuter: 30000}", "nStepOuter: 3000}").replace("outStepInterval: 500", "outStepInterval: 250")
cfg_p = re.sub(r"^space: \{", f"space: {{pRef: {PREF}, ", base, count=1, flags=re.M)
if cfg_p == base: und("基準の設定に 'space: {' が無く pRef の期待値を作れない")
for k, runs in ARMS.items():
    for r in runs:
        d = HERE / r
        rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
        if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
        want = cfg_p if k.endswith("p") else base
        if not (d / "solverConfig.yaml").exists() or (d / "solverConfig.yaml").read_text() != want: und(f"{r}: 設定が登録と同一でない")
        rows = list(csv.DictReader(open(d / "residual_history.csv")))
        miss = [c for c in NEED if c not in rows[0]]
        if miss: und(f"{r}: 残差の必須列が無い {miss}")
        if any(not math.isfinite(float(x[c])) for x in rows for c in rows[0] if c.startswith("rms_")): und(f"{r}: 残差に非有限")
        oe = [x for x in rows if x.get("phase", "outer_end") == "outer_end"]
        if int(oe[-1]["step"]) + 1 < 3000: und(f"{r}: 3,000 step に届かない")
        OUT["records"][f"{r} rms_ro の窓平均 (2,000〜2,500 / 2,500〜3,000)"] = [float(np.mean([float(x["rms_ro"]) for x in oe if a <= int(x["step"]) <= b])) for a, b in ((2000, 2500), (2500, 3000))]
        for s in SAVED:
            p = d / f"res_{s}.h5"
            if not p.exists(): und(f"{r}: res_{s} が無い"); break
            with h5py.File(p, "r") as h:
                if "VALUE/T" not in h: und(f"{r}: res_{s} に T が無い"); break
                bad = [q for q in h["VALUE"] if h["VALUE"][q].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][q][:]))]
                if bad: und(f"{r}: res_{s} に非有限 {bad[:3]}"); break
        m = re.findall(rf"{r} 起動 .* 入力 ([0-9a-f]{{16}})", log)
        if not m or m[-1] != "3baafa3d23d31e5e": und(f"{r}: 入力の sha が継続場と違う")
def q9(run):
    with h5py.File(HERE / run / "res_0.h5", "r") as h: return {q: np.asarray(h["VALUE/" + q][:]) for q in Q9 if "VALUE/" + q in h}
def coords(run):
    with h5py.File(HERE / run / "nozzle.h5", "r") as h: return np.asarray(h["MESH/COORD"][:])
if not OUT["undecidable"]:
    for prec, grp in (("FP64", ARMS["C"] + ARMS["Cp"]), ("float", ARMS["F"] + ARMS["Fp"])):
        r0, x0 = q9(grp[0]), coords(grp[0])
        if len(r0) != 9: und(f"{prec}: res_0 の 9 保存量がそろわない")
        for r in grp[1:]:
            v = q9(r)
            if len(v) != 9 or any(not np.array_equal(v[q], r0[q]) for q in Q9): und(f"{r}: res_0 の 9 保存量が {grp[0]} と違う")
            if not np.array_equal(coords(r), x0): und(f"{r}: 座標が {grp[0]} と違う")
verdict = None
if not OUT["undecidable"]:
    xy = coords(ARMS["C"][0]).reshape(-1, 3)
    S = float(np.load(XC.OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")["scale"])
    ni = len(xy) // NJ; X = (xy[:, 0] / S).reshape(ni, NJ); JW = np.tile(NJ - 1 - np.arange(NJ), (ni, 1))
    REG = ((X >= -9) & (X <= -6) & (JW >= 40) & (JW <= 90)).ravel()
    cache = {}
    def T(run, s):
        if (run, s) not in cache:
            with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: cache[(run, s)] = np.asarray(h["VALUE/T"][:], float)[REG]
        return cache[(run, s)]
    def rms(f, ss): return float(np.sqrt(np.mean([np.mean(f(s) ** 2) for s in ss])))
    pre, ra, rb = [], [], []
    for w, ss in WIN.items():
        mC = lambda s: 0.5 * (T(ARMS["C"][0], s) + T(ARMS["C"][1], s))  # noqa: E731
        mCp = lambda s: 0.5 * (T(ARMS["Cp"][0], s) + T(ARMS["Cp"][1], s))  # noqa: E731
        E0 = [rms(lambda s, r=r: T(r, s) - mC(s), ss) for r in ARMS["F"]]
        Ep = [rms(lambda s, r=r: T(r, s) - mCp(s), ss) for r in ARMS["Fp"]]
        N = max(rms(lambda s: T(ARMS["C"][0], s) - T(ARMS["C"][1], s), ss), rms(lambda s: T(ARMS["Cp"][0], s) - T(ARMS["Cp"][1], s), ss))
        G = max(rms(lambda s, a=a, b=b: T(a, s) - T(b, s), ss) for a in ARMS["Cp"] for b in ARMS["C"])
        mx = {r: float(max(np.max(np.abs(T(r, s) - (mCp(s) if r in ARMS["Fp"] else mC(s)))) for s in ss)) for r in ARMS["F"] + ARMS["Fp"]}
        OUT["E"][w] = dict(E0=E0, Ep=Ep, N=N, G=G, max_abs_dT=mx)
        print(f"{w}: E0 {E0[0]:.3f}/{E0[1]:.3f} K | Ep {Ep[0]:.3f}/{Ep[1]:.3f} K | N {N:.4f} K | G {G:.4f} K | max|ΔT| F {mx[ARMS['F'][0]]:.1f}/{mx[ARMS['F'][1]]:.1f} Fp {mx[ARMS['Fp'][0]]:.1f}/{mx[ARMS['Fp'][1]]:.1f} K")
        pre.append(min(E0) >= 10 * N and G <= 0.05 * min(E0))
        ra.append(max(Ep) + N + G <= 0.5 * min(E0)); rb.append(min(Ep) - N - G > 0.5 * max(E0))
    OUT["checks"]["前提 (両窓で min E0 ≥ 10 N かつ G ≤ 0.05 min E0)"] = pre
    if not all(pre): und("前提: 陽性対照 (min E0 ≥ 10 N) か FP64 の pRef の影響 (G ≤ 0.05 min E0) が窓によって成り立たない")
    elif all(ra): verdict = "結果 A: pRef だけで登録窓の float の温度差が半分以下になった (第 1 仮説の有限期間の予測を支持)"
    elif all(rb): verdict = "結果 B: pRef を入れても float の温度差は半分を超えて残った (第 1 仮説を退ける。圧力に関わる誤差全般の除外とはしない)"
    else: verdict = "判別不能"
    # 記録: check_convergence と評価域の平均温度の準定常
    for r in ALLR:
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
        sc = HERE / r / "rpc_series.csv"
        with open(sc, "w", newline="") as fh:
            wcsv = csv.writer(fh); wcsv.writerow(["step", "T_region_mean"])
            for s in SAVED[1:]: wcsv.writerow([s, float(np.mean(T(r, s)))])
        qs = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r), "--series-csv", str(sc), "--series-cols", "T_region_mean", "--tail", "0.5"], capture_output=True, text=True)
        OUT["records"][f"{r} check_quasisteady (評価域の平均温度、後半)"] = [l.strip() for l in qs.stdout.splitlines() if "csv" in l][:1]
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.55: {verdict}")
(HERE / "_band_ab/cold_pair/rpc_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rpc_judge.json")
