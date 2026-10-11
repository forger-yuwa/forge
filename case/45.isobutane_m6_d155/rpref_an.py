"""plan architecture-float-state-double-geometry §6.44 (事前登録、codex diagnose 2026-10-11) の判定: pRef (対流流束の圧力の基準値) の短い試験。A で回す。
腕: Cp = run_0700〜0702 (FP64 + space.pRef 5.499e6)、Fp = run_0703〜0705 (float + 同じ pRef)。既存: C = run_0595〜0597 (FP64)、F = run_0621〜0623 (float)。
前段の関門 (cp_gate.py、FP64 の同等性) を通ったことが前提。ほかの前提: Cp・Fp の設定の差は space.pRef だけ (run_0595・run_0621 との差)、
バイナリ・入力が §6.40 と同じ、res_0 の 9 保存量が同じ精度の既存の腕と一致、全残差・保存した場が有限。
指標 a は §6.40 の定義のまま。e_Fp = mean(a_Fp) − mean(a_Cp)、d = [mean(a_F) − mean(a_C)] − e_Fp (4 腕の対比)。
Welch (d は SE² = Σ v_g、ν = (Σ v_g)² / Σ v_g²/(n_g − 1))、2 比較で各 t の分位 1 − 0.05/4。
結果 A: CI(e_Fp) ⊂ [−0.5, +0.5] かつ L_d ≥ 0.5。結果 B: L_Fp > 0.5、U_Fp < −0.5、または U_d < 0.5。それ以外は判別不能。→ _band_ab/cold_pair/rpref_judge.json"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
from scipy import stats
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
C = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3"); F = ("run_0621_rq_f1", "run_0622_rq_f2", "run_0623_rq_f3")
CP = ("run_0700_rp_cp1", "run_0701_rp_cp2", "run_0702_rp_cp3"); FP = ("run_0703_rp_fp1", "run_0704_rp_fp2", "run_0705_rp_fp3")
REF = {"A": "run_0487_tr_fp64", "B": "run_0488_tr_f32", "Bp": "run_0489_tr_f32b"}
S = (1000, 1500, 2000); Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "pre": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
g = HERE / "_band_ab/cold_pair/cp_gate.json"
gate = json.loads(g.read_text()) if g.exists() else {}
OUT["pre"]["前段の関門"] = gate.get("verdict")
if gate.get("verdict") != "関門を通った": und(f"前段の関門を通っていない: {gate.get('verdict')}")
def heads(log):
    t = (HERE / log).read_text(errors="replace"); return (re.findall(r"F64 ([0-9a-f]{16})", t) or [None])[-1], (re.findall(r"F32 ([0-9a-f]{16})", t) or [None])[-1]
OUT["pre"]["バイナリ (rq.log / rpref.log)"] = [heads("rq.log"), heads("rpref.log")]
if heads("rq.log") != heads("rpref.log"): und("バイナリが §6.40 と同じと確かめられない")
for r in CP + FP:
    d = HERE / r
    if (d / "RUN_RC").read_text().strip() != "0": und(f"{r}: RUN_RC が 0 でない"); continue
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 2000: und(f"{r}: 2,000 step に届かない")
    for s in (0, 500, 1000, 1500, 2000):
        with h5py.File(d / f"res_{s}.h5", "r") as h:
            if any(h["VALUE"][k].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][k][:])) for k in h["VALUE"]): und(f"{r}: res_{s} に非有限")
    base = "run_0595_rq_c1" if r in CP else "run_0621_rq_f1"
    import difflib
    a_, b_ = (d / "solverConfig.yaml").read_text().splitlines(), (HERE / base / "solverConfig.yaml").read_text().splitlines()
    dl = [l for l in difflib.unified_diff(b_, a_, lineterm="", n=0) if l[:1] in "+-" and not l.startswith(("+++", "---"))]
    ok = len(dl) == 2 and dl[0].startswith("-space: {") and dl[1].startswith("+space: {pRef: 5499000.0, ") and dl[1].replace("pRef: 5499000.0, ", "")[1:] == dl[0][1:]
    if not ok: und(f"{r}: {base} との設定の差が space.pRef だけでない: {dl}")
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(d)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
def q9(run):
    with h5py.File(HERE / run / "res_0.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in Q9 if k in h["VALUE"]}
if not OUT["undecidable"]:
    for runs, base in ((CP, C[0]), (FP, F[0])):
        r0 = q9(base)
        for r in runs:
            bad = [k for k, v in q9(r).items() if not np.array_equal(v, r0[k])]
            if bad: und(f"{r}: res_0 の保存量が {base} と違う {bad}")
verdict = None
if not OUT["undecidable"]:
    with h5py.File(HERE / REF["A"] / "nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
    ni, nj, Sc = XC.mesh_info(xy0.size // 3); X = (xy0.reshape(-1, 3)[:, 0] / Sc).reshape(ni, nj)
    I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
    DOM = ((X < -2) & (I >= 1) & (J >= 1) & (J <= nj - 2)).ravel()
    def T(run, s):
        with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: return np.asarray(h["VALUE/T"][:], float)
    TA0 = T(REF["A"], 0); z = lambda run, s: (T(run, s) - T(run, 0)) / TA0  # noqa: E731
    V = {s: 0.5 * (z(REF["B"], s) + z(REF["Bp"], s)) - z(REF["A"], s) for s in S}
    a = lambda run: float(np.mean([np.dot(V[s][DOM], z(run, s)[DOM]) / np.dot(V[s][DOM], V[s][DOM]) for s in S]))  # noqa: E731
    A_ = {k: np.array([a(r) for r in runs]) for k, runs in (("C", C), ("F", F), ("Cp", CP), ("Fp", FP))}
    OUT["a"] = {k: v.tolist() for k, v in A_.items()}
    v = {k: x.var(ddof=1) / len(x) for k, x in A_.items()}
    tq = 1 - 0.05 / 4
    def ci(e, vs):
        se = math.sqrt(sum(vs))
        if se == 0: return e, e, se, float("inf"), 0.0
        nu = sum(vs) ** 2 / sum(x ** 2 / 2 for x in vs if x); t = float(stats.t.ppf(tq, nu)); return e - t * se, e + t * se, se, nu, t
    eFp = A_["Fp"].mean() - A_["Cp"].mean(); LFp, UFp, se1, nu1, t1 = ci(eFp, [v["Fp"], v["Cp"]])
    eF = A_["F"].mean() - A_["C"].mean(); d = eF - eFp; Ld, Ud, se2, nu2, t2 = ci(d, [v["F"], v["C"], v["Fp"], v["Cp"]])
    OUT["compare"] = dict(e_Fp=dict(e=eFp, L=LFp, U=UFp, se=se1, df=nu1, t=t1), d=dict(e=d, L=Ld, U=Ud, se=se2, df=nu2, t=t2), e_F_ref=eF)
    print(f"a: C {A_['C'].round(4).tolist()} F {A_['F'].round(4).tolist()} Cp {A_['Cp'].round(4).tolist()} Fp {A_['Fp'].round(4).tolist()}")
    print(f"e_Fp = {eFp:+.4f} 区間 [{LFp:+.4f}, {UFp:+.4f}] (df {nu1:.1f})、d = e_F − e_Fp = {d:+.4f} 区間 [{Ld:+.4f}, {Ud:+.4f}] (df {nu2:.1f})、参考 e_F = {eF:+.4f}")
    if -0.5 <= LFp and UFp <= 0.5 and Ld >= 0.5: verdict = "結果 A: この pRef の変更で、登録した成分を参照の尺度で 0.5 以上減らし、残りも ±0.5 以内 (圧力ゲージの変更の十分性。残差の加算が真因とは言わない)"
    elif LFp > 0.5 or UFp < -0.5 or Ud < 0.5: verdict = "結果 B: この pRef の変更だけでは抑制の目標を満たさない (圧力の丸めが無関係とは言わない)"
    else: verdict = "判別不能"
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.44: {verdict}")
(HERE / "_band_ab/cold_pair/rpref_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rpref_judge.json")
