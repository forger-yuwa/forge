"""plan architecture-float-state-double-geometry §6.54 (事前登録、codex diagnose `2026-10-12-ic-result-next-diagnose.md`) の判定:
G0c・B0 (上限 50)・継続場 `_ic/cont_q32.h5` から、同じ FP64 のビルド (`~/forge-fgeom11-fp64`、commit b0303547) で step 末尾の保存量 9 つの量子化を OFF (C) / ON (Q)、
陽性対照は同じソースの float のビルド (`~/forge-fgeom11-f32`) (F)。各 2 本・3,000 step・250 ごと。B で回す。
評価域: x/r_t ∈ [−9, −6] かつ壁から 40〜90 節点 (G0c の 121 節点)。窓: W1 = step 2,000・2,250・2,500、W2 = 2,500・2,750・3,000。
E_run(W) = 評価域・窓の中の (T_run − T̄_C)² の平均の平方根 (T̄_C = C 2 本の同じ時点の平均) [K]。N(W) = 同じ域・窓の (T_C1 − T_C2)² の平均の平方根 [K] (C の再実行の差、全差で保守側)。
前提 (欠ければ判別不能): 陽性対照 両窓で min E_F ≥ 10 N。バイナリ・入力の sha、6 本の設定の一致、完走・有限、C と Q の res_0 の 9 保存量と座標のビット一致、
  Q の自己点検 (step 1・2,750・3,000 の Q・Q_N・Q_M で表せない値 0、step 1 の化学種の履歴は全 CV か入口のピンを除く内点で 0)、C と F に量子化の起動の記録が無いこと、必要なファイルがそろうこと。
結果 A: 両窓で min E_Q − N ≥ 0.5 × max E_F → 「末尾の量子化だけでは不十分」(第 1 仮説) を退ける (この介入の十分性まで。float の主因・同じモードの証明とはしない)。
結果 B: 両窓で max E_Q + N ≤ 0.1 × min E_F → 第 2 仮説 (継続場では末尾の量子化だけで十分) をこの起点・3,000 step に限って退ける。
それ以外は判別不能。→ _band_ab/cold_pair/rqs_judge.json"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
ARMS = {"C": ("run_0737_qs_c1", "run_0740_qs_c2"), "Q": ("run_0738_qs_q1", "run_0741_qs_q2"), "F": ("run_0739_qs_f1", "run_0742_qs_f2")}
ALLR = [r for g in ARMS.values() for r in g]
WIN = {"W1": (2000, 2250, 2500), "W2": (2500, 2750, 3000)}
STEPS = sorted({0} | {s for w in WIN.values() for s in w})
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
QA = {"Q": list(Q9), "Q_N": [q + "N" for q in Q9[:7]], "Q_M": [q + "M" for q in Q9[:7]]}
QB = ["roY0N", "roY1N", "roY0M", "roY1M"]
NJ = 121
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "E": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
log = (HERE / "rqs.log").read_text(errors="replace") if (HERE / "rqs.log").exists() else ""
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 継続場 ([0-9a-f]{16})", log)
OUT["checks"]["バイナリ・入力"] = hd[-1] if hd else None
if not hd or hd[-1][1:] != ("74954e7be5ef38a8", "3baafa3d23d31e5e"): und(f"FP64 のバイナリ・継続場の sha が登録と違う ({hd[-1] if hd else None})")
bl = (HERE / "rqc_build.log").read_text(errors="replace") if (HERE / "rqc_build.log").exists() else ""
m = re.findall(r"make rc=0 .* ([0-9a-f]{64})  ", bl)
OUT["checks"]["float のビルド (rqc_build.log)"] = m[-1][:16] if m else None
if not (m and hd and m[-1][:16] == hd[-1][0]): und("float のビルドの sha が rqc_build.log と同じと確かめられない")
cref = (HERE / "run_0717_c50_d1/solverConfig.yaml").read_text().replace("nStepOuter: 30000}", "nStepOuter: 3000}").replace("outStepInterval: 500", "outStepInterval: 250")
for r in ALLR:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
    if not (d / "solverConfig.yaml").exists() or (d / "solverConfig.yaml").read_text() != cref: und(f"{r}: 設定が登録 (run_0717 から 2 行だけ変更) と同一でない")
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 3000: und(f"{r}: 3,000 step に届かない")
    OUT["records"][f"{r} rms_ro の窓平均 (2,000〜2,500 / 2,500〜3,000)"] = [float(np.mean([float(x["rms_ro"]) for x in rows if a <= int(x["step"]) <= b])) for a, b in ((2000, 2500), (2500, 3000))]
    for s in STEPS:
        p = d / f"res_{s}.h5"
        if not p.exists(): und(f"{r}: res_{s} が無い"); break
        with h5py.File(p, "r") as h:
            if "VALUE/T" not in h or not np.all(np.isfinite(h["VALUE/T"][:])): und(f"{r}: res_{s} の T が無いか非有限"); break
    m = re.findall(rf"{r} 起動 .* 入力 ([0-9a-f]{{16}})", log)
    if not m or m[-1] != "3baafa3d23d31e5e": und(f"{r}: 入力の sha が継続場と違う ({m[-1] if m else None})")
    lines = [l for l in (d / "forge_run.log").read_text(errors="replace").splitlines() if "[quantizeState]" in l] if (d / "forge_run.log").exists() else []
    started = any("診断用: step の末尾で保存量" in l for l in lines)
    if r in ARMS["Q"]:
        if not any("保存量 9 本" in l for l in lines): und(f"{r}: 9 本を丸める起動の記録が無い")
    elif started: und(f"{r}: 量子化の起動の記録がある (C・F では無いはず)")
def q9(run):
    p = HERE / run / "res_0.h5"
    if not p.exists(): return None
    with h5py.File(p, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in Q9 if "VALUE/" + k in h}
def coords(run):
    with h5py.File(HERE / run / "nozzle.h5", "r") as h: return np.asarray(h["MESH/COORD"][:])
if not OUT["undecidable"]:
    r0 = q9(ARMS["C"][0]); x0 = coords(ARMS["C"][0])
    if r0 is None or len(r0) != 9: und("C1 の res_0 の 9 保存量がそろわない")
    else:
        for r in ARMS["C"][1:] + ARMS["Q"]:
            v = q9(r)
            if v is None or len(v) != 9 or any(not np.array_equal(v[k], r0[k]) for k in Q9): und(f"{r}: res_0 の 9 保存量が {ARMS['C'][0]} と違う")
            if not np.array_equal(coords(r), x0): und(f"{r}: 座標が {ARMS['C'][0]} と違う")
        vf = q9(ARMS["F"][0])
        if vf is not None: OUT["records"]["C1 と F1 の res_0 の最大相対差"] = {k: float(np.max(np.abs(vf[k].astype(float) - r0[k].astype(float)) / np.maximum(np.abs(r0[k].astype(float)), 1e-300))) for k in Q9 if k in vf}
def pairs(seg): return dict((k, int(v)) for k, v in re.findall(r"(ro[A-Za-z0-9]*) (\d+)", seg))
for r in ARMS["Q"]:
    p = HERE / r / "forge_run.log"
    lines = [l for l in p.read_text(errors="replace").splitlines() if "[quantizeState]" in l] if p.exists() else []
    OUT["checks"][f"{r} 自己点検"] = [l[:300] for l in lines if re.search(r"step \d+:", l)][:10]
    for st in (1, 2750, 3000):
        la = [l for l in lines if f"step {st}:" in l and "表せない値" in l and "(Q)" in l]
        if not la: und(f"{r}: step {st} の (a) の行が無い"); continue
        segs = dict(zip(("Q", "Q_N", "Q_M"), re.split(r"\((?:Q|Q_N|Q_M)\)", la[0].split("表せない値", 1)[1])[:3]))
        for g, names in QA.items():
            pv = pairs(segs.get(g, ""))
            miss = [n for n in names if n not in pv]; nz = {n: pv[n] for n in names if pv.get(n)}
            if miss or nz: und(f"{r}: step {st} の (a) {g} で 名前が無い {miss}・表せない値 {nz}")
    lb = [l for l in lines if "step 1:" in l and "(roY_N)" in l]
    if not lb: und(f"{r}: step 1 の (b) の行が無い")
    else:
        pv = pairs(lb[0].split("表せない値", 1)[1].split("合計")[0])
        miss = [n for n in QB if n not in pv]
        if miss: und(f"{r}: (b) で名前が無い {miss}")
        elif any(pv[n] for n in QB):
            note = [l for l in lines if "step 1:" in l and "入口のピン" in l and "を除く内点では" in l]
            ok = bool(note) and re.search(r"を除く内点では.*合計 0(?!\d)", note[-1])
            OUT["checks"][f"{r} (b) 0 でない (全 CV)"] = {n: pv[n] for n in QB}
            if not ok: und(f"{r}: (b) で入口のピンを除く内点にも表せない値がある: {note[-1][:200] if note else '注意の行が無い'}")
verdict = None
if not OUT["undecidable"]:
    xy = coords(ARMS["C"][0]).reshape(-1, 3)
    S = float(np.load(XC.OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")["scale"])
    ni = len(xy) // NJ; X = (xy[:, 0] / S).reshape(ni, NJ); JW = np.tile(NJ - 1 - np.arange(NJ), (ni, 1))
    REG = ((X >= -9) & (X <= -6) & (JW >= 40) & (JW <= 90)).ravel()
    OUT["records"]["評価域の節点数"] = int(REG.sum())
    def T(run, s):
        with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: return np.asarray(h["VALUE/T"][:], float)[REG]
    Tc = {s: (T(ARMS["C"][0], s), T(ARMS["C"][1], s)) for s in STEPS}
    ok_pc, ra, rb = [], [], []
    for w, ss in WIN.items():
        N = float(np.sqrt(np.mean([np.mean((Tc[s][0] - Tc[s][1]) ** 2) for s in ss])))
        E = {r: float(np.sqrt(np.mean([np.mean((T(r, s) - 0.5 * (Tc[s][0] + Tc[s][1])) ** 2) for s in ss]))) for r in ALLR}
        mx = {r: float(max(np.max(np.abs(T(r, s) - 0.5 * (Tc[s][0] + Tc[s][1]))) for s in ss)) for r in ALLR}
        eQ = [E[r] for r in ARMS["Q"]]; eF = [E[r] for r in ARMS["F"]]
        OUT["E"][w] = dict(N=N, E=E, max_abs_dT=mx)
        print(f"{w}: N {N:.4f} K | E_C {E[ARMS['C'][0]]:.4f}/{E[ARMS['C'][1]]:.4f} | E_Q {eQ[0]:.4f}/{eQ[1]:.4f} | E_F {eF[0]:.4f}/{eF[1]:.4f} K | max|ΔT| Q {mx[ARMS['Q'][0]]:.2f}/{mx[ARMS['Q'][1]]:.2f} F {mx[ARMS['F'][0]]:.2f}/{mx[ARMS['F'][1]]:.2f} K")
        ok_pc.append(min(eF) >= 10 * N); ra.append(min(eQ) - N >= 0.5 * max(eF)); rb.append(max(eQ) + N <= 0.1 * min(eF))
    OUT["checks"]["陽性対照 (両窓で min E_F ≥ 10 N)"] = ok_pc
    if not all(ok_pc): und("陽性対照が成立しない (min E_F < 10 N の窓がある)")
    elif all(ra): verdict = "結果 A: 継続場では末尾の 9 量の量子化だけで float の差の半分以上を再現した (「量子化だけでは不十分」を退ける。この介入の十分性まで)"
    elif all(rb): verdict = "結果 B: 末尾の 9 量の量子化は float の差の 1 割以下 (「継続場では量子化だけで十分」をこの起点・3,000 step に限って退ける)"
    else: verdict = "判別不能"
    for r in ALLR:
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.54: {verdict}")
(HERE / "_band_ab/cold_pair/rqs_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rqs_judge.json")
