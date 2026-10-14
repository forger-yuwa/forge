"""plan architecture-float-state-double-geometry §6.52 (事前登録、codex diagnose 2026-10-11 `2026-10-11-g0c-precision-test-diagnose.md`) の判定:
縮流部を粗くした格子 G0c (ni 4467) で、B0 のまま float と FP64 を比べる。B で回す。
D = run_0725_g0c_d1・run_0726_g0c_d2 (FP64)、F = run_0727_g0c_f1・run_0728_g0c_f2 (float)。
S_c = run_0725 の step 0 の |θ_r| (断面ごと)。Δ = 100 × (窓の平均 − 自分の step 0) / S_c、W1 = 20,000〜25,000・W2 = 25,000〜30,000。
D_c = D_old × S_old / S_c (D_old = (+1.892, +1.647, +1.431) %、§6.25 から持ち込んだ尺度。S_old = 旧 A run_0487 の step 0 の |θ_r|)。
e_ij = ΔF_i − ΔD_j、N = max(|ΔD1 − ΔD2|, |ΔF1 − ΔF2|) (2 反復で観測した幅。信頼区間ではない)。
前提: (1) §6.51 の上限 50 の陽性対照 (run_0717・0721・0719・0723) で両窓・3 断面とも min sign(D)(ΔF − ΔD) − N ≥ 0.5|D| (rcap_judge.json)。
  (2) G0c の準備の照合 (g0c_prep.log の終了、g0c_check.json の VERDICT が照合 OK、入力の sha)。
  (3) 抽出の検査: run_0725 と run_0717 (G0、同じ FP64 のバイナリ・同じ設定) の step 0 の θ_r(40/70/94) の相対差 ≤ 1e-9 (区間の外の断面はビット一致の入力から作る)。
  (4) バイナリと入力の sha256、各精度の中の res_0 の 9 保存量のビット一致、全期間の残差と使う出力の有限性、4 本の設定が run_0717 と同一。
結果 A (大差が残る): 両窓・3 断面で min|e| − N ≥ 0.5|D_c| (符号は別に記録)。
結果 B (差が小さい): 両窓・3 断面で max|e| + N ≤ 0.1|D_c|。結論は「登録した窓で θ_r の増分の差が小さい」まで。
  窓の平均の値そのものの差 V_ij = 100 × (F_i の窓の平均 − D_j の窓の平均) / S_c が、N_v (精度の中の 2 本の窓の平均の差) と合わせて max|V| + N_v ≤ 0.1|D_c| を満たすときだけ「値の差が小さい」と書く。
それ以外は判別不能。Q_w は副 (記録だけ)。→ _band_ab/cold_pair/rg0c_judge.json"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
ARMS = {"D": ("run_0725_g0c_d1", "run_0726_g0c_d2"), "F": ("run_0727_g0c_f1", "run_0728_g0c_f2")}
ALLR = [r for g in ARMS.values() for r in g]
G0REF = "run_0717_c50_d1"   # G0 の FP64・上限 50 (§6.51) — 抽出の検査と設定の照合の基準
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
STEPS = [0] + list(range(20000, 30001, 500))
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
DOLD = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
NJ = 121
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "delta": {}, "value": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
S_RT = float(np.load(XC.OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")["scale"])
yb_x, yb = XC.common_yb()
def snap(run, st):
    """cold_series.snapshot と同じ量。ni は節点数 / 121 (XC.mesh_info は G0 の ni に固定なので使わない)。"""
    xy, ro, ux, uy, T, tke = XC.load_forge(run, st)
    if len(ro) % NJ: raise SystemExit(f"{run.name}: 節点数 {len(ro)} が {NJ} で割り切れない")
    o = XC.reduce_fields(xy, ro, ux, uy, T, tke, False, len(ro) // NJ, NJ, S_RT, yb_x, yb)
    rec = {"nonfinite": int(sum(np.count_nonzero(~np.isfinite(a)) for a in (ro, ux, uy, T, tke))), "Q_w": float(o["Q_w"])}
    for xs in XC.XS:
        i = int(np.argmin(np.abs(o["x"] - xs))); rec[f"theta_r_{int(xs)}"] = float(o["theta_r"][i])
    return rec
# 前提 (1): 陽性対照。§6.51 の上限 50 の 4 本から直接計算する (rcap_judge.json は §6.51 の別の前提で止まると Δ を残さないため)。式は §6.51 の前提 1 と同じ
PC = {"D": ("run_0717_c50_d1", "run_0721_c50_d2"), "F": ("run_0719_c50_f1", "run_0723_c50_f2")}
def pos_control():
    import cold_series as CS
    A0 = {int(k): v for k, v in json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())["series"]["A"].items()}[0]
    sr = {r: {s: CS.snapshot(HERE / r, s, yb_x, yb) for s in STEPS} for g in PC.values() for r in g}
    pc = {}
    for w, (w0, w1) in WIN.items():
        for q in KEYS:
            dl = lambda r: 100 * (np.mean([sr[r][s][q] for s in STEPS if w0 <= s <= w1]) - sr[r][0][q]) / abs(A0[q])
            dD = [dl(r) for r in PC["D"]]; dF = [dl(r) for r in PC["F"]]; N = max(abs(dD[0] - dD[1]), abs(dF[0] - dF[1]))
            pc[f"{w}|{q}"] = float(min(np.sign(DOLD[q]) * (f - d) for f in dF for d in dD) - N)
    return pc, all(pc[f"{w}|{q}"] >= 0.5 * abs(DOLD[q]) for w in WIN for q in KEYS)
if "--pc-only" in sys.argv:
    pc, ok = pos_control()
    for k, v in pc.items(): print(f"陽性対照 {k}: min sgn(D)(ΔF−ΔD)−N = {v:+.3f} (0.5|D| = {0.5 * abs(DOLD[k.split('|')[1]]):.3f})")
    print("陽性対照:", "成立" if ok else "不成立"); sys.exit(0 if ok else 1)
# 前提 (4): sha
log = (HERE / "rg0c.log").read_text(errors="replace")
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 共通の入力 ([0-9a-f]{16})", log)
OUT["checks"]["バイナリ・入力"] = hd[-1] if hd else None
if not hd or hd[-1] != ("793215917d156680", "7e95e5bbd8f50c9e", "6f8fcda6d97482e3"): und("バイナリ・G0c の入力の sha256 が登録と同じと確かめられない")
# 前提 (2): G0c の準備
pl = (HERE / "g0c_prep.log").read_text(errors="replace")
gc = json.loads((HERE / "_band_ab/cold_pair/g0c_check.json").read_text())
OUT["checks"]["G0c の準備"] = dict(prep_end="== 終了" in pl and "中止" not in pl, check=gc.get("verdict"), seg_same_j=gc.get("seg_same_j_fraction"), seg_dist_max_rt=gc["seg_nearest_dist_rt"]["max"])
if not ("== 終了" in pl and "中止" not in pl and "共通の入力 _g0c/q32_init.h5 sha 6f8fcda6d97482e3" in pl): und("G0c の準備の記録が登録と違う")
if not str(gc.get("verdict", "")).startswith("照合 OK"): und("G0c と実際の G0 の入力の照合が OK でない")
# 前提 (1): 陽性対照
try:
    pc, pc_ok = pos_control()
    OUT["checks"]["陽性対照 (§6.51 上限 50) min sgn(D)(ΔF−ΔD)−N"] = pc
    if not pc_ok: und("§6.51 の上限 50 の陽性対照で大差が再現しない組がある")
except (FileNotFoundError, KeyError, OSError) as e:
    und(f"§6.51 の陽性対照を計算できない ({e})")
# 前提 (4): 完走・有限性・設定・res_0
def q9(path):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"] if k in Q9}
cref = (HERE / G0REF / "solverConfig.yaml").read_text()
for r in ALLR:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
    if (d / "solverConfig.yaml").read_text() != cref: und(f"{r}: 設定が {G0REF} と同一でない")
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 30000: und(f"{r}: 30,000 step に届かない")
    for s in STEPS:
        p = d / f"res_{s}.h5"
        if not p.exists(): und(f"{r}: res_{s} が無い"); break
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][q].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][q][:])) for q in h["VALUE"]): und(f"{r}: res_{s} に非有限"); break
for k, (r1, r2) in ARMS.items():
    if (HERE / r1 / "res_0.h5").exists() and (HERE / r2 / "res_0.h5").exists():
        a, b = q9(HERE / r1 / "res_0.h5"), q9(HERE / r2 / "res_0.h5")
        if set(a) != set(Q9) or set(b) != set(Q9) or any(not np.array_equal(a[q], b[q]) for q in Q9): und(f"{k}: res_0 の 9 保存量が精度の中で一致しない")
if (HERE / ARMS["D"][0] / "res_0.h5").exists() and (HERE / ARMS["F"][0] / "res_0.h5").exists():
    a, b = q9(HERE / ARMS["D"][0] / "res_0.h5"), q9(HERE / ARMS["F"][0] / "res_0.h5")
    OUT["records"]["精度の間の res_0 の最大相対差 (記録だけ)"] = {q: float(np.max(np.abs(a[q] - b[q]) / np.maximum(np.abs(a[q]), 1e-300))) for q in Q9 if q in a and q in b}
# 前提 (3): 抽出の検査
verdict = None
if not OUT["undecidable"]:
    g0 = snap(HERE / G0REF, 0); c0 = snap(HERE / ARMS["D"][0], 0)
    rel = {q: abs(c0[q] - g0[q]) / abs(g0[q]) for q in KEYS}
    OUT["checks"]["抽出の検査: step 0 の θ_r の G0c 対 G0 の相対差"] = rel
    OUT["records"]["step 0 の Q_w: G0c / G0 (区間を含むので違ってよい)"] = [c0[QW], g0[QW]]
    print("抽出の検査 (step 0 の θ_r、G0c 対 G0 の相対差):", {q: f"{v:.2e}" for q, v in rel.items()})
    if not all(v <= 1e-9 for v in rel.values()): und("抽出の検査: 区間の外の断面の θ_r が G0 と合わない")
if not OUT["undecidable"]:
    ser = {}
    for r in ALLR:
        ser[r] = {}
        for s in STEPS:
            x = snap(HERE / r, s)
            if x["nonfinite"] or any(not np.isfinite(x[q]) for q in KEYS + (QW,)): und(f"{r}: step {s} の系列に非有限"); break
            ser[r][s] = {q: x[q] for q in KEYS + (QW,)}
if not OUT["undecidable"]:
    tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
    A0 = {int(k): v for k, v in tj["series"]["A"].items()}[0]
    Sc = {q: abs(ser[ARMS["D"][0]][0][q]) for q in KEYS + (QW,)}
    Dc = {q: DOLD[q] * abs(A0[q]) / Sc[q] for q in KEYS + (QW,)}
    OUT["records"]["S_c / S_old / D_c"] = {q: [Sc[q], abs(A0[q]), Dc[q]] for q in KEYS + (QW,)}
    def wm(r, q, w): w0, w1 = WIN[w]; return float(np.mean([ser[r][s][q] for s in ser[r] if w0 <= s <= w1]))
    def dl(r, q, w): return 100 * (wm(r, q, w) - ser[r][0][q]) / Sc[q]
    ra, rb, rv = [], [], []
    print("窓  断面         ΔD1      ΔD2      ΔF1      ΔF2      N      min|e|−N  max|e|+N  符号つき平均  0.5|D_c| 0.1|D_c|  | 値: max|V|+N_v")
    for w in WIN:
        for q in KEYS + (QW,):
            dD = [dl(r, q, w) for r in ARMS["D"]]; dF = [dl(r, q, w) for r in ARMS["F"]]
            N = max(abs(dD[0] - dD[1]), abs(dF[0] - dF[1])); e = [f - d for f in dF for d in dD]
            lo_abs = min(abs(x) for x in e) - N; hi = max(abs(x) for x in e) + N
            mD = [wm(r, q, w) for r in ARMS["D"]]; mF = [wm(r, q, w) for r in ARMS["F"]]
            V = [100 * (f - d) / Sc[q] for f in mF for d in mD]; Nv = 100 * max(abs(mD[0] - mD[1]), abs(mF[0] - mF[1])) / Sc[q]
            hv = max(abs(x) for x in V) + Nv; Dk = abs(Dc[q])
            OUT["delta"][f"{w}|{q}"] = dict(D=dD, F=dF, N=N, e=e, lo_abs=lo_abs, hi=hi, sign_mean=float(np.mean(dF) - np.mean(dD)), Dc=Dc[q])
            OUT["value"][f"{w}|{q}"] = dict(V=V, Nv=Nv, hi=hv)
            print(f"{w}  {q:11s} {dD[0]:+7.3f}  {dD[1]:+7.3f}  {dF[0]:+7.3f}  {dF[1]:+7.3f}  {N:.3f}  {lo_abs:+.3f}    {hi:.3f}     {np.mean(dF) - np.mean(dD):+.3f}       {0.5 * Dk:.3f}    {0.1 * Dk:.3f}   | {hv:.3f}{' (副)' if q == QW else ''}")
            if q == QW: continue
            ra.append(lo_abs >= 0.5 * Dk); rb.append(hi <= 0.1 * Dk); rv.append(hv <= 0.1 * Dk)
    for r in ALLR:
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
        sc = HERE / r / "rg0c_series.csv"
        with open(sc, "w", newline="") as fh:
            wcsv = csv.writer(fh); wcsv.writerow(["step"] + list(KEYS) + [QW])
            for s in sorted(t for t in ser[r] if t > 0): wcsv.writerow([s] + [ser[r][s][q] for q in KEYS + (QW,)])
        qs = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + (QW,)),
                             "--drift", "0.0001", "--osc", "0.0001", "--tail", "1.0"], capture_output=True, text=True)
        OUT["records"][f"{r} check_quasisteady"] = [l.strip() for l in qs.stdout.splitlines() if "csv" in l][:1]
    if all(ra): verdict = "結果 A: G0c でも B0 の float と FP64 の θ_r の増分の差が大きく残る (第 1 仮説の有限期間の予測を退ける)"
    elif all(rb): verdict = ("結果 B: G0c では登録した窓で θ_r の増分の差が小さい" + ("。窓の平均の値の差も小さい" if all(rv) else "。窓の平均の値の差は小差の条件を満たさない (値の差が小さいとは書かない)")
                             + " (遅いモードの同定・double の不要とはしない)")
    else: verdict = "判別不能"
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.52: {verdict}")
(HERE / "_band_ab/cold_pair/rg0c_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rg0c_judge.json")
