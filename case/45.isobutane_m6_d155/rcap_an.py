"""plan architecture-float-state-double-geometry §6.51 (事前登録、codex diagnose 2026-10-11) の判定: 方向別 dt の上限だけを 50 → 1 にした float と FP64 の比較 (案 C1)。B で回す。
上限 50: D = run_0717・0721 (FP64)、F = run_0719・0723 (float)。上限 1: D = run_0718・0722、F = run_0720・0724。
Δ = 100 × (窓の平均 − 自分の step 0) / |旧 A (run_0487) の step 0|、W1 = 20,000〜25,000・W2 = 25,000〜30,000、D = (+1.892, +1.647, +1.431) %。N は上限ごと。
前提: (1) 上限 50 で両窓・3 断面とも min sign(D)(ΔF − ΔD) − N ≥ 0.5|D|。(2) step 20,000 の dt_local の上限 1 / 上限 50 (FP64 の 1 本目どうし) の比の中央値が、
  縮流部 (x/r_t ∈ [−8, −2)) の壁側 (j ≥ 100) と全域で < 0.5。(3) 入力・バイナリの sha256、各精度の中の res_0 の 9 保存量のビット一致、有限性、設定の差は上限だけ。
結果 A: 上限 1 で両窓・3 断面とも max|ΔF − ΔD| + N ≤ 0.1|D|。結果 B: 上限 1 でも両窓・3 断面とも min|ΔF − ΔD| − N ≥ 0.5|D|。それ以外は判別不能。→ _band_ab/cold_pair/rcap_judge.json"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ARMS = {50: {"D": ("run_0717_c50_d1", "run_0721_c50_d2"), "F": ("run_0719_c50_f1", "run_0723_c50_f2")},
        1: {"D": ("run_0718_c1_d1", "run_0722_c1_d2"), "F": ("run_0720_c1_f1", "run_0724_c1_f2")}}
ALLR = [r for c in ARMS.values() for g in c.values() for r in g]
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
STEPS = [0] + list(range(20000, 30001, 500))
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "delta": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
log = (HERE / "rcap.log").read_text(errors="replace")
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 共通の入力 ([0-9a-f]{16})", log)
OUT["checks"]["バイナリ・入力"] = hd[-1] if hd else None
if not hd or hd[-1] != ("793215917d156680", "7e95e5bbd8f50c9e", "ba859e02f3ec236a"): und("バイナリ・入力の sha256 が §6.26 と同じと確かめられない")
def q9(path):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"] if k in Q9}
for r in ALLR:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 30000: und(f"{r}: 30,000 step に届かない")
    for s in STEPS:
        p = d / f"res_{s}.h5"
        if not p.exists(): und(f"{r}: res_{s} が無い"); break
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][q].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][q][:])) for q in h["VALUE"]): und(f"{r}: res_{s} に非有限"); break
for cap, g in ARMS.items():
    for k, (r1, r2) in g.items():
        a, b = q9(HERE / r1 / "res_0.h5"), q9(HERE / r2 / "res_0.h5")
        if set(a) != set(Q9) or set(b) != set(Q9) or any(not np.array_equal(a[q], b[q]) for q in Q9): und(f"上限 {cap} {k}: res_0 の 9 保存量が精度の中で一致しない")
    for k in ("D", "F"):
        c50, c1 = (HERE / ARMS[50][k][0] / "solverConfig.yaml").read_text().splitlines(), (HERE / ARMS[1][k][0] / "solverConfig.yaml").read_text().splitlines()
        df = [(x, y) for x, y in zip(c50, c1) if x != y]
        ok = len(c50) == len(c1) and len(df) == 1 and df[0][0].replace("lineDtDirectionalCap: 50.0", "lineDtDirectionalCap: 1.0") == df[0][1]
        if not ok: und(f"{k}: 上限 50 と 1 の設定の差が上限だけでない: {df[:2]}")
# 前提 2: dt の縮小
with h5py.File(HERE / "_tr/q32_init.h5", "r") as h: xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
ni, nj, S = XC.mesh_info(len(xy)); X = (xy[:, 0] / S).reshape(ni, nj); J = np.tile(np.arange(nj), (ni, 1))
def dtl(r):
    with h5py.File(HERE / r / "res_20000.h5", "r") as h: return np.asarray(h["VALUE/dt_local"][:], float).reshape(ni, nj) if "VALUE/dt_local" in h else None
if not OUT["undecidable"]:
    t50, t1 = dtl(ARMS[50]["D"][0]), dtl(ARMS[1]["D"][0])
    if t50 is None or t1 is None: und("dt_local が出力に無い")
    else:
        ratio = t1 / np.where(t50 > 0, t50, np.nan)
        cw = (X >= -8) & (X < -2) & (J >= 100) & (J <= nj - 2)
        m_cw, m_all = float(np.nanmedian(ratio[cw])), float(np.nanmedian(ratio))
        OUT["checks"]["dt_local の比 (上限 1 / 上限 50、FP64、step 20,000) 中央値: 縮流部の壁側 / 全域"] = [m_cw, m_all]
        OUT["records"]["dt_local の比の 10 % 点・90 % 点 (縮流部の壁側)"] = [float(np.nanpercentile(ratio[cw], 10)), float(np.nanpercentile(ratio[cw], 90))]
        print(f"dt_local の比 (上限 1 / 50) 中央値: 縮流部の壁側 {m_cw:.3f}、全域 {m_all:.3f}")
        if not (m_cw < 0.5 and m_all < 0.5): und(f"実効の dt が縮んでいない (中央値 {m_cw:.3f}・{m_all:.3f})")
tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
A0 = {int(k): v for k, v in tj["series"]["A"].items()}[0]
yb_x, yb = XC.common_yb()
verdict = None
if not OUT["undecidable"]:
    ser = {}
    for r in ALLR:
        ser[r] = {}
        for s in STEPS:
            x = CS.snapshot(HERE / r, s, yb_x, yb)
            if x.get("nonfinite") or any(not np.isfinite(x[q]) for q in KEYS + (QW,)): und(f"{r}: step {s} の系列に非有限"); break
            ser[r][s] = {q: float(x[q]) for q in KEYS + (QW,)}
if not OUT["undecidable"]:
    def dl(r, q, w):
        w0, w1 = WIN[w]; return 100 * (np.mean([ser[r][s][q] for s in ser[r] if w0 <= s <= w1]) - ser[r][0][q]) / abs(A0[q])
    pos, ra, rb = [], [], []
    print("上限  窓  断面         ΔD1      ΔD2      ΔF1      ΔF2      N       min sgn(F−D)−N  max|F−D|+N  min|F−D|−N  0.5|D|  0.1|D|")
    for cap in (50, 1):
        for w in WIN:
            for q in KEYS + (QW,):
                dD = [dl(r, q, w) for r in ARMS[cap]["D"]]; dF = [dl(r, q, w) for r in ARMS[cap]["F"]]
                sg, Dk = np.sign(D[q]), abs(D[q]); N = max(abs(dD[0] - dD[1]), abs(dF[0] - dF[1]))
                prs = [(f, d) for f in dF for d in dD]
                lo = min(sg * (f - d) for f, d in prs) - N; hi = max(abs(f - d) for f, d in prs) + N; lo_abs = min(abs(f - d) for f, d in prs) - N
                OUT["delta"][f"{cap}|{w}|{q}"] = dict(D=dD, F=dF, N=N, lo=lo, hi=hi, lo_abs=lo_abs, sign_mean=float(np.mean(dF) - np.mean(dD)))
                print(f"{cap:4d}  {w}  {q:11s} {dD[0]:+7.3f}  {dD[1]:+7.3f}  {dF[0]:+7.3f}  {dF[1]:+7.3f}  {N:.3f}  {lo:+.3f}          {hi:.3f}       {lo_abs:+.3f}      {0.5 * Dk:.3f}   {0.1 * Dk:.3f}{' (副)' if q == QW else ''}")
                if q == QW: continue
                if cap == 50: pos.append(lo >= 0.5 * Dk)
                else: ra.append(hi <= 0.1 * Dk); rb.append(lo_abs >= 0.5 * Dk)
    for r in ALLR:
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
        sc = HERE / r / "rcap_series.csv"
        with open(sc, "w", newline="") as fh:
            wcsv = csv.writer(fh); wcsv.writerow(["step"] + list(KEYS) + [QW])
            for s in sorted(t for t in ser[r] if t > 0): wcsv.writerow([s] + [ser[r][s][q] for q in KEYS + (QW,)])
        qs = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + (QW,)),
                             "--drift", "0.0001", "--osc", "0.0001", "--tail", "1.0"], capture_output=True, text=True)
        OUT["records"][f"{r} check_quasisteady"] = [l.strip() for l in qs.stdout.splitlines() if "csv" in l][:1]
    if not all(pos): verdict = "判別不能 (前提 1: 上限 50 の陽性対照で大差が再現しない組がある)"
    elif all(ra): verdict = "結果 A: 上限 1 でこの期間の float と FP64 の差が小さい (方向別 dt の拡大が関わる。長期の一致・遅いモードの同定とはしない)"
    elif all(rb): verdict = "結果 B: 上限 1 でも大差が残る (上限の変更だけでは小さくならない)"
    else: verdict = "判別不能"
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.51: {verdict}")
(HERE / "_band_ab/cold_pair/rcap_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rcap_judge.json")
