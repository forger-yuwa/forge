"""plan architecture-float-state-double-geometry §6.53 (事前登録、codex diagnose `2026-10-12-cap-g0c-result-diagnose.md`) の判定:
G0c・B0 (上限 50) で、初期場だけを変える A/B。旧起点 = `_g0c/q32_init.h5`、継続場 = `run_0725_g0c_d1/res_30000.h5` の 9 保存量を float32 に丸めた `_ic/cont_q32.h5`。
各起点で FP64 ×2・float ×2、各 6,000 step・250 ごと。B で回す。
窓 V1 = 4,000〜5,000、V2 = 5,000〜6,000 (residual_history.csv の outer_end の rms_ro の窓平均)。
前提 (欠ければ判別不能):
  (a) 旧起点で再上昇が再現: 旧起点の float 2 本が両窓 ≥ 1.5e-5、旧起点の FP64 2 本が両窓 ≤ 3e-6。
  (b) 継続場は起動直後の残差が小さい: 継続場の各 run の step 1〜250 の rms_ro の平均 ≤ 0.1 × 旧起点の 4 本の同じ平均の最小。
  (c) 継続場で精度を切り替えた直後に大差が無い: 継続場の各 float の step 1〜250 の平均 ≤ 3 × 継続場の FP64 2 本の同じ平均の最大。
  (d) バイナリ・入力の sha、8 本の設定が run_0717_c50_d1 から nStepOuter・outStepInterval だけを変えたものとバイト一致、各起点・各精度の中の res_0 の 9 保存量のビット一致、完走、残差の有限性。
結果 A: 継続場の float 2 本が両窓 ≥ 1.5e-5 かつ継続場の FP64 2 本が両窓 ≤ 3e-6 → 移植直後の大きな過渡を必要条件とする読み (H3 の限定版) を退ける。
結果 B: 継続場の 4 本がすべて両窓 ≤ 3e-6 → 初期場によらず 6,000 step で戻るという第 1 仮説を退け、初期場への依存を優先する。
それ以外は判別不能。rms_ro は既知の再上昇の検出に使うだけで、収束の判定ではない。θ_r の小差や double の要否はこの試験で判定しない。
記録: 全残差列の窓平均、各起点の F1 対 D1 の温度差の最大と位置 (250 step ごと)、精度の間の res_0 の差、check_convergence。→ _band_ab/cold_pair/ric_judge.json"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
ARMS = {"old": {"D": ("run_0729_ic_old_d1", "run_0733_ic_old_d2"), "F": ("run_0731_ic_old_f1", "run_0735_ic_old_f2")},
        "cont": {"D": ("run_0730_ic_cont_d1", "run_0734_ic_cont_d2"), "F": ("run_0732_ic_cont_f1", "run_0736_ic_cont_f2")}}
ALLR = [r for a in ARMS.values() for g in a.values() for r in g]
WIN = {"V1": (4000, 5000), "V2": (5000, 6000)}
HI, LO = 1.5e-5, 3e-6
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
NJ = 121; SEG = list(range(0, 6001, 250))
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "rms_ro": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
log = (HERE / "ric.log").read_text(errors="replace")
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 旧起点 ([0-9a-f]{16}) 継続場 ([0-9a-f]{16})", log)
OUT["checks"]["バイナリ・入力"] = hd[-1] if hd else None
if not hd or hd[-1][:3] != ("793215917d156680", "7e95e5bbd8f50c9e", "6f8fcda6d97482e3"): und("バイナリ・旧起点の sha が登録と違う")
cref = (HERE / "run_0717_c50_d1/solverConfig.yaml").read_text().replace("nStepOuter: 30000}", "nStepOuter: 6000}").replace("outStepInterval: 500", "outStepInterval: 250")
rows = {}
for r in ALLR:
    d = HERE / r
    if (d / "RUN_RC").read_text().strip() != "0" if (d / "RUN_RC").exists() else True: und(f"{r}: RUN_RC が 0 でない"); continue
    if (d / "solverConfig.yaml").read_text() != cref: und(f"{r}: 設定が登録 (run_0717 から 2 行だけ変更) と同一でない")
    rr = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rr[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rr for c in cols): und(f"{r}: 残差に非有限")
    if int(rr[-1]["step"]) + 1 < 6000: und(f"{r}: 6,000 step に届かない")
    rows[r] = (rr, cols)
# 各起点の入力の sha (ric.sh が run を作るときに記録した値と照合)
for a, g in ARMS.items():
    for k, runs in g.items():
        for r in runs:
            m = re.findall(rf"{r} 起動 .* 入力 ([0-9a-f]{{16}})", log)
            want = hd[-1][2] if (hd and a == "old") else (hd[-1][3] if hd else None)
            if not m or m[-1] != want: und(f"{r}: 入力の sha が起点の登録と違う ({m[-1] if m else None} 対 {want})")
def q9(p):
    with h5py.File(p, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in Q9 if "VALUE/" + k in h}
for a, g in ARMS.items():
    for k, (r1, r2) in g.items():
        if all((HERE / r / "res_0.h5").exists() for r in (r1, r2)):
            x, y = q9(HERE / r1 / "res_0.h5"), q9(HERE / r2 / "res_0.h5")
            if set(x) != set(Q9) or any(not np.array_equal(x[q], y[q]) for q in Q9): und(f"{a} {k}: res_0 の 9 保存量が精度の中で一致しない")
    if all((HERE / g[k][0] / "res_0.h5").exists() for k in ("D", "F")):
        x, y = q9(HERE / g["D"][0] / "res_0.h5"), q9(HERE / g["F"][0] / "res_0.h5")
        OUT["records"][f"{a}: 精度の間の res_0 の最大相対差"] = {q: float(np.max(np.abs(x[q] - y[q]) / np.maximum(np.abs(x[q]), 1e-300))) for q in Q9}
verdict = None
if not OUT["undecidable"]:
    def wmean(r, c, w0, w1): return float(np.mean([float(x[c]) for x in rows[r][0] if w0 <= int(x["step"]) <= w1]))
    for r in ALLR:
        OUT["rms_ro"][r] = {w: wmean(r, "rms_ro", *WIN[w]) for w in WIN} | {"s1_250": wmean(r, "rms_ro", 1, 250)}
        OUT["records"][f"{r} 全残差列の窓平均 (V1, V2)"] = {c: [wmean(r, c, *WIN["V1"]), wmean(r, c, *WIN["V2"])] for c in rows[r][1]}
    R = OUT["rms_ro"]
    print("run                    s1〜250     V1          V2")
    for r in ALLR: print(f"{r:22s} {R[r]['s1_250']:.3e}  {R[r]['V1']:.3e}  {R[r]['V2']:.3e}")
    pa = all(R[r][w] >= HI for r in ARMS["old"]["F"] for w in WIN) and all(R[r][w] <= LO for r in ARMS["old"]["D"] for w in WIN)
    pb = max(R[r]["s1_250"] for g in ARMS["cont"].values() for r in g) <= 0.1 * min(R[r]["s1_250"] for g in ARMS["old"].values() for r in g)
    pc = max(R[r]["s1_250"] for r in ARMS["cont"]["F"]) <= 3 * max(R[r]["s1_250"] for r in ARMS["cont"]["D"])
    OUT["checks"]["前提 (a) 旧起点で再上昇の再現 / (b) 継続場の起動直後が小さい / (c) 継続場の精度の切替直後に大差が無い"] = [pa, pb, pc]
    if not pa: und("前提 (a): 旧起点で再上昇が再現しない")
    if not pb: und("前提 (b): 継続場の起動直後の残差が旧起点の 0.1 倍以下でない")
    if not pc: und("前提 (c): 継続場で精度を切り替えた直後に大差がある")
    ra = all(R[r][w] >= HI for r in ARMS["cont"]["F"] for w in WIN) and all(R[r][w] <= LO for r in ARMS["cont"]["D"] for w in WIN)
    rb = all(R[r][w] <= LO for g in ARMS["cont"].values() for r in g for w in WIN)
    if not OUT["undecidable"]:
        verdict = ("結果 A: 継続場からも float だけが 6,000 step 内に再上昇した (移植直後の大きな過渡を必要条件とする読みを退ける)" if ra else
                   "結果 B: 継続場からは両精度とも低い残差を保った (初期場によらず戻るという第 1 仮説を退け、初期場への依存を優先する)" if rb else "判別不能")
# 記録: F1 対 D1 の温度差の最大と位置 (250 step ごと)
if all(r in rows for r in ALLR):
    with h5py.File(HERE / ARMS["old"]["D"][0] / "nozzle.h5", "r") as h: xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
    S = float(np.load(XC.OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")["scale"])
    ni = len(xy) // NJ; X = (xy[:, 0] / S).reshape(ni, NJ)
    for a, g in ARMS.items():
        tr = []
        for s in SEG:
            pd, pf = HERE / g["D"][0] / f"res_{s}.h5", HERE / g["F"][0] / f"res_{s}.h5"
            if not (pd.exists() and pf.exists()): continue
            with h5py.File(pd, "r") as hd_, h5py.File(pf, "r") as hf_:
                dT = (np.asarray(hf_["VALUE/T"][:], float) - np.asarray(hd_["VALUE/T"][:], float)).reshape(ni, NJ)
            i, j = np.unravel_index(np.argmax(np.abs(dT)), dT.shape)
            tr.append(dict(step=s, max_abs_dT_K=float(abs(dT[i, j])), x_rt=float(X[i, j]), j_from_wall=int(NJ - 1 - j), rms_dT_K=float(np.sqrt(np.mean(dT ** 2)))))
        OUT["records"][f"{a}: F1 − D1 の温度差の推移"] = tr
        if tr: print(f"{a}: F1−D1 の |ΔT| 最大 (step: K @ x/r_t, 壁から j):", " ".join(f"{t['step']}:{t['max_abs_dT_K']:.2f}@{t['x_rt']:.1f},{t['j_from_wall']}" for t in tr[::4]))
    for r in ALLR:
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / r)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4]) if OUT["undecidable"] else "判別不能"
OUT["verdict"] = verdict
print(f"== VERDICT §6.53: {verdict}")
(HERE / "_band_ab/cold_pair/ric_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/ric_judge.json")
