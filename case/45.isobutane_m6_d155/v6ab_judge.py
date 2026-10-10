"""plan architecture-float-state-double-geometry §6.18 (事前登録) の判定。v6ab.sh の後に AWS の case/45 で回す。結果は _band_ab/cold_pair/v6ab_judge.json。
D = V4 の到達時の差 (float − FP64)/FP64 [%]、Δ = 30,000 step 後 − 起点 (FP64 の到達時の値で割った %)。
支持: 4 量すべてで B が FP64 の向きに 0.5|D| 以上動き、|Δ_A| ≤ 0.1|D|、|Δ_B − Δ_A| > 10·|Δ_A2 − Δ_A|。
棄却: 4 量すべてで |Δ_B − Δ_A| ≤ 0.1|D|、かつ A・B がともに同じ向きに 0.5|D| 以上動く。それ以外は判別不能。"""
import json, sys, subprocess, numpy as np, h5py
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
K = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
ARMS = {"A": HERE / "run_0460_ab_f32cont", "A2": HERE / "run_0461_ab_f32cont2", "B": HERE / "run_0462_ab_fp64switch"}
STEPS = list(range(2500, 30001, 2500))
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
v4f = {r["step"]: r for r in json.load(open(HERE / "run_0451_v4_f32/m9_watch.json"))["rows"]}
v4d = {r["step"]: r for r in json.load(open(HERE / "run_0452_v4_fp64/m9_watch.json"))["rows"]}
start = v4f[140000]; ref = v4d[125000]
D = {k: 100 * (start[k] - ref[k]) / abs(ref[k]) for k in K}
yb_x, yb = XC.common_yb()
ser = {}
for a, run in ARMS.items():
    ser[a] = {}
    for s in STEPS:
        r = CS.snapshot(run, s, yb_x, yb)
        if r["nonfinite"] or any(not np.isfinite(r[k]) for k in K): raise SystemExit(f"判定不能: {a} の step {s} に非有限")
        ser[a][s] = r
print("step   " + "  ".join(f"{k}: A / A2 / B [% 起点から]" for k in K))
for s in STEPS:
    print(f"{s:6d} " + "  ".join("%+.3f/%+.3f/%+.3f" % tuple(100 * (ser[a][s][k] - start[k]) / abs(ref[k]) for a in ("A", "A2", "B")) for k in K))
dA = {k: 100 * (ser["A"][30000][k] - start[k]) / abs(ref[k]) for k in K}
dA2 = {k: 100 * (ser["A2"][30000][k] - start[k]) / abs(ref[k]) for k in K}
dB = {k: 100 * (ser["B"][30000][k] - start[k]) / abs(ref[k]) for k in K}
sup = []; rej = []
for k in K:
    toward = -np.sign(D[k]) * dB[k] >= 0.5 * abs(D[k])
    s_ = toward and abs(dA[k]) <= 0.1 * abs(D[k]) and abs(dB[k] - dA[k]) > 10 * abs(dA2[k] - dA[k])
    r_ = abs(dB[k] - dA[k]) <= 0.1 * abs(D[k]) and np.sign(dA[k]) == np.sign(dB[k]) and min(abs(dA[k]), abs(dB[k])) >= 0.5 * abs(D[k])
    sup.append(s_); rej.append(r_)
    print(f"  {k:11s} D {D[k]:+.3f} %  Δ_A {dA[k]:+.3f}  Δ_A2 {dA2[k]:+.3f}  Δ_B {dB[k]:+.3f}  → 支持の条件 {s_}・棄却の条件 {r_}")
verdict = "支持 (精度に依存して差が保たれている)" if all(sup) else "棄却 (この期間で精度の変更は支配的でない)" if all(rej) else "判別不能"
print(f"== VERDICT §6.18: {verdict}")
# 記録: 収縮部の領域の A と B の差
RT, NJ = 0.076807, 121
def fld(run, s):
    with h5py.File(run / f"res_{s}.h5", "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in ("Ux", "Uy", "T", "vis_turb")}
with h5py.File(HERE / "run_0460_ab_f32cont/nozzle.h5", "r") as h: X = np.asarray(h["MESH/COORD"][:], np.float64).reshape(-1, 3)
n = X.shape[0]; jj = np.arange(n) % NJ; xr = X[:, 0] / RT
m = (xr >= -5) & (xr < -1) & (jj >= 20) & (jj <= 60)
print("== 記録: 収縮部 (x/r_t ∈ [−5, −1)、j 20〜60) の A と B の差の平均 (|U| は |B| ≥ 1 m/s で正規化)")
rows = []
for s in (2500, 10000, 20000, 30000):
    a, b = fld(ARMS["A"], s), fld(ARMS["B"], s)
    ua, ub = np.hypot(a["Ux"], a["Uy"]), np.hypot(b["Ux"], b["Uy"])
    du = np.mean(np.abs(ua[m] - ub[m]) / np.maximum(ub[m], 1.0)); dt = np.mean(np.abs(a["T"][m] - b["T"][m]) / b["T"][m])
    dmu = np.mean(np.abs(a["vis_turb"][m] - b["vis_turb"][m]) / np.maximum(b["vis_turb"][m], 1e-12))
    print(f"  step {s}: |U| {du:.3e}  T {dt:.3e}  μt {dmu:.3e}"); rows.append((s, du, dt, dmu))
print("== 記録: 判定ツール (この 30,000 step の区間)")
tools = {}
for a, run in ARMS.items():
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    tools[a] = ov[-1] if ov else f"(行なし) rc={p.returncode}"; print(f"  {a} check_convergence: {tools[a]}")
out = dict(D=D, dA=dA, dA2=dA2, dB=dB, support=sup, reject=rej, verdict=verdict, contraction=rows, tools=tools,
           series={a: {s: {k: ser[a][s][k] for k in K} for s in STEPS} for a in ARMS})
(HERE / "_band_ab/cold_pair/v6ab_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
