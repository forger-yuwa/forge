"""plan time_integration-line-implicit-speed §6.23 (事前登録) の集計: hoopAreaFromClosure 0/1 の 2000 step の短い比較。
AWS の case/45 で回す (fg3.sh の後)。結果を標準出力と _band_ab/cold_pair/hc_ab.json に書く。"""
import csv, json, subprocess, sys, h5py, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
RT = 0.076807; NJ = 121; OUTS = (1600, 1700, 1800, 1900, 2000)
A, B, A2 = HERE / "run_0409_hc_A", HERE / "run_0410_hc_B", HERE / "run_0411_hc_A2"
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
R = {}
def coords(run):
    for f in sorted(run.glob("res_*.h5")):
        with h5py.File(f, "r") as h:
            if "MESH/COORD" in h: return np.asarray(h["MESH/COORD"][:], np.float64).reshape(-1, 3)
    raise SystemExit(f"{run}: 座標の入った res が無い")
X = coords(A); n = X.shape[0]; jj = np.arange(n) % NJ; xr = X[:, 0] / RT
OM = (xr >= -5.0) & (xr < 0.0); J24 = OM & (jj >= 2) & (jj <= 4); J2059 = OM & (jj >= 20) & (jj <= 59)
def res_ro(run, s):
    with h5py.File(run / f"res_{s}.h5", "r") as h:
        a = np.asarray(h["VALUE/res_ro"][:], np.float64)
        bad = sum(int(np.count_nonzero(~np.isfinite(np.asarray(h["VALUE/" + k][:])))) for k in h["VALUE"].keys())
    return a, bad
print("== §6.23 の判定 (Ω = x/r_t ∈ [−5, 0)、S = Ω の Σ|res_ro|)")
S = {}; bad_total = 0
for tag, run in (("A", A), ("B", B), ("A2", A2)):
    S[tag] = {}
    for s in OUTS:
        a, bad = res_ro(run, s); bad_total += bad
        S[tag][s] = {"S": float(np.sum(np.abs(a[OM]))), "S_j2_4": float(np.sum(np.abs(a[J24]))), "S_j20_59": float(np.sum(np.abs(a[J2059]))),
                     "max": float(np.max(np.abs(a[OM]))), "argmax": int(np.flatnonzero(OM)[np.argmax(np.abs(a[OM]))])}
rB = [S["B"][s]["S"] / S["A"][s]["S"] for s in OUTS]; rA2 = [S["A2"][s]["S"] / S["A"][s]["S"] for s in OUTS]
for s, b, a2 in zip(OUTS, rB, rA2):
    k = S["A"][s]["argmax"]; kb = S["B"][s]["argmax"]
    print(f"  step {s}: S(B)/S(A) {b:.3f}  S(A2)/S(A) {a2:.3f} | j 2〜4: B/A {S['B'][s]['S_j2_4'] / S['A'][s]['S_j2_4']:.3f} | j 20〜59: B/A {S['B'][s]['S_j20_59'] / S['A'][s]['S_j20_59']:.3f}"
          f" | 最大 A {S['A'][s]['max']:.3e} (x/r_t {xr[k]:.2f}・j {jj[k]})  B {S['B'][s]['max']:.3e} (x/r_t {xr[kb]:.2f}・j {jj[kb]})")
if bad_total: verdict = "判定不能 (非有限あり)"
elif not all(0.8 <= r <= 1.25 for r in rA2): verdict = "判定不能 (再実行の比が [0.8, 1.25] の外)"
elif all(r <= 0.5 for r in rB): verdict = "支持 (closure の補正が床に強く効く)"
elif all(r >= 0.9 for r in rB): verdict = "棄却 (短い時間で床を大きく下げる期待を棄却、延長しない)"
else: verdict = "判別不能"
print(f"  VERDICT §6.23: {verdict}")
R["judge"] = {"ratio_B_over_A": rB, "ratio_A2_over_A": rA2, "verdict": verdict, "nonfinite": bad_total, "S": S}
print("== 記録: step 1500〜2000 の rms_* の中央値")
def tail(run):
    rows = [x for x in csv.DictReader(open(run / "residual_history.csv")) if x["phase"] == "outer_end" and 1500 <= int(x["step"]) <= 2000]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
    return {c: float(np.median([float(x[c]) for x in rows])) for c in cols}
T = {t: tail(r) for t, r in (("A", A), ("B", B), ("A2", A2))}
def ratio(a, b): return f"{a / b:.3f}" if b != 0 else ("—" if a == 0 else "inf")
for c in T["A"]:
    print(f"  {c:14s} A {T['A'][c]:.3e}  B {T['B'][c]:.3e} (B/A {ratio(T['B'][c], T['A'][c])})  A2 {T['A2'][c]:.3e} (A2/A {ratio(T['A2'][c], T['A'][c])})")
R["rms_tail_median"] = T
print("== 記録: step 2000 の θ_r・Q_w (B と A の相対差、フラグ 0.05 %・0.1 %)")
yb_x, yb = XC.common_yb()
snap = {t: CS.snapshot(r, 2000, yb_x, yb) for t, r in (("A", A), ("B", B), ("A2", A2))}
for k in ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w", "deficit"):
    a, b, a2 = snap["A"][k], snap["B"][k], snap["A2"][k]
    lim = 0.1 if k == "Q_w" else 0.05
    flag = "" if k == "deficit" else (" ← フラグ" if abs(100 * (b / a - 1)) > lim else "")
    print(f"  {k:11s} A {a:.6g}  B {b:.6g} (B/A−1 {100 * (b / a - 1):+.4f} %)  A2 {a2:.6g} (A2/A−1 {100 * (a2 / a - 1):+.4f} %){flag}")
R["snapshot_2000"] = snap
print("== 記録: 判定ツール (この 2000 step の区間)")
for t, r in (("A", A), ("B", B), ("A2", A2)):
    for tool in ("check_convergence.py", "check_quasisteady.py"):
        p = subprocess.run([sys.executable, str(TOOLS / tool), str(r)], capture_output=True, text=True)
        v = [l for l in (p.stdout + p.stderr).splitlines() if "VERDICT" in l]
        print(f"  {t} {tool}: {v[-1] if v else '(VERDICT の行なし) rc=' + str(p.returncode)}")
        R.setdefault("tools", {})[f"{t}:{tool}"] = v[-1] if v else None
print("== 記録: 同じ状態での照合 (1 step の run、ΔR = R_B − R_A)")
a1, b1, a21 = HERE / "run_0412_hc1_A", HERE / "run_0413_hc1_B", HERE / "run_0414_hc1_A2"
dump = h5py.File(HERE / "run_0416_fg2_dumphc64/stage2.h5", "r")
Acy, Acx, Ap = (np.asarray(dump[f"/axisym/{k}"][:n], np.float64) for k in ("A_closure_y", "A_closure_x", "A_planar"))
def v1(run, keys):
    with h5py.File(run / "res_1.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in keys}
keys = ("res_ro", "res_roUx", "res_roUy", "res_roe", "P", "Uy", "vis_lam", "vis_turb", "axisym_divU", "volume")
FA, FB, FA2 = v1(a1, keys), v1(b1, keys), v1(a21, keys)
mu = FA["vis_lam"] + FA["vis_turb"]; reff = np.where(Ap > 0, FA["volume"] / np.where(Ap > 0, Ap, 1), 0.0)
uor = np.where(reff > 0, FA["Uy"] / np.where(reff > 0, reff, 1), 0.0); tau = 2 * mu * uor - (2.0 / 3.0) * mu * FA["axisym_divU"]
na = jj > 0
predP = FA["P"] * (Acy - Ap); predV = -tau * (Acy - Ap); predX = FA["P"] * Acx
for name, d, noise, preds in (("res_roUy", FB["res_roUy"] - FA["res_roUy"], FA2["res_roUy"] - FA["res_roUy"], (("圧力", predP), ("圧力+粘性", predP + predV))),
                              ("res_roUx", FB["res_roUx"] - FA["res_roUx"], FA2["res_roUx"] - FA["res_roUx"], (("圧力", predX),))):
    nd = np.linalg.norm(d[na]); print(f"  Δ{name}: ‖ΔR‖ {nd:.3e} (再実行 {np.linalg.norm(noise[na]):.3e})")
    for pn, p in preds:
        print(f"     予測 ({pn}) との差 ‖ΔR − 予測‖/‖ΔR‖ = {np.linalg.norm(d[na] - p[na]) / nd if nd > 0 else float('nan'):.3e}")
for k in ("res_ro", "res_roe"):
    print(f"  Δ{k}: 最大 |B − A| {np.max(np.abs(FB[k] - FA[k])):.3e}、再実行 {np.max(np.abs(FA2[k] - FA[k])):.3e}")
o = HERE / "_band_ab/cold_pair/hc_ab.json"; o.write_text(json.dumps(R, ensure_ascii=False, indent=1, default=float)); print("→", o)
