"""plan architecture-float-state-double-geometry §6.26 (事前登録) の判定: 共通の Q32 の起点からの float と FP64 の軌跡の A/B。tr.sh の後に B で回す。
A = run_0487_tr_fp64、B = run_0488_tr_f32、B' = run_0489_tr_f32b。各出力 (0〜30,000 step、500 ごと) から cold_series.snapshot で θ_r(40/70/94)・Q_w を作る。
主判定 (θ_r の 3 断面、窓 W1 = 20,000〜25,000・W2 = 25,000〜30,000):
  D = (+1.892, +1.647, +1.431) % (§6.25)。Δ = (窓の平均 − 自分の起点) / A の起点 × 100。
  精度依存を支持: 両窓・3 断面とも ΔB と ΔB' が ΔA から D の向きへ 0.5|D| 以上離れ、|mean(ΔB,ΔB') − ΔA| > 10·|ΔB − ΔB'|、かつ |ΔA| ≤ 0.1|D|。
  共通の過渡を支持: 両窓・3 断面とも |mean(ΔB,ΔB') − ΔA| ≤ 0.1|D|、かつ ΔA と mean(ΔB,ΔB') が同じ向きに 0.5|D| 以上。
  それ以外・非有限・欠損は判別不能。Q_w は副判定として同じ式で記録する (主判定に混ぜない)。
場の記録 (判定に使わない): 量ごと・領域ごとに、起点 (A の step 0) の領域の RMS で正規化した差の L2。立ち上がり = ‖B−A‖ > 10‖B−B'‖ が 3 出力続いた最初の step。
結果は標準出力と _band_ab/cold_pair/tr_judge.json。"""
import csv, json, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ARMS = {"A": HERE / "run_0487_tr_fp64", "B": HERE / "run_0488_tr_f32", "Bp": HERE / "run_0489_tr_f32b"}
STEPS = list(range(0, 30001, 500))
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
NJ, RT = 121, 0.076807
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "series": {}, "delta": {}, "fields": {}}
yb_x, yb = XC.common_yb()
for a, run in ARMS.items():
    ser = {}
    for s in STEPS:
        if not (run / f"res_{s}.h5").exists():
            if s == 0: continue
            OUT["undecidable"].append(f"{a}: res_{s} が無い"); break
        r = CS.snapshot(run, s, yb_x, yb)
        if r.get("nonfinite") or any(not np.isfinite(r[k]) for k in KEYS + (QW,)): OUT["undecidable"].append(f"{a}: step {s} に非有限"); break
        ser[s] = {k: float(r[k]) for k in KEYS + (QW,)}
    OUT["series"][a] = ser
def start(a): return OUT["series"][a][min(OUT["series"][a])]
verdict = None
if not OUT["undecidable"]:
    A0 = start("A")
    OUT["start_values"] = {a: start(a) for a in ARMS}
    for wn, (w0, w1) in WIN.items():
        for k in KEYS + (QW,):
            dl = {}
            for a in ARMS:
                v = [OUT["series"][a][s][k] for s in OUT["series"][a] if w0 <= s <= w1]
                if len(v) != 11: OUT["undecidable"].append(f"{a} {wn} {k}: {len(v)} 点"); continue
                dl[a] = 100 * (np.mean(v) - start(a)[k]) / abs(A0[k])
            OUT["delta"][f"{wn}|{k}"] = dl
    if not OUT["undecidable"]:
        sup, com = [], []
        print("窓  断面         ΔA [%]    ΔB [%]   ΔB' [%]   |mB−A|   10|B−B'|   0.5|D|   0.1|D|")
        for wn in WIN:
            for k in KEYS + (QW,):
                dl = OUT["delta"][f"{wn}|{k}"]; dA, dB, dBp = dl["A"], dl["B"], dl["Bp"]; mB = 0.5 * (dB + dBp); sg = np.sign(D[k]); Dk = abs(D[k])
                s_ = (sg * (dB - dA) >= 0.5 * Dk) and (sg * (dBp - dA) >= 0.5 * Dk) and abs(mB - dA) > 10 * abs(dB - dBp) and abs(dA) <= 0.1 * Dk
                c_ = abs(mB - dA) <= 0.1 * Dk and np.sign(dA) == np.sign(mB) and min(abs(dA), abs(mB)) >= 0.5 * Dk
                print(f"{wn}  {k:11s} {dA:+8.3f}  {dB:+8.3f}  {dBp:+8.3f}  {abs(mB - dA):7.3f}  {10 * abs(dB - dBp):8.3f}  {0.5 * Dk:6.3f}  {0.1 * Dk:6.3f}  {'支持' if s_ else ''}{'共通' if c_ else ''}{' (副)' if k == QW else ''}")
                if k != QW: sup.append(s_); com.append(c_)
        verdict = ("精度依存を支持 (float に依存する作用素・更新の差が偏りを作る。SST の commit の証明ではない)" if all(sup)
                   else "共通の過渡を支持 (この期間では精度の差は支配的でない)" if all(com) else "判別不能")
if verdict is None: verdict = "判別不能: " + "; ".join(OUT["undecidable"][:5])
OUT["verdict"] = verdict
print(f"== VERDICT §6.26: {verdict}")
# ---- 場の記録 ----
QF = ("ro", "roK", "roOmega", "k", "omega", "vis_turb", "P", "T", "U", "omg_prod", "omg_dest", "omg_cross", "omg_trans")
with h5py.File(ARMS["A"] / "nozzle.h5", "r") as h: X = np.asarray(h["MESH/COORD"][:], np.float64).reshape(-1, 3)
n = X.shape[0]; jj = np.arange(n) % NJ; xr = X[:, 0] / RT
REG = {"収縮部全体": (xr >= -15) & (xr < -1), "収縮部の内部": (xr >= -5) & (xr < -1) & (jj >= 20) & (jj <= 60),
       "軸の近く": jj <= 8, "壁際": jj >= 117, "下流": xr >= 1}
def fields(run, s):
    with h5py.File(run / f"res_{s}.h5", "r") as h:
        V = h["VALUE"]; f = {}
        for q in QF:
            if q == "U": f[q] = np.hypot(np.asarray(V["Ux"][:], np.float64), np.asarray(V["Uy"][:], np.float64))
            elif q in V: f[q] = np.asarray(V[q][:], np.float64)
        return f
try:
    s0 = min(OUT["series"]["A"]) if OUT["series"].get("A") else 0
    F0 = fields(ARMS["A"], s0)
    scale = {(q, rn): float(np.sqrt(np.mean(F0[q][m] ** 2))) for q in F0 for rn, m in REG.items()}
    hist = {}
    for s in [x for x in STEPS if x > 0]:
        Fa, Fb, Fp = fields(ARMS["A"], s), fields(ARMS["B"], s), fields(ARMS["Bp"], s)
        for q in F0:
            if q not in Fb or q not in Fp: continue
            for rn, m in REG.items():
                sc = scale[(q, rn)]
                if not sc > 0: continue
                ba = np.sqrt(np.mean((Fb[q][m] - Fa[q][m]) ** 2)) / sc; bp = np.sqrt(np.mean((Fb[q][m] - Fp[q][m]) ** 2)) / sc
                hist.setdefault(f"{q}|{rn}", []).append((s, float(ba), float(bp)))
    onset = {}
    for key, rows in hist.items():
        run_ = 0
        for s, ba, bp in rows:
            run_ = run_ + 1 if ba > 10 * bp else 0
            if run_ >= 3: onset[key] = s - 1000; break
    OUT["fields"] = {"hist": hist, "onset_step": onset}
    print("== 記録: 立ち上がり (‖B−A‖ > 10‖B−B'‖ が 3 出力続いた最初の step、早い順)")
    for key, s in sorted(onset.items(), key=lambda kv: kv[1])[:25]: print(f"  {s:6d}  {key}  (step 30000 の ‖B−A‖ {hist[key][-1][1]:.2e}, ‖B−B'‖ {hist[key][-1][2]:.2e})")
    print("  立ち上がりの無い組:", sorted(set(hist) - set(onset))[:20])
except Exception as e:
    print("場の記録に失敗:", e); OUT["fields_error"] = str(e)
print("== 記録: 判定ツール")
for a, run in ARMS.items():
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT.setdefault("tools", {})[f"{a} check_convergence"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
    sc = run / "tr_series_20k_30k.csv"
    with open(sc, "w") as fh:
        fh.write("step," + ",".join(KEYS + (QW,)) + "\n")
        for s in sorted(OUT["series"].get(a, {})):
            if 20000 <= s <= 30000: fh.write(f"{s}," + ",".join(f"{OUT['series'][a][s][k]:.10g}" for k in KEYS + (QW,)) + "\n")
    p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + (QW,)), "--tail", "1.0"], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("OVERALL")]
    OUT["tools"][f"{a} check_quasisteady (20k〜30k)"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
    cl = run / "commit_loss.csv"
    if cl.exists():
        rows = list(csv.DictReader(open(cl))); last = max(int(x["step"]) for x in rows) if rows else None
        OUT["tools"][f"{a} commit_loss (最後の step {last}, interior)"] = {x["qty"]: x["S_lostfrac"] for x in rows if int(x["step"]) == last and x["region"] == "interior"}
for k, v in OUT.get("tools", {}).items(): print(f"  {k}: {v}")
o = HERE / "_band_ab/cold_pair/tr_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
