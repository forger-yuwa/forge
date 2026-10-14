"""plan architecture-float-state-double-geometry §6.46 (事前登録、codex diagnose 2026-10-11) の判定: pRef の長い試験 (30,000 step、下流の θ_r・Q_w)。A で回す。
新: F0 = run_0709_lp_f0 (float・pRef なし)、Bp = run_0710_lp_bp1・run_0711_lp_bp2 (float + pRef)、Ap = run_0712_lp_ap (FP64 + pRef)。
旧: A = run_0487 (FP64)・B = run_0488・B' = run_0489 (float) の系列は _band_ab/cold_pair/tr_judge.json。
照合 (欠ければ判定不能、結果を記録): バイナリ・入力の sha256 (lp.log の見出し)、設定の差 (F0 は差 0 行、Bp・Ap は space.pRef の 1 行だけ)、
  初期化後の保存量 9 つのキーの完全一致と値のビット一致 (F0・Bp は旧 B の res_0、Ap は旧 A の res_0)、座標の一致、全残差と使う出力の有限性。
Δ = 100 × (窓の平均 − 自分の step 0) / |旧 A の step 0|。W1 = 20,000〜25,000、W2 = 25,000〜30,000、D = (+1.892, +1.647, +1.431) %。
前提: 全 6 組で |ΔAp − ΔA| ≤ 0.1|D|。F0 ∈ [min(ΔB, ΔB') − N_B, max(ΔB, ΔB') + N_B] かつ sign(D)(ΔF0 − ΔAp) ≥ 0.5|D| (6 組すべて)。
結果 A: 全 6 組で Bp の 2 本がそれぞれ |ΔBp − ΔAp| ≤ 0.1|D|。結果 B: 全 6 組で 2 本がそれぞれ sign(D)(ΔBp − ΔAp) ≥ 0.5|D|。それ以外は判別不能。
Q_w は独立した副判定 (D = 2.114 %)。→ _band_ab/cold_pair/lp_judge.json"""
import csv, gzip, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
NEW = {"F0": "run_0709_lp_f0", "Bp1": "run_0710_lp_bp1", "Bp2": "run_0711_lp_bp2", "Ap": "run_0712_lp_ap"}
BASE = {"F0": "run_0488_tr_f32", "Bp1": "run_0488_tr_f32", "Bp2": "run_0488_tr_f32", "Ap": "run_0487_tr_fp64"}
STEPS = [0] + list(range(20000, 30001, 500)); KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "checks": {}, "records": {}, "delta": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
log = (HERE / "lp.log").read_text(errors="replace")
hd = re.findall(r"F32 ([0-9a-f]{16}) F64 ([0-9a-f]{16}) 共通の入力 ([0-9a-f]{16})", log)
OUT["checks"]["バイナリ・入力 (lp.log)"] = hd[-1] if hd else None
if not hd or hd[-1] != ("793215917d156680", "7e95e5bbd8f50c9e", "ba859e02f3ec236a"): und(f"バイナリ・入力の sha256 が §6.26 と同じと確かめられない: {hd[-1] if hd else None}")
def q9(path):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"] if k in Q9}
with h5py.File(HERE / "run_0487_tr_fp64/nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
for k, r in NEW.items():
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
    with h5py.File(d / "nozzle.h5", "r") as h: same_xy = bool(np.array_equal(np.asarray(h["MESH/COORD"][:], float), xy0))
    if not same_xy: und(f"{r}: 座標が旧 A と違う")
    dl = [l for l in (d / "config_vs_base.diff").read_text().splitlines() if l[:1] in "<>"]
    if k == "F0": ok = (dl == [])
    else: ok = len(dl) == 2 and dl[0].startswith("< space: {pRef: 5499000.0, ") and dl[1].startswith("> space: {") and dl[0].replace("pRef: 5499000.0, ", "")[2:] == dl[1][2:]
    if not ok: und(f"{r}: {BASE[k]} との設定の差が想定と違う: {dl}")
    a, b = q9(d / "res_0.h5"), q9(HERE / BASE[k] / "res_0.h5")
    keys_ok = (set(a) == set(Q9) and set(b) == set(Q9)); bit_ok = keys_ok and all(np.array_equal(a[q], b[q]) for q in Q9)
    if not keys_ok: und(f"{r}: res_0 の保存量のキーが 9 つそろわない (新 {sorted(a)}、旧 {sorted(b)})")
    elif not bit_ok: und(f"{r}: res_0 の保存量が {BASE[k]} とビット一致しない {[q for q in Q9 if not np.array_equal(a[q], b[q])]}")
    OUT["checks"][r] = dict(RUN_RC=rc, 座標一致=same_xy, 設定の差=dl, res0_9キー=keys_ok, res0_ビット一致=bit_ok)
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(d)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
ser = {a: {int(k): v for k, v in tj["series"][a].items()} for a in ("A", "B", "Bp")}
yb_x, yb = XC.common_yb()
verdict = None
if not OUT["undecidable"]:
    for k, r in NEW.items():
        ser[k] = {}
        for s in STEPS:
            x = CS.snapshot(HERE / r, s, yb_x, yb)
            if x.get("nonfinite") or any(not np.isfinite(x[q]) for q in KEYS + (QW,)): und(f"{r}: step {s} の系列に非有限"); break
            ser[k][s] = {q: float(x[q]) for q in KEYS + (QW,)}
        sc = HERE / r / "lp_series.csv"
        with open(sc, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["step"] + list(KEYS) + [QW])
            for s in sorted(t for t in ser[k] if t >= 20000): w.writerow([s] + [ser[k][s][q] for q in KEYS + (QW,)])
        q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + (QW,)),
                            "--drift", "0.0001", "--osc", "0.0001", "--tail", "1.0"], capture_output=True, text=True)
        OUT["records"][f"{r} check_quasisteady (θ_r・Q_w、20,000〜30,000、1e-4)"] = [l.strip() for l in q.stdout.splitlines() if any(t in l for t in KEYS + (QW, "csv"))][:6]
if not OUT["undecidable"]:
    A0 = ser["A"][0]
    def dl_(a, k, w):
        w0, w1 = WIN[w]; return 100 * (np.mean([ser[a][s][k] for s in ser[a] if w0 <= s <= w1]) - ser[a][0][k]) / abs(A0[k])
    pre_ok, ra, rb = True, [], []
    print("窓  断面         ΔA       ΔB       ΔB'      ΔF0      ΔAp      ΔBp1     ΔBp2     0.1|D|  0.5|D|  判定")
    for w in WIN:
        for k in KEYS + (QW,):
            v = {a: dl_(a, k, w) for a in ("A", "B", "Bp", "F0", "Ap", "Bp1", "Bp2")}
            sg, Dk = np.sign(D[k]), abs(D[k]); NB = abs(v["B"] - v["Bp"])
            p1 = abs(v["Ap"] - v["A"]) <= 0.1 * Dk
            p2 = (min(v["B"], v["Bp"]) - NB <= v["F0"] <= max(v["B"], v["Bp"]) + NB) and (sg * (v["F0"] - v["Ap"]) >= 0.5 * Dk)
            A_ = all(abs(v[b] - v["Ap"]) <= 0.1 * Dk for b in ("Bp1", "Bp2")); B_ = all(sg * (v[b] - v["Ap"]) >= 0.5 * Dk for b in ("Bp1", "Bp2"))
            OUT["delta"][f"{w}|{k}"] = dict(**v, pre_Ap=p1, pre_F0=p2, resultA=A_, resultB=B_, resid_Bp1=v["Bp1"] - v["Ap"], resid_Bp2=v["Bp2"] - v["Ap"], vs_F0_Bp1=v["Bp1"] - v["F0"], vs_F0_Bp2=v["Bp2"] - v["F0"])
            print(f"{w}  {k:11s} {v['A']:+7.3f}  {v['B']:+7.3f}  {v['Bp']:+7.3f}  {v['F0']:+7.3f}  {v['Ap']:+7.3f}  {v['Bp1']:+7.3f}  {v['Bp2']:+7.3f}  {0.1 * Dk:.3f}  {0.5 * Dk:.3f}  "
                  f"{'' if p1 else '前提Ap外 '}{'' if p2 else '前提F0外 '}{'A' if A_ else ''}{'B' if B_ else ''}{' (副)' if k == QW else ''}")
            if k != QW: pre_ok &= (p1 and p2); ra.append(A_); rb.append(B_)
    if not pre_ok: verdict = "判別不能 (前提: Ap が旧 A と同じ、または F0 が陽性対照として成立、を満たさない組がある)"
    else: verdict = ("結果 A: 今回の 2 本・登録の窓では、pRef で下流の θ_r の差が 0.1|D| 以内 (90 % 削減と同じ意味ではない。残差の加算が真因・double を廃止できる、とは言わない)" if all(ra)
                     else "結果 B: pRef の変更だけでは下流の θ_r の差を抑えられない" if all(rb) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.46: {verdict}")
(HERE / "_band_ab/cold_pair/lp_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/lp_judge.json")
