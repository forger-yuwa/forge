"""plan architecture-float-state-double-geometry §6.32 (事前登録) の判定: block の系の組み立ての精度の感度試験 (implicitSolvePrecision 0 → 1)。ds.sh の後に B で回す。
対照: A = run_0487_tr_fp64 (FP64、参照の軌跡)、B = run_0488_tr_f32・B' = run_0489_tr_f32b (float、§6.27)。追加: C = run_0590_ds_f32・C' = run_0591_ds_f32b (float + implicitSolvePrecision 1)。
前提 (総合の判定の前に確かめる。欠ければ判定不能):
  C・C' の全期間の残差・全出力が有限で 61 出力、check_convergence が実行できる。C と B の設定の差は implicitSolvePrecision の 1 行だけ。
  バイナリの sha256 が §6.26 の float と同じ (tr.log と ds.log)。初期化後の保存量 (res_0 の ro..roe・roK・roOmega・roY*) が B・B'・C・C' でビット一致。
判定 (θ_r の 3 断面、W1 = 20,000〜25,000・W2 = 25,000〜30,000、Δ は §6.26 と同じ、D = (+1.892, +1.647, +1.431) %):
  mB = mean(ΔB, ΔB')、mC = mean(ΔC, ΔC')、N_B = |ΔB − ΔB'|、N_C = |ΔC − ΔC'|。
  結果 A (抑制): max(|ΔC − ΔA|, |ΔC' − ΔA|) + 10 N_C ≤ 0.1|D|、かつ sign(D)(mB − mC) ≥ 0.5|D|、かつ |mB − mC| > 10 max(N_B, N_C)。
  結果 B (抑制不十分): C・C' がそれぞれ A から D の向きへ 0.5|D| 以上、かつ |mC − ΔA| > 10 N_C。
  両窓・3 断面で揃わなければ判別不能。Q_w は副判定として記録。結果は _band_ab/cold_pair/ds_judge.json。"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
NEW = {"C": HERE / "run_0590_ds_f32", "Cp": HERE / "run_0591_ds_f32b"}
OLD = {"A": HERE / "run_0487_tr_fp64", "B": HERE / "run_0488_tr_f32", "Bp": HERE / "run_0489_tr_f32b"}
STEPS = list(range(0, 30001, 500)); KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "pre": {}, "delta": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
# ---- 前提 ----
for a, run in NEW.items():
    rows = [x for x in csv.DictReader(open(run / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{a}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 30000: und(f"{a}: 30,000 step に届かない")
    n_out = 0
    for s in STEPS:
        p = run / f"res_{s}.h5"
        if not p.exists(): und(f"{a}: res_{s} が無い"); break
        n_out += 1
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][k].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][k][:])) for k in h["VALUE"]): und(f"{a}: res_{s} に非有限"); break
    OUT["pre"][f"{a} 出力の数"] = n_out
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    if not ov: und(f"{a}: check_convergence が実行できない")
    OUT["pre"][f"{a} check_convergence"] = ov[-1] if ov else None
    diff = [l for l in (run / "config_vs_run_0488.diff").read_text().splitlines() if l[:1] in "<>"]
    if len(diff) != 2 or not all("implicitSolvePrecision" in l or "deltaT" in l for l in diff): und(f"{a}: run_0488 との設定の差が想定と違う: {diff}")
    elif diff[0].replace("implicitSolvePrecision: 1, ", "")[2:] != diff[1][2:]: und(f"{a}: 設定の差が implicitSolvePrecision だけでない")
sha = lambda f, pat: (re.findall(pat, (HERE / f).read_text()) or [None])[0]
s_tr = sha("tr.log", r"F32 ([0-9a-f]{16})"); s_ds = sha("ds.log", r"F32 ([0-9a-f]{64})")
OUT["pre"]["バイナリ sha256 (tr.log 先頭 16 / ds.log)"] = [s_tr, s_ds]
if not (s_tr and s_ds and s_ds.startswith(s_tr)): und("バイナリが §6.26 の float と同じと確かめられない")
with h5py.File(OLD["B"] / "res_0.h5", "r") as h0:
    keys = [k for k in h0["VALUE"] if k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega") or re.fullmatch(r"roY\d+", k)]
    ref = {k: np.asarray(h0["VALUE/" + k][:]) for k in keys}
for a, run in {"Bp": OLD["Bp"], **NEW}.items():
    with h5py.File(run / "res_0.h5", "r") as h:
        bad = [k for k in keys if not np.array_equal(np.asarray(h["VALUE/" + k][:]), ref[k])]
    if bad: und(f"{a}: 初期化後の保存量が B と違う {bad}")
OUT["pre"]["初期化後の保存量の照合 (res_0)"] = keys
print("前提:", json.dumps(OUT["pre"], ensure_ascii=False, default=str))
# ---- 系列 ----
tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
ser = {a: {int(k): v for k, v in tj["series"][a].items()} for a in ("A", "B", "Bp")}
yb_x, yb = XC.common_yb()
if not OUT["undecidable"]:
    for a, run in NEW.items():
        ser[a] = {}
        for s in STEPS:
            r = CS.snapshot(run, s, yb_x, yb)
            if r.get("nonfinite") or any(not np.isfinite(r[k]) for k in KEYS + (QW,)): und(f"{a}: step {s} の系列に非有限"); break
            ser[a][s] = {k: float(r[k]) for k in KEYS + (QW,)}
verdict = None
if not OUT["undecidable"]:
    st = lambda a: ser[a][min(ser[a])]; A0 = st("A"); ra, rb = [], []
    print("窓  断面         ΔA       ΔB       ΔB'      ΔC       ΔC'      mB−mC   N_B    N_C    0.5|D|  0.1|D|")
    for wn, (w0, w1) in WIN.items():
        for k in KEYS + (QW,):
            dl = {a: 100 * (np.mean([ser[a][s][k] for s in ser[a] if w0 <= s <= w1]) - st(a)[k]) / abs(A0[k]) for a in ("A", "B", "Bp", "C", "Cp")}
            mB, mC = 0.5 * (dl["B"] + dl["Bp"]), 0.5 * (dl["C"] + dl["Cp"]); NB, NC = abs(dl["B"] - dl["Bp"]), abs(dl["C"] - dl["Cp"])
            sg = np.sign(D[k]); Dk = abs(D[k])
            a_ = (max(abs(dl["C"] - dl["A"]), abs(dl["Cp"] - dl["A"])) + 10 * NC <= 0.1 * Dk) and (sg * (mB - mC) >= 0.5 * Dk) and (abs(mB - mC) > 10 * max(NB, NC))
            b_ = (sg * (dl["C"] - dl["A"]) >= 0.5 * Dk) and (sg * (dl["Cp"] - dl["A"]) >= 0.5 * Dk) and (abs(mC - dl["A"]) > 10 * NC)
            OUT["delta"][f"{wn}|{k}"] = dict(**dl, mB=mB, mC=mC, NB=NB, NC=NC, resultA=bool(a_), resultB=bool(b_))
            print(f"{wn}  {k:11s} {dl['A']:+7.3f}  {dl['B']:+7.3f}  {dl['Bp']:+7.3f}  {dl['C']:+7.3f}  {dl['Cp']:+7.3f}  {mB - mC:+6.3f}  {NB:.3f}  {NC:.3f}  {0.5 * Dk:.3f}  {0.1 * Dk:.3f}  {'A' if a_ else ''}{'B' if b_ else ''}{' (副)' if k == QW else ''}")
            if k != QW: ra.append(a_); rb.append(b_)
    verdict = ("結果 A: 抑制 (block の系の組み立ての double 化で、登録の窓の偏りを抑えられた。発生源・必要性の証明ではない)" if all(ra)
               else "結果 B: 抑制不十分 (この変更だけでは登録の窓の偏りを抑えられない。残差・commit だけが原因とは言わない)" if all(rb) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.32: {verdict}")
OUT["series_new"] = {a: ser.get(a, {}) for a in NEW}
(HERE / "_band_ab/cold_pair/ds_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/ds_judge.json")
