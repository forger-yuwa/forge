"""plan architecture-float-state-double-geometry §6.34 (事前登録) の判定: 平均流の状態を FP64 で積む float の軌跡 (qAccumulatorFP64 の感度試験)。qa.sh の後に B で回す。
対照: A = run_0487_tr_fp64 (FP64、保存済みの参照の軌跡)、B = run_0488_tr_f32・B' = run_0489_tr_f32b (float、§6.27)。
追加: E = run_0592_qa_f32・E' = run_0593_qa_f32b (float + qAccumulatorFP64 1)、F = run_0594_qa_off (同じ新しいバイナリ、qAccumulatorFP64 0)。
前提 (一つでも欠ければ判定不能。結果を見て対照を差し替えない):
  1. 単体試験・SASS の照合は qa.sh の外 (implementer の記録) で確かめ、ここでは qa_pre.json の記録を読む。
  2. E・E'・F の全期間の残差・全出力が有限で 61 出力、check_convergence が実行できる。初期化後の保存量 9 つ (res_0) が B・E・E'・F でビット一致。
     設定の差: E・E' は run_0488 と qAccumulatorFP64 の 1 語だけ、F は差なし。
  3. E・E' の両方で、採用の mask (step 1 と 30,000) が軸の節点・内点で 0、採用は壁の節点だけ。monitorInterval ごとの採用数 ≤ 等温壁の節点数 × monitorInterval。
  4. F が B・B' の範囲: 6 組すべてで min(ΔB, ΔB') − N_B ≤ ΔF ≤ max(ΔB, ΔB') + N_B。
判定 (θ_r の 3 断面、W1 = 20,000〜25,000・W2 = 25,000〜30,000、Δ は §6.26 と同じ、D = (+1.892, +1.647, +1.431) %、k = 10):
  mB = mean(ΔB, ΔB')、mE = mean(ΔE, ΔE')、N_E = |ΔE − ΔE'|、N_BF = max(ΔB, ΔB', ΔF) − min(ΔB, ΔB', ΔF)。
  結果 S (抑制): max(|ΔE − ΔA|, |ΔE' − ΔA|) + 10 N_E ≤ 0.1|D|、かつ sign(D)(mB − mE) ≥ 0.5|D|、かつ |mB − mE| > 10 max(N_BF, N_E)。
  結果 R (抑制不十分): E・E' がそれぞれ A から D の向きへ 0.5|D| 以上、かつ |mE − ΔA| > 10 N_E。
  両窓・3 断面で揃わなければ判別不能。Q_w は副判定として記録。結果は _band_ab/cold_pair/qa_judge.json。"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
NEW = {"E": HERE / "run_0592_qa_f32", "Ep": HERE / "run_0593_qa_f32b", "F": HERE / "run_0594_qa_off"}
QA_ON = ("E", "Ep")
OLD = {"A": HERE / "run_0487_tr_fp64", "B": HERE / "run_0488_tr_f32", "Bp": HERE / "run_0489_tr_f32b"}
STEPS = list(range(0, 30001, 500)); KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
K = 10.0
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "pre": {}, "delta": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
# ---- 前提 1: 単体試験・SASS (implementer の記録を qa_pre.json に写したもの) ----
pre1 = HERE / "_band_ab/cold_pair/qa_pre.json"
if not pre1.exists(): und("qa_pre.json (単体試験・SASS の照合の記録) が無い")
else:
    pj = json.loads(pre1.read_text()); OUT["pre"]["単体試験・SASS"] = pj
    if not (pj.get("unit_test_pass") is True and pj.get("sass_mismatch_except_axis") == []): und(f"単体試験・SASS の照合が通っていない: {pj}")
# ---- 前提 2 ----
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
    if n_out != len(STEPS): und(f"{a}: 出力が {n_out} (61 のはず)")
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    if not ov: und(f"{a}: check_convergence が実行できない")
    OUT["pre"][f"{a} check_convergence"] = ov[-1] if ov else None
    diff = [l for l in (run / "config_vs_run_0488.diff").read_text().splitlines() if l[:1] in "<>"]
    if a in QA_ON:
        if len(diff) != 2 or diff[0].replace("qAccumulatorFP64: 1, ", "")[2:] != diff[1][2:]: und(f"{a}: run_0488 との設定の差が qAccumulatorFP64 の 1 語だけでない: {diff}")
    elif diff: und(f"{a}: run_0488 と設定が違う: {diff}")
with h5py.File(OLD["B"] / "res_0.h5", "r") as h0:
    keys = [k for k in h0["VALUE"] if k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega") or re.fullmatch(r"roY\d+", k)]
    ref = {k: np.asarray(h0["VALUE/" + k][:]) for k in keys}
if len(keys) != 9: und(f"初期化後の照合の量が 9 でない: {keys}")
for a, run in NEW.items():
    with h5py.File(run / "res_0.h5", "r") as h:
        bad = [k for k in keys if not np.array_equal(np.asarray(h["VALUE/" + k][:]), ref[k])]
    if bad: und(f"{a}: 初期化後の保存量が B と違う {bad}")
OUT["pre"]["初期化後の保存量の照合 (res_0)"] = keys
# ---- 前提 3: 採用の mask と採用数 ----
for a in QA_ON:
    run = NEW[a]; log = (run / "forge_run.log").read_text(errors="replace")
    if "FORGE_QACC_AXISYM_DIAG" not in log and "診断用" not in log: und(f"{a}: 軸対称の診断用の経路で起動した記録がログに無い")
    cnt = [int(m) for m in re.findall(r"\[qAccumulatorFP64\] step \d+: reconcile 採用 (\d+)", log)]
    OUT["pre"][f"{a} 採用数 (monitorInterval ごと) 最小・最大・個数"] = [min(cnt), max(cnt), len(cnt)] if cnt else None
    if not cnt: und(f"{a}: 採用数の記録が無い")
    for st in (1, 30000):
        f = run / f"qacc_adopt_{st}.h5"
        if not f.exists(): und(f"{a}: step {st} の採用の mask {f.name} が無い"); continue
        with h5py.File(f, "r") as h:
            ax = np.asarray(h["cells/axis_flag"][:]) == 1; wl = np.asarray(h["cells/wall_flag"][:]) == 1
            iso = np.asarray(h["cells/iso_wall_flag"][:]) == 1; bn = np.asarray(h["cells/bnode_flag"][:]) == 1
            rec = {}
            for q in ("ro", "roUx", "roUy", "roUz", "roe"):
                m = np.asarray(h[f"{q}/mask"][:]).astype(bool)
                rec[q] = dict(total=int(m.sum()), axis=int((m & ax).sum()), interior=int((m & ~bn & ~ax).sum()), wall=int((m & wl).sum()),
                              iso_wall=int((m & iso).sum()), outside_wall=int((m & ~wl).sum()))
                if rec[q]["axis"] or rec[q]["interior"] or rec[q]["outside_wall"]:
                    und(f"{a}: step {st} の {q} で軸 {rec[q]['axis']}・内点 {rec[q]['interior']}・壁の外 {rec[q]['outside_wall']} 節点が採用された")
            rec["n_axis"] = int(ax.sum()); rec["n_wall"] = int(wl.sum()); rec["n_iso_wall"] = int(iso.sum())
            rec["axis_and_wall_overlap"] = int((ax & wl).sum()); rec["axis_and_iso_wall_overlap"] = int((ax & iso).sum())
        OUT["pre"][f"{a} 採用の mask step {st}"] = rec
    niso = OUT["pre"].get(f"{a} 採用の mask step 1", {}).get("n_iso_wall")
    stp = [int(m) for m in re.findall(r"\[qAccumulatorFP64\] step (\d+): reconcile 採用", log)]
    mi = min(np.diff(stp)) if len(stp) > 1 else None   # monitorInterval (採用数を出す間隔)
    OUT["pre"][f"{a} monitorInterval (ログの間隔)"] = int(mi) if mi else None
    if not mi: und(f"{a}: monitorInterval をログから決められない")
    elif cnt and niso is not None and max(cnt) > niso * mi: und(f"{a}: 採用数 {max(cnt)} が 等温壁の節点数 × monitorInterval ({niso} × {mi}) を超えた")
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
    st = lambda a: ser[a][min(ser[a])]; A0 = st("A"); rs, rr, fok = [], [], []
    print("窓  断面         ΔA       ΔB       ΔB'      ΔF       ΔE       ΔE'      mB−mE   N_B    N_BF   N_E    0.5|D|  0.1|D|")
    for wn, (w0, w1) in WIN.items():
        for k in KEYS + (QW,):
            dl = {a: 100 * (np.mean([ser[a][s][k] for s in ser[a] if w0 <= s <= w1]) - st(a)[k]) / abs(A0[k]) for a in ("A", "B", "Bp", "F", "E", "Ep")}
            mB, mE = 0.5 * (dl["B"] + dl["Bp"]), 0.5 * (dl["E"] + dl["Ep"])
            NB, NE = abs(dl["B"] - dl["Bp"]), abs(dl["E"] - dl["Ep"])
            NBF = max(dl["B"], dl["Bp"], dl["F"]) - min(dl["B"], dl["Bp"], dl["F"])
            sg = np.sign(D[k]); Dk = abs(D[k])
            f_ = (min(dl["B"], dl["Bp"]) - NB <= dl["F"] <= max(dl["B"], dl["Bp"]) + NB)
            s_ = (max(abs(dl["E"] - dl["A"]), abs(dl["Ep"] - dl["A"])) + K * NE <= 0.1 * Dk) and (sg * (mB - mE) >= 0.5 * Dk) and (abs(mB - mE) > K * max(NBF, NE))
            r_ = (sg * (dl["E"] - dl["A"]) >= 0.5 * Dk) and (sg * (dl["Ep"] - dl["A"]) >= 0.5 * Dk) and (abs(mE - dl["A"]) > K * NE)
            OUT["delta"][f"{wn}|{k}"] = dict(**dl, mB=mB, mE=mE, NB=NB, NBF=NBF, NE=NE, F_in_range=bool(f_), resultS=bool(s_), resultR=bool(r_))
            print(f"{wn}  {k:11s} {dl['A']:+7.3f}  {dl['B']:+7.3f}  {dl['Bp']:+7.3f}  {dl['F']:+7.3f}  {dl['E']:+7.3f}  {dl['Ep']:+7.3f}  {mB - mE:+6.3f}  {NB:.3f}  {NBF:.3f}  {NE:.3f}  {0.5 * Dk:.3f}  {0.1 * Dk:.3f}  {'F外 ' if not f_ else ''}{'S' if s_ else ''}{'R' if r_ else ''}{' (副)' if k == QW else ''}")
            if k != QW: rs.append(s_); rr.append(r_); fok.append(f_)
    if not all(fok):
        verdict = "判定不能: F が B・B' の範囲に入らない組がある (前提 4。対照を差し替えない)"
    else:
        verdict = ("結果 S: 抑制 (平均流の保存量の累積の方法を変えると、登録の窓の参照 A からの偏りを抑えられた。吸収が原因とは言わない)" if all(rs)
                   else "結果 R: 抑制不十分 (この変更だけでは抑えられない。残差の評価・増分の生成・float のままの k/ω・ρY のどれかは決まらない)" if all(rr) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.34: {verdict}")
OUT["series_new"] = {a: ser.get(a, {}) for a in NEW}
(HERE / "_band_ab/cold_pair/qa_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/qa_judge.json")
