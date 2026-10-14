"""plan tooling-nozzle-core-grid §4.100 P-hix の判定 (事前登録、codex 諮問 2026-10-14 hi-next): 入口半径を半分にした形 H の継続 対 固定した今の形 C。
H = run_0916_hi_h_cont (run_0915_hi_h_main の res_60000 からの継続、1000 ごと、step は継続の中の番号 1000〜40000 = 本段の通算 61k〜100k)。
C = run_0912_hi_c_main の本段の窓 (40k, 60k] の平均を固定の対照にする (C は延長しない。同じ計算時間での速度の比ではなく、固定した停滞の水準との比)。
主指標: 縮流部 x/r_t ∈ [−8, −2] の q = √(ΣV(R/V)²/ΣV) の ρ・ρE。H の窓 W3 = 継続の (0, 20k]・W4 = (20k, 40k] (本段の通算 (60k, 80k]・(80k, 100k]) の平均 q̄。
  比 = q̄_H(Wk) / q̄_C(40k, 60k]、H の d = q̄_H(W4) / q̄_H(W3)。
結果 A: 両窓・ρ・ρE で比 ≤ 1/3、H の d ≤ 1.1 (ρ・ρE)、かつ両窓で 縮流部の ρUx・ρUy・ρk・ρω、縮流部の外の ρ・ρUx・ρUy・ρE、throatwall の ρ・ρE の比 ≤ 1.5。
結果 B: 両窓・ρ・ρE で比 ≥ 2/3、かつ H の d ≥ 0.8 (ρ・ρE)。それ以外は中間 (記録する)。残差が下がり続けることは資格違反にしない。
資格 (測定の成立条件だけ): (1) restart.log が「VERDICT: OK (9 量を移した、SRC とビット一致)」。(2) cg_hix.log に設定の一致の行 (H の本段と nStepOuter だけが違い bcond 同一)。
  (3) RUN_RC 0・継続の各窓に 20 出力・C の (40k, 60k] に 20 出力。(4) 全出力の残差・物理量が有限、ρ・P・T > 0。資格不足は判定不能。
別の記録 (判定 A/B を変えない): θ_r[65,75]・Q_w の継続の 40 出力を check_quasisteady (CSV 系列モード、既定の閾値) で判定 → 両方 STEADY なら A は「改善がその二量の過渡中だけだった」を退け、
  B は第 1 仮説 (改善が持続する) を退ける。STEADY でなければその切り分けは未完了と書く。確認できる範囲はこの二量に限る。
  全保存量の check_convergence: 継続だけ (再開直後の跳ねから測る) と、本段 + 継続を連結したもの (同じ実効設定) の 2 つの VERDICT を残す。A でも「収束した」とは呼ばない。
usage (AWS の case/45 で): python3 cg_hix_judge.py → _band_ab/core_grid/hix_judge.json"""
import csv, json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
TL = HERE / ".." / ".." / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TL))
from check_quasisteady import classify_series  # noqa: E402
C, HM, H = "run_0912_hi_c_main", "run_0915_hi_h_main", "run_0916_hi_h_cont"
WC = (40000, 60000)
WH = ((0, 20000), (20000, 40000))
OFF = 60000
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
R = {r: jl(f"region_res_{r}.jsonl") for r in (C, H)}; PH = {r: jl(f"phys_{r}.jsonl") for r in (C, H)}
def wm(run, key, lo, hi):
    v = [x[key] for x in R[run] if lo < x["step"] <= hi and key in x]; return float(np.mean(v)) if len(v) == 20 else None
qual = {}
rl = (HERE / H / "restart.log"); qual["restart_bitwise_9"] = rl.exists() and "VERDICT: OK (9 量を移した、SRC とビット一致)" in rl.read_text()
lg = HERE / "cg_hix.log"; qual["settings_same_as_H_main"] = lg.exists() and f"{H}: 設定は H の本段と nStepOuter だけが違う (bcond 同一)" in lg.read_text()
rc = HERE / H / "RUN_RC"
main = ("res_ro", "res_roe")
side = [("contraction", c) for c in ("res_roUx", "res_roUy", "res_roK", "res_roOmega")] + [("outside", c) for c in ("res_ro", "res_roUx", "res_roUy", "res_roe")] \
    + [("throatwall", c) for c in main]
keys = [f"contraction_{c}" for c in main] + [f"{rg}_{c}" for rg, c in side]
qual["completed_and_20_outputs"] = rc.exists() and rc.read_text().strip() == "0" and all(wm(H, k, *w) is not None for k in keys for w in WH) \
    and all(wm(C, k, *WC) is not None for k in keys)
qual["finite_physical"] = not any(x.get("nonfinite") for r in (H,) for x in R[r] + PH[r]) and len(PH[H]) >= 40 \
    and all(x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for x in PH[H])
ratio, d = {}, {}
for wi, w in enumerate(WH):
    ratio[f"W{wi + 3}"] = {k: (wm(H, k, *w) / wm(C, k, *WC)) if (wm(H, k, *w) is not None and wm(C, k, *WC)) else None for k in keys}
for c in main:
    a, b = wm(H, f"contraction_{c}", *WH[0]), wm(H, f"contraction_{c}", *WH[1])
    d[c] = (b / a) if (a and b is not None) else None
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
else:
    mr = [ratio[w][f"contraction_{c}"] for w in ratio for c in main]
    sr = [ratio[w][f"{rg}_{c}"] for w in ratio for rg, c in side]
    if all(v <= 1 / 3 for v in mr) and all(v <= 1.1 for v in d.values()) and all(v <= 1.5 for v in sr):
        res = "A: H の縮流部の残差の改善は、継続の 4 万 step (本段の通算 60k〜100k) でも固定した C の停滞の水準の 1/3 以下に保たれた"
    elif all(v >= 2 / 3 for v in mr) and all(v >= 0.8 for v in d.values()):
        res = "B: H の改善は継続の窓で C の停滞の水準の 2/3 以上に戻った (改善が持続する、を退ける)"
    else:
        res = "中間 (A・B のどちらの条件も満たさない。記録する)"
# 物理量の準定常 (別の記録)
csvp = OUTD / "hix_phys_series.csv"; rows = []
for k in range(1000, 40001, 1000):
    p = OUTD / f"gref_series_{H}_{k}.csv"
    if p.exists():
        r = list(csv.DictReader(open(p)))
        if r: rows.append({"step": OFF + k, "theta_r_w65_75": r[-1]["theta_r_w65_75"], "Q_w": r[-1]["Q_w"]})
with open(csvp, "w", newline="") as f:
    w_ = csv.DictWriter(f, fieldnames=["step", "theta_r_w65_75", "Q_w"]); w_.writeheader(); w_.writerows(rows)
qs = subprocess.run([sys.executable, str(TL / "check_quasisteady.py"), "--series-csv", str(csvp), "--series-cols", "theta_r_w65_75,Q_w"], capture_output=True, text=True)
qs_v = {}
for q in ("theta_r_w65_75", "Q_w"):
    st = [r["step"] for r in rows]; vv = [float(r[q]) for r in rows]
    v, det, ma = classify_series(st, vv, 0.4, 0.05, 0.10, 4) if len(vv) == 40 else ("points", f"{len(vv)} 点", None)
    qs_v[q] = {"verdict": v, "detail": det, "n": len(vv), "first": vv[0] if vv else None, "last": vv[-1] if vv else None,
               "w3_mean": float(np.mean(vv[:20])) if len(vv) == 40 else None, "w4_mean": float(np.mean(vv[20:])) if len(vv) == 40 else None}
steady = all(qs_v[q]["verdict"] == "STEADY" for q in qs_v)
if res.startswith("A"):
    reading = "A かつ θ_r・Q_w が STEADY: 観測した改善がこの二量の過渡中だけだった、を退ける" if steady else "A だが θ_r・Q_w の少なくとも一方が STEADY でない: 過渡中だけだったかの切り分けは未完了 (物理量は未定常)"
elif res.startswith("B"):
    reading = "B かつ θ_r・Q_w が STEADY: 第 1 仮説 (改善が持続する) を退ける" if steady else "B だが物理量は未定常: 切り分けは未完了"
else:
    reading = "中間 または 判定不能: 物理量の判定は記録だけ"
# 全保存量の収束 (継続だけ と 本段 + 継続の連結)
cc = OUTD / "hix_concat_main_cont"; cc.mkdir(exist_ok=True)
with open(cc / "residual_history.csv", "w") as fo:
    for i, rr in enumerate((HERE / HM / "residual_history.csv", HERE / H / "residual_history.csv")):
        with open(rr) as fi:
            hdr = fi.readline()
            if i == 0: fo.write(hdr)
            for ln in fi: fo.write(ln)
conv = {}
for lab, rd in (("cont_only", HERE / H), ("main_plus_cont", cc)):
    o = subprocess.run([sys.executable, str(TL / "check_convergence.py"), str(rd)], capture_output=True, text=True)
    conv[lab] = o.stdout.strip().splitlines()[-25:]
(cc / "residual_history.csv").unlink()
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.100", "C_fixed_window": [C, WC], "H": H, "H_windows_cont_steps": WH, "offset": OFF,
       "qualification": qual, "ratio_H_over_Cfixed": ratio, "d_H": d, "result": res, "phys_quasisteady": qs_v, "phys_tool_stdout": qs.stdout.strip().splitlines()[-8:],
       "reading": reading, "convergence": conv,
       "series_contraction_ro_H": [(OFF + x["step"], x["contraction_res_ro"]) for x in R[H] if x["step"] % 5000 == 0],
       "note": "継続の窓 W3 の始めには再開直後の残差の跳ねが入る (C の固定窓には無い)。A に不利な向き。"}
(OUTD / "hix_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
fmt = lambda v: "-" if v is None else f"{v:.3g}"
for w, rr in ratio.items(): print(w, "  ".join(f"{k} {fmt(v)}" for k, v in rr.items()))
print("H の d:", {c: fmt(v) for c, v in d.items()}); print("資格:", qual); print("結果:", res)
print("物理量:", {q: (v["verdict"], v["detail"]) for q, v in qs_v.items()}); print("読み:", reading)
for lab, ls in conv.items(): print(lab, *ls, sep="\n  ")
