"""plan tooling-nozzle-core-grid §4.65 の判定 (事前登録、codex 2026-10-12 twin-grefc)。比べる組は「格子とその背景軌道」(格子間隔だけの効果とは呼ばない)。
入力 (明示): Gref の twin_main_cap50.jsonl・twin_noise_cap50.jsonl・twin_step01_cap50.jsonl (既存、run_0667/0668)、
  Grefc の twin_main_grefc.jsonl・twin_noise_grefc.jsonl・twin_step01_grefc.jsonl (run_0677/0678)、Grefc の摂動 run_0678_twc50_pert/perturb.json。
測定限界 L (各格子): §4.63 と同じ (ρE の雑音比を T に換算)。主比較の窓は Gref と同じ step 1000〜2750 (250 ごとの 8 点) に固定。
資格 (Grefc): 窓の全点で島の D_T ≥ 100 L・max|ΔT/T| < 1e−3・非有限なし (基準・摂動とも)、摂動の実現振幅 max|δρE/ρE| ≤ 1.0e−7 (+丸め)。
λ_c = 窓での ln D_T の step への傾き。λ_ref = Gref の同じ窓の傾き (3.83e−3)。持続増幅の条件は §4.62 と同じ (R² ≥ 0.95、前半・後半とも正で比 2 倍以内、ρ・ρE も正)。
結果 A: 資格を満たし λ_c ≤ λ_ref/3 → 第 1 仮説 (Grefc の背景軌道では有限時間の摂動増幅が弱い) を支持し、「両背景で同程度に増幅する」を退ける。
結果 B: 持続増幅を満たし λ_c ≥ 2λ_ref/3 → 第 1 仮説を退ける。中間・資格不足・境界域 (§4.62 と同じ) は判定不能。
窓より後に強い増幅が現れた場合は「増幅が弱い」でなく「立ち上がりが遅い」と書く (後の窓の傾きも報告する)。
場所: 島・縮流部・全域・上流 (x < −8)・遷移区間の D_T の系列も報告し、島だけ弱く別の場所で育つ場合は全域の増幅の抑制と判定しない (注記)。
usage: python3 cg_twin_grefc_judge.py → _band_ab/core_grid/twin_grefc_judge.json"""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
rd = lambda n: sorted({r["step"]: r for r in (json.loads(l) for l in open(OUT / f"twin_{n}.jsonl"))}.values(), key=lambda r: r["step"])
def fit(st, y):
    A = np.vstack([st, np.ones_like(st)]).T; c, *_ = np.linalg.lstsq(A, y, rcond=None)
    yh = A @ c; ss = np.sum((y - np.mean(y)) ** 2)
    return float(c[0]), float(1 - np.sum((y - yh) ** 2) / ss) if ss > 0 else float("nan")
def floor(tag):
    noise = rd(f"noise_{tag}")[-1]; p1 = [r for r in rd(f"step01_{tag}") if r["step"] == 1][0]
    return p1["island_D_T"] * noise["island_D_roe"] / p1["island_D_roe"]
def window(M, L, lo, hi):
    pts = [r for r in M if lo <= r["step"] <= hi]
    st = np.array([r["step"] for r in pts], float); DT = np.array([r["island_D_T"] for r in pts])
    qual = len(pts) == int((hi - lo) / 250) + 1 and all(r["island_D_T"] >= 100 * L and r["max_dT_rel"] < 1e-3 and r["nonfinite"] == 0 for r in pts)
    if len(pts) < 5: return {"qual": False, "n": len(pts)}
    h = len(pts) // 2
    lam, r2 = fit(st, np.log(DT)); l1, _ = fit(st[:h + 1], np.log(DT[:h + 1])); l2, _ = fit(st[h:], np.log(DT[h:]))
    lro, _ = fit(st, np.log([r["island_D_ro"] for r in pts])); lre, _ = fit(st, np.log([r["island_D_roe"] for r in pts]))
    pers = lam >= 1e-3 and r2 >= 0.95 and l1 > 0 and l2 > 0 and max(l1, l2) / min(l1, l2) <= 2 and lro > 0 and lre > 0
    border = (0.8e-3 <= lam <= 1.2e-3) or (0.93 <= r2 <= 0.97)
    return {"qual": bool(qual), "n": len(pts), "lambda": lam, "R2": r2, "halves": [l1, l2], "lambda_ro": lro, "lambda_roe": lre, "persistent": bool(pers), "borderline": bool(border)}
Mg, Mc = rd("main_cap50"), rd("main_grefc")
Lg, Lc = floor("cap50"), floor("grefc")
wg, wc = window(Mg, Lg, 1000, 2750), window(Mc, Lc, 1000, 2750)
later = {f"{lo}-{lo + 1750}": window(Mc, Lc, lo, lo + 1750) for lo in (2000, 3000, 4000)}
pj = json.loads((HERE / "run_0678_twc50_pert" / "perturb.json").read_text())
amp_ok = pj["max_rel"] <= 1.0e-7 * (1 + 1e-6)
lr = wg["lambda"]
if not (wc.get("qual") and amp_ok):
    res = "判定不能 (資格不足: " + ("窓の信号・振幅・非有限" if not wc.get("qual") else "摂動の振幅") + ")"
elif wc["borderline"] or wg["borderline"]:
    res = "判定不能 (境界域)"
elif wc["lambda"] <= lr / 3:
    res = "A: Grefc の背景軌道では有限時間の摂動増幅が弱い (λ_c ≤ λ_ref/3)"
    if any(v.get("persistent") and v.get("lambda", 0) >= 2 * lr / 3 for v in later.values()): res += "。ただし後の窓で強い増幅 → 「立ち上がりが遅い」"
elif wc["persistent"] and wc["lambda"] >= 2 * lr / 3:
    res = "B: Grefc でも同程度に持続増幅 (λ_c ≥ 2λ_ref/3)"
else:
    res = "判定不能 (中間)"
reg = ("island", "contraction", "all", "upstream", "transition")
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.65", "L": {"Gref": Lg, "Grefc": Lc}, "perturb_grefc": pj,
       "window_1000_2750": {"Gref": wg, "Grefc": wc}, "later_windows_grefc": later,
       "series_grefc": [{"step": r["step"], **{f"{k}_D_T": r.get(f"{k}_D_T") for k in reg}, "max_dT": r["max_dT"], "max_x": r["max_x"], "max_jw": r["max_jw"], "max_dx": r["max_dx"]} for r in Mc],
       "result": res}
(OUT / "twin_grefc_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("L:", {k: f"{v:.2e}" for k, v in out["L"].items()}, " 摂動:", pj["max_rel"], pj["nodes"], pj["mode"])
print("窓 1000〜2750:", {k: {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "halves"} for k, v in out["window_1000_2750"].items()})
print("Grefc の後の窓:", {k: (round(v["lambda"], 5) if "lambda" in v else None, v.get("persistent")) for k, v in later.items()})
for r in out["series_grefc"]:
    print(f"   {r['step']:6d}  島 {r['island_D_T']:.3e}  縮流部 {r['contraction_D_T']:.3e}  上流 {r['upstream_D_T']:.3e}  遷移 {r['transition_D_T']:.3e}  max {r['max_dT']:.2e} @ x {r['max_x']:.2f} jw {r['max_jw']} Δx {r['max_dx']:.4f}")
print("結果:", res)
