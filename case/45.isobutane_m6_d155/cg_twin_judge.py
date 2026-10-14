"""plan tooling-nozzle-core-grid §4.62 の判定 (事前登録、codex 2026-10-12 twin)。入力: _band_ab/core_grid/twin_{noise,step01,main}_cap{50,1}.jsonl と各 pert の perturb.json。
主指標は島の D_T(n) = √[Σ V (T⁺ − T)² / Σ V]。測定限界 L = 雑音の床。**訂正 (2026-10-12、本番の結果を見る前、§4.63)**: 出力の T は更新前の状態の値 (1 step 遅れ) で、
同じ入力の 2 本の 1 step 後の T の差は 0 になるので、L は保存量 ρE の雑音比から T に換算する: L = D_T,pert(1) × D_ρE,noise(1) / D_ρE,pert(1) (島)。
適格な点: step ≥ 1000、島の D_T ≥ 100 L、max|ΔT/T| < 1e−3。適格な点が連続する区間 (250 ごと) のうち 5 点以上・1000 step 以上の最大の区間をすべて報告する。
各区間で ln D_T を step に当てはめて λ (1/step)・R²・前半と後半の傾きを出し、ρ・ρE (島) の λ も出す。
持続増幅: λ ≥ 1e−3、R² ≥ 0.95、前半・後半の傾きがともに正で比が 2 倍以内、ρ・ρE の λ も正。境界付近 (λ が 0.8e−3〜1.2e−3、R² が 0.93〜0.97) は判定不能扱い。
結果 A: 上限 50 のどれかの区間が持続増幅を満たし、上限 1 は同じ step の範囲で信号を保ったまま (適格な点で) 減衰する (λ < 0)。
結果 B: 上限 50 に適格な区間があり、どれも持続増幅を満たさない (境界付近でない)。または上限 1 も持続増幅を満たし λ が上限 50 の 1/2 以上。
それ以外 (適格な区間が無い・境界付近・成分の不一致・上限 1 に比較区間の信号が無い) → 判定不能。λ = 0 や λ < 0 を無信号の区間に割り当てない。
usage: python3 cg_twin_judge.py → _band_ab/core_grid/twin_judge.json"""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
def rd(name):
    return sorted((json.loads(l) for l in open(OUT / f"twin_{name}.jsonl")), key=lambda r: r["step"])
def fit(st, y):
    A = np.vstack([st, np.ones_like(st)]).T; c, *_ = np.linalg.lstsq(A, y, rcond=None)
    yh = A @ c; ss = np.sum((y - np.mean(y)) ** 2)
    return float(c[0]), float(1 - np.sum((y - yh) ** 2) / ss) if ss > 0 else float("nan")
res = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.62", "caps": {}}
sha = {}
for cap in ("50", "1"):
    for k in ("pert1", "pert"):
        pj = next(HERE.glob(f"run_*_tw{cap}_{k}/perturb.json"), None)
        sha[f"{cap}_{k}"] = json.loads(pj.read_text())["sha256_delta"] if pj else None
res["perturb_same"] = len(set(sha.values())) == 1 and None not in sha.values()
for cap in ("50", "1"):
    noise = rd(f"noise_cap{cap}")[-1]
    s01 = rd(f"step01_cap{cap}"); M = rd(f"main_cap{cap}")
    p1 = [r for r in s01 if r["step"] == 1][0]
    L = p1["island_D_T"] * noise["island_D_roe"] / p1["island_D_roe"]   # §4.63 の訂正 (T は 1 step 遅れ)
    st = np.array([r["step"] for r in M], float); D = {k: np.array([r[f"island_D_{k}"] for r in M]) for k in ("T", "ro", "roe")}
    rel = np.array([r["max_dT_rel"] for r in M]); nf = sum(r["nonfinite"] for r in M)
    q = (st >= 1000) & (D["T"] >= 100 * max(L, 1e-300)) & (rel < 1e-3)
    runs, cur = [], []
    for k in range(len(st)):
        if q[k] and (not cur or st[k] - st[cur[-1]] == 250): cur.append(k)
        else:
            if len(cur) >= 5 and st[cur[-1]] - st[cur[0]] >= 1000: runs.append(cur)
            cur = [k] if q[k] else []
    if len(cur) >= 5 and st[cur[-1]] - st[cur[0]] >= 1000: runs.append(cur)
    segs = []
    for r in runs:
        s = st[r]; h = len(r) // 2
        lam, r2 = fit(s, np.log(D["T"][r])); l1, _ = fit(s[:h + 1], np.log(D["T"][r][:h + 1])); l2, _ = fit(s[h:], np.log(D["T"][r][h:]))
        lro, _ = fit(s, np.log(D["ro"][r])); lre, _ = fit(s, np.log(D["roe"][r]))
        border = (0.8e-3 <= lam <= 1.2e-3) or (0.93 <= r2 <= 0.97)
        pers = lam >= 1e-3 and r2 >= 0.95 and l1 > 0 and l2 > 0 and max(l1, l2) / min(l1, l2) <= 2 and lro > 0 and lre > 0
        segs.append({"steps": [float(s[0]), float(s[-1])], "n": len(r), "lambda": lam, "R2": r2, "lambda_halves": [l1, l2], "lambda_ro": lro, "lambda_roe": lre,
                     "persistent": bool(pers), "borderline": bool(border)})
    res["caps"][cap] = {"noise_island_D_T": L, "noise_island_D_T_raw": noise["island_D_T"], "noise_island_D_roe": noise["island_D_roe"], "step0_1": s01, "nonfinite": nf, "segments": segs,
                        "series": [{"step": r["step"], "D_T": r["island_D_T"], "D_ro": r["island_D_ro"], "D_roe": r["island_D_roe"], "max_dT": r["max_dT"], "max_dT_rel": r["max_dT_rel"],
                                    "max_x": r["max_x"], "max_jw": r["max_jw"], "max_dx": r["max_dx"], "dt_ratio": r.get("island_dt_ratio_median")} for r in M]}
c50, c1 = res["caps"]["50"], res["caps"]["1"]
pers50 = [s for s in c50["segments"] if s["persistent"]]
border50 = any(s["borderline"] for s in c50["segments"])
def cap1_decays(seg):
    """上限 1 の、seg と重なる step の適格な点で λ < 0 か (信号を保ったまま減衰)"""
    S = c1["series"]; L1 = c1["noise_island_D_T"]
    pts = [r for r in S if seg["steps"][0] <= r["step"] <= seg["steps"][1] and r["D_T"] >= 100 * L1 and r["max_dT_rel"] < 1e-3]
    if len(pts) < 5: return None
    lam, _ = fit(np.array([r["step"] for r in pts], float), np.log([r["D_T"] for r in pts]))
    return lam
if not res["perturb_same"] or c50["nonfinite"] or c1["nonfinite"]:
    out = "判定不能 (摂動が両上限で同じでない、または非有限)"
elif pers50 and not border50:
    lams = [cap1_decays(s) for s in pers50]
    if all(l is not None and l < 0 for l in lams): out = "A: 上限 50 だけが持続増幅、上限 1 は同じ範囲で信号を保ったまま減衰"
    elif any(s["persistent"] and s["lambda"] >= 0.5 * max(p["lambda"] for p in pers50) for s in c1["segments"]): out = "B: 上限 1 も同程度に持続増幅"
    else: out = "判定不能 (上限 1 の比較区間に信号が無い、または減衰でない)"
elif c50["segments"] and not border50:
    out = "B: 上限 50 に十分な信号の区間があるが持続増幅を満たさない"
else:
    out = "判定不能 (適格な区間が無い、または境界付近)"
res["result"] = out
(OUT / "twin_judge.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
for cap, c in res["caps"].items():
    print(f"上限 {cap}: 雑音の床 {c['noise_island_D_T']:.2e} K、区間 {[(s['steps'], round(s['lambda'], 5), round(s['R2'], 3), s['persistent']) for s in c['segments']]}")
    for r in c["series"]:
        print(f"   {r['step']:6d}  D_T {r['D_T']:.3e}  max|ΔT| {r['max_dT']:.3e} @ x {r['max_x']:.2f} jw {r['max_jw']}  Δx {r['max_dx']:.4f}")
print("摂動の一致:", res["perturb_same"], " 結果:", out)
