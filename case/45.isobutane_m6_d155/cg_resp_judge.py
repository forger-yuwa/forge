"""plan tooling-nozzle-core-grid §4.46 M-r の判定 (事前登録、codex 2026-10-11 response)。入力は cg_resp_an.py の resp_<tag50>.json・resp_<tag1>.json。
資格 (両腕): q1・q2 ≤ 1e−6 (平均流の増分 = δ、R₀ = b₀)、q3 SST の凍結の実効、q4 組成の変化 ≤ 0.01·max|δρ/ρ|、q5 (ρ・ρE) ≤ 0.1ε、
  q6 島の ρ・ρE で ‖g‖ ≤ 0.01‖ΔR_lhs‖、q7 島の ρ・ρE で ‖ΔR_lhs‖/‖R₀‖ ≥ 100ε かつ |d_lhs| ≥ 100ε。
結果 A: 上限 50 の島の ρ・ρE とも m ≥ 0.5 かつ d_lhs > 0 かつ d_act ≤ 0.5·d_lhs、さらに上限 1 の不足 (d_lhs − d_act) ≤ 上限 50 の不足の 1/3
  → 上限に敏感な応答の不整合を支持し、「時間項・緩和項だけで実際の応答を説明できる」を退ける。
結果 B: 両腕の島の ρ・ρE とも m ≤ 0.1 かつ |d_act − d_lhs| ≤ 0.1·|d_lhs| → 大きな応答の不整合を退ける
  (上限 50 の遅さまで予測が再現すれば第 2 仮説 (時間項と緩和が作る実効の前処理そのものが遅い) を支持し、Tδ と cLδ の射影から内訳を読む)。
それ以外・資格不足・信号不足 → 判定不能。m が大きいだけでは結果 A にしない。
usage: python3 cg_resp_judge.py <tag50> <tag1> → _band_ab/core_grid/resp_judge.json"""
import json, sys
from pathlib import Path
OUT = Path(__file__).resolve().parent / "_band_ab" / "core_grid"
t50, t1 = sys.argv[1:3]
J = {t: json.loads((OUT / f"resp_{t}.json").read_text()) for t in (t50, t1)}
MAIN = ("ro", "roe")
qual = {}
for t, X in J.items():
    q = X["qual_raw"]; I = X["region"]["island"]
    qual[t] = {"q1_q2": all(q[f"q1_{c}"] <= 1e-6 and q[f"q2_{c}"] <= 1e-6 for c in ("ro", "roUx", "roUy", "roe")),
               "q3": q["q3_sst_frozen"], "q4": q["q4_ok"],
               "q5": all(q[f"q5_{c}"] <= 0.1 * q[f"eps_{c}"] for c in MAIN),
               "q6": all(I[c]["g_over_dRlhs"] <= 0.01 for c in MAIN),
               "q7": all(I[c]["dRlhs_over_R0"] >= 100 * q[f"eps_{c}"] and abs(I[c]["d_lhs"]) >= 100 * q[f"eps_{c}"] for c in MAIN)}
ok = all(all(v.values()) for v in qual.values())
I50, I1 = J[t50]["region"]["island"], J[t1]["region"]["island"]
summary = {t: {c: {k: J[t]["region"]["island"][c][k] for k in ("m", "rho_act", "rho_lhs", "d_act", "d_lhs", "deficit", "dRlhs_over_R0")} for c in MAIN} for t in J}
if not ok:
    res = "判定不能 (資格不足: " + ", ".join(f"{t}:{k}" for t, v in qual.items() for k, b in v.items() if not b) + ")"
elif all(I50[c]["m"] >= 0.5 and I50[c]["d_lhs"] > 0 and I50[c]["d_act"] <= 0.5 * I50[c]["d_lhs"] and I1[c]["deficit"] <= I50[c]["deficit"] / 3 for c in MAIN):
    res = "A: 上限に敏感な応答の不整合を支持 (時間項・緩和項だけでは実際の応答を説明できない)"
elif all(X[c]["m"] <= 0.1 and abs(X[c]["d_act"] - X[c]["d_lhs"]) <= 0.1 * abs(X[c]["d_lhs"]) for X in (I50, I1) for c in MAIN):
    res = "B: 大きな応答の不整合を退ける (内訳は Tδ・cLδ の射影で読む)"
else:
    res = "判定不能 (中間・成分の不一致)"
j = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.46", "cap50": t50, "cap1": t1, "qualification": qual, "island": summary, "result": res}
(OUT / "resp_judge.json").write_text(json.dumps(j, ensure_ascii=False, indent=1))
print(json.dumps(qual, ensure_ascii=False))
for t, s in summary.items():
    print(t, {c: {k: round(v, 4) for k, v in d.items()} for c, d in s.items()})
print("結果:", res)
