"""plan tooling-nozzle-core-grid §4.34 の D1 の判定 (2 本の比較)。
usage: python3 cg_b0_diag_cmp.py <上限 50 の run> <上限 1 の run> → _band_ab/core_grid/d1_judge.json
結果 A: 上限 50 の ρ・ρE のそれぞれで、持続する交互更新の節点が残差二乗和の過半を担い、上限 1 では同じ節点の更新の振幅が 1/3 以下
       → 「上限に敏感な step ごとの交互更新が主因」を支持 (平均のずれの全量を説明したとはしない)。
結果 B: 残差を担う節点 (上位 0.1 %) に交互更新が無く (lag-1 の中央値 > −0.8)、提案の増分が実際の増分で失われる (ずれ ≥ 10 %) → 第 3 仮説 (床・上書き) を支持。
それ以外 → 判定不能 (混在・閾値の中間・窓で結果が変わる)。"""
import json, sys
from pathlib import Path
import numpy as np
OUTD = Path(__file__).resolve().parent / "_band_ab" / "core_grid"
A, B = sys.argv[1:3]
ja = json.loads((OUTD / f"d1_{A}.json").read_text()); jb = json.loads((OUTD / f"d1_{B}.json").read_text())
na = np.load(OUTD / f"d1_{A}.npz"); nb = np.load(OUTD / f"d1_{B}.npz")
out = {"A": A, "B": B, "vars": {}}
for q in ("ro", "roUx", "roUy", "roe", "roK", "roOmega"):
    alt = na[f"{q}_alt"]
    ampA = np.sqrt((na[f"{q}_amp1"].astype(float) ** 2 + na[f"{q}_amp2"].astype(float) ** 2) / 2)
    ampB = np.sqrt((nb[f"{q}_amp1"].astype(float) ** 2 + nb[f"{q}_amp2"].astype(float) ** 2) / 2)
    r = float(np.sqrt((ampB[alt] ** 2).mean() / (ampA[alt] ** 2).mean())) if alt.any() and (ampA[alt] ** 2).mean() > 0 else None
    out["vars"][q] = {"A_alt_nodes": int(alt.sum()), "A_alt_residual_share": ja["vars"][q]["alt_residual_share"],
                      "B_over_A_amp_at_A_alt": r, "A_top_lag1": ja["vars"][q]["lag1_median_top_residual"], "B_top_lag1": jb["vars"][q]["lag1_median_top_residual"],
                      "A_lost": ja["vars"][q]["proposed_vs_actual"], "B_lost": jb["vars"][q]["proposed_vs_actual"]}
def resA(q):
    v = out["vars"][q]; return (v["A_alt_residual_share"] or 0) > 0.5 and v["B_over_A_amp_at_A_alt"] is not None and v["B_over_A_amp_at_A_alt"] <= 1 / 3
def resB(q):
    v = out["vars"][q]; top = v["A_top_lag1"]; lost = v["A_lost"]["interior"]["lost_frac"] or 0
    return top is not None and min(top) > -0.8 and lost >= 0.1
out["result"] = ("A: 上限に敏感な step ごとの交互更新が主因を支持 (ρ・ρE)" if resA("ro") and resA("roe") else
                 "B: 交互更新は無く、提案の増分が失われる (床・上書きの仮説を支持)" if resB("ro") and resB("roe") else "判定不能 (混在・中間)")
(OUTD / "d1_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for q, v in out["vars"].items():
    print(f"{q:7s} A の交互 {v['A_alt_nodes']} 残差の割合 {v['A_alt_residual_share']} B/A の振幅 {v['B_over_A_amp_at_A_alt']} 上位の lag1 A {v['A_top_lag1']} B {v['B_top_lag1']} "
          f"提案のずれ A 内部 {v['A_lost']['interior']['lost_frac']} B 内部 {v['B_lost']['interior']['lost_frac']}")
print("結果:", out["result"])
