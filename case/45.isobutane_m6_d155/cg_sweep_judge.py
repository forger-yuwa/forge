"""plan tooling-nozzle-core-grid §4.42 M-d の判定 (事前登録、codex 2026-10-11 sweep)。
入力: cg_sweep_an.py の sweep_<tag50>.json (上限 50) と sweep_<tag1>.json (上限 1)。同じ停滞場から 1 外反復・41 sweep。
資格: (1) 両腕の rhs_s0・状態・節点が一致し非有限 0、(2) 各 sweep のライン solve の後退誤差の最大 ≤ 1e−8、
      (3) 両腕・島とライン全体・ρ/ρE/全 5 成分で ‖g₄₀‖/‖b₀‖ ≤ 1e−3 かつ ‖δ₄₀ − δ₃₀‖/‖δ₄₀‖ ≤ 1e−3。
結果 A: 資格を満たし、島の ρ・ρE の両方で 上限 50 の u5 ≥ 0.3 かつ 上限 1 の u5 ≤ 上限 50 の 1/3
        → 「5 sweep の不足が大きく、上限操作で減る」(第 1 仮説の中間機構) を支持。外反復の停滞原因・格子粗化の機構は確定しない。
結果 B: 資格を満たし、上限 50 の島で ρ・ρE・全 5 成分の u5 ≤ 0.1 → 「大きな 5→40 sweep の解き残し」を退ける。H2′ の確定ではない。
それ以外・資格不足 → 判定不能。
usage: python3 cg_sweep_judge.py <tag50> <tag1> → _band_ab/core_grid/sweep_judge.json"""
import json, sys
from pathlib import Path
OUT = Path(__file__).resolve().parent / "_band_ab" / "core_grid"
t50, t1 = sys.argv[1], sys.argv[2]
A, B = (json.loads((OUT / f"sweep_{t}.json").read_text()) for t in (t50, t1))
q = {}
q["same_input"] = all(A[k] == B[k] for k in ("rhs_s0_sha256", "state_sha256", "nodes_sha256")) and A["nonfinite"] == 0 and B["nonfinite"] == 0
q["backward_error"] = max(A["backward_error_max"], B["backward_error_max"]) <= 1e-8
q["solved_40"] = all(X["region"][r][c]["g40_over_b0"] <= 1e-3 and X["region"][r][c]["conv_40_30"] <= 1e-3
                     for X in (A, B) for r in ("island", "lines") for c in ("ro", "roe", "all5"))
qual = all(q.values())
u50 = {c: A["region"]["island"][c]["u5"] for c in ("ro", "roe", "all5")}
u1 = {c: B["region"]["island"][c]["u5"] for c in ("ro", "roe", "all5")}
if not qual:
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in q.items() if not v) + ")"
elif all(u50[c] >= 0.3 and u1[c] <= u50[c] / 3 for c in ("ro", "roe")):
    res = "A: 上限 50 の 5 sweep の解き残しが大きく、上限 1 で減る (第 1 仮説の中間機構を支持)"
elif all(u50[c] <= 0.1 for c in ("ro", "roe", "all5")):
    res = "B: 上限 50 でも 5 sweep の解き残しは小さい (大きな解き残しを退ける)"
else:
    res = "判定不能 (中間・成分の不一致)"
    if all(u50[c] >= 0.3 for c in ("ro", "roe")):
        res += "。上限 50 の解き残しは大きいが上限 1 で 1/3 まで減らない → 上限の効果を内部反復不足の改善で説明する接続は支持されない"
j = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.42", "cap50": t50, "cap1": t1, "qualification": q,
     "island_u5_cap50": u50, "island_u5_cap1": u1, "result": res}
(OUT / "sweep_judge.json").write_text(json.dumps(j, ensure_ascii=False, indent=1))
print(json.dumps(q, ensure_ascii=False)); print("島の u5 上限 50:", {k: round(v, 3) for k, v in u50.items()}, " 上限 1:", {k: round(v, 3) for k, v in u1.items()}); print("結果:", res)
