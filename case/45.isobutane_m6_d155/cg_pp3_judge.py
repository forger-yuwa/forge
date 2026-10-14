"""plan tooling-nozzle-core-grid §4.82 P-m3 の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 point-precond): point の lowMachPrecond 0 (C) 対 3 (P)。
主指標: 縮流部 x/r_t ∈ [−8, −2] の q = √(ΣV(R/V)²/ΣV) の P/C、窓 W1 = 10k〜15k・W2 = 15k〜20k の出力の平均 (各 3 点)。
資格: cg_pp3_check.py の 4 項 (前処理の効き・同じ起点と右辺・境界節点の凍結なし・組成の一様) + 全出力で非有限なし・ρ・P・T > 0。C の停滞は資格にしない。
結果 A: 両窓で P/C が ρ・ρE とも ≤ 1/3、かつ両窓で縮流部の ρUx・ρUy・ρk・ρω と縮流部の外の ρ・ρUx・ρUy・ρE・ρk・ρω の P/C ≤ 1.5 (C の窓平均が 0 の成分は除く)。
結果 B: 資格を満たし、両窓で P/C が ρ・ρE とも ≥ 0.8。P だけが発散 → 「この構成が不安定」(加速は未判定)。それ以外は判定不能。
usage: python3 cg_pp3_judge.py <C の run> <P の run> → _band_ab/core_grid/pp3_judge.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
C, P = sys.argv[1:3]; W = ((10000, 15000), (15000, 20000))
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
RC, RP = jl(f"region_res_{C}.jsonl"), jl(f"region_res_{P}.jsonl"); PC, PP = jl(f"phys_{C}.jsonl"), jl(f"phys_{P}.jsonl")
def wm(rows, key, lo, hi):
    v = [r[key] for r in rows if lo <= r["step"] <= hi and key in r]
    return float(np.mean(v)) if len(v) == 3 else None
rc_ = lambda r: (HERE / r / "RUN_RC").read_text().strip() if (HERE / r / "RUN_RC").exists() else None
chk = json.loads((OUTD / "pp3_check.json").read_text()) if (OUTD / "pp3_check.json").exists() else {"qualification": {"check_missing": False}}
fin = {r: (not any(x["nonfinite"] for x in R)) and (not any(x["nonfinite"] for x in Ph)) and rc_(r) == "0" for r, R, Ph in ((C, RC, PC), (P, RP, PP))}
phys = {r: all(x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for x in Ph) for r, Ph in ((C, PC), (P, PP))}
qual = dict(chk["qualification"]); qual["control_finite"] = fin[C]; qual["physical"] = phys[C] and phys[P]
main = ("res_ro", "res_roe")
side = [("contraction", c) for c in ("res_roUx", "res_roUy", "res_roK", "res_roOmega")] + [("outside", c) for c in ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega")]
ratios = {}
for lo, hi in W:
    rr = {}
    for reg, c in [("contraction", m) for m in main] + side:
        a, b = wm(RP, f"{reg}_{c}", lo, hi), wm(RC, f"{reg}_{c}", lo, hi)
        rr[f"{reg}_{c}"] = None if (a is None or b is None) else ("C=0" if b == 0 else a / b)
    ratios[f"{lo}-{hi}"] = rr
def first_last(R, key):
    return (R[0][key], R[-1][key]) if R else (None, None)
drop = {r: {c: first_last(R, f"contraction_{c}") for c in main} for r, R in ((C, RC), (P, RP))}
if not fin[P] and fin[C]:
    res = "P だけが発散: この構成が不安定 (安定なときの加速は未判定)"
elif not all(qual.values()) or any(v is None for rr in ratios.values() for v in rr.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + (", 窓の点が足りない" if any(v is None for rr in ratios.values() for v in rr.values()) else "") + ")"
else:
    mainr = [rr[f"contraction_{c}"] for rr in ratios.values() for c in main]
    sider = [v for rr in ratios.values() for k, v in rr.items() if k not in ("contraction_res_ro", "contraction_res_roe") and v != "C=0"]
    if all(v <= 1 / 3 for v in mainr) and all(v <= 1.5 for v in sider):
        res = "A: point で lowMachPrecond 3 が縮流部の平均流の残差を固定予算内で大きく下げ、悪化を他へ移していない"
    elif all(v >= 0.8 for v in mainr):
        res = "B: この設定・予算では lowMachPrecond 3 による大きな改善はない"
    else:
        res = "判定不能 (中間・成分で割れる・悪化の移動)"
series = {r: [{"step": x["step"], "ro": x.get("contraction_res_ro"), "roe": x.get("contraction_res_roe"), "roOmega_all": x.get("all_res_roOmega")} for x in R] for r, R in ((C, RC), (P, RP))}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.82", "C": C, "P": P, "qualification": qual, "ratios_P_over_C": ratios,
       "first_last_contraction": drop, "series": series, "check": chk, "result": res}
(OUTD / "pp3_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for w, rr in ratios.items():
    print(f"{w}: " + "  ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}" for k, v in rr.items()))
for r in (C, P):
    print(r, "縮流部の最初と最後:", {c: (f"{a:.3e}" if a is not None else None, f"{b:.3e}" if b is not None else None) for c, (a, b) in drop[r].items()})
print("資格:", qual); print("結果:", res)
