"""plan tooling-nozzle-core-grid §4.84 P-m2 の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 point-precond2): point の lowMachPrecond 2 (P2) 対 既存の C (0)。
縮流部 x/r_t ∈ [−8, −2] の q (各腕の自分の演算子)。窓 W1 = 10k〜15k・W2 = 15k〜20k の出力の平均 q̄ (各 3 点)。d = q̄(W2)/q̄(W1)、H = d_2/d_C。
結果 A: ρ・ρE の両方で d_2 ≤ 0.5 かつ H ≤ 0.5、P2 の log q の傾き (各窓 3 点の最小二乗) が両窓で負、
  縮流部の ρUx・ρUy・ρk・ρω と縮流部の外の ρ・ρUx・ρUy・ρE・ρk・ρω で d_2 ≤ 1 かつ d_2/d_C ≤ 1.5。
結果 B: ρ・ρE の両方で H ≥ 0.8、または両方で d_2 ≥ 1。
P2 に非有限 → この構成が不安定。RUN_RC ≠ 0 で非有限なし → 実行失敗。窓の点の不足 → 証拠不足。それ以外・資格不足 → 判定不能。
水準比 q̄_2/q̄_C は記述だけ (主判定に使わない)。資格は cg_pp2_check.py の 5 項 + 全出力の非有限なし・ρ・P・T > 0。
usage: python3 cg_pp2_judge.py <P2 の run> [C の run (既定 run_0688_pp_c)] → _band_ab/core_grid/pp2_judge.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
P2 = sys.argv[1]; C = sys.argv[2] if len(sys.argv) > 2 else "run_0688_pp_c"; W = ((10000, 15000), (15000, 20000))
def jl(name):
    p = OUTD / name
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"]) if p.exists() else []
RC, RP = jl(f"region_res_{C}.jsonl"), jl(f"region_res_{P2}.jsonl"); PC, PP = jl(f"phys_{C}.jsonl"), jl(f"phys_{P2}.jsonl")
def pts(rows, key, lo, hi): return [(r["step"], r[key]) for r in rows if lo <= r["step"] <= hi and key in r]
def wm(rows, key, lo, hi):
    v = pts(rows, key, lo, hi); return float(np.mean([x for _, x in v])) if len(v) == 3 else None
def slope(rows, key, lo, hi):
    v = pts(rows, key, lo, hi)
    if len(v) != 3 or any(x <= 0 for _, x in v): return None
    s = np.array([a for a, _ in v], float); y = np.log([b for _, b in v]); return float(np.polyfit(s, y, 1)[0])
rc_txt = (HERE / P2 / "RUN_RC").read_text().strip() if (HERE / P2 / "RUN_RC").exists() else None
nonfinite_P2 = any(r["nonfinite"] for r in RP) or any(r["nonfinite"] for r in PP)
chk = json.loads((OUTD / "pp2_check.json").read_text()) if (OUTD / "pp2_check.json").exists() else {"qualification": {"check_missing": False}}
qual = dict(chk["qualification"])
qual["finite_and_physical"] = (not nonfinite_P2) and (not any(r["nonfinite"] for r in RC + PC)) and all(
    x.get("min_ro", 1) > 0 and x.get("min_P", 1) > 0 and x.get("min_T", 1) > 0 for x in PC + PP)
main = [("contraction", "res_ro"), ("contraction", "res_roe")]
side = [("contraction", c) for c in ("res_roUx", "res_roUy", "res_roK", "res_roOmega")] + [("outside", c) for c in ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega")]
tab = {}
for reg, c in main + side:
    k = f"{reg}_{c}"
    w = {arm: [wm(R, k, *W[0]), wm(R, k, *W[1])] for arm, R in (("C", RC), ("P2", RP))}
    d = {arm: (w[arm][1] / w[arm][0] if None not in w[arm] and w[arm][0] > 0 else None) for arm in w}
    tab[k] = {"qbar_W1_W2": w, "d": d, "H": (d["P2"] / d["C"] if d["P2"] is not None and d["C"] else None),
              "level_ratio_P2_over_C_W2": (w["P2"][1] / w["C"][1] if w["P2"][1] is not None and w["C"][1] else None),
              "slope_P2": [slope(RP, k, *W[0]), slope(RP, k, *W[1])]}
missing = any(v is None for k in tab for v in (tab[k]["d"]["P2"], tab[k]["d"]["C"], tab[k]["H"]))
if nonfinite_P2:
    res = "P2 に非有限: この構成が不安定 (安定なときの減衰率は未判定)"
elif rc_txt != "0":
    res = f"判定不能 (実行失敗: RUN_RC {rc_txt}、非有限なし)"
elif missing:
    res = "判定不能 (証拠不足: 窓の点が足りない、またはゼロ分母)"
elif not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
else:
    m = [tab[f"{r}_{c}"] for r, c in main]
    slopes_ok = all(s is not None and s < 0 for t in m for s in t["slope_P2"])
    side_ok = all(tab[f"{r}_{c}"]["d"]["P2"] <= 1 and tab[f"{r}_{c}"]["d"]["P2"] / tab[f"{r}_{c}"]["d"]["C"] <= 1.5 for r, c in side)
    if all(t["d"]["P2"] <= 0.5 and t["H"] <= 0.5 for t in m) and slopes_ok and side_ok:
        res = "A: 変更後の方程式に対する有限区間の残差の減衰が C より速い (「=2 でも有限予算内の改善が得られない」を退ける)"
    elif all(t["H"] >= 0.8 for t in m) or all(t["d"]["P2"] >= 1 for t in m):
        res = "B: この設定・予算では =2 による減衰の改善はない (第 1 仮説を退ける)"
    else:
        res = "判定不能 (成分で割れる・中間・悪化の移動)"
series = {arm: [{"step": r["step"], "ro": r.get("contraction_res_ro"), "roe": r.get("contraction_res_roe")} for r in R] for arm, R in (("C", RC), ("P2", RP))}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.84", "C": C, "P2": P2, "qualification": qual, "table": tab, "series": series, "check": chk, "result": res}
(OUTD / "pp2_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
fmt = lambda v: "-" if v is None else f"{v:.3g}"
for k, t in tab.items():
    print(f"{k:24s} d_C {fmt(t['d']['C'])} d_P2 {fmt(t['d']['P2'])} H {fmt(t['H'])} 水準比(W2) {fmt(t['level_ratio_P2_over_C_W2'])} 傾き P2 {[fmt(s) for s in t['slope_P2']]}")
for arm, s in series.items():
    print(arm, " ".join(f"{x['step']}:{fmt(x['ro'])}" for x in s))
print("資格:", qual); print("結果:", res)
