"""plan tooling-nozzle-core-grid §4.40 M-a の判定 (一般形: 変更側の run 対 対照 run_0619 の 0〜20k、縮流部の残差)。
窓 W1 = 10k〜15k・W2 = 15k〜20k の出力の平均で、縮流部 x∈[−8, −2] の q = sqrt(ΣV(R/V)²/ΣV) (ρ・ρE) の変更側/対照。
対照の停滞 = ρ・ρE とも W2/W1 ≥ 1/2。結果 A: 対照が停滞し両窓・両成分で ≤ 1/3。結果 B: 両窓・両成分で ≥ 2/3 かつ縮流部が全域の残差を担う割合の比 0.5〜2 (集中の場所が維持)。
それ以外 (中間・成分で割れる・対照の停滞が不成立) → 判定不能。非有限 → 発散。SST (ρk・ρω の縮流部の残差) が変更側で大きければ「平均流のみ改善、連成系は未解決」と付ける。延長しない。
usage: python3 cg_ab_region_judge.py <変更側の run> <対照の run> <出力名> [登録の節] → _band_ab/core_grid/<出力名>.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
X, C, name = sys.argv[1:4]; W = ((10000, 15000), (15000, 20000))
def region(run):
    rows = [json.loads(l) for l in open(OUTD / f"region_res_{run}.jsonl")]
    return sorted({r["step"]: r for r in rows}.values(), key=lambda r: r["step"])
RX, RC = region(X), region(C)
def wm(R, key, lo, hi): return float(np.mean([r[key] for r in R if lo <= r["step"] <= hi]))
out = {"registration": "plans/active/tooling-nozzle-core-grid.md " + (sys.argv[4] if len(sys.argv) > 4 else "§4.40"), "changed": X, "control": C, "region": {}}
nonfinite = any(r["nonfinite"] for r in RX)
stall = all(wm(RC, f"contraction_{c}", *W[1]) / wm(RC, f"contraction_{c}", *W[0]) >= 0.5 for c in ("res_ro", "res_roe"))
A, B, sst_worse = True, True, False
for lo, hi in W:
    rr = {}
    for c in ("res_ro", "res_roe", "res_roUx", "res_roUy", "res_roK", "res_roOmega"):
        r = wm(RX, f"contraction_{c}", lo, hi) / wm(RC, f"contraction_{c}", lo, hi)
        sh = (wm(RX, f"contraction_{c}", lo, hi) / wm(RX, f"all_{c}", lo, hi)) / (wm(RC, f"contraction_{c}", lo, hi) / wm(RC, f"all_{c}", lo, hi))
        rr[c] = {"ratio": r, "share_ratio": sh}
        if c in ("res_ro", "res_roe"):
            A &= r <= 1 / 3; B &= r >= 2 / 3 and 0.5 <= sh <= 2
        if c in ("res_roK", "res_roOmega") and r > 1: sst_worse = True
    out["region"][f"{lo}-{hi}"] = rr
out["control_stall"] = stall
out["result"] = ("発散 (非有限)" if nonfinite else "判定不能 (対照が停滞していない)" if not stall else
                 ("A: 変更で縮流部の平均流の停滞が大きく緩和" + ("。ただし SST の縮流部の残差は変更側が大きい (平均流のみ改善、連成系は未解決)" if sst_worse else "")) if A else
                 "B: 変更しても縮流部の停滞が続き、集中の場所も維持" if B else "判定不能 (中間・成分で割れる)")
(OUTD / f"{name}.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for w, rr in out["region"].items():
    print(f"縮流部 {w}: " + "  ".join(f"{k[4:]} {v['ratio']:.3g} (割合 {v['share_ratio']:.2f})" for k, v in rr.items()))
print("対照の停滞:", stall); print("結果:", out["result"])
