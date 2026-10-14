"""plan tooling-nozzle-core-grid §4.38 F1 の判定: SST を凍結した B0 (run_0626) 対 凍結しない B0 (run_0619 の 0〜20k)。
窓 W1 = 10k〜15k・W2 = 15k〜20k。全域: residual_history (outer_end) の各列の中央値 (rms_ro・roUx・roUy・roe・roY0・roY1)。
領域: 縮流部 x∈[−8, −2] の q = sqrt(ΣV(R/V)²/ΣV) の窓の出力の平均 (ρ・ρE)。
凍結の実効性: 内部の節点で ρk・ρω の変化の RMS 比 < 1e-10 (無ければ判定不能)。μt の変化と dt の比は記録。
結果 A: 両窓で全域の全列と領域の ρ・ρE がすべて凍結/対照 ≤ 1/3 → SST の更新への強い依存を支持 (μt の遅れや x の格子依存の機構は確定しない)。
結果 B: 両窓で領域の ρ・ρE がともに ≥ 2/3、かつ縮流部が全域の残差を担う割合が対照と同程度 (比 0.5〜2) → SST の更新なしでも停滞が続く (第 2 仮説を優先)。
それ以外 → 判定不能。非有限は別扱い。延長しない。
usage: python3 cg_freeze_judge.py <凍結の run> <対照の run> → _band_ab/core_grid/freeze_judge.json"""
import csv, json, math, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
F, C = sys.argv[1:3]; W = ((10000, 15000), (15000, 20000))
COLS = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roY0", "rms_roY1"); SSTC = ("rms_roK", "rms_roOmega")
def hist(run):
    rows = [x for x in csv.DictReader(open(HERE / run / "residual_history.csv")) if x["phase"] == "outer_end"]
    return {k: np.array([float(x[k]) for x in rows]) for k in ("step",) + COLS + SSTC}
def region(run):
    rows = [json.loads(l) for l in open(OUTD / f"region_res_{run}.jsonl")]
    return sorted({r["step"]: r for r in rows}.values(), key=lambda r: r["step"])
HF, HC = hist(F), hist(C); RF, RC_ = region(F), region(C)
nonfinite = any(not np.all(np.isfinite(H[k])) for H in (HF, HC) for k in COLS) or any(r["nonfinite"] for r in RF)
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.38 F1", "frozen": F, "control": C, "global": {}, "region": {}}
allA = True
for lo, hi in W:
    g = {}
    for k in COLS + SSTC:
        mf = float(np.median(HF[k][(HF["step"] >= lo) & (HF["step"] < hi)])); mc = float(np.median(HC[k][(HC["step"] >= lo) & (HC["step"] < hi)]))
        g[k] = {"frozen": mf, "control": mc, "ratio": mf / mc if mc > 0 else None}
        if k in COLS: allA &= g[k]["ratio"] is not None and g[k]["ratio"] <= 1 / 3
    out["global"][f"{lo}-{hi}"] = g
regA, regB, share = True, True, {}
for lo, hi in W:
    rr = {}
    for comp in ("res_ro", "res_roe"):
        f = [r[f"contraction_{comp}"] for r in RF if lo <= r["step"] <= hi]; c = [r[f"contraction_{comp}"] for r in RC_ if lo <= r["step"] <= hi]
        fa = [r[f"all_{comp}"] for r in RF if lo <= r["step"] <= hi]; ca = [r[f"all_{comp}"] for r in RC_ if lo <= r["step"] <= hi]
        ratio = float(np.mean(f) / np.mean(c)); sh = float((np.mean(f) / np.mean(fa)) / (np.mean(c) / np.mean(ca)))
        rr[comp] = {"frozen": float(np.mean(f)), "control": float(np.mean(c)), "ratio": ratio, "n": [len(f), len(c)], "contraction_share_ratio": sh}
        regA &= ratio <= 1 / 3; regB &= ratio >= 2 / 3 and 0.5 <= sh <= 2
    out["region"][f"{lo}-{hi}"] = rr
fc = [json.loads(l) for l in open(OUTD / f"freeze_check_{F}.jsonl")]
eff = all(x["roK"]["rms_rel"] < 1e-10 and x["roOmega"]["rms_rel"] < 1e-10 for x in fc)
out["freeze_check"] = {"effective": eff, "roK_rms_rel_max": max(x["roK"]["rms_rel"] for x in fc), "roOmega_rms_rel_max": max(x["roOmega"]["rms_rel"] for x in fc),
                       "vis_turb_rms_rel": [x["vis_turb"]["rms_rel"] for x in fc], "dt_ratio_median": [x.get("dt_ratio_median") for x in fc]}
out["result"] = ("発散 (非有限)" if nonfinite else "判定不能 (凍結が効いていない)" if not eff else
                 "A: SST の更新を止めると停滞が大幅に緩和 (SST の更新への強い依存を支持。μt の遅れ・x の格子依存の機構は確定しない)" if allA and regA else
                 "B: SST の更新を止めても停滞が続く (平均流の側の第 2 仮説を優先)" if regB else "判定不能 (中間・成分で割れる)")
(OUTD / "freeze_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for w, g in out["global"].items():
    print(f"全域 {w}: " + "  ".join(f"{k[4:]} {v['ratio']:.3g}" for k, v in g.items()))
for w, rr in out["region"].items():
    print(f"縮流部 {w}: " + "  ".join(f"{k} 比 {v['ratio']:.3g} (縮流部の割合の比 {v['contraction_share_ratio']:.2f})" for k, v in rr.items()))
print("凍結の実効性:", {k: v for k, v in out["freeze_check"].items() if k != "dt_ratio_median"})
print("結果:", out["result"])
