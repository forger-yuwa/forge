"""plan tooling-nozzle-core-grid §4.76 の判定 (事前登録、codex 2026-10-13 after-small-amp)。
入力: _band_ab/core_grid/twin_rzp_noise.jsonl (R 対 Z = ビット一致の複製)・twin_rzp_sig.jsonl (R 対 P = 摂動 1e−7)。
主指標: 固定領域 x/r_t ∈ [−9.5, −2] (内部の節点) の体積重み付きの温度差の RMS。N(n) = ‖Z − R‖、S(n) = ‖P − R‖ (ノルムどうしを引かない)。
主判定区間 step 1500〜2500 の 5 点。資格: 非有限なし、両差の max|ΔT/T| < 1e−3。
結果 A: 全点で N/S ≥ 0.2 かつ N が区間内で増大 (N(2500) > N(1500)) → 「初期摂動への応答を十分分離できている」を退ける。
結果 B: 全点で N/S ≤ 0.05 → 1e−7 の測定に大きな実行間差が混ざるという仮説を退ける (1e−10 側の雑音支配までは否定しない)。
それ以外・資格不足は判定不能。ρ・ρE の比も併記する。
usage: python3 cg_twin_rzp_judge.py → _band_ab/core_grid/twin_rzp_judge.json"""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parent / "_band_ab" / "core_grid"
rd = lambda t: {r["step"]: r for r in (json.loads(l) for l in open(OUT / f"twin_{t}.jsonl"))}
Nn, Ss = rd("rzp_noise"), rd("rzp_sig")
steps = sorted(set(Nn) & set(Ss))
win = [s for s in steps if 1500 <= s <= 2500]
ser = [{"step": s, **{f"{c}_{k}": r[f"fixed95_D_{k}"] for c, r in (("N", Nn[s]), ("S", Ss[s])) for k in ("T", "ro", "roe")},
        "N_over_S_T": Nn[s]["fixed95_D_T"] / Ss[s]["fixed95_D_T"] if Ss[s]["fixed95_D_T"] > 0 else None,
        "N_over_S_ro": Nn[s]["fixed95_D_ro"] / Ss[s]["fixed95_D_ro"] if Ss[s]["fixed95_D_ro"] > 0 else None,
        "N_over_S_roe": Nn[s]["fixed95_D_roe"] / Ss[s]["fixed95_D_roe"] if Ss[s]["fixed95_D_roe"] > 0 else None,
        "N_max_rel": Nn[s]["max_dT_rel"], "S_max_rel": Ss[s]["max_dT_rel"], "nonfinite": Nn[s]["nonfinite"] + Ss[s]["nonfinite"]} for s in steps]
W = [r for r in ser if r["step"] in win]
qual = {"points": len(W) == 5, "nonfinite": all(r["nonfinite"] == 0 for r in ser), "linear": all(r["N_max_rel"] < 1e-3 and r["S_max_rel"] < 1e-3 for r in W)}
rat = [r["N_over_S_T"] for r in W]
if not all(qual.values()) or any(v is None for v in rat):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif all(v >= 0.2 for v in rat) and W[-1]["N_T"] > W[0]["N_T"]:
    res = "A: 実行間差 N が摂動への応答 S の 0.2 倍以上で、N も育つ (初期摂動への応答を十分分離できていない)"
elif all(v <= 0.05 for v in rat):
    res = "B: 実行間差 N は摂動への応答 S の 5 % 以下 (1e−7 の測定は実行間差に大きく汚されていない)"
else:
    res = "判定不能 (中間)"
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.76", "qualification": qual, "series": ser, "result": res}
(OUT / "twin_rzp_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for r in ser:
    print(f"  {r['step']:5d}  N_T {r['N_T']:.3e}  S_T {r['S_T']:.3e}  N/S T {r['N_over_S_T'] if r['N_over_S_T'] is None else round(r['N_over_S_T'], 4)}  ρ {r['N_over_S_ro'] if r['N_over_S_ro'] is None else round(r['N_over_S_ro'], 4)}  ρE {r['N_over_S_roe'] if r['N_over_S_roe'] is None else round(r['N_over_S_roe'], 4)}")
print("資格:", qual, " 結果:", res)
