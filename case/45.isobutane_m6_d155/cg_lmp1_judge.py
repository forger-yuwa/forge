"""plan tooling-nozzle-core-grid §4.80 L-m1 の判定 (事前登録、codex 諮問 2026-10-13 lowmach): 低マッハの散逸 (lowMachPrecond 1) だけを変えた A/B。
窓 W1 = 10k〜15k・W2 = 15k〜20k の出力の平均 (各 3 点)。主指標: 縮流部 x/r_t ∈ [−8, −2] の q = √(ΣV(R/V)²/ΣV) の L/R (ρ・ρE)。
資格: (1) R の ρ・ρE とも W2/W1 ≥ 1/2 (対照が停滞)。(2) R/run_0619 が ρ・ρE・両窓とも 0.5〜2 (対照が再現)。
  (3) forge_run.log の 'lowMachPrecond' in 'time.deltaT' が L で 1、R で 0。(4) R に非有限なし。
結果 A: 両窓で L/R が ρ・ρE とも ≤ 1/3、かつ両窓で縮流部の ρUx・ρUy・ρk・ρω と縮流部の外の ρ・ρUx・ρUy・ρE の L/R がすべて ≤ 1.5。
結果 B: 両窓で L/R が ρ・ρE とも ≥ 0.8、または L が発散 (非有限)。それ以外・資格不足は判定不能。
記述: P の a/q・η・f (cg_lmp1_saw.py、10k・15k・20k)、8 量の L − R (準定常の判定なし)。
usage: python3 cg_lmp1_judge.py <R の run> <L の run> → _band_ab/core_grid/lmp1_judge.json"""
import csv, json, re, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
R, L = sys.argv[1:3]; CTL = "run_0619_cg_gref_b0_ctl"; W = ((10000, 15000), (15000, 20000))
def region(run):
    p = OUTD / f"region_res_{run}.jsonl"
    if not p.exists(): return []
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"])
def wm(rows, key, lo, hi):
    v = [r[key] for r in rows if lo <= r["step"] <= hi and key in r]
    return float(np.mean(v)) if len(v) == 3 else None
def lmp(run):
    m = re.search(r"'lowMachPrecond' in 'time\.deltaT': (\d+)", (HERE / run / "forge_run.log").read_text(errors="replace"))
    return int(m.group(1)) if m else None
RR, RL, RC = region(R), region(L), region(CTL)
nonfinite_L = any(r["nonfinite"] for r in RL) or (HERE / L / "RUN_RC").read_text().strip() != "0"
main = ("res_ro", "res_roe")
qual = {"control_stall": all((wm(RR, f"contraction_{c}", *W[1]) or 0) / (wm(RR, f"contraction_{c}", *W[0]) or np.inf) >= 0.5 for c in main),
        "control_reproduces_0619": all(wm(RR, f"contraction_{c}", *w) is not None and wm(RC, f"contraction_{c}", *w) is not None
                                       and 0.5 <= wm(RR, f"contraction_{c}", *w) / wm(RC, f"contraction_{c}", *w) <= 2 for c in main for w in W),
        "lmp_effective": lmp(L) == 1 and lmp(R) == 0, "control_finite": not any(r["nonfinite"] for r in RR)}
ratios = {}
for lo, hi in W:
    rr = {}
    for reg, cs in (("contraction", ("res_ro", "res_roe", "res_roUx", "res_roUy", "res_roK", "res_roOmega")), ("outside", ("res_ro", "res_roUx", "res_roUy", "res_roe"))):
        for c in cs:
            a, b = wm(RL, f"{reg}_{c}", lo, hi), wm(RR, f"{reg}_{c}", lo, hi)
            rr[f"{reg}_{c}"] = a / b if (a is not None and b) else None
    ratios[f"{lo}-{hi}"] = rr
if nonfinite_L:
    res = "B: L が発散 (非有限)。「`=1` だけで緩和できる」を退ける (低マッハの機構そのものは退けない)"
elif not all(qual.values()) or any(v is None for rr in ratios.values() for v in rr.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + (", 窓の点が足りない" if any(v is None for rr in ratios.values() for v in rr.values()) else "") + ")"
else:
    mainr = [rr[f"contraction_{c}"] for rr in ratios.values() for c in main]
    side = [v for rr in ratios.values() for k, v in rr.items() if k not in ("contraction_res_ro", "contraction_res_roe")]
    if all(v <= 1 / 3 for v in mainr) and all(v <= 1.5 for v in side):
        res = "A: `=1` で縮流部の停滞が大きく緩和し、悪化を他へ移していない (第 1 仮説を支持)"
    elif all(v >= 0.8 for v in mainr):
        res = "B: `=1` を足しても縮流部の停滞が続く (「`=1` だけで緩和できる」を退ける)"
    else:
        res = "判定不能 (中間・指標で割れる)"
def saw(run):
    p = OUTD / f"lmp1_saw_{run}.jsonl"
    return {r["step"]: r for r in (json.loads(l) for l in open(p))} if p.exists() else {}
SR, SL = saw(R), saw(L)
sawtab = {str(n): {d: {"R": {k: SR[n][d][k] for k in ("aP_over_q", "eta_P", "f_P")}, "L": {k: SL[n][d][k] for k in ("aP_over_q", "eta_P", "f_P")}}
                   for d in ("wallnormal", "streamwise")} for n in sorted(set(SR) & set(SL))}
def series(run, n):
    p = OUTD / f"gref_series_{run}_{n}.csv"
    if not p.exists(): return None
    rows = list(csv.DictReader(open(p)))
    return rows[-1] if rows else None
qty = {}
for n in range(2500, 20001, 2500):
    a, b = series(L, n), series(R, n)
    if a and b:
        qty[str(n)] = {k: (float(a[k]) - float(b[k])) / abs(float(b[k])) if float(b[k]) != 0 else None
                       for k in a if k not in ("step", "nonfinite", "run", "tag") and re.fullmatch(r"-?[0-9.eE+-]+", a[k] or "x") and re.fullmatch(r"-?[0-9.eE+-]+", b[k] or "x")}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.80", "R": R, "L": L, "qualification": qual, "ratios_L_over_R": ratios,
       "result": res, "saw_P": sawtab, "quantities_rel_diff_L_minus_R": qty, "lowMachPrecond": {"R": lmp(R), "L": lmp(L)}}
(OUTD / "lmp1_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for w, rr in ratios.items():
    print(f"{w}: " + "  ".join(f"{k} {v:.3g}" if v is not None else f"{k} -" for k, v in rr.items()))
for n, t in sawtab.items():
    print(f"P の上下 {n}: " + "  ".join(f"{d} a/q R {t[d]['R']['aP_over_q']:.2e} L {t[d]['L']['aP_over_q']:.2e} η R {t[d]['R']['eta_P']:.3g} L {t[d]['L']['eta_P']:.3g}" for d in t))
print("資格:", qual); print("結果:", res)
